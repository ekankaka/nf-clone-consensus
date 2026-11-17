"""Tests for read grouping, reference selection, PCR weighting, and reports.

Alignment records are mocked; these tests do not assess MAFFT alignments.
"""
import csv
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'tests' / 'data'


def load_module(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / 'bin' / f'{name}.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


pmod = load_module('prepare_patient')
fmod = load_module('finalize_patient')


def csvrows(path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f))


def test_select_reference_with_different_timepoints():
    records = [
        {'sample_id':'PT1-1','timepoint':'1','n_reads':3,'ref_sha256':'abc','reference':'rA'},
        {'sample_id':'PT1-2','timepoint':'2','n_reads':9,'ref_sha256':'def','reference':'rB'},
    ]
    selected, audit = pmod.select_reference('PT1', records)
    assert selected['sample_id'] == 'PT1-2'
    assert audit['reference_selection_rule'] == 'HIGHEST_READ_TIMEPOINT_THEN_SAMPLE'


def test_reference_ties_raise():
    records = [
        {'sample_id':'PT1-1','timepoint':'1','n_reads':3,'ref_sha256':'abc','reference':'rA'},
        {'sample_id':'PT1-2','timepoint':'2','n_reads':3,'ref_sha256':'def','reference':'rB'},
    ]
    with pytest.raises(ValueError, match='tie'):
        pmod.select_reference('PT1', records)


def test_prepare_and_audit(tmp_path):
    samples = []
    for sid, tp, ref in [('PT1-1','1','refA.fa'),('PT1-2','2','refB.fa')]:
        samples.append({'sample_id':sid,'timepoint':tp,
                        'fasta':str(DATA/'fastas'/f'{sid}.fa'),
                        'integration_site':str(DATA/'sites'/f'{sid}.csv'),
                        'reference':str(DATA/'references'/ref),
                        'original_best_hiv_ref_filepath':str(DATA/'references'/ref)})
    staged = tmp_path/'staged.json'
    staged.write_text(json.dumps(samples))
    subprocess.run([sys.executable, str(ROOT/'bin/prepare_patient.py'),
                    '--patient','PT1','--samples',str(staged)],cwd=tmp_path,check=True)
    clone_info=json.loads((tmp_path/'clone_info.json').read_text())
    clone_one=next(c for c in clone_info if c['chromosome']=='chr1')
    assert clone_one['n_reads'] == 5 and clone_one['n_unique_reads'] == 4
    assert clone_one['consensus_header'] == 'PT1_chr1_100_5_4'
    assert len(clone_one['timepoints']) == 2
    assert 'NN' not in (tmp_path/'clone_inputs'/f"{clone_one['clone_id']}.fasta").read_text()
    assert 'NNAAAACCCCNN' in (tmp_path/'PT1_merged.fasta').read_text()
    assert len(csvrows(tmp_path/'PT1.unmatched.csv')) == 1
    assert len(csvrows(tmp_path/'PT1.reference_selection.csv')) == 1
    assert csvrows(tmp_path/'PT1.reference_selection.csv')[0]['selected_sample_id'] == 'PT1-2'
    assert '>PT1-1_readNS' in (tmp_path/'PT1_unintegrated.fasta').read_text()
    assert '>PT1-2_readNA' in (tmp_path/'PT1_unintegrated.fasta').read_text()
    assert len(csvrows(tmp_path/'PT1.unused_sites.csv')) == 1


def make_alignment_file(tmp_path, name, lines):
    folder=tmp_path/'aligned'
    folder.mkdir(exist_ok=True)
    path=folder/(name+'.fasta')
    with path.open('w') as f:
        for h, seq in lines:
            f.write(f'>{h}\n{seq}\n')
    return path


def test_pcr_duplicates_one_vote(tmp_path):
    info = {'clone_id':'PT1_clone000001','n_reads':5,'n_unique_reads':4,
            'reads':[
                {'internal_id':'R0000001','sample_id':'PT1-1','PCR_ID':'pcr01'},
                {'internal_id':'R0000002','sample_id':'PT1-1','PCR_ID':'pcr01'},
                {'internal_id':'R0000003','sample_id':'PT1-1','PCR_ID':'pcr02'},
                {'internal_id':'R0000004','sample_id':'PT1-2','PCR_ID':'pcr01'},
                {'internal_id':'R0000005','sample_id':'PT1-2','PCR_ID':'pcr02'},
            ]}
    aln = make_alignment_file(tmp_path, info['clone_id'], [
        ('NF_REF', 'AAAACCCCGGGGTTTT'),
        ('R0000001','AAAACCCC--------'),
        ('R0000002','AAAATCCC--------'),
        ('R0000003','AAAACCCC--------'),
        ('R0000004','----CCCCGGGGTTTT'),
        ('R0000005','----CCCCGGGGTTTT'),
    ])
    seq, qc, aligned_reads, aligned_consensus=fmod.make_clone(info, aln, .60, .60, 2)
    assert seq == 'AAAACCCCGGGGTTTT'
    assert qc['n_low_support_positions'] == 0
    assert qc['consensus_type'] == 'clone_consensus'
    assert qc['percent_positions_multi_family'] == 100.0
    assert aligned_consensus == 'AAAACCCCGGGGTTTT'
    assert len(aligned_reads) == 5


def test_uncovered_reference_is_not_filled(tmp_path):
    info={'clone_id':'PT1_clone000002','n_reads':2,'n_unique_reads':2,'reads':[
        {'internal_id':'R0000001','sample_id':'PT1-1','PCR_ID':'A'},
        {'internal_id':'R0000002','sample_id':'PT1-2','PCR_ID':'B'}]}
    aln=make_alignment_file(tmp_path,info['clone_id'],[
        ('NF_REF','AAAACCCCGGGGTTTT'),
        ('R0000001','AAAA------------'),
        ('R0000002','------------TTTT')])
    seq,qc,aligned_reads,aligned_consensus=fmod.make_clone(info,aln,.60,.60,2)
    assert seq == 'N'*16
    assert qc['n_internal_uncovered_positions'] == 8
    assert qc['qc_status'] == 'REVIEW'
    assert len(aligned_consensus) == len(aligned_reads[0][1])


def test_supported_internal_deletion(tmp_path):
    info={'clone_id':'PT1_clone000003','n_reads':2,'n_unique_reads':2,'reads':[
        {'internal_id':'R0000001','sample_id':'PT1-1','PCR_ID':'A'},
        {'internal_id':'R0000002','sample_id':'PT1-2','PCR_ID':'B'}]}
    aln=make_alignment_file(tmp_path, info['clone_id'],[
        ('NF_REF','AACC'),('R0000001','AA-C'),('R0000002','AA-C')])
    seq,qc,aligned_reads,aligned_consensus=fmod.make_clone(info,aln,.60,.60,2)
    assert seq == 'AAC'
    assert aligned_consensus == 'AA-C'
    assert qc['n_called_deletions'] == 1
    assert qc['qc_status']=='REVIEW'


def test_single_read_consensus(tmp_path):
    info={'clone_id':'PT1_clone000004','n_reads':1,'n_unique_reads':1,'reads':[
        {'internal_id':'R0000001','sample_id':'PT1-1','PCR_ID':'A'}]}
    path=make_alignment_file(tmp_path,info['clone_id'],[('R0000001','GGTACC')])
    seq,qc,aligned_reads,aligned_consensus=fmod.make_clone(info,path,.60,.60,2)
    assert seq=='GGTACC' and qc['consensus_type']=='single_read_consensus'
    assert aligned_consensus == seq
    assert aligned_reads == [('R0000001', seq)]


def test_component_reports_end_to_end_without_mafft(tmp_path):
    """Exercise Python components with a known-good mocked fragment alignment.

    This is NOT an execution test of actual MAFFT, container, or Nextflow.
    """
    source_refs = {'PT1-1':'refA.fa', 'PT1-2':'refB.fa', 'PT2-3':'refA.fa'}
    for pid, sample_ids in [('PT1',['PT1-1','PT1-2']), ('PT2',['PT2-3'])]:
        folder = tmp_path / pid
        folder.mkdir()
        samples = []
        for sid in sample_ids:
            ref = DATA/'references'/source_refs[sid]
            samples.append({'sample_id':sid,'timepoint':sid.split('-')[1],
                            'fasta':str(DATA/'fastas'/f'{sid}.fa'),
                            'integration_site':str(DATA/'sites'/f'{sid}.csv'),
                            'reference':str(ref),
                            'original_best_hiv_ref_filepath':str(ref)})
        staged = folder/'staged.json'
        staged.write_text(json.dumps(samples))
        subprocess.run([sys.executable, str(ROOT/'bin/prepare_patient.py'),
                        '--patient',pid,'--samples',str(staged)],cwd=folder,check=True)
        info=json.loads((folder/'clone_info.json').read_text())
        (folder/'aligned_clones').mkdir()
        for clone in info:
            path=folder/'aligned_clones'/(clone['clone_id']+'.fasta')
            if clone['n_reads']==1:
                path.write_text((folder/'clone_inputs'/(clone['clone_id']+'.fasta')).read_text())
            else:
                assert pid=='PT1' and clone['chromosome']=='chr1'
                aligned=['AAAACCCC--------', 'AAAATCCC--------', 'AAAACCCC--------',
                         '----CCCCGGGGTTTT','----CCCCGGGGTTTT']
                records=[('NF_REF','AAAACCCCGGGGTTTT')]
                records += [(r['internal_id'],seq) for r, seq in zip(clone['reads'],aligned)]
                with path.open('w') as output:
                    for header, sequence in records:
                        output.write(f'>{header}\n{sequence}\n')
        subprocess.run([sys.executable, str(ROOT/'bin/finalize_patient.py'),
                        '--patient',pid,'--aligned-clones','aligned_clones',
                        '--clone-info','clone_info.json',
                        '--counts',pid+'.counts.csv',
                        '--unmatched',pid+'.unmatched.csv',
                        '--unintegrated',pid+'_unintegrated.fasta',
                        '--reference-selection',pid+'.reference_selection.csv'],
                       cwd=folder, check=True)
        if pid=='PT1':
            clone_fa=list(fmod.fasta_iter(folder/'PT1_consensus.fasta'))
            assert clone_fa[0] == ('PT1_chr1_100_5_4','AAAACCCCGGGGTTTT')
            assert clone_fa[1] == ('PT1_chr2_200_1_1','GGGTTT')
            assert len(list(fmod.fasta_iter(folder/'PT1_unintegrated.fasta'))) == 2
            a = list(fmod.fasta_iter(folder/'PT1'/'PT1_clone000001_alignment.fasta'))
            assert a[0] == ('PT1_chr1_100_5_4', 'AAAACCCCGGGGTTTT')
            assert len(a) == 6
            assert a[1][0] == 'PT1-1_readA_more'
            assert a[1][1] == 'AAAACCCC--------'
            assert all(len(sequence) == len(a[0][1]) for _, sequence in a)
            b = list(fmod.fasta_iter(folder/'PT1'/'PT1_clone000002_alignment.fasta'))
            assert len(b) == 2 and b[0][1] == b[1][1] == 'GGGTTT'
            assert len(csvrows(folder/'PT1.unmatched.csv')) == 1
            assert csvrows(folder/'PT1.patient_stats.csv')[0]['n_input_sequences'] == '9'

    summary=tmp_path/'summary'
    summary.mkdir()
    from shutil import copy
    for pid in ('PT1','PT2'):
        src=tmp_path/pid
        for suffix in ('counts','clone_summary','patient_stats','unmatched',
                       'assignments','unused_sites','reference_selection'):
            copy(src/(pid+'.'+suffix+'.csv'),summary/(pid+'.'+suffix+'.csv'))
    subprocess.run([sys.executable,str(ROOT/'bin/make_summaries.py'),
                    '--counts','*.counts.csv',
                    '--clone-summaries','*.clone_summary.csv',
                    '--stats','*.patient_stats.csv',
                    '--unmatched','*.unmatched.csv',
                    '--assignments','*.assignments.csv',
                    '--unused-sites','*.unused_sites.csv',
                    '--references','*.reference_selection.csv'],cwd=summary,check=True)
    assert len(csvrows(summary/'clone_summary.csv'))==3
    assert len(csvrows(summary/'reference_selection.csv'))==2
    assert len(csvrows(summary/'read_assignments.csv'))==11
    assert len(csvrows(summary/'unmatched_fasta_reads.csv'))==1


def test_merged_fasta_is_declared_as_nextflow_output():
    workflow = (ROOT / 'main.nf').read_text()
    assert 'path "${patient_id}_merged.fasta", emit: merged_fasta' in workflow


def test_reverse_complemented_read_label(tmp_path):
    clone = {'clone_id':'PT1_clone000005','n_reads':2,'n_unique_reads':2,'reads':[
        {'internal_id':'R0000001','sample_id':'PT1-1','PCR_ID':'A', 'fasta_header':'PT1-1_readA'},
        {'internal_id':'R0000002','sample_id':'PT1-2','PCR_ID':'B', 'fasta_header':'PT1-2_readB'}]}
    aln = make_alignment_file(tmp_path,clone['clone_id'],[
        ('NF_REF','AAAA'), ('R0000001','AAAA'),('_R_R0000002','AAAA')])
    seq, qc, reads, aligned = fmod.make_clone(clone,aln,.60,.60,2)
    assert seq == aligned == 'AAAA'
    assert reads == [('PT1-1_readA','AAAA'),('PT1-2_readB [reverse_complemented]','AAAA')]
    assert qc['n_reverse_complemented_reads'] == 1
