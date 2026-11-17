#!/usr/bin/env python3
"""Prepare patient FASTAs, integration-site groups, and the reference for alignment."""
import argparse
import csv
import gzip
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import quote

NULLS = {'', 'na', 'n/a', 'null', 'none', 'nan'}
DNA = set('ACGTRYSWKMBDHVN')
ASSIGN_COLUMNS = ['patient_id', 'sample_id', 'timepoint', 'fasta_header', 'READ',
                  'PCR_ID', 'pcr_family_key', 'chromosome', 'INTEGRATION_SITE',
                  'classification', 'clone_id', 'reason', 'original_length', 'trimmed_length']
UNUSED_COLUMNS = ['patient_id', 'sample_id', 'READ', 'PCR_ID', 'chromosome',
                  'INTEGRATION_SITE', 'reason']
REF_COLUMNS = ['patient_id', 'selected_sample_id', 'selected_timepoint',
               'selected_best_hiv_ref_filepath', 'reference_sha256',
               'n_distinct_reference_sequences', 'reference_selection_rule',
               'selected_timepoint_n_reads', 'selected_sample_n_reads']


def missing(value):
    return value is None or str(value).strip().lower() in NULLS


def fasta_iter(path):
    opener = gzip.open if str(path).lower().endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8-sig') as inp:
        head, chunks = None, []
        for line in inp:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if head is not None:
                    if not chunks:
                        raise ValueError(f'Empty FASTA sequence: {path} >{head}')
                    yield head, ''.join(chunks)
                head, chunks = line[1:], []
                if not head or not head.split()[0]:
                    raise ValueError(f'Empty FASTA header in {path}')
            else:
                if head is None:
                    raise ValueError(f'Sequence precedes FASTA header in {path}')
                chunks.append(line)
        if head is not None:
            if not chunks:
                raise ValueError(f'Empty FASTA sequence: {path} >{head}')
            yield head, ''.join(chunks)


def fasta_write(stream, header, sequence):
    stream.write('>' + header + '\n')
    for pos in range(0, len(sequence), 80):
        stream.write(sequence[pos:pos + 80] + '\n')


def write_csv(path, columns, rows):
    with open(path, 'w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction='ignore', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def load_sites(path):
    opener = gzip.open if str(path).lower().endswith('.gz') else open
    with opener(path, 'rt', newline='', encoding='utf-8-sig') as handle:
        first = handle.readline()
        if not first:
            raise ValueError(f'Empty integration site table: {path}')
        handle.seek(0)
        delim = '\t' if '\t' in first else ','
        reader = csv.DictReader(handle, delimiter=delim)
        headers = {str(x).strip() for x in (reader.fieldnames or [])}
        required = {'READ', 'chromosome', 'INTEGRATION_SITE', 'PCR_ID'}
        if not required.issubset(headers):
            raise ValueError(f'{path} missing required columns: {sorted(required - headers)}')
        for num, row in enumerate(reader, 2):
            if not any(v and str(v).strip() for v in row.values()):
                continue
            rec = {k.strip(): (v or '').strip() for k, v in row.items() if k}
            if missing(rec['READ']):
                raise ValueError(f'{path}, line {num}: READ missing')
            yield dict(rec, line_number=num)


def sha_reference(path):
    records = list(fasta_iter(path))
    if len(records) != 1:
        raise ValueError(f'Reference FASTA must contain exactly one sequence: {path}; found {len(records)}')
    seq = records[0][1].upper()
    if set(seq) - DNA:
        raise ValueError(f'Reference contains invalid nucleotide symbols: {path}')
    return seq, hashlib.sha256(seq.encode('ascii')).hexdigest()


def select_reference(patient, samples):
    """Use a shared reference if available; otherwise select by read count."""
    all_refs = {s['ref_sha256'] for s in samples}
    if len(all_refs) == 1:
        selected = sorted(samples, key=lambda x: (x['timepoint'], x['sample_id']))[0]
        reason = 'IDENTICAL_REFERENCE_SEQUENCES'
    else:
        tp_totals = Counter()
        for s in samples:
            tp_totals[s['timepoint']] += s['n_reads']
        best_total = max(tp_totals.values())
        winners = sorted(tp for tp, total in tp_totals.items() if total == best_total)
        if len(winners) > 1:
            winning_refs = {s['ref_sha256'] for s in samples if s['timepoint'] in winners}
            if len(winning_refs) > 1:
                raise ValueError(f'{patient}: different HIV references and a tie among maximum-read '
                                 f'timepoints {winners} ({best_total} reads each). Resolve before running.')
        winner_tp = winners[0]
        contenders = [s for s in samples if s['timepoint'] == winner_tp]
        largest = max(s['n_reads'] for s in contenders)
        finalists = [s for s in contenders if s['n_reads'] == largest]
        if len({s['ref_sha256'] for s in finalists}) > 1:
            raise ValueError(f'{patient}: multiple distinct references among equally high-read '
                             f'samples at timepoint {winner_tp}. Resolve before running.')
        selected = sorted(finalists, key=lambda x: x['sample_id'])[0]
        reason = 'HIGHEST_READ_TIMEPOINT_THEN_SAMPLE'
    total_selected_tp = sum(s['n_reads'] for s in samples if s['timepoint'] == selected['timepoint'])
    metadata = {'patient_id': patient, 'selected_sample_id': selected['sample_id'],
                'selected_timepoint': selected['timepoint'],
                'selected_best_hiv_ref_filepath': selected.get('original_best_hiv_ref_filepath', selected['reference']),
                'reference_sha256': selected['ref_sha256'],
                'n_distinct_reference_sequences': len(all_refs),
                'reference_selection_rule': reason,
                'selected_timepoint_n_reads': total_selected_tp,
                'selected_sample_n_reads': selected['n_reads']}
    return selected, metadata


def longest_read_prefix(key, indexed):
    """Match a FASTA header to the longest READ prefix in the site table."""
    prefixes = [read for read in indexed if key.startswith(read)]
    if not prefixes:
        return None
    return max(prefixes, key=len)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--patient', required=True)
    parser.add_argument('--samples', required=True)
    args = parser.parse_args()
    patient = args.patient
    samples = json.loads(Path(args.samples).read_text(encoding='utf-8'))
    Path('clone_inputs').mkdir(exist_ok=True)

    records = []
    site_index = {}
    families_by_sample = {}
    selected_samples = []
    for sample in samples:
        sid, timepoint = sample['sample_id'], sample['timepoint']
        rs = list(fasta_iter(sample['fasta']))
        if not rs:
            raise ValueError(f'{patient}/{sid}: empty input FASTA')
        sequence, digest = sha_reference(sample['reference'])
        selected_samples.append(dict(sample, n_reads=len(rs), ref_sequence=sequence, ref_sha256=digest))
        per_read = defaultdict(list)
        family_sites = defaultdict(set)
        for r in load_sites(sample['integration_site']):
            per_read[r['READ']].append(r)
            if not missing(r['PCR_ID']):
                family_sites[r['PCR_ID']].add((r['chromosome'], r['INTEGRATION_SITE']))
        site_index[sid] = per_read
        families_by_sample[sid] = family_sites
        for head, seq in rs:
            records.append({'sample_id': sid, 'timepoint': timepoint,
                            'header': sid + '_' + head, 'sequence': seq,
                            'original_length': len(seq)})

    # Select one reference for all clones from this patient.
    selected, ref_metadata = select_reference(patient, selected_samples)
    with open('patient_reference.fasta', 'w', encoding='utf-8') as handle:
        fasta_write(handle, 'NF_REF', selected['ref_sequence'])
    write_csv(f'{patient}.reference_selection.csv', REF_COLUMNS, [ref_metadata])

    records.sort(key=lambda r: (r['sample_id'], r['header']))
    ids = [r['header'].split()[0] for r in records]
    if len(set(ids)) != len(ids):
        raise ValueError(f'{patient}: duplicate sample-prefixed FASTA header identifiers')
    with open(f'{patient}_merged.fasta', 'w', encoding='utf-8') as output:
        for r in records:
            fasta_write(output, r['header'], r['sequence'])

    # Keep unmatched or conflicting reads separate; they should not enter a consensus.
    assignments, unresolved, unintegrated, clone_groups = [], [], [], defaultdict(list)
    matched_records = set()
    for rec in records:
        sid, head, sequence = rec['sample_id'], rec['header'], rec['sequence']
        unprefixed = head[len(sid) + 1:].split()[0]
        key = longest_read_prefix(unprefixed, site_index[sid])
        row = {'patient_id': patient, 'sample_id': sid, 'timepoint': rec['timepoint'],
               'fasta_header': head, 'READ': key or '', 'PCR_ID': '', 'pcr_family_key': '',
               'chromosome': '', 'INTEGRATION_SITE': '', 'classification': '',
               'clone_id': '', 'reason': '', 'original_length': rec['original_length'],
               'trimmed_length': ''}
        if key is None:
            row.update(classification='REVIEW', reason='NO_MATCHING_INTEGRATION_RECORD')
        else:
            matched_records.add((sid, key))
            matching_rows = site_index[sid][key]
            assignments_for_read = {(r['chromosome'], r['INTEGRATION_SITE'], r['PCR_ID'])
                                    for r in matching_rows}
            if len(assignments_for_read) != 1:
                row.update(classification='REVIEW', reason='CONFLICTING_READ_ASSIGNMENTS')
            else:
                chrom, site, pcr = next(iter(assignments_for_read))
                row.update(chromosome=chrom, INTEGRATION_SITE=site, PCR_ID=pcr)
                if chrom.upper() == 'HIV' or missing(site):
                    row.update(classification='NO_SITE')
                    unintegrated.append(rec)
                elif missing(chrom):
                    row.update(classification='REVIEW', reason='MISSING_CHROMOSOME_WITH_SITE')
                elif missing(pcr):
                    row.update(classification='REVIEW', reason='MISSING_PCR_ID')
                elif len(families_by_sample[sid][pcr]) > 1:
                    row.update(classification='REVIEW', reason='PCR_ID_CONFLICTING_SITES')
                else:
                    trimmed = sequence.upper().strip('N')
                    if not trimmed:
                        row.update(classification='REVIEW', reason='ALL_TERMINAL_NS')
                    elif set(trimmed) - DNA:
                        row.update(classification='REVIEW', reason='INVALID_NUCLEOTIDE_SYMBOL')
                    else:
                        row.update(classification='CLONE', trimmed_length=len(trimmed),
                                   pcr_family_key=sid + ':' + pcr)
                        clone_groups[(chrom, site)].append(dict(rec, trimmed=trimmed, READ=key,
                                                                PCR_ID=pcr, assignment=row))
        assignments.append(row)
        if row['classification'] == 'REVIEW':
            unresolved.append(rec)

    with open(f'{patient}_unintegrated.fasta', 'w', encoding='utf-8') as f:
        for rec in unintegrated:
            fasta_write(f, rec['header'], rec['sequence'])
    with open(f'{patient}.unmatched.fasta', 'w', encoding='utf-8') as f:
        for rec in unresolved:
            fasta_write(f, rec['header'], rec['sequence'])

    clone_metadata = []
    header_seen = set()
    for idx, ((chrom, site), reads) in enumerate(sorted(clone_groups.items()), 1):
        clone_id = f'{patient}_clone{idx:06d}'
        # Encode unusual characters so clone headers remain valid FASTA identifiers.
        printable_chrom = quote(chrom, safe='.-')
        printable_site = quote(site, safe='.-')
        header = f'{patient}_{printable_chrom}_{printable_site}_{len(reads)}_' \
                 f'{len({(r["sample_id"], r["PCR_ID"]) for r in reads})}'
        if header in header_seen:
            raise ValueError(f'Clone FASTA header collision: {header}')
        header_seen.add(header)
        metadata_reads = []
        with open(f'clone_inputs/{clone_id}.fasta', 'w', encoding='utf-8') as f:
            for r_idx, r in enumerate(reads, 1):
                internal_id = f'R{r_idx:07d}'
                fasta_write(f, internal_id, r['trimmed'])
                r['assignment']['clone_id'] = clone_id
                metadata_reads.append({'internal_id': internal_id, 'sample_id': r['sample_id'],
                                       'timepoint': r['timepoint'], 'READ': r['READ'],
                                       'PCR_ID': r['PCR_ID'], 'fasta_header': r['header'],
                                       'trimmed_length': len(r['trimmed'])})
        clone_metadata.append({'clone_id': clone_id, 'patient_id': patient,
                               'chromosome': chrom, 'INTEGRATION_SITE': site,
                               'consensus_header': header, 'n_reads': len(reads),
                               'n_unique_reads': len({(r['sample_id'], r['PCR_ID']) for r in reads}),
                               'sample_ids': sorted({r['sample_id'] for r in reads}),
                               'timepoints': sorted({r['timepoint'] for r in reads}),
                               'reads': metadata_reads})

    unused_rows = []
    for sample in samples:
        sid = sample['sample_id']
        for read, variants in site_index[sid].items():
            if (sid, read) not in matched_records:
                for r in variants:
                    unused_rows.append({'patient_id': patient, 'sample_id': sid, 'READ': read,
                                        'PCR_ID': r['PCR_ID'], 'chromosome': r['chromosome'],
                                        'INTEGRATION_SITE': r['INTEGRATION_SITE'],
                                        'reason': 'NO_MATCHED_FASTA_HEADER'})

    Path('clone_info.json').write_text(json.dumps(clone_metadata, indent=2), encoding='utf-8')
    write_csv(f'{patient}.counts.csv', ['patient_id','n_fasta_files','sequence_count'],
              [{'patient_id': patient, 'n_fasta_files': len(samples), 'sequence_count': len(records)}])
    write_csv(f'{patient}.unmatched.csv', ASSIGN_COLUMNS,
              [r for r in assignments if r['classification'] == 'REVIEW'])
    write_csv(f'{patient}.assignments.csv', ASSIGN_COLUMNS, assignments)
    write_csv(f'{patient}.unused_sites.csv', UNUSED_COLUMNS, unused_rows)
    # Each input read must be assigned to exactly one output category.
    if len(records) != sum([len(unintegrated), len(unresolved),
                            sum(len(group) for group in clone_groups.values())]):
        raise ValueError(f'{patient}: not all FASTA records classified exactly once')


if __name__ == '__main__':
    main()
