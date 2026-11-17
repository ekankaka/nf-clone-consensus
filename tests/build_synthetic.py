#!/usr/bin/env python3
"""Build the small synthetic FASTA and integration-site test dataset."""
import csv
from pathlib import Path

base = Path(__file__).resolve().parent / 'data'
for folder in ('fastas','sites','references'):
    (base / folder).mkdir(parents=True, exist_ok=True)


def fasta(path, records):
    with path.open('w') as f:
        for name, seq in records:
            f.write(f'>{name}\n{seq}\n')


def site(path, rows):
    with path.open('w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['READ','chromosome','INTEGRATION_SITE','PCR_ID'])
        writer.writerows(rows)

fasta(base / 'references/refA.fa', [('RefA','AAAACCCCGGGGTTTT')])
fasta(base / 'references/refB.fa', [('RefB','AAAACCCCGGGGTTTA')])
fasta(base / 'fastas/PT1-1.fa', [
    ('readA_more','NNAAAACCCCNN'), ('readB','AAAATCCC'),
    ('readC','AAAACCCC'), ('readNS','GGGTT'),
])
site(base / 'sites/PT1-1.csv', [
    ('readA','chr1','100','pcr01'), ('readB','chr1','100','pcr01'),
    ('readC','chr1','100','pcr02'), ('readNS','HIV','NA','pcr99'),
])
fasta(base / 'fastas/PT1-2.fa', [
    ('readD','CCCCGGGGTTTT'), ('readE','CCCCGGGGTTTT'),
    ('readS','NNGGGTTTNN'), ('readNA','TCTCT'),
    ('readUNMATCHED','ATATAT'),
])
site(base / 'sites/PT1-2.csv', [
    ('readD','chr1','100','pcr01'), ('readE','chr1','100','pcr02'),
    ('readS','chr2','200','pcr03'), ('readNA','chr4','NA','pcr04'),
    ('unused_id','chr5','1000','pcr05'),
])
fasta(base / 'fastas/PT2-3.fa', [('readOnly','NNGGTACCNN'),('readHIV','AAGGT')])
site(base / 'sites/PT2-3.csv', [('readOnly','chr8','789','pcr1'),('readHIV','HIV','NA','pcr2')])
with (base / 'annotation.csv').open('w', newline='') as f:
    w=csv.writer(f)
    w.writerow(['sample_id','patient_id','timepoint','hiv_fasta_filepath','integration_site_filepath','best_hiv_ref_filepath'])
    w.writerows([
        ['PT1-1','PT1','1','fastas/PT1-1.fa','sites/PT1-1.csv','references/refA.fa'],
        ['PT1-2','PT1','2','fastas/PT1-2.fa','sites/PT1-2.csv','references/refB.fa'],
        ['PT2-3','PT2','3','fastas/PT2-3.fa','sites/PT2-3.csv','references/refA.fa'],
    ])
print(base)
