#!/usr/bin/env python3
"""Combine patient-level CSVs into the final summary files."""
import argparse
import csv
import glob

REPORTS = [
    ('counts','patient_sequence_counts.csv'),
    ('clone_summaries','clone_summary.csv'),
    ('stats','patient_clone_statistics.csv'),
    ('unmatched','unmatched_fasta_reads.csv'),
    ('assignments','read_assignments.csv'),
    ('unused_sites','unused_integration_records.csv'),
    ('references','reference_selection.csv')
]


def combine(pattern, output):
    files = sorted(glob.glob(pattern))
    if not files:
        raise ValueError(f'{output}: no CSV files matched {pattern}')
    header = None
    rows = []
    for path in files:
        with open(path, newline='', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            if not reader.fieldnames:
                raise ValueError(f'{path}: missing CSV header')
            if header is None:
                header = reader.fieldnames
            elif reader.fieldnames != header:
                raise ValueError(f'{path}: CSV columns inconsistent with other patients')
            rows.extend(reader)
    with open(output, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=header, lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    for key, _ in REPORTS:
        ap.add_argument('--' + key.replace('_','-'), required=True)
    args = ap.parse_args()
    for key, name in REPORTS:
        combine(getattr(args, key), name)


if __name__ == '__main__':
    main()
