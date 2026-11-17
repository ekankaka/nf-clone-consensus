#!/usr/bin/env python3
"""Call a consensus for each integrated clone from its aligned reads.

The patient reference sets the alignment coordinates but does not supply bases.
Each sample_id + PCR_ID family contributes at most one vote per position.
"""
import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

BASES = 'ACGT'
STATES = set(BASES + '-')


def fasta_iter(path):
    head, parts = None, []
    with open(path, encoding='utf-8-sig') as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if head is not None:
                    yield head, ''.join(parts)
                head, parts = line[1:], []
            else:
                if head is None:
                    raise ValueError(f'Not a FASTA file: {path}')
                parts.append(line)
    if head is not None:
        yield head, ''.join(parts)


def fasta_write(handle, header, seq):
    handle.write('>' + header + '\n')
    for pos in range(0, len(seq), 80):
        handle.write(seq[pos:pos + 80] + '\n')


def read_csv(path):
    with open(path, encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle))


def write_csv(path, fields, rows):
    with open(path, 'w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator='\n', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def supported_call(votes, threshold):
    """Return the majority base or gap; use N if support is insufficient."""
    counts = Counter(v for v in votes if v in STATES)
    if not counts:
        return 'N', 0
    top = max(counts.values())
    winners = [v for v, n in counts.items() if n == top]
    if len(winners) != 1 or top / sum(counts.values()) < threshold:
        return 'N', sum(counts.values())
    return winners[0], sum(counts.values())


def column_states(seq):
    """Mark terminal gaps as missing coverage, keeping internal gaps as calls."""
    seq = seq.upper().replace('.', '-')
    first = next((i for i, v in enumerate(seq) if v != '-'), None)
    if first is None:
        raise ValueError('Alignment contains an entirely gapped read')
    last = len(seq) - 1 - next(i for i, v in enumerate(reversed(seq)) if v != '-')
    return [v if first <= i <= last else None for i, v in enumerate(seq)]


def make_clone(clone, aln_path, min_fraction, min_family_fraction, min_families):
    reads = clone['reads']
    if clone['n_reads'] == 1:
        # A site with one read keeps that read after terminal-N trimming.
        obs = list(fasta_iter(aln_path))
        if len(obs) != 1 or obs[0][0].split()[0] != reads[0]['internal_id']:
            raise ValueError(f"{clone['clone_id']}: failed to verify single-read FASTA")
        seq = obs[0][1].upper()
        return seq, {'consensus_type': 'single_read_consensus',
                     'n_ambiguous_bases': sum(c not in BASES for c in seq),
                     'n_evaluated_alignment_columns': len(seq),
                     'n_called_deletions': 0, 'n_internal_uncovered_positions': 0,
                     'max_uncovered_run': 0, 'n_low_support_positions': 0,
                     'n_multi_family_supported_positions': 0,
                     'percent_positions_multi_family': 0.0,
                     'n_reverse_complemented_reads': 0,
                     'qc_status': 'LIMITED_SUPPORT', 'qc_flags': 'SINGLE_READ'}

    observed = list(fasta_iter(aln_path))
    ids = {r['internal_id'] for r in reads}
    aligned = {}
    reverse_count = 0
    ref = None
    for head, seq in observed:
        label = head.split()[0]
        if label == 'NF_REF':
            if ref is not None:
                raise ValueError(f"{clone['clone_id']}: duplicate reference in MAFFT output")
            ref = seq.upper().replace('.', '-')
            continue
        reversed_orientation = label.startswith('_R_')
        if reversed_orientation:
            label = label[3:]
            reverse_count += 1
        if label not in ids or label in aligned:
            raise ValueError(f"{clone['clone_id']}: unexpected/duplicate aligned read {label}")
        aligned[label] = seq.upper().replace('.', '-')
    if ref is None or set(aligned) != ids:
        raise ValueError(f"{clone['clone_id']}: missing reference or missing reads in alignment")
    width = len(ref)
    if any(len(seq) != width for seq in aligned.values()):
        raise ValueError(f"{clone['clone_id']}: inconsistent MAFFT alignment widths")
    states = {label: column_states(seq) for label, seq in aligned.items()}
    families = defaultdict(list)
    for read in reads:
        families[(read['sample_id'], read['PCR_ID'])].append(states[read['internal_id']])
    if len(families) != clone['n_unique_reads']:
        raise ValueError(f"{clone['clone_id']}: PCR family count mismatch")

    # First make one call per PCR family at each aligned position.
    # Terminal gaps do not count as deletions.
    family_columns = []
    for family in sorted(families):
        original = families[family]
        votes = []
        for i in range(width):
            local_votes = [read[i] for read in original if read[i] is not None]
            if not local_votes:
                votes.append(None)
            else:
                call, _ = supported_call(local_votes, min_family_fraction)
                votes.append(call)
        family_columns.append(votes)

    aggregate_coverage = [sum(f[i] is not None for f in family_columns) for i in range(width)]
    covered = [i for i, n in enumerate(aggregate_coverage) if n > 0]
    if not covered:
        raise ValueError(f"{clone['clone_id']}: no supported read positions")
    start, stop = min(covered), max(covered)
    answer = []
    internal_uncovered = 0
    max_uncovered_run = 0
    run = 0
    low_support = 0
    multi_supported = 0
    called_deletions = 0
    for i in range(start, stop + 1):
        if aggregate_coverage[i] == 0:
            # Uncovered positions stay as N, not reference bases.
            # This also keeps the distance between non-overlapping fragments.
            if ref[i] != '-':
                answer.append('N')
                internal_uncovered += 1
                run += 1
                max_uncovered_run = max(max_uncovered_run, run)
            continue
        run = 0
        votes = [family[i] for family in family_columns if family[i] is not None]
        call, informative = supported_call(votes, min_fraction)
        if informative >= min_families:
            multi_supported += 1
        elif len(families) >= min_families:
            low_support += 1
            call = 'N'
        if call == '-':
            # Remove independently supported deletions from the FASTA sequence.
            called_deletions += 1
            continue
        answer.append(call)

    seq = ''.join(answer)
    n_evaluated = stop - start + 1
    if not seq:
        raise ValueError(f"{clone['clone_id']}: consensus has zero length")
    flags = []
    if len(families) == 1:
        flags.append('SINGLE_PCR_FAMILY')
    if internal_uncovered:
        flags.append('INTERNAL_NO_COVERAGE')
    if low_support:
        flags.append('LOW_INDEPENDENT_COVERAGE')
    if called_deletions:
        flags.append('INDEL_CALLS_REVIEW')
    if 'N' in seq:
        flags.append('AMBIGUOUS_CONSENSUS_BASES')
    # Keep flagged clones in the FASTA so their alignments can be reviewed.
    qc_status = ('REVIEW' if internal_uncovered or called_deletions or 'N' in seq else
                 'LIMITED_SUPPORT' if len(families) == 1 or low_support else 'PASS')
    return seq, {'consensus_type': 'single_family_consensus' if len(families) == 1 else 'clone_consensus',
                 'n_ambiguous_bases': sum(c not in BASES for c in seq),
                 'n_evaluated_alignment_columns': n_evaluated,
                 'n_called_deletions': called_deletions,
                 'n_internal_uncovered_positions': internal_uncovered,
                 'max_uncovered_run': max_uncovered_run,
                 'n_low_support_positions': low_support,
                 'n_multi_family_supported_positions': multi_supported,
                 'percent_positions_multi_family': round(100 * multi_supported / n_evaluated, 2),
                 'n_reverse_complemented_reads': reverse_count,
                 'qc_status': qc_status,
                 'qc_flags': ';'.join(flags)}


CLONE_COLUMNS = ['patient_id', 'clone_id', 'consensus_header', 'chromosome', 'INTEGRATION_SITE',
                 'n_reads', 'n_unique_reads', 'n_samples', 'n_timepoints', 'sample_ids',
                 'timepoints', 'consensus_length', 'n_evaluated_alignment_columns', 'n_ambiguous_bases', 'percent_N',
                 'n_multi_family_supported_positions', 'percent_positions_multi_family',
                 'n_low_support_positions', 'n_internal_uncovered_positions',
                 'max_uncovered_run', 'n_called_deletions', 'n_reverse_complemented_reads',
                 'consensus_type', 'qc_status', 'qc_flags', 'reference_sha256']
STAT_COLUMNS = ['patient_id', 'n_input_fastas', 'n_input_sequences', 'n_integrated_reads_used',
                'n_unique_integration_sites', 'n_multifamily_clones', 'n_single_family_clones',
                'n_single_read_clones', 'n_unintegrated_reads', 'n_review_reads',
                'n_clones_flagged_for_review', 'n_final_integrated_sequences']


def main():
    p = argparse.ArgumentParser()
    for key in ['patient', 'aligned-clones', 'clone-info', 'counts', 'unmatched',
                'unintegrated', 'reference-selection']:
        p.add_argument('--' + key, required=True)
    p.add_argument('--min-base-fraction', type=float, default=0.60)
    p.add_argument('--min-family-fraction', type=float, default=0.60)
    p.add_argument('--min-families-per-call', type=int, default=2)
    args = p.parse_args()
    if not (0.5 < args.min_base_fraction <= 1.0 and 0.5 < args.min_family_fraction <= 1.0):
        raise ValueError('Both fractional support thresholds must be >0.5 and <=1')
    if args.min_families_per_call < 2:
        raise ValueError('min-families-per-call must be at least 2')

    clones = json.loads(Path(args.clone_info).read_text(encoding='utf-8'))
    count = read_csv(args.counts)[0]
    ref = read_csv(args.reference_selection)[0]
    unmatched_count = len(read_csv(args.unmatched))
    no_site_count = len(list(fasta_iter(args.unintegrated)))
    rows = []
    with open(f'{args.patient}_consensus.fasta', 'w', encoding='utf-8') as output:
        for clone in clones:
            aligned_file = Path(args.aligned_clones) / (clone['clone_id'] + '.fasta')
            seq, qc = make_clone(clone, aligned_file, args.min_base_fraction,
                                 args.min_family_fraction, args.min_families_per_call)
            fasta_write(output, clone['consensus_header'], seq)
            rows.append({
                'patient_id': args.patient, 'clone_id': clone['clone_id'],
                'consensus_header': clone['consensus_header'],
                'chromosome': clone['chromosome'], 'INTEGRATION_SITE': clone['INTEGRATION_SITE'],
                'n_reads': clone['n_reads'], 'n_unique_reads': clone['n_unique_reads'],
                'n_samples': len(clone['sample_ids']), 'n_timepoints': len(clone['timepoints']),
                'sample_ids': ';'.join(clone['sample_ids']),
                'timepoints': ';'.join(clone['timepoints']),
                'consensus_length': len(seq),
                'percent_N': round(100 * seq.count('N') / len(seq), 2),
                'reference_sha256': ref['reference_sha256'], **qc
            })

    all_reads = int(count['sequence_count'])
    used = sum(int(c['n_reads']) for c in clones)
    if all_reads != used + unmatched_count + no_site_count:
        raise ValueError(f'{args.patient}: accounting failed; {all_reads} input vs '
                         f'{used} integrated + {no_site_count} unintegrated + {unmatched_count} review')
    write_csv(f'{args.patient}.clone_summary.csv', CLONE_COLUMNS, rows)
    write_csv(f'{args.patient}.patient_stats.csv', STAT_COLUMNS, [{
        'patient_id': args.patient,
        'n_input_fastas': int(count['n_fasta_files']),
        'n_input_sequences': all_reads,
        'n_integrated_reads_used': used,
        'n_unique_integration_sites': len(clones),
        'n_multifamily_clones': sum(c['n_unique_reads'] > 1 for c in clones),
        'n_single_family_clones': sum(c['n_unique_reads'] == 1 and c['n_reads'] > 1 for c in clones),
        'n_single_read_clones': sum(c['n_reads'] == 1 for c in clones),
        'n_unintegrated_reads': no_site_count,
        'n_review_reads': unmatched_count,
        'n_clones_flagged_for_review': sum(r['qc_status'] == 'REVIEW' for r in rows),
        'n_final_integrated_sequences': len(clones)
    }])


if __name__ == '__main__':
    main()
