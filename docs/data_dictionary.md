# Data dictionary

This dictionary describes the inputs and final outputs of `nf-clone-consensus`. Field names are case-sensitive. Counts are per patient unless stated otherwise.

## 1. Input manifest (`--input`)

One row represents one sample. File paths can be absolute or relative to the directory containing the manifest.

| Column | Meaning |
|---|---|
| `sample_id` | Unique sample identifier within a patient; used to prefix FASTA headers and define PCR families |
| `patient_id` | Patient identifier; reads from the same patient are grouped across timepoints |
| `timepoint` | Sampling timepoint; the synthetic examples use 1, 2, and 3 |
| `hiv_fasta_filepath` | FASTA containing HIV reads for the sample |
| `integration_site_filepath` | CSV/TSV containing integration-site assignments for the sample |
| `best_hiv_ref_filepath` | FASTA containing a single best HIV reference sequence for the sample |

The integration-site file requires these columns:

| Column | Meaning |
|---|---|
| `READ` | Read identifier used to match an input FASTA header; longest matching prefix is selected |
| `PCR_ID` | PCR-family identifier, unique **within** a sample; duplicates share `(sample_id, PCR_ID)` |
| `chromosome` | Host chromosome of integration; `HIV` indicates no identified host integration site |
| `INTEGRATION_SITE` | Integration coordinate or site label; missing/`NA` means no identified integration site |

A valid integration-site clone needs a host chromosome, an integration-site value, and a PCR ID. Reads without a matched integration-site record, or with conflicting/invalid annotations, go to `review/`.

## 2. FASTA outputs

| Path | Meaning |
|---|---|
| `merged_fastas/<patient_id>_merged.fasta` | All sample FASTA sequences combined for a patient, with `sample_id_` prefixed to each original header |
| `clone_consensus/<patient_id>_consensus.fasta` | One sequence for each host integration site, including single-read and single-PCR-family consensuses |
| `unintegrated_sequences/<patient_id>_unintegrated.fasta` | Reads marked `chromosome=HIV` or with a missing integration site; preserves sample-prefixed original headers and sequences |
| `review/unmatched_fastas/<patient_id>.unmatched.fasta` | Reads excluded from clone construction because their integration-site assignment could not be used |

**Clone FASTA header**

```text
>patient_id_chromosome_INTEGRATION_SITE_n_reads_n_unique_reads
```

`n_reads` counts accepted FASTA reads assigned to the clone. `n_unique_reads` counts distinct `(sample_id, PCR_ID)` families, **not** unique nucleotide strings. For example, `>PT1_chr1_100_5_4` has five reads from four PCR families. Unusual characters in chromosome/site labels may be URL-encoded in FASTA headers; the original strings appear in `clone_summary.csv`.

## 3. `patient_sequence_counts.csv`

| Column | Meaning |
|---|---|
| `patient_id` | Patient identifier |
| `n_fasta_files` | Number of input FASTA files listed for this patient |
| `sequence_count` | Total FASTA sequence records across those files, including integrated, unintegrated, and review reads |

## 4. `clone_summary.csv`

One row per constructed clone consensus.

| Column | Meaning |
|---|---|
| `patient_id` | Patient identifier |
| `clone_id` | Internal sequential clone identifier within patient, e.g. `PT1_clone000001` |
| `consensus_header` | Full FASTA header without `>` |
| `chromosome` | Host chromosome of the integration site |
| `INTEGRATION_SITE` | Original integration-site value from the annotation |
| `n_reads` | Number of accepted input FASTA reads assigned to this clone |
| `n_unique_reads` | Distinct `(sample_id, PCR_ID)` families assigned to the clone |
| `n_samples` | Number of distinct sample IDs contributing reads |
| `n_timepoints` | Number of distinct timepoints contributing reads |
| `sample_ids` | Semicolon-separated contributing sample IDs |
| `timepoints` | Semicolon-separated contributing timepoints |
| `consensus_length` | Length in nucleotides of the final consensus FASTA sequence, excluding called deletions |
| `n_evaluated_alignment_columns` | Number of alignment columns between the first and last covered columns, including deleted columns |
| `n_ambiguous_bases` | Number of consensus characters not equal to A/C/G/T |
| `percent_N` | Percentage of the final consensus sequence represented by `N` |
| `n_multi_family_supported_positions` | Alignment positions with at least two informative PCR-family calls; may include deletion columns |
| `percent_positions_multi_family` | `100 × n_multi_family_supported_positions / n_evaluated_alignment_columns` |
| `n_low_support_positions` | Covered alignment positions with too few informative PCR families for a multi-family call |
| `n_internal_uncovered_positions` | Internal reference positions between covered segments with no PCR-family coverage; represented by `N` |
| `max_uncovered_run` | Longest consecutive span of internally uncovered reference positions |
| `n_called_deletions` | Alignment columns called as deletions and omitted from the FASTA sequence |
| `n_reverse_complemented_reads` | Reads reverse-complemented by MAFFT's orientation adjustment |
| `consensus_type` | `clone_consensus`, `single_family_consensus`, or `single_read_consensus` |
| `qc_status` | `PASS`, `LIMITED_SUPPORT`, or `REVIEW` |
| `qc_flags` | Semicolon-separated reasons for caution; empty when none |
| `reference_sha256` | SHA-256 hash of the selected reference sequence (uppercase nucleotide string) |

**Consensus rules (defaults):** Each PCR family contributes at most one vote per aligned position. A base or deletion must have a unique majority with at least 60% support. Multi-family clones require at least two informative families to call a position; insufficiently supported positions are `N`. Terminal gaps are ignored as missing coverage. A single original read retains its trimmed sequence.

**QC flags**

| Value | Meaning |
|---|---|
| `SINGLE_READ` | Consensus contains a single original read |
| `SINGLE_PCR_FAMILY` | Multiple reads, but all belong to one PCR family |
| `INTERNAL_NO_COVERAGE` | Internal span with no read coverage |
| `LOW_INDEPENDENT_COVERAGE` | Some covered positions lack sufficient independent PCR-family support |
| `INDEL_CALLS_REVIEW` | One or more deletions were called; inspect alignment |
| `AMBIGUOUS_CONSENSUS_BASES` | Consensus contains one or more `N` calls |

`PASS` means no flagged issue under the current rules. `LIMITED_SUPPORT` means single-family or low independent support without another reason requiring `REVIEW`. `REVIEW` indicates ambiguous bases, missing internal coverage, or called deletions. These labels are QC aids, not a measure of biological truth.

## 5. `patient_clone_statistics.csv`

One row per patient.

| Column | Meaning |
|---|---|
| `patient_id` | Patient identifier |
| `n_input_fastas` | Number of input FASTA files |
| `n_input_sequences` | Total input FASTA sequences |
| `n_integrated_reads_used` | Accepted reads used to construct integrated-site clones |
| `n_unique_integration_sites` | Number of distinct host integration sites assigned to clones |
| `n_multifamily_clones` | Clones supported by two or more distinct PCR families |
| `n_single_family_clones` | Clones with multiple reads but only one PCR family |
| `n_single_read_clones` | Clones represented by exactly one read |
| `n_unintegrated_reads` | Reads classified as having no identified host integration site |
| `n_review_reads` | Reads excluded from consensus and sent for review |
| `n_clones_flagged_for_review` | Clones with `qc_status=REVIEW` |
| `n_final_integrated_sequences` | Number of integrated clone sequences written (one per clone) |

For each patient, `n_input_sequences = n_integrated_reads_used + n_unintegrated_reads + n_review_reads`.

## 6. `reference_selection.csv`

One row per patient.

| Column | Meaning |
|---|---|
| `patient_id` | Patient identifier |
| `selected_sample_id` | Sample whose reference FASTA was selected |
| `selected_timepoint` | Timepoint of the selected sample |
| `selected_best_hiv_ref_filepath` | Original path to the selected reference FASTA |
| `reference_sha256` | SHA-256 of the selected reference sequence |
| `n_distinct_reference_sequences` | Number of distinct reference sequences among samples for the patient |
| `reference_selection_rule` | `IDENTICAL_REFERENCE_SEQUENCES` or `HIGHEST_READ_TIMEPOINT_THEN_SAMPLE` |
| `selected_timepoint_n_reads` | Total reads across samples at the selected timepoint |
| `selected_sample_n_reads` | Input reads in the selected sample |

If all reference sequences are identical, the pipeline selects one deterministically. Otherwise it chooses the timepoint with the most reads, then its highest-read sample if needed. An unresolved tie among different references stops the run.

## 7. `review/read_assignments.csv` and `review/unmatched_fasta_reads.csv`

`read_assignments.csv` contains one row for each input FASTA read. `unmatched_fasta_reads.csv` contains the subset classified `REVIEW`.

| Column | Meaning |
|---|---|
| `patient_id` | Patient identifier |
| `sample_id` | Source sample |
| `timepoint` | Source timepoint |
| `fasta_header` | Header of the sample-prefixed merged FASTA sequence, without `>` |
| `READ` | Matched READ identifier; empty if no match |
| `PCR_ID` | PCR identifier from the matched annotation, if available |
| `pcr_family_key` | `sample_id:PCR_ID` for accepted clone reads; otherwise empty |
| `chromosome` | Matched chromosome, if available |
| `INTEGRATION_SITE` | Matched integration-site value, if available |
| `classification` | `CLONE`, `NO_SITE`, or `REVIEW` |
| `clone_id` | Clone identifier for accepted integrated reads; empty otherwise |
| `reason` | Reason for review, when classification is `REVIEW` |
| `original_length` | Input FASTA sequence length before terminal-N trimming |
| `trimmed_length` | Sequence length after trimming, for accepted clone reads |

**Review reasons:** `NO_MATCHING_INTEGRATION_RECORD`, `CONFLICTING_READ_ASSIGNMENTS`, `MISSING_CHROMOSOME_WITH_SITE`, `MISSING_PCR_ID`, `PCR_ID_CONFLICTING_SITES`, `ALL_TERMINAL_NS`, and `INVALID_NUCLEOTIDE_SYMBOL`.

## 8. `review/unused_integration_records.csv`

Integration-site records that could not be linked to any input FASTA header.

| Column | Meaning |
|---|---|
| `patient_id` | Patient identifier |
| `sample_id` | Source sample |
| `READ` | Unmatched READ value in the integration-site table |
| `PCR_ID` | Its PCR-family identifier |
| `chromosome` | Its chromosome value |
| `INTEGRATION_SITE` | Its integration-site value |
| `reason` | Currently `NO_MATCHED_FASTA_HEADER` |

## Interpreting the files

The merged FASTAs contain all input sequences. The consensus FASTAs contain only integrated-site clones; unintegrated sequences and reads requiring review are kept separately. Check `read_assignments.csv` and `clone_summary.csv` before downstream analyses. The selected reference guides alignment, but does not contribute bases to the consensus. Clones with disjoint fragments or low coverage may not represent a fully reconstructed HIV genome.
