# nf-clone-consensus

Nextflow pipeline for HIV clone consensus sequences from integration-site-linked reads, with PCR duplicate handling and reference-guided alignment.

The pipeline produces one consensus per integration site. Reads with no identified integration site are kept separately. Unmatched or conflicting annotations are reported for review.

## Requirements

- Nextflow and a compatible Java installation
- One container runtime: Docker (PC), Singularity, or Apptainer
- SLURM, only when running on an HPC with that scheduler
- Python 3 and `pytest`, only for the Python tests

Python and MAFFT are provided by the process containers. The workflow was developed with Nextflow 25.10.4. Container tags are specified in `main.nf`; pin image digests for validated production releases.

## Input

Supply a CSV with these six columns:

```csv
sample_id,patient_id,timepoint,hiv_fasta_filepath,integration_site_filepath,best_hiv_ref_filepath
PT1-1,PT1,1,fastas/PT1-1.fa,sites/PT1-1.csv,references/refA.fa
PT1-2,PT1,2,fastas/PT1-2.fa,sites/PT1-2.csv,references/refB.fa
PT2-3,PT2,3,fastas/PT2-3.fa,sites/PT2-3.csv,references/refA.fa
```

Each row is one sample. Paths may be absolute or relative to the CSV's directory. The integration-site table must contain `READ`, `PCR_ID`, `chromosome`, and `INTEGRATION_SITE` (case-sensitive). FASTA headers must begin with the corresponding `READ` value, which may be followed by additional characters. The pipeline matches the longest `READ` prefix and adds `sample_id_` to merged FASTA headers.

The synthetic files in `tests/data/` illustrate the required format. See [data dictionary](docs/data_dictionary.md) for fields and outputs.

## Run the test data

Clone the repository, then run from its root:

```bash
git clone https://github.com/ekankaka/nf-clone-consensus.git
cd nf-clone-consensus

nextflow run main.nf \
    -c configs/local-singularity.config \
    --input tests/data/annotation.csv \
    --outdir results_test
```

If you use Apptainer or Docker, replace the config file with `configs/local-apptainer.config` or `configs/local-docker.config`. Only one runtime is needed. The default Nextflow working directory is `work/`.

The example contains PT1 and PT2, timepoints 1, 2, and 3, PCR duplicates, reads without integration sites, and an unmatched read.

## Run your own data

**On a PC (Docker example)**

```bash
nextflow run main.nf \
    -c configs/local-docker.config \
    --input /path/to/annotation.csv \
    --outdir /path/to/results
```

**On an HPC with SLURM (Singularity example)**

```bash
export NXF_SINGULARITY_CACHEDIR=/shared/container-cache

nextflow run main.nf \
    -c configs/slurm-singularity.config \
    -c /path/to/my_hpc.config \
    --input /research/study/annotation.csv \
    --outdir /research/study/results \
    -w /shared/scratch/study/work
```

Use `configs/slurm-apptainer.config` for Apptainer. Set the partition, account, time limits, and site-specific paths in your own `my_hpc.config`, outside the repository. The work directory and container cache must be accessible from compute nodes. Start Nextflow in an allowed HPC allocation if your login node prohibits long-running jobs.

For a repeated run, add `-resume`.

## Method

1. Merge sample FASTAs by patient and prefix headers with `sample_id_`.
2. Match reads to integration-site annotations. Define each clone by `(patient_id, chromosome, INTEGRATION_SITE)` across all timepoints.
3. Select one HIV reference per patient. If references differ, use the reference from the timepoint with the most reads; ties that cannot be resolved unambiguously stop the run.
4. Trim terminal Ns from integrated reads and align the fragments to the selected reference with MAFFT `--addfragments`. The reference sets alignment coordinates, but does not supply consensus bases.
5. Group PCR duplicates by `(sample_id, PCR_ID)`. Form a call for each PCR family and give each family one vote at each aligned position.
6. Call one clone consensus. With the default settings, a base needs at least 60% support and at least two informative PCR families for a multi-family clone. Positions without enough support are `N`; terminal gaps are missing coverage rather than deletions.

Single-read clones retain the trimmed read. Multiple reads from one PCR family are labeled as single-family consensuses and flagged for limited support. Internal deletion calls and ambiguous positions require review. These thresholds are provisional and may be changed in `nextflow.config`.

## Outputs

```text
results/
├── merged_fastas/
│   └── PT1_merged.fasta
├── clone_consensus/
│   └── PT1_consensus.fasta
├── clone_alignments/
│   └── PT1/
│       ├── PT1_clone000001_alignment.fasta
│       └── PT1_clone000002_alignment.fasta
├── unintegrated_sequences/
│   └── PT1_unintegrated.fasta
├── patient_sequence_counts.csv
├── clone_summary.csv
├── patient_clone_statistics.csv
├── reference_selection.csv
└── review/
    ├── read_assignments.csv
    ├── unmatched_fasta_reads.csv
    ├── unused_integration_records.csv
    └── unmatched_fastas/
```

`clone_alignments/` contains one FASTA per clone, with the aligned consensus followed by its terminal-N-trimmed contributing reads. Gaps (`-`) mark positions used to maintain alignment. Reads reversed by MAFFT are labeled `[reverse_complemented]`.

`clone_consensus/` contains integrated-site consensuses only. `unintegrated_sequences/` retains the original sample-prefixed sequence headers and nucleotide sequences. Ambiguous or unmatched reads are kept in `review/` and excluded from clone consensus construction.

Clone FASTA headers are:

```text
>patient_id_chromosome_INTEGRATION_SITE_n_reads_n_unique_reads
```

For example, `>PT1_chr1_100_5_4` is based on five accepted reads representing four distinct `(sample_id, PCR_ID)` families. Review the QC fields in `clone_summary.csv` before interpreting these sequences.

## Tests and limitations

Run the component tests:

```bash
python3 -m pytest -q tests/test_consensus.py
```

The tests cover read assignment, reference selection, PCR-family voting, uncovered positions, deletion calls, and summary reports. They use simulated alignments, **not** a full Nextflow/MAFFT execution. A full run of the synthetic data on a configured machine remains necessary before using the pipeline for a study.

These consensuses are computational summaries, not proof that a complete proviral genome has been reconstructed. Reads covering separate regions may yield an `N`-filled consensus. Review low-coverage clones, single-family consensuses, insertions/deletions, and reference-guided alignments. FASTA input does not provide per-base quality scores. Integration-site coordinates and chromosome names must be consistent across samples.

## Project layout

- `main.nf` — workflow
- `bin/` — Python analysis scripts
- `nextflow.config` — analysis defaults
- `configs/` — example local and SLURM configurations
- `tests/` — synthetic data and Python tests
- `docs/data_dictionary.md` — inputs, outputs, and QC field definitions

License: MIT. See [LICENSE](LICENSE).
