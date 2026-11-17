# Changelog

## Next update

- Publish patient-level merged FASTA files.
- Add consensus-versus-read alignments for each clone.
- Simplify repository description and comments.

## v3.1 documentation revision

- Edited comments, documentation, and setup instructions for clarity and consistency.
- No changes to the analytical code or output files.

## v3.1 — Output file naming

- Merged FASTA: `<patient_id>_merged.fasta`.
- Integrated clone FASTA: `<patient_id>_consensus.fasta`.
- Unintegrated FASTA: `<patient_id>_unintegrated.fasta`.
- No change to clone grouping, sequence processing, or consensus calls.

## v3 — Clone consensus and reference selection

- Added `best_hiv_ref_filepath` to the sample manifest and `PCR_ID` to integration-site records.
- Selected a reference for each patient, using the highest-read timepoint when references differ; unresolved ties are reported as errors.
- Grouped PCR duplicates using `(sample_id, PCR_ID)` and weighted each family once per aligned position.
- Trimmed terminal Ns from integrated reads and used MAFFT fragment alignment against the patient reference.
- Added coverage-aware consensus calls, clone QC statistics, and patient-level read accounting.
- Separated integrated consensuses, unintegrated reads, and reads requiring review.
- Added synthetic tests for preparation, consensus logic, and reporting.


## Public repository preparation

- Moved the Nextflow entry point to the repository root.
- Added portable local and SLURM configuration examples.
- Replaced synthetic identifiers with PT1/PT2 and timepoints 1, 2, and 3.
- Added output data dictionary, license, and GitHub installation instructions.
- No change to the scientific analysis logic.
