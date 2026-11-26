# Changelog

## Latest update

- Retired legacy local/SLURM execution config files; use composable `-profile` options instead.
- Added optional orientation checking of individual reads against each patient's selected reference, including single-read clones, with QC flags.
- Published merged FASTA and clone-to-read alignment outputs.
- Simplified the README and code comments.
- Added the corresponding output definitions to the data dictionary.

## Public GitHub version

- Added local and SLURM execution configurations.
- Added PT1/PT2 synthetic examples, tests, a data dictionary and a license.

## Earlier development (v3–v3.1)

- Grouped reads by patient and integration site across timepoints.
- Added reference-guided fragment alignment, PCR-family weighting and consensus QC.
- Separated unintegrated reads and unresolved assignments.
- Standardized output filenames.
