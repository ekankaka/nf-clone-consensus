#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"

nextflow run main.nf \
    -profile local,singularity \
    --input tests/data/annotation.csv \
    --outdir results_test \
    "$@"
