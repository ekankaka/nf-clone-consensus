#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"

nextflow run main.nf \
    -c configs/local-singularity.config \
    --input tests/data/annotation.csv \
    --outdir results_test \
    "$@"
