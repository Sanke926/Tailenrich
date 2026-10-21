#!/usr/bin/env sh
set -eu
cd "$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)"
if [ "$#" -lt 2 ]; then
    echo "Usage: sh scripts/run_example.sh expression.csv groups.tsv [output_dir] [cores]" >&2
    exit 1
fi
exec Rscript scripts/run_example.R "$@"
