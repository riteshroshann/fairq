#!/usr/bin/env sh
# Reproduces every result and figure of the manuscript (about 60 minutes on one CPU core).
set -e
mkdir -p results
python verify.py
python circuit_spec.py
python gaps.py udg er reg3
python tracking.py udg er reg3
python dynamic_run.py
python figures.py
