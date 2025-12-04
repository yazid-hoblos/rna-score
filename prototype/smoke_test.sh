#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

echo "Preparing data (smoke set)..."
python3 src/data_preparation.py \
  --ids-file data/examples/smoke_ids.txt \
  --output data/examples/prepared_smoke.csv \
  --max-distance 20 \
  --min-separation 4 \
  --atom "C3'"

echo "Training scoring tables..."
python3 src/training.py \
  --input data/examples/prepared_smoke.csv \
  --output-dir data/examples/training_output_smoke

echo "Scoring test structure (1EHZ)..."
python3 src/scoring.py \
  --structure data/examples/test/1EHZ.cif \
  --scores-dir data/examples/training_output_smoke

echo "Smoke test completed."
