#!/usr/bin/env bash
set -e
mkdir -p data figures
python3 scripts/make_bcs_table.py
python3 tests/test_abs.py
python3 tests/test_noise.py
python3 scripts/exp_universal.py
python3 scripts/exp_design.py
python3 scripts/exp_matched_points.py
python3 scripts/exp_calorimetry.py
python3 scripts/fig1.py
python3 scripts/fig2.py
python3 scripts/fig3.py
python3 scripts/fig4.py
python3 scripts/figS1.py
python3 scripts/make_numbers.py
echo "ALL DONE"
