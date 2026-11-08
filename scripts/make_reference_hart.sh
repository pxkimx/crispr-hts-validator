#!/bin/bash
# Reference outputs from the OFFICIAL Hart lab tools, run once (agreed in SPEC.md, "Validation"):
#   BAGEL2  https://github.com/hart-lab/bagel  (MIT)
#   drugZ   https://github.com/hart-lab/drugz  (MIT)
# Their own virtual environment (numpy<2, pandas 2.1.4: the versions this 2021 code was written for) lives in
#   ~/Library/Application Support/CRISPR HTS Validator/reference-venv-hart
# Outputs go to test-data/reference/bagel and test-data/reference/drugz, run from those folders with relative paths.
set -euo pipefail
TOOL="$(cd "$(dirname "$0")/.." && pwd)"
HOME_DIR="$HOME/Library/Application Support/CRISPR HTS Validator"
DL="$HOME_DIR/downloads"
RV="$HOME_DIR/reference-venv-hart"
PY="$RV/bin/python"
for r in bagel drugz; do [ -d "$DL/$r" ] || git clone -q --depth 1 "https://github.com/hart-lab/$r.git" "$DL/$r"; done
if [ ! -x "$PY" ]; then
  /usr/local/bin/python3.12 -m venv "$RV"
  "$RV/bin/pip" install -q --upgrade pip
  "$RV/bin/pip" install -q "numpy<2" "pandas==2.1.4" scipy click scikit-learn six
fi
K=../../hart-bagel/tkov3-example/HAP1-TKOv3-EXAMPLE.txt
E=../../hart-bagel/CEGv2.txt
N=../../hart-bagel/NEGv1.txt

mkdir -p "$TOOL/test-data/reference/bagel"
cd "$TOOL/test-data/reference/bagel"
{
  echo "# Reference outputs made $(date '+%Y-%m-%d %H:%M'); BAGEL commit $(git -C "$DL/bagel" rev-parse --short HEAD)"
  echo "# $("$PY" -c 'import numpy,pandas,scipy;print("numpy",numpy.__version__,"pandas",pandas.__version__,"scipy",scipy.__version__)')"
  echo "# run from test-data/reference/bagel; paths are relative to it"
} > COMMANDS.txt
run() { echo "$*" | sed "s#$DL/##g" >> COMMANDS.txt; "$@" > /dev/null 2>> bagel.log; }
: > bagel.log
run "$PY" "$DL/bagel/BAGEL.py" fc -i "$K" -o tkov3 -c T0_-
for t in T3 T18 T18_Starved; do
  case $t in
    T3) c=T3_A,T3_B,T3_C ;;
    T18) c=T18_A,T18_B,T18_C ;;
    T18_Starved) c=T18_A_Starved,T18_B_Starved,T18_C_Starved ;;
  esac
  run "$PY" "$DL/bagel/BAGEL.py" bf -i tkov3.foldchange -o "tkov3_$t.bf" -e "$E" -n "$N" -c "$c" -s 1
  run "$PY" "$DL/bagel/BAGEL.py" pr -i "tkov3_$t.bf" -o "tkov3_$t.pr" -e "$E" -n "$N"
done

mkdir -p "$TOOL/test-data/reference/drugz"
cd "$TOOL/test-data/reference/drugz"
{
  echo "# Reference outputs made $(date '+%Y-%m-%d %H:%M'); drugZ commit $(git -C "$DL/drugz" rev-parse --short HEAD)"
  echo "# run from test-data/reference/drugz; paths are relative to it"
  echo "# TKOv3 T18 (untreated) vs T18 starved (treated): replicates paired A-A, B-B, C-C; and unpaired"
} > COMMANDS.txt
: > drugz.log
run2() { echo "$*" | sed "s#$DL/##g" >> COMMANDS.txt; "$@" 2>> drugz.log; }
run2 "$PY" "$DL/drugz/drugz.py" -i "$K" -o tkov3_starved.drugz.txt -f tkov3_starved.drugz_fc.txt \
     -c T18_A,T18_B,T18_C -x T18_A_Starved,T18_B_Starved,T18_C_Starved -q
run2 "$PY" "$DL/drugz/drugz.py" -i "$K" -o tkov3_starved_unpaired.drugz.txt \
     -c T18_A,T18_B,T18_C -x T18_A_Starved,T18_B_Starved,T18_C_Starved -unpaired -q
# no personal paths in anything kept; repeated warnings collapsed to one line each
for f in "$TOOL/test-data/reference/bagel/bagel.log" "$TOOL/test-data/reference/drugz/drugz.log"; do
  awk '!seen[$0]++' "$f" > "$f.tmp" && mv "$f.tmp" "$f"
done
sed -i '' "s#$HOME#~#g" "$TOOL/test-data/reference/bagel/"{bagel.log,COMMANDS.txt} "$TOOL/test-data/reference/drugz/"{drugz.log,COMMANDS.txt}
echo "BAGEL2 and drugZ reference outputs written."
