#!/bin/bash
# Make the reference outputs the validation tests compare against, by running the OFFICIAL tools once.
# Agreed with Paul H. Kim (SPEC.md, "Validation"). Everything runs in its own virtual environment under
#   ~/Library/Application Support/CRISPR HTS Validator/reference-venv
# and the downloaded sources stay in .../downloads. Outputs go to test-data/reference/<tool>/.
#
#   MAGeCK 0.5.9.5  https://sourceforge.net/projects/mageck/files/0.5/mageck-0.5.9.5.tar.gz  (BSD)
#   BAGEL2          https://github.com/hart-lab/bagel (MIT)            -- used from version 0.3
set -euo pipefail
TOOL="$(cd "$(dirname "$0")/.." && pwd)"
HOME_DIR="$HOME/Library/Application Support/CRISPR HTS Validator"
DL="$HOME_DIR/downloads"
RV="$HOME_DIR/reference-venv"
OUT="$TOOL/test-data/reference/mageck"
TD="$TOOL/test-data"
mkdir -p "$DL" "$OUT"

if [ ! -x "$RV/bin/mageck" ]; then
  [ -f "$DL/mageck-0.5.9.5.tar.gz" ] || curl -sSL -o "$DL/mageck-0.5.9.5.tar.gz" \
      "https://sourceforge.net/projects/mageck/files/0.5/mageck-0.5.9.5.tar.gz/download"
  rm -rf "$DL/mageck-src" && mkdir -p "$DL/mageck-src" && tar -xzf "$DL/mageck-0.5.9.5.tar.gz" -C "$DL/mageck-src"
  SRC="$(ls -d "$DL"/mageck-src/*/ | head -1)"
  [ -x "$RV/bin/python" ] || /usr/local/bin/python3.12 -m venv "$RV"
  "$RV/bin/pip" install -q --upgrade pip setuptools numpy scipy
  (cd "$SRC" && "$RV/bin/python" setup.py install) > "$DL/mageck-install.log" 2>&1
fi
export PATH="$RV/bin:$PATH"
# call through python: the "#!" line of the installed script breaks on the space in "Application Support"
mageck() { "$RV/bin/python" "$RV/bin/mageck" "$@"; }
cd "$OUT"
{
  echo "# Reference outputs made $(date '+%Y-%m-%d %H:%M') on $(sw_vers -productName) $(sw_vers -productVersion)"
  echo "# mageck $(mageck --version 2>/dev/null)  python $("$RV/bin/python" -V 2>&1 | cut -d' ' -f2)"
  echo "# run from test-data/reference/mageck; paths are relative to it"
} > COMMANDS.txt
run() { echo "$*" >> COMMANDS.txt; "$@" > /dev/null 2>> mageck.log; }
: > mageck.log

# 1. read counting (demo 2, raw FASTQ)
run mageck count -l "../../mageck-demo/demo2/library.txt" -n demo2 --sample-label L1,CTRL \
    --fastq "../../mageck-demo/demo2/test1.fastq" "../../mageck-demo/demo2/test2.fastq"
run mageck test -k demo2.count.txt -t L1 -c CTRL -n demo2
# 2. count tables (demo 1 and demo 3), as in the demos' run.sh
run mageck test -k "../../mageck-demo/demo1/sample.txt" -t HL60.final,KBM7.final -c HL60.initial,KBM7.initial -n demo1
run mageck test -k "../../mageck-demo/demo3/leukemia.new.csv" -t HL60.final,KBM7.final -c HL60.initial,KBM7.initial -n demo3
# 3. HAP1 TKOv3 (Hart lab), each timepoint against T0
K="../../hart-bagel/tkov3-example/HAP1-TKOv3-EXAMPLE.txt"
run mageck test -k "$K" -t T3_A,T3_B,T3_C -c T0_- -n tkov3_T3
run mageck test -k "$K" -t T18_A,T18_B,T18_C -c T0_- -n tkov3_T18
run mageck test -k "$K" -t T18_A_Starved,T18_B_Starved,T18_C_Starved -c T0_- -n tkov3_T18_Starved
# MAGeCK's R report templates (*.Rmd) are left in place and kept out of git by .gitignore
# keep the logs free of local folder names
sed -i '' "s#$HOME#~#g" ./*.log COMMANDS.txt 2>/dev/null || true
echo "Reference outputs written to $OUT"
