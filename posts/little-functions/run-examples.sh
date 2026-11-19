#!/usr/bin/env bash
set -euo pipefail
if [[ $# -ne 3 || "$1" != /* || "$2" != /* || "$3" != /* ]]; then
  echo 'usage: run-examples.sh /absolute/REPO /absolute/ROOT /absolute/OUTPUT' >&2
  exit 2
fi
REPO="$1"; ROOT="$2"; OUTPUT="$3"
if [[ ! -d "$REPO" || ! -d "$ROOT" || -e "$OUTPUT" ]]; then
  echo 'REPO and ROOT must exist; OUTPUT must be a new path' >&2
  exit 2
fi
REPO="$(realpath -e -- "$REPO")"
ROOT="$(realpath -e -- "$ROOT")"
OUTPUT="$(realpath -m -- "$OUTPUT")"
case "$OUTPUT/" in
  "$REPO/"*|"$ROOT/"*)
    echo 'OUTPUT must be outside REPO and ROOT' >&2
    exit 2 ;;
esac
mkdir -p "$OUTPUT/r" "$OUTPUT/python"
/usr/bin/Rscript --vanilla "$REPO/articles/little-functions/examples/read_saved.R" "$REPO" "$ROOT" "$OUTPUT/r"
"$ROOT/venv/bin/python" -B "$REPO/articles/little-functions/examples/read_saved.py" "$REPO" "$ROOT" "$OUTPUT/python"
