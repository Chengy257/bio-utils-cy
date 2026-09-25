#!/usr/bin/env bash
# ==============================================================================
# Run the functional test suite (tests/test_*.py).
# Usage:  bash tests/run_tests.sh [pattern]
#   pattern  optional unittest pattern, e.g. "test_multi_file_join"
# Uses only the Python standard library (unittest); every bin/ script is
# exercised end-to-end through subprocess calls with synthetic data.
# ==============================================================================
set -euo pipefail
cd "$(dirname "$0")/.."
PATTERN="${1:-test_*.py}"
if [[ "$PATTERN" != *"*"* && "$PATTERN" != *"?"* ]]; then
    PATTERN="${PATTERN}*.py"   # bare name -> glob
fi
exec python3 -m unittest discover -s tests -p "$PATTERN" -v
