#!/bin/sh
# Format HemoCell's own C++ sources with clang-format (style: .clang-format).
#
# Usage (from anywhere in the repository):
#   scripts/format.sh                 # format all HemoCell C++ files in place
#   scripts/format.sh core/hemoCell.cpp helper/
#   scripts/format.sh --check [paths] # only report files that need formatting
#
# For day-to-day work prefer formatting only what you changed:
#   git clang-format            # staged changes
#   git clang-format HEAD~1     # changes since the previous commit
#
# Set CLANG_FORMAT to choose the binary (default: clang-format, >= 14).
# Vendored (external/, palabos/) and historical (legacy/) code is skipped.
set -eu

CLANG_FORMAT=${CLANG_FORMAT:-clang-format}
if ! command -v "$CLANG_FORMAT" >/dev/null 2>&1; then
  echo "error: $CLANG_FORMAT not found (install clang-format >= 14 or set CLANG_FORMAT)" >&2
  exit 2
fi

cd "$(git rev-parse --show-toplevel)"

check=0
if [ "${1:-}" = "--check" ]; then
  check=1
  shift
fi

if [ "$#" -eq 0 ]; then
  set -- config core helper io mechanics tests examples cases tools hemocell.h
fi

# Tracked C/C++ sources under the given paths, excluding vendored/legacy code.
files=$(git ls-files -- "$@" \
  | grep -E '\.(c|cc|cpp|h|hh|hpp)$' \
  | grep -vE '^(external|legacy|palabos)/' || true)

if [ -z "$files" ]; then
  echo "no C/C++ files to format"
  exit 0
fi

if [ "$check" -eq 1 ]; then
  # Lists every file that would change; non-zero exit if any.
  echo "$files" | xargs "$CLANG_FORMAT" --dry-run --Werror
else
  echo "$files" | xargs "$CLANG_FORMAT" -i
  echo "formatted $(echo "$files" | wc -l) files"
fi
