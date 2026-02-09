#!/usr/bin/env bash
set -euo pipefail

# Run tests (produces .gcda files)
fpm @test-gfortran-coverage test

# Collect only directories that actually contain .gcda files
COV_DIRS=()
while IFS= read -r -d '' d; do
  COV_DIRS+=("$d")
done < <(find build -path '*/copa' -type d -print0)

# Filter to those with .gcda files
VALID_DIRS=()
for d in "${COV_DIRS[@]}"; do
  if compgen -G "$d/*.gcda" > /dev/null; then
    VALID_DIRS+=("$d")
  fi
done

if [ "${#VALID_DIRS[@]}" -eq 0 ]; then
  echo "ERROR: No .gcda files found anywhere under build/"
  exit 1
fi

# Build lcov command safely
LCOV_ARGS=()
for d in "${VALID_DIRS[@]}"; do
  LCOV_ARGS+=("-d" "$d")
done

# Capture coverage
lcov -c \
     "${LCOV_ARGS[@]}" \
     -b "$(pwd)" \
     -o coverage.info

# Keep only project sources
lcov --extract coverage.info "$(pwd)/src/*.f90" -o coverage.filtered.info

# Generate HTML, ignoring format errors
genhtml --ignore-errors format coverage.filtered.info -o coverage-html

# Print summary
echo
lcov --list coverage.filtered.info

