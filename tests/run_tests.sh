#!/bin/bash
# Regression tests for TCRMatch: run the binary on small fixtures and compare
# with the expected outputs in tests/expected/.
#
#   make test                       # or: bash tests/run_tests.sh
#   TCRMATCH=/path/to/tcrmatch bash tests/run_tests.sh
#   UPDATE=1 bash tests/run_tests.sh   # regenerate tests/expected/ (review the diff!)
#
# Fixtures (tests/data/):
#   legacy_db.tsv   rows of the legacy IEDB_data.tsv (6 columns, trailing tab)
#   schema2_db.tsv  rows of IEDB_data.schema2.tsv (7 columns, '|' lists,
#                   receptor_species); the last row (receptor group 999999999)
#                   is synthetic: empty species and empty list placeholders
#   plain_db.txt    headerless list of CDR3b sequences
#   query*.txt/tsv  text, AIRR and TRUST4 inputs (query.txt has an invalid
#                   sequence, which is skipped with a warning)
# per_group_*.tsv and collapsed_*.tsv in tests/expected/ are the -p (one line
# per receptor group) and -c (collapsed by epitope) files.
set -u
cd "$(dirname "$0")/.."
TCRMATCH="${TCRMATCH:-./tcrmatch}"
D=tests/data
E=tests/expected
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
pass=0
fail=0

check() {  # name, then tcrmatch arguments; -t 1 keeps output order deterministic
    local name="$1"
    shift
    "$TCRMATCH" "$@" -t 1 > "$TMP/$name.tsv" 2> "$TMP/$name.err"
    local status=$?
    if [[ "${UPDATE:-0}" == "1" ]]; then
        cp "$TMP/$name.tsv" "$E/$name.tsv"
        echo "updated $name"
    elif [[ $status -eq 0 ]] && cmp -s "$TMP/$name.tsv" "$E/$name.tsv"; then
        pass=$((pass + 1))
        echo "ok      $name"
    else
        fail=$((fail + 1))
        echo "FAIL    $name (exit $status)"
        diff "$E/$name.tsv" "$TMP/$name.tsv" | head -10
    fi
}

for db in legacy schema2 plain; do
    f="$D/${db}_db.tsv"
    [[ $db == plain ]] && f="$D/plain_db.txt"
    check "text_$db" -i "$D/query.txt" -d "$f"
    check "airr_$db" -i "$D/query_airr.tsv" -a -d "$f"
    check "trust4_$db" -i "$D/query_trust4.tsv" -r -d "$f"
done
check text_schema2_untrimmed -i "$D/query.txt" -k -d "$D/schema2_db.tsv"

check_view() {  # name, view flag (-p or -c), then tcrmatch arguments; compares that file
    local name="$1" flag="$2"
    shift 2
    "$TCRMATCH" "$@" -t 1 "$flag" "$TMP/$name.tsv" > "$TMP/$name.stdout" 2> "$TMP/$name.err"
    local status=$?
    if [[ "${UPDATE:-0}" == "1" ]]; then
        cp "$TMP/$name.tsv" "$E/$name.tsv"
        echo "updated $name"
    elif [[ $status -eq 0 ]] && cmp -s "$TMP/$name.tsv" "$E/$name.tsv"; then
        pass=$((pass + 1))
        echo "ok      $name"
    else
        fail=$((fail + 1))
        echo "FAIL    $name (exit $status)"
        diff "$E/$name.tsv" "$TMP/$name.tsv" | head -10
    fi
}

check_view collapsed_text_schema2 -c -i "$D/query.txt" -d "$D/schema2_db.tsv"
check_view collapsed_trust4_schema2 -c -i "$D/query_trust4.tsv" -r -d "$D/schema2_db.tsv"
check_view per_group_text_schema2 -p -i "$D/query.txt" -d "$D/schema2_db.tsv"
check_view per_group_trust4_schema2 -p -i "$D/query_trust4.tsv" -r -d "$D/schema2_db.tsv"

[[ "${UPDATE:-0}" == "1" ]] && exit 0

# -p and -c write separate files; the main output must not change.
unchanged=1
for v in collapsed per_group; do
    cmp -s "$TMP/${v}_text_schema2.stdout" "$E/text_schema2.tsv" || unchanged=0
    cmp -s "$TMP/${v}_trust4_schema2.stdout" "$E/trust4_schema2.tsv" || unchanged=0
done
if [[ $unchanged -eq 1 ]]; then
    pass=$((pass + 1)); echo "ok      views_stdout_unchanged"
else
    fail=$((fail + 1)); echo "FAIL    views_stdout_unchanged"
fi

# -p and -c refuse a legacy database (pooled, misaligned lists) before scoring.
for flag in -p -c; do
    rm -f "$TMP/legacy_view.tsv"
    "$TCRMATCH" -i "$D/query.txt" -t 1 -d "$D/legacy_db.tsv" "$flag" "$TMP/legacy_view.tsv" \
        > /dev/null 2> "$TMP/legacy_view.err"
    status=$?
    if [[ $status -ne 0 ]] && [[ ! -e "$TMP/legacy_view.tsv" ]] &&
            grep -q "needs a schema2 database" "$TMP/legacy_view.err"; then
        pass=$((pass + 1)); echo "ok      view_rejects_legacy $flag"
    else
        fail=$((fail + 1)); echo "FAIL    view_rejects_legacy $flag (exit $status)"
    fi
done

# Columns are read by header name: reordering them, adding an unknown column
# and using Windows line endings must not change the output.
awk -F'\t' -v OFS='\t' '{print $1, $7, (NR == 1 ? "some_new_column" : "x"), $3, $2, $6, $5, $4}' \
    "$D/schema2_db.tsv" | sed $'s/$/\r/' > "$TMP/reordered_db.tsv"
if "$TCRMATCH" -i "$D/query.txt" -t 1 -d "$TMP/reordered_db.tsv" 2>/dev/null \
        | cmp -s - "$E/text_schema2.tsv"; then
    pass=$((pass + 1)); echo "ok      reordered_columns_crlf"
else
    fail=$((fail + 1)); echo "FAIL    reordered_columns_crlf"
fi

# A header without a required column is a clear error, not a crash.
cut -f1-5 "$D/schema2_db.tsv" > "$TMP/missing_db.tsv"
"$TCRMATCH" -i "$D/query.txt" -t 1 -d "$TMP/missing_db.tsv" > /dev/null 2> "$TMP/missing.err"
status=$?
if [[ $status -ne 0 ]] && grep -q "missing the required column 'source_antigens'" "$TMP/missing.err"; then
    pass=$((pass + 1)); echo "ok      missing_required_column"
else
    fail=$((fail + 1)); echo "FAIL    missing_required_column (exit $status)"
fi

echo "$pass passed, $fail failed"
[[ $fail -eq 0 ]]
