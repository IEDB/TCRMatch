#!/bin/bash
# Regression tests for TCRMatch: run the binary on small fixtures and compare
# with the expected outputs in tests/expected/.
#
#   make test                       # or: bash tests/run_tests.sh
#   TCRMATCH=/path/to/tcrmatch bash tests/run_tests.sh
#   UPDATE=1 bash tests/run_tests.sh   # regenerate tests/expected/ (review the diff!)
#
# Fixtures (tests/data/):
#   schema2_db.tsv  rows of IEDB_data.schema2.tsv (7 columns, '|' lists,
#                   receptor_species); the last row (receptor group 999999999)
#                   is synthetic: empty species and empty list placeholders
#   cedar_db.tsv    rows of CEDAR_data.schema2.tsv (8 columns: the 7 above
#                   plus antigen_types) matched by query.txt
#   legacy_db.tsv   rows of the legacy IEDB_data.tsv (TCRMatch <= 1.3), which
#                   must be rejected with a clear error
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

check text_schema2 -i "$D/query.txt" -d "$D/schema2_db.tsv"
check airr_schema2 -i "$D/query_airr.tsv" -a -d "$D/schema2_db.tsv"
check trust4_schema2 -i "$D/query_trust4.tsv" -r -d "$D/schema2_db.tsv"
check text_schema2_untrimmed -i "$D/query.txt" -k -d "$D/schema2_db.tsv"
# --species keeps only the receptor groups of the named species (any case).
check text_schema2_mouse -i "$D/query.txt" -d "$D/schema2_db.tsv" --species "mus MUSCULUS"

# CEDAR: no --antigen-type means no filtering; with it, only epitopes of
# those types are kept (any case), and groups left without one are dropped.
check text_cedar -i "$D/query.txt" -d "$D/cedar_db.tsv"
check text_cedar_viral -i "$D/query.txt" -d "$D/cedar_db.tsv" --antigen-type viral
check text_cedar_neoantigen_germline -i "$D/query.txt" -d "$D/cedar_db.tsv" --antigen-type "Neoantigen, GERMLINE"

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
# The per-group view shows each group's remaining epitopes, lists in step.
check_view per_group_cedar_viral -p -i "$D/query.txt" -d "$D/cedar_db.tsv" --antigen-type viral

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

expect_error() {  # name, expected stderr text, then tcrmatch arguments
    local name="$1" message="$2"
    shift 2
    "$TCRMATCH" "$@" > "$TMP/$name.out" 2> "$TMP/$name.err"
    local status=$?
    if [[ $status -ne 0 ]] && grep -qF "$message" "$TMP/$name.err"; then
        pass=$((pass + 1)); echo "ok      $name"
    else
        fail=$((fail + 1)); echo "FAIL    $name (exit $status)"; head -3 "$TMP/$name.err"
    fi
}

# Databases of TCRMatch <= 1.3 are rejected with a pointer to the new file,
# before scoring and before -p/-c files are created.
expect_error rejects_legacy_db "legacy format" -i "$D/query.txt" -t 1 -d "$D/legacy_db.tsv" -c "$TMP/legacy_view.tsv"
if [[ -e "$TMP/legacy_view.tsv" ]]; then
    fail=$((fail + 1)); echo "FAIL    rejects_legacy_db created the -c file"
fi
printf 'ASSSRSSYEQY\nASSPGLGTYEQY\n' > "$TMP/plain_db.txt"
expect_error rejects_headerless_db "has no header line" -i "$D/query.txt" -t 1 -d "$TMP/plain_db.txt"
expect_error missing_db "Cannot read the database file" -i "$D/query.txt" -t 1 -d "$TMP/no_such_file.tsv"
expect_error antigen_type_needs_cedar_db "needs a database with an antigen_types column" -i "$D/query.txt" -t 1 -d "$D/schema2_db.tsv" --antigen-type viral
expect_error unknown_antigen_type "Unknown antigen type 'cancer'" -i "$D/query.txt" -t 1 -d "$D/cedar_db.tsv" --antigen-type viral,cancer
expect_error no_antigen_type_matches "match the requested antigen types and species" -i "$D/query.txt" -t 1 -d "$D/cedar_db.tsv" --antigen-type germline --species "Mus musculus"
expect_error unknown_species "No receptors of the requested species" -i "$D/query.txt" -t 1 -d "$D/schema2_db.tsv" --species "Danio rerio"

# --list-species needs no input file.
if "$TCRMATCH" -l -d "$D/schema2_db.tsv" 2>/dev/null | cmp -s - "$E/list_species.tsv"; then
    pass=$((pass + 1)); echo "ok      list_species"
else
    fail=$((fail + 1)); echo "FAIL    list_species"
fi

# A species name also matches its strains.
awk -F'\t' -v OFS='\t' '$3 == "21200" {$7 = "Mus musculus C57BL/6"} {print}' \
    "$D/schema2_db.tsv" > "$TMP/strain_db.tsv"
if "$TCRMATCH" -i "$D/query.txt" -t 1 -d "$TMP/strain_db.tsv" --species "Mus musculus" 2>/dev/null \
        | sed 's/Mus musculus C57BL\/6/Mus musculus/g' | cmp -s - "$E/text_schema2_mouse.tsv"; then
    pass=$((pass + 1)); echo "ok      species_matches_strains"
else
    fail=$((fail + 1)); echo "FAIL    species_matches_strains"
fi

# Sequences over 30 residues would overflow k3_sum's buffer: input ones are
# skipped with a warning, a database one is an error.
long=ASSSRSSYEQYASSSRSSYEQYASSSRSSYEQ
printf 'ASSPGLGTYEQY\n%s\n' "$long" > "$TMP/long_query.txt"
"$TCRMATCH" -i "$TMP/long_query.txt" -t 1 -d "$D/schema2_db.tsv" > "$TMP/long.tsv" 2> "$TMP/long.err"
status=$?
if [[ $status -eq 0 ]] && grep -q "longer than 30 residues. Skipping." "$TMP/long.err" &&
        grep -q "^ASSPGLGTYEQY" "$TMP/long.tsv" && ! grep -q "$long" "$TMP/long.tsv"; then
    pass=$((pass + 1)); echo "ok      skips_long_input"
else
    fail=$((fail + 1)); echo "FAIL    skips_long_input (exit $status)"
fi
{ cat "$D/schema2_db.tsv"; printf '%s\t%s\t1\tGILGFVFTL\tInfluenza A virus\tMatrix protein 1\tHomo sapiens\n' "$long" "$long"; } > "$TMP/long_db.tsv"
expect_error rejects_long_db_sequence "longer than 30 residues" -i "$D/query.txt" -t 1 -d "$TMP/long_db.tsv"

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
