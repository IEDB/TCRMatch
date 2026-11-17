# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

TCRMatch is a single-binary C++11/OpenMP CLI that scores input TCR CDR3β amino-acid sequences against a reference database of characterized receptors (IEDB or CEDAR). It reports every pair whose normalized kernel similarity (arXiv:1205.6031) exceeds a threshold, defaulting to 0.97. Reference: Chronister et al., *Frontiers in Immunology* 2021.

The IEDB database files are built weekly by a separate repo, `update-repspec` (gitlab.lji.org: `iedb/tools/djangotools-deps/update-repspec`), and published at https://downloads.iedb.org/misc/TCRMatch/.

## Build, test and run

```shell
make                    # Linux: g++ with -fopenmp
make CC=g++-16          # macOS: Apple clang rejects -fopenmp, so use Homebrew GCC (any version)
make test               # (add the same CC=... on macOS) runs tests/run_tests.sh
make clean
```

- **Single translation unit:** `src/main.cpp` does `#include "tcrmatch.cpp"`, and the Makefile compiles only `main.cpp`. Don't add `tcrmatch.cpp` to the compile line, because that produces duplicate symbols. New functions go in `tcrmatch.cpp` (or a new file included the same way).
- **Databases aren't in git.** `data/IEDB_data*.tsv` is gitignored.
  - `cd scripts && ./update.sh` downloads `data/IEDB_data.schema2.tsv`, the default `-d` path (resolved relative to the CWD), and `data/CEDAR_data.schema2.tsv` (CEDAR references only; `-d data/CEDAR_data.schema2.tsv`).
  - For each file, the script downloads first, checks the header has `receptor_species`, then moves the old copy to `data_archive/`.
  - `INSTALL.sh` and the Dockerfile run it.
- **Tests:** `tests/run_tests.sh` runs the binary on fixtures in `tests/data/` (IEDB and CEDAR schema2 databases; text, AIRR and TRUST4 inputs) and compares the output byte for byte with `tests/expected/`.
  - It runs with `-t 1` so the output order is deterministic.
  - It also checks reordered columns with Windows line endings, a missing required column, the legacy and headerless databases being rejected, `--species`, `--antigen-type` and `--list-species`, and sequences over 30 residues.
  - `TCRMATCH=path` tests another binary. `UPDATE=1` regenerates the expected files; review that diff, because it is the regression baseline.
- There is no linter. Format code with the repo's `.clang-format` (LLVM-based): `clang-format -i src/*.cpp`.

## CLI flags (from `getopt` in `main.cpp`)

Parsed with `getopt_long`. `-i` input file and `-t` thread count are required. Optional flags: `-s` threshold, `-d` database path, `-m` memory cap in GB (default 2), `-a` AIRR TSV input (reads the `cdr3_aa` column), `-r` TRUST4 report input, `-k` disable trimming, `-p FILE` also write one line per receptor group to FILE, `-c FILE` also write the matches collapsed by epitope to FILE, `--species "A,B"` keep only receptor groups of those species, `--antigen-type a,b` keep only epitopes of those CEDAR antigen types (no filtering without it), `-v`/`--version` print version.

`-l`/`--list-species` prints each `receptor_species` in the database with its receptor group and CDR3b counts, then exits; it needs neither `-i` nor `-t`.

## Database format

Only the schema2 format is supported: `IEDB_data.schema2.tsv` and `CEDAR_data.schema2.tsv`, both built weekly by `update-repspec` from the IEDB and CEDAR query APIs. The header is `trimmed_seq, original_seq, receptor_group, epitopes, source_organisms, source_antigens, receptor_species`; the CEDAR file adds `antigen_types`.

- One row per receptor group, alpha-beta receptors only.
- The CDR3b is IEDB's calculated CDR3 (IMGT, without the conserved C and F/W), the same form TCRMatch's `trim()` gives input. `original_seq` equals `trimmed_seq` in the IEDB file. Groups without a calculated CDR3 are left out.
- Lists are `|`-separated and stay in step, with empty placeholders for missing values.
- `antigen_types` (CEDAR only) is a `|` list in step with `epitopes`; each entry is a `,` list of `neoantigen`, `viral`, `germline`, `other` (CEDAR's per-epitope flags; may be empty). It is used only for filtering and isn't printed.
- TCRMatch ≤ 1.3 segfaults on this file, because the old parser used a fixed `values[6]` array.
- The legacy `IEDB_data.tsv` (6 columns, pooled `,` lists) and headerless files are rejected by `create_IEDB_map` with an error pointing to `scripts/update.sh`. Don't add support back: the legacy file attributes epitopes to the wrong receptor groups.

## Pipeline architecture

1. **Load DB:** `create_IEDB_map` reads the file once.
   - It builds a map from trimmed sequence to metadata rows, and fills a vector of the **distinct** trimmed sequences to score, in file order (this order sets the output order).
   - Columns are located **by header name**, and unknown columns are ignored. A missing required column, no header, or no `receptor_species` column (the legacy file) is a hard error. Missing fields become empty, so it's bounds-safe, and line endings are normalized with `erase_newline`.
   - `--species` filters rows here (`species_match`: case-insensitive, and a name also matches `name + " " + strain`), so other species' sequences are never scored. Names that match nothing get a warning; if nothing is left, the run stops.
   - `--antigen-type` filters epitopes here too (`keep_antigen_types`): each row keeps only the epitopes with one of the requested types, with the epitope, organism and antigen lists cut in step, and a row with none left is dropped. Unknown type names and databases without an `antigen_types` column (the IEDB file) are errors.
   - DB sequences containing characters outside `ARNDCQEGHILKMFPSTWYV`, or longer than `MAX_SEQ_LEN` (30), cause a hard failure.
2. **Parse input** into `peptide` structs (`peplist1`). There are three modes: plain text, AIRR (`-a`) and TRUST4 (`-r`).
   - TRUST4 requires an exact header and keeps only rows with a `TRB*` V/D/J gene.
   - Input sequences with invalid characters, or longer than `MAX_SEQ_LEN` after trimming, are skipped with a warning.
3. **Trimming:** `trim()` removes a leading `C` together with a trailing `F` or `W`, so input matches the pre-trimmed DB sequences.
   - It applies to text and TRUST4 input but **not to AIRR input**. `-k` turns it off.
   - When trimming is on, the output header's first column becomes `trimmed_input_sequence`.
4. **Scoring:** `fmatrix_k1()` computes the global `k1` matrix from the hardcoded BLOSUM62 `blm_qij` frequencies with `P_BETA`.
   - `k3_sum` sums the k-mer products for k=1..30 using a stack array `float[31][31][31]`, which is why sequences are capped at `MAX_SEQ_LEN` = 30 residues.
   - Each peptide's self-score (`aff`) is precomputed, and the final score is `k3(a,b)/sqrt(aff_a*aff_b)`.
5. **Matching:** `find_matches` runs an OpenMP all-vs-all loop that writes into per-thread result vectors. It increments a global atomic `current_memory_usage` and calls `exit(0)` once `MEMORY_LIMIT` is exceeded. Next, `remove_duplicate_matches` merges and dedups the vectors. `multi_calc_k3` then builds the outputs from them; each view reports an (input, match) pair once (`distinct_matches`), even if the input repeats a sequence.
6. **Standard output** (`group_by_cdr3`) is TSV with one line per input, match and untrimmed CDR3β (the DB's `original_seq`; with the IEDB file that is one line per input and match). Its columns are `input_sequence, match_sequence, score, untrimmed_match_sequence, n_receptor_groups, receptor_groups, n_epitopes, epitopes, n_antigens, antigens, n_organisms, organisms, n_receptor_species, receptor_species`.
   - Each list column has one `;` entry per receptor group, in step. Each entry is that group's `|` list, unchanged.
   - Every `n_*` column comes just before its list and counts distinct non-empty values (`count_distinct`).
   - `without_semicolons` writes a `;` inside a name as `,`. Only the antigen "UL123; IE1" has one in current data.
   - In TRUST4 mode, each TRUST4 row is prepended to its matches. Diagnostics go to stderr.
7. **Per-receptor-group output (`-p FILE`)** (`per_receptor_group`): one line per input, match and receptor group, which was the stdout layout before. Its columns are `input_sequence, match_sequence, score, untrimmed_match_sequence, receptor_group, n_epitopes, epitopes, n_antigens, antigens, n_organisms, organisms, receptor_species`.
8. **Collapsed output (`-c FILE`)** (`collapse_by_epitope`): one line per input, match and epitope. Its columns are `input_sequence, match_sequence, score, epitope, antigen, organism, n_receptor_groups, receptor_groups, n_receptor_species, receptor_species`.
   - Each `|` list is split with empties kept, and each distinct receptor group counts once per epitope.
   - Antigen and organism come from the label given by the most groups, preferring a non-empty one. Epitopes are sorted by group count (descending), then by name.
9. **Rules shared by `-p` and `-c`.** Both are written in addition to stdout and never change it, through `open_summary_file` and `write_summary_file`.
   - `open_summary_file` opens the files before scoring, so a bad path fails fast.
   - In TRUST4 mode neither file has TRUST4 columns.

## Release

- The version is the `VERSION` macro at the top of `src/main.cpp`. Keep the help text and README in sync with it, as earlier commits did.
- Publishing a GitHub release triggers `.github/workflows/harbor.yml`, which builds and pushes `harbor.lji.org/iedb-public/tcrmatch`. The `Dockerfile` `git clone`s `IEDB/TCRMatch` master instead of using the local build context. The image therefore reflects pushed master, not the release tag or local changes.
