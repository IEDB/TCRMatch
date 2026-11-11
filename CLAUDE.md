# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

TCRMatch is a single-binary C++11/OpenMP CLI that scores input TCR CDR3β amino-acid sequences against a reference database of characterized receptors (IEDB, or the bundled CEDAR set). It reports every pair whose normalized kernel similarity (arXiv:1205.6031) exceeds a threshold, defaulting to 0.97. Reference: Chronister et al., *Frontiers in Immunology* 2021.

The IEDB database files are built weekly by a separate repo, `update-repspec` (gitlab.lji.org: `iedb/tools/djangotools-deps/update-repspec`), and published at https://downloads.iedb.org/misc/TCRMatch/.

## Build, test and run

```shell
make                    # Linux: g++ with -fopenmp
make CC=g++-16          # macOS: Apple clang rejects -fopenmp, so use Homebrew GCC
make test               # (add the same CC=... on macOS) runs tests/run_tests.sh
make clean
```

- **Single translation unit:** `src/main.cpp` does `#include "tcrmatch.cpp"`, and the Makefile compiles only `main.cpp`. Don't add `tcrmatch.cpp` to the compile line, because that produces duplicate symbols. New functions go in `tcrmatch.cpp` (or a new file included the same way).
- **Databases aren't in git.** `data/IEDB_data*.tsv` is gitignored.
  - `cd scripts && ./update.sh` downloads `data/IEDB_data.tsv`, and `./update.sh --schema2` downloads `data/IEDB_data.schema2.tsv`.
  - The script downloads first, then moves the old copy to `data_archive/`.
  - `INSTALL.sh` and the Dockerfile run it. `data/IEDB_data.tsv` is the default `-d` path, resolved relative to the CWD.
  - `data/CEDAR_data.tsv` (legacy format) is committed and works offline.
- **Tests:** `tests/run_tests.sh` runs the binary on fixtures in `tests/data/` (legacy, schema2 and headerless databases; text, AIRR and TRUST4 inputs) and compares the output byte for byte with `tests/expected/`.
  - It runs with `-t 1` so the output order is deterministic.
  - It also checks reordered columns with Windows line endings, and that a missing required column produces an error.
  - `TCRMATCH=path` tests another binary. `UPDATE=1` regenerates the expected files; review that diff, because it is the regression baseline.
- There is no linter. Format code with the repo's `.clang-format` (LLVM-based): `clang-format -i src/*.cpp`.

## CLI flags (from `getopt` in `main.cpp`)

`-i` input file and `-t` thread count are required. Optional flags: `-s` threshold, `-d` database path, `-m` memory cap in GB (default 2), `-a` AIRR TSV input (reads the `cdr3_aa` column), `-r` TRUST4 report input, `-k` disable trimming, `-v` print version.

## Database formats

Both formats start with the header `trimmed_seq, original_seq, receptor_group, epitopes, source_organisms, source_antigens`.

- **Legacy `IEDB_data.tsv`.** One row per CDR3β. When receptor groups share a CDR3β, epitopes are pooled across them under one `receptor_group`. Lists are `,`-separated and may not line up. Data lines end with a stray tab.
- **`IEDB_data.schema2.tsv`.** One row per receptor group, alpha-beta receptors only.
  - Lists are `|`-separated and stay in step, with empty placeholders for missing values.
  - A 7th column, `receptor_species`, holds an organism name.
  - TCRMatch ≤ 1.3 segfaults on this file, because the old parser used a fixed `values[6]` array.
- **Headerless.** A plain list of CDR3β sequences, or rows in the legacy column order without a header. Also supported.

## Pipeline architecture

1. **Load DB:** the DB file is parsed twice.
   - `create_IEDB_map` builds a map from trimmed sequence to metadata rows. If the first field of line 1 is `trimmed_seq`, columns are located **by header name**, unknown columns are ignored, and a missing required column is a hard error. Otherwise columns are read positionally in the legacy order. Missing fields become empty, so it's bounds-safe, and line endings are normalized with `erase_newline`.
   - `read_IEDB_data` returns the **distinct** trimmed sequences (the first whitespace token of each line) to score against. `trimmed_seq` must therefore stay the first column.
   - Sequences containing characters outside `ARNDCQEGHILKMFPSTWYV` cause a hard failure for the DB.
2. **Parse input** into `peptide` structs (`peplist1`). There are three modes: plain text, AIRR (`-a`) and TRUST4 (`-r`).
   - TRUST4 requires an exact header and keeps only rows with a `TRB*` V/D/J gene.
   - Input sequences with invalid characters are skipped with a warning.
3. **Trimming:** `trim()` removes a leading `C` together with a trailing `F` or `W`, so input matches the pre-trimmed DB sequences.
   - It applies to text and TRUST4 input but **not to AIRR input**. `-k` turns it off.
   - When trimming is on, the output header's first column becomes `trimmed_input_sequence`.
4. **Scoring:** `fmatrix_k1()` computes the global `k1` matrix from the hardcoded BLOSUM62 `blm_qij` frequencies with `P_BETA`.
   - `k3_sum` sums the k-mer products for k=1..30 using a stack array `float[31][31][31]`, so it assumes sequences are ≤30 residues.
   - Each peptide's self-score (`aff`) is precomputed, and the final score is `k3(a,b)/sqrt(aff_a*aff_b)`.
5. **Matching:** `find_matches` runs an OpenMP all-vs-all loop that writes into per-thread result vectors. It increments a global atomic `current_memory_usage` and calls `exit(0)` once `MEMORY_LIMIT` is exceeded. Next, `remove_duplicate_matches` merges and dedups the vectors, and `multi_calc_k3` expands each match into one output line per DB metadata row (one per receptor group with schema2).
6. **Output** is TSV on stdout: `input_sequence, match_sequence, score, receptor_group, epitope, antigen, organism, receptor_species`.
   - `receptor_species` is empty with legacy or headerless databases.
   - In TRUST4 mode, each TRUST4 row is prepended to its matches.
   - Diagnostics go to stderr.

## Release

- The version is the `VERSION` macro at the top of `src/main.cpp`. Keep the help text and README in sync with it, as earlier commits did.
- Publishing a GitHub release triggers `.github/workflows/harbor.yml`, which builds and pushes `harbor.lji.org/iedb-public/tcrmatch`. The `Dockerfile` `git clone`s `IEDB/TCRMatch` master instead of using the local build context. The image therefore reflects pushed master, not the release tag or local changes.
