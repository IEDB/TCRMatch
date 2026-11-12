# TCRMatch
#### Please open any issues on Github


## Requirements:
- Linux OS and the g++ compiler (on macOS, use Homebrew GCC: `make CC=g++-14`, or whichever version you have; Apple clang does not support OpenMP)
- Python 3 (to generate certain output)
- curl (to download the database)

## Installation:
Please download TCRMatch from the latest [release](https://github.com/IEDB/TCRMatch/releases/) and compile from source by running:
```shell
make
```
Then download the IEDB database, which is not included in the repository:
```shell
cd scripts && ./update.sh
```
`INSTALL.sh` runs both steps.

To run the regression tests:
```shell
make test
```


## Usage:
### Input  
-  A text file containing newline separated CDR3beta amino acid sequences
  ```
ASSLRWGDEQY
SADRVGNTLY
ASSFGREQY
ASSQDQVQDTQY
ASSPGQGPYEQY
ATSSGTGAYEQY
ASEAGGSYEQY
ASGDRGADTQY
ASSQAGAYEQY
  ```

Alternatively, you may upload files using the AIRR tsv format specified [here](https://docs.airr-community.org/en/stable/datarep/rearrangements.html) (the `cdr3_aa` column is used), or a TRUST4 report.

To generate the initial scores, please run
```shell
./tcrmatch -i input_file -t num_threads [-s threshold] [-d /path/to/database] [-a | -r] [-k] [-m max_ram_gb] [-p per_group_file] [-c collapsed_file] > output_file
```
- -i specifies the path to the input file
- -t specifies the number of threads to use if operating on more than 1 core
- -s optional parameter to specify the threshold (default is .97, in alignment with manuscript)
- -d optional parameter to specify the database (default: `data/IEDB_data.tsv`). Use `-d data/IEDB_data.schema2.tsv` for the [schema2 database](#databases). You can also point this to a file of newline separated CDR3b sequences to test your own private set
- -a specifies input file is in AIRR TSV format
- -r specifies input file is a TRUST4 report (`*_report.tsv`); only rows with a TRB V, D or J gene are used, and the TRUST4 columns are kept in the output
- -k keeps input sequences as they are. By default, a leading C together with a trailing F or W is trimmed from text and TRUST4 input, to match the trimmed database sequences (AIRR input is never trimmed)
- -m sets the maximum amount of RAM that TCRMatch can use (in GB. E.g.: -m 4 means 4GB. Default: 2GB)
- -p also writes the matches with [one line per receptor group](#per-receptor-group-output--p) to the named file (schema2 database only)
- -c also writes the matches [collapsed by epitope](#collapsed-output--c) to the named file (schema2 database only)
- -v prints the version

Input sequences containing characters other than the 20 standard amino acids are skipped with a warning.

### Databases
To download the latest database from the IEDB, run
```shell
cd scripts && ./update.sh             # data/IEDB_data.tsv (legacy format, the default)
cd scripts && ./update.sh --schema2   # data/IEDB_data.schema2.tsv
```
The previous copy is moved to the `data_archive` folder in the top directory and renamed with the current day_month_year (e.g. `10_11_2026_IEDB_data.tsv`). If the download fails, the current database is left unchanged.

Two database formats are available:
- **`IEDB_data.tsv` (legacy).** One line per CDR3b sequence. When several receptor groups share a CDR3b, the line shows one receptor group, but its epitopes are pooled from all of them. Missing antigens and organisms are left out of the comma-separated lists, so those lists may not line up with the epitopes.
- **`IEDB_data.schema2.tsv`.** One line per receptor group, with that group's own epitopes. Lists are separated by `|`, and entry *n* of the epitope, antigen and organism lists always describe the same epitope (missing values are empty). Only alpha-beta receptors are included, and a `receptor_species` column is added. Requires TCRMatch 1.4 or later. **Recommended.**

`data/CEDAR_data.tsv` (cancer-related receptors, legacy format) is also included.

TCRMatch reads database columns by header name, so it also accepts custom databases with the header `trimmed_seq, original_seq, receptor_group, epitopes, source_organisms, source_antigens` (and optionally `receptor_species`), or a plain list of CDR3b sequences with no header.


### Output
TCRMatch writes TSV to standard output, with **one line per input sequence, matched sequence and untrimmed matched sequence**. With the schema2 database, each line lists every receptor group with that CDR3b, keeping each group's epitopes separate. Two optional files give other views of the same matches: [one line per receptor group](#per-receptor-group-output--p) (`-p`) and [one line per epitope](#collapsed-output--c) (`-c`).

Columns:
-  `input_sequence`: the input sequence (`trimmed_input_sequence` when trimming is on)
-  `match_sequence`: a CDR3b sequence from the IEDB scoring above the threshold, as trimmed for matching
-  `score`: the kernel similarity score (from the threshold up to 1, an exact match)
-  `untrimmed_match_sequence`: the matched CDR3b as recorded in the IEDB, before trimming. Most are identical to `match_sequence`; others include the conserved C and F/W or extra residues. When IEDB recorded a sequence in several forms, each form gets its own line
-  `receptor_groups`: the receptor groups with this CDR3b, `;`-separated. They can be looked up in the IEDB: receptor group 12843 is https://iedb.org/receptor/12843
-  `epitopes`, `antigens`, `organisms`: the epitopes each receptor group recognizes, with their source antigens and organisms
-  `receptor_species`: the species of each receptor group
-  `n_receptor_groups`, `n_epitopes`, `n_antigens`, `n_organisms`, `n_receptor_species`: placed before each list, the number of distinct non-empty values in it

The list columns have one `;`-separated entry per receptor group, in the same order as `receptor_groups`. Each entry is that group's `|`-separated list, so entry *n* of the epitope, antigen and organism lists describe the same epitope. To read a line, split on `;` to get the receptor groups, then on `|` to get the epitopes. A `;` inside a name (e.g. `UL123; IE1`) is written as `,`. Antigen, organism and species information is not available for every receptor (empty entries).

Example with the schema2 database (one line, shown as a list):

| Column | Value |
|---|---|
| trimmed_input_sequence | ASSIRAAETQY |
| match_sequence | ASSIRAAETQY |
| score | 1.000000 |
| untrimmed_match_sequence | ASSIRAAETQY |
| n_receptor_groups | 6 |
| receptor_groups | 1517;1525;8679;40437;40438;40439 |
| n_epitopes | 2 |
| epitopes | GILGFVFTL;GILGFVFTL;GILGFVFTL\|LLLDRLNQL;GILGFVFTL;GILGFVFTL;GILGFVFTL |
| n_antigens | 2 |
| antigens | Matrix protein 1;Matrix protein 1;Matrix protein 1\|Nucleoprotein;Matrix protein 1;Matrix protein 1;Matrix protein 1 |
| n_organisms | 2 |
| organisms | Influenza A virus;Influenza A virus;Influenza A virus\|SARS coronavirus BJ01;Influenza A virus;Influenza A virus;Influenza A virus |
| n_receptor_species | 1 |
| receptor_species | Homo sapiens;Homo sapiens;Homo sapiens;Homo sapiens;Homo sapiens;Homo sapiens |

With the legacy database (or the CEDAR file), each line has the single receptor group the database row names, the lists are copied as they are in the file (comma-separated, with the caveats described under [Databases](#databases)), and the `n_*` columns are empty because those lists cannot be counted reliably. In TRUST4 mode, the TRUST4 columns come first. An input sequence that appears more than once in the input is reported once.

### Per-receptor-group output (-p)
With `-p per_group.tsv` and the schema2 database, TCRMatch also writes one line per input sequence, matched sequence and receptor group. A CDR3b shared by many receptor groups gives many lines (194 for the exact match to ASSIRSSYEQY).

-  `input_sequence`, `match_sequence`, `score`, `untrimmed_match_sequence`: as in the standard output
-  `receptor_group`: one receptor group
-  `epitopes`, `antigens`, `organisms`: that group's `|`-separated lists, aligned by position
-  `n_epitopes`, `n_antigens`, `n_organisms`: the number of distinct non-empty values in each list
-  `receptor_species`: the group's species

| trimmed_input_sequence | match_sequence | score | untrimmed_match_sequence | receptor_group | n_epitopes | epitopes | n_antigens | antigens | n_organisms | organisms | receptor_species |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ASSIRAAETQY | ASSIRAAETQY | 1.000000 | ASSIRAAETQY | 1517 | 1 | GILGFVFTL | 1 | Matrix protein 1 | 1 | Influenza A virus | Homo sapiens |
| ASSIRAAETQY | ASSIRAAETQY | 1.000000 | ASSIRAAETQY | 8679 | 2 | GILGFVFTL\|LLLDRLNQL | 2 | Matrix protein 1\|Nucleoprotein | 2 | Influenza A virus\|SARS coronavirus BJ01 | Homo sapiens |

### Collapsed output (-c)
With `-c collapsed.tsv` and the schema2 database, TCRMatch also writes a summary with one line per input sequence, matched sequence and epitope. Each receptor group of the matched sequence counts once for every epitope it recognizes.

-  `input_sequence`, `match_sequence`, `score`: as in the standard output
-  `epitope`: a single epitope
-  `antigen`, `organism`: its source antigen and organism (when receptor groups disagree, the label given by the most groups, preferring a non-empty one)
-  `n_receptor_groups`, `receptor_groups`: how many receptor groups with the matched CDR3b recognize this epitope, and which (`|`-separated)
-  `n_receptor_species`, `receptor_species`: their distinct species (`|`-separated)

For each match, epitopes are listed by `n_receptor_groups`, highest first.

| trimmed_input_sequence | match_sequence | score | epitope | antigen | organism | n_receptor_groups | receptor_groups | n_receptor_species | receptor_species |
|---|---|---|---|---|---|---|---|---|---|
| ASSIRAAETQY | ASSIRAAETQY | 1.000000 | GILGFVFTL | Matrix protein 1 | Influenza A virus | 6 | 1517\|1525\|8679\|40437\|40438\|40439 | 1 | Homo sapiens |
| ASSIRAAETQY | ASSIRAAETQY | 1.000000 | LLLDRLNQL | Nucleoprotein | SARS coronavirus BJ01 | 1 | 8679 | 1 | Homo sapiens |

In TRUST4 mode, the `-p` and `-c` files are keyed on the CDR3 sequence and do not repeat the TRUST4 columns. Both stop with an error if the database is in the legacy format, because the legacy file pools epitopes across receptor groups.

## How it works:
- TCRMatch implements a similarity metric as defined in [arXiv:1205.6031v2](https://arxiv.org/abs/1205.6031v2)
- Each input sequence is compared to all TCR CDR3beta sequences characterized in the IEDB
- Sequences that have a  score greater than the threshold (.97 by default) are returned as matches

## References:
[Chronister, William D. et al. "Tcrmatch: Predicting T-Cell Receptor Specificity Based On Sequence Similarity To Previously Characterized Receptors". Frontiers In Immunology, vol 12, 2021. Frontiers Media SA, doi:10.3389/fimmu.2021.640725.](https://www.frontiersin.org/articles/10.3389/fimmu.2021.640725/full)

["Towards a Mathematical Foundation of Immunology and Amino Acid Chains" arXiv:1205.6031v2](https://arxiv.org/abs/1205.6031v2)

