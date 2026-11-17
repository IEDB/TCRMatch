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
Then download the databases, which are not included in the repository:
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
./tcrmatch -i input_file -t num_threads [-s threshold] [-d /path/to/database] [-a | -r] [-k] [-m max_ram_gb] [-p per_group_file] [-c collapsed_file] [--species "Homo sapiens"] [--antigen-type neoantigen,viral,germline,other] > output_file
```
- -i specifies the path to the input file
- -t specifies the number of threads to use if operating on more than 1 core
- -s optional parameter to specify the threshold (default is .97, in alignment with manuscript)
- -d optional parameter to specify the [database](#databases) (default: `data/IEDB_data.schema2.tsv`)
- -a specifies input file is in AIRR TSV format
- -r specifies input file is a TRUST4 report (`*_report.tsv`); only rows with a TRB V, D or J gene are used, and the TRUST4 columns are kept in the output
- -k keeps input sequences as they are. By default, a leading C together with a trailing F or W is trimmed from text and TRUST4 input, to match the trimmed database sequences (AIRR input is never trimmed)
- -m sets the maximum amount of RAM that TCRMatch can use (in GB. E.g.: -m 4 means 4GB. Default: 2GB)
- -p also writes the matches with [one line per receptor group](#per-receptor-group-output--p) to the named file
- -c also writes the matches [collapsed by epitope](#collapsed-output--c) to the named file
- --species only matches receptors of the given species, e.g. `--species "Homo sapiens"`. Separate several species with commas. Names are not case-sensitive, and a species also matches its strains (`Mus musculus` includes `Mus musculus C57BL/6`)
- --antigen-type only matches epitopes of the given CEDAR antigen types: `neoantigen`, `viral`, `germline` or `other`, separated by commas (e.g. `--antigen-type neoantigen,germline`). Epitopes of other types are left out of each receptor group, and receptor groups left with none are not matched. Without it, nothing is filtered. It needs the CEDAR database (`-d data/CEDAR_data.schema2.tsv`), the only one with antigen types
- -l, --list-species lists the receptor species in the database, with the number of receptor groups and CDR3b sequences of each, and exits (no `-i` or `-t` needed)
- -v prints the version

Input sequences containing characters other than the 20 standard amino acids, or longer than 30 residues, are skipped with a warning.

### Databases
Two databases are built weekly from the IEDB and CEDAR query APIs:
- `data/IEDB_data.schema2.tsv`: all receptors in the IEDB (the default)
- `data/CEDAR_data.schema2.tsv`: all receptors in [CEDAR](https://cedar.iedb.org), the Cancer Epitope Database and Analysis Resource. Use it with `-d data/CEDAR_data.schema2.tsv`. CEDAR includes every epitope tested in its references, not only cancer epitopes, so this database also has an `antigen_types` column with CEDAR's classification of each epitope, which `--antigen-type` filters on:
  - `neoantigen`: epitopes from mutated (tumor-specific) proteins
  - `viral`: epitopes from cancer-associated viruses such as EBV, CMV, HBV and HCV
  - `germline`: epitopes from unmutated self proteins, such as cancer-testis antigens
  - `other`: everything else, e.g. influenza and SARS-CoV-2 epitopes
  
  An epitope can have more than one type.

To download the latest versions, run
```shell
cd scripts && ./update.sh
```
The previous copies are moved to the `data_archive` folder in the top directory and renamed with the current day_month_year (e.g. `10_11_2026_IEDB_data.schema2.tsv`). If a download fails, that database is left unchanged.

Each database has one line per receptor group, with that group's own epitopes. Lists are separated by `|`, and entry *n* of the epitope, antigen and organism lists always describe the same epitope (missing values are empty). Only alpha-beta receptors are included. The CDR3b is the IEDB's calculated CDR3 (IMGT definition, without the conserved C and F/W); receptors for which the IEDB could not calculate one are left out. A `receptor_species` column gives each receptor group's species.

TCRMatch 1.4 no longer reads the legacy `IEDB_data.tsv` and `CEDAR_data.tsv` used by version 1.3 and earlier, which pooled epitopes across receptor groups that share a CDR3b. It stops with an error that points to `scripts/update.sh`.

Columns are read by header name. A custom database needs the header columns `trimmed_seq, original_seq, receptor_group, epitopes, source_organisms, source_antigens, receptor_species` (in any order; other columns are ignored), and `trimmed_seq` values of at most 30 residues.


### Output
TCRMatch writes TSV to standard output, with **one line per input sequence and matched sequence**. Each line lists every receptor group with that CDR3b, keeping each group's epitopes separate. Two optional files give other views of the same matches: [one line per receptor group](#per-receptor-group-output--p) (`-p`) and [one line per epitope](#collapsed-output--c) (`-c`).

Columns:
-  `input_sequence`: the input sequence (`trimmed_input_sequence` when trimming is on)
-  `match_sequence`: a CDR3b sequence from the IEDB scoring above the threshold, as trimmed for matching
-  `score`: the kernel similarity score (from the threshold up to 1, an exact match)
-  `untrimmed_match_sequence`: the database's `original_seq` for the matched CDR3b. In the IEDB database it is the same as `match_sequence`; a custom database can record another form, and each form gets its own line
-  `receptor_groups`: the receptor groups with this CDR3b, `;`-separated. They can be looked up in the IEDB: receptor group 12843 is https://iedb.org/receptor/12843
-  `epitopes`, `antigens`, `organisms`: the epitopes each receptor group recognizes, with their source antigens and organisms
-  `receptor_species`: the species of each receptor group
-  `n_receptor_groups`, `n_epitopes`, `n_antigens`, `n_organisms`, `n_receptor_species`: placed before each list, the number of distinct non-empty values in it

The list columns have one `;`-separated entry per receptor group, in the same order as `receptor_groups`. Each entry is that group's `|`-separated list, so entry *n* of the epitope, antigen and organism lists describe the same epitope. To read a line, split on `;` to get the receptor groups, then on `|` to get the epitopes. A `;` inside a name (e.g. `UL123; IE1`) is written as `,`. Antigen, organism and species information is not available for every receptor (empty entries).

Example (one line, shown as a list):

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

In TRUST4 mode, the TRUST4 columns come first. An input sequence that appears more than once in the input is reported once.

### Per-receptor-group output (-p)
With `-p per_group.tsv`, TCRMatch also writes one line per input sequence, matched sequence and receptor group. A CDR3b shared by many receptor groups gives many lines (193 for the exact match to ASSIRSSYEQY).

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
With `-c collapsed.tsv`, TCRMatch also writes a summary with one line per input sequence, matched sequence and epitope. Each receptor group of the matched sequence counts once for every epitope it recognizes.

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

In TRUST4 mode, the `-p` and `-c` files are keyed on the CDR3 sequence and do not repeat the TRUST4 columns.

## How it works:
- TCRMatch implements a similarity metric as defined in [arXiv:1205.6031v2](https://arxiv.org/abs/1205.6031v2)
- Each input sequence is compared to all TCR CDR3beta sequences characterized in the IEDB
- Sequences that have a  score greater than the threshold (.97 by default) are returned as matches

## References:
[Chronister, William D. et al. "Tcrmatch: Predicting T-Cell Receptor Specificity Based On Sequence Similarity To Previously Characterized Receptors". Frontiers In Immunology, vol 12, 2021. Frontiers Media SA, doi:10.3389/fimmu.2021.640725.](https://www.frontiersin.org/articles/10.3389/fimmu.2021.640725/full)

["Towards a Mathematical Foundation of Immunology and Amino Acid Chains" arXiv:1205.6031v2](https://arxiv.org/abs/1205.6031v2)

