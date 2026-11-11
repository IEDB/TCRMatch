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
./tcrmatch -i input_file -t num_threads [-s threshold] [-d /path/to/database] [-a | -r] [-k] [-m max_ram_gb] > output_file
```
- -i specifies the path to the input file
- -t specifies the number of threads to use if operating on more than 1 core
- -s optional parameter to specify the threshold (default is .97, in alignment with manuscript)
- -d optional parameter to specify the database (default: `data/IEDB_data.tsv`). Use `-d data/IEDB_data.schema2.tsv` for the [schema2 database](#databases). You can also point this to a file of newline separated CDR3b sequences to test your own private set
- -a specifies input file is in AIRR TSV format
- -r specifies input file is a TRUST4 report (`*_report.tsv`); only rows with a TRB V, D or J gene are used, and the TRUST4 columns are kept in the output
- -k keeps input sequences as they are. By default, a leading C together with a trailing F or W is trimmed from text and TRUST4 input, to match the trimmed database sequences (AIRR input is never trimmed)
- -m sets the maximum amount of RAM that TCRMatch can use (in GB. E.g.: -m 4 means 4GB. Default: 2GB)
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
-  Output file has 8 columns in TSV format.
-  First column is the user provided input sequence (`trimmed_input_sequence` when trimming is on).
-  Second column contains a CDR3beta sequence from the IEDB scoring higher than the threshold
-  The third column contains the score based on the kernel similarity metric (ranging from threshold to 1, which is an exact match)
-  The fourth column contains the receptor group (these can be used to reference the sequences back to the IEDB. If the receptor group is 12843, you can find all the receptors associated by visiting https://iedb.org/receptor/12843)
-  The fifth column contains the epitopes associated with the matched receptor group
-  The sixth column contains the name of the source antigen for each epitope
-  The seventh column contains the name of the source organism for each epitope
-  The eighth column contains the species of the matched receptor (schema2 database only; empty otherwise)

With the schema2 database, a matched sequence shared by several receptor groups gives one line per receptor group, and the epitope, antigen and organism lists are `|`-separated and aligned. With the legacy database, the lists are comma-separated, with the caveats described under [Databases](#databases). Antigen, organism and species information is not available for every receptor.

Example with the schema2 database:

| trimmed_input_sequence | match_sequence | score | receptor_group | epitope | antigen | organism | receptor_species |
|---|---|---|---|---|---|---|---|
| ASSSRSSYEQY | ASSSRSSYEQY | 1.000000 | 222969 | GILGFVFTL | Matrix protein 1 | Influenza A virus | Homo sapiens |
| ASSSRSSYEQY | ASSSRSSYEQY | 1.000000 | 223097 | GILGFVFTL | Matrix protein 1 | Influenza A virus | Homo sapiens |
| ASSSRSSYEQY | ASSARSSYEQY | 0.982682 | 39732 | GILGFVFTL\|RAKFKQLL | Matrix protein 1\|Lytic switch protein BZLF1 | Influenza A virus\|Human herpesvirus 4 (Epstein Barr virus) | Homo sapiens |
| ASSSRSSYEQY | ASSSSRSSYEQY | 0.978749 | 9409 | NLVPMVATV | HCMVUL83 | Human herpesvirus 5 (Human cytomegalovirus) | Homo sapiens |
| ASSSRSSYEQY | ASSRSSYEQY | 0.975811 | 45711 | AVFDRKSDAK | Epstein-Barr nuclear antigen 4 | Human herpesvirus 4 (Epstein Barr virus) | Homo sapiens |
| ASSSRSSYEQY | ASSERSSYEQY | 0.972149 | 21200 | AMKRHGLDNYREYSLGN | | | Mus musculus |

## How it works:
- TCRMatch implements a similarity metric as defined in [arXiv:1205.6031v2](https://arxiv.org/abs/1205.6031v2)
- Each input sequence is compared to all TCR CDR3beta sequences characterized in the IEDB
- Sequences that have a  score greater than the threshold (.97 by default) are returned as matches

## References:
[Chronister, William D. et al. "Tcrmatch: Predicting T-Cell Receptor Specificity Based On Sequence Similarity To Previously Characterized Receptors". Frontiers In Immunology, vol 12, 2021. Frontiers Media SA, doi:10.3389/fimmu.2021.640725.](https://www.frontiersin.org/articles/10.3389/fimmu.2021.640725/full)

["Towards a Mathematical Foundation of Immunology and Amino Acid Chains" arXiv:1205.6031v2](https://arxiv.org/abs/1205.6031v2)

