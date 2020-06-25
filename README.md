# TCRMatch
### Authors: Austin Crinklaw, Will Chronister

## Requirements:
- Linux OS
- Python 3+
  - Python packages: numpy

## Installation:
TCRMatch can be downloaded through PyPI using the following pip command.
```shell
pip install TCRMatch
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

- Alternatively, you can upload a TSV file in the AIRR format [(see specifications here)](https://docs.airr-community.org/en/stable/datarep/rearrangements.html).
This can be specified by passing in the command line parameter ```-f airr```


### Commands
-  Using a .txt file as an input:
```shell
python -m TCRMatch match -i /path/to/input.txt -o /path/to/output.txt
```

There are a few other command line parameters:
- -t specifies the threshold (should be a value between 0 and 1)
- -p specifies the number of threads to use if operating on more than 1 core
- -f specifies the input file format, this is optional but needs to be specified as ```airr``` if using AIRR tsv format files

### API
The TCRMatch module also exposes a simple method for calculating the distance between two CDR3b sequences. This can be called by importing the module like so:
```Python
from tcrmatch_c import tcrmatch

res = tcrmatch("ASSQDRDTQY", "ASGDAGGGYEQY")
```
The tcrmatch method returns a tuple containing the following format (seq1, seq2, tcrmatch_score)
```Python
("ASSQDRDTQY", "ASGDAGGGYEQY", .74)
```
### Output  
-  Output file has 5 columns in TSV format. 
-  First column is the user provided input sequence.  
-  Second column contains a CDR3beta sequence from the IEDB scoring higher than the threshold
-  The third column contains the score based on the kernel similarity metric (ranging from threshold to .99)
- The fourth column contains the epitopes associated with the matched IEDB sequence. If there are multiple, they will be comma separated.
- The fifth column contains the receptor groups (these can be used to reference the sequences back to the IEDB. If the receptor group is 12843, you can find all the receptors associated by visiting https://iedb.org/receptor/12843)*

*Note that due to the trimming of IEDB CDR3b sequences, there may be multiple receptor groups that map to the same trimmed sequence. In these cases, the receptor groups will be comma separated and each link to their own unique receptor on the IEDB.

| input_sequence | match_sequence | score | epitopes      | receptor_group |
|----------------|----------------|-------|---------------|----------------|
| ASSQDRDTQY     | ASSGSGPLRGYT   | 0.77  | PKYVKQNTLKLAT | 26382          |
| ASSQDRDTQY     | ASSLVASNYGYT   | 0.78  | KLGGALQAK     | 54210          |
| ASSQDRDTQY     | AWSVPPAASYGYT  | 0.71  | KLGGALQAK     | 51332          |
| ASSQDRDTQY     | GASWGNTGQLY    | 0.77  | HGIRNASFI     | 20067          |
| ASSQDRDTQY     | ASSTIAGGYNEQF  | 0.77  | NEGVKAAW      | 27704          |
| ASSQDRDTQY     | ASRTRLDGYT     | 0.84  | ELAGIGILTV    | 29587          |

## How it works:
- TCRMatch implements a similarity metric as defined in [arXiv:1205.6031v2](https://arxiv.org/abs/1205.6031v2)
- Each input sequence is compared to all TCR CDR3beta sequences characterized in the IEDB
- Sequences that have a  score greater than the threshold are returned as matches

## References:
["Towards a Mathematical Foundation of Immunology and Amino Acid Chains" arXiv:1205.6031v2](https://arxiv.org/abs/1205.6031v2)
