# TCRMatch
### Authors: Austin Crinklaw, Will Chronister

## Requirements:
- Linux OS
- Python 3+
  - Python packages: Pandas, Cython

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


### Commands
-  Using a .txt file as an input:
```shell
python -m TCRMatch match -i /path/to/input.txt -o /path/to/output.txt
```
### Output  
-  Output file has 3 columns in TSV format. 
-  First column is the user provided input sequence.  
-  Second column contains a CDR3beta sequence from the IEDB scoring higher than the threshold
-  The third column contains the score based on the kernel similarity metric (ranging from threshold to .99)

| input_sequence | match_sequence | score              |
|----------------|----------------|--------------------|
| ASSQDRDTQY     | ASGDAGGGYEQY   | 0.7474520339781123 |
| ASSQDRDTQY     | ASGDASGAETLY   | 0.8084889120179228 |
| ASSQDRDTQY     | ASGDASGGNTLY   | 0.7825980249961477 |
| ASSQDRDTQY     | ASGDFWGDTLY    | 0.7614741870953899 |
| ASSQDRDTQY     | ASRYRDDSYNEQF  | 0.8165516456455025 |

## How it works:
- TCRMatch implements a kernel based similarity metric as defined in [arXiv:1205.6031v2](https://arxiv.org/abs/1205.6031v2)
- Each input sequence is compared to all TCR CDR3beta sequences characterized in the IEDB
- Sequences that have a kernel score greater than the threshold are returned as matches

## References:
["Towards a Mathematical Foundation of Immunology and Amino Acid Chains" arXiv:1205.6031v2](https://arxiv.org/abs/1205.6031v2)
