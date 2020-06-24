#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

//Peptide struct to make life easier
typedef struct peptide
{
    char *seq;
    int len;
    int *i;
    float aff;

} PEPTIDE;

//global variable declarations
float k1[20][20];
int p_kmin = 1;
int p_kmax = 30;
float p_beta = 0.11387;
//Hardcoded to remove file I/O and annoying file parsing
float blm_qij[20][20] = {
    {0.0215, 0.0023, 0.0019, 0.0022, 0.0016, 0.0019, 0.003, 0.0058, 0.0011, 0.0032, 0.0044, 0.0033, 0.0013, 0.0016, 0.0022, 0.0063, 0.0037, 0.0004, 0.0013, 0.0051},
    {0.0023, 0.0178, 0.002, 0.0016, 0.0004, 0.0025, 0.0027, 0.0017, 0.0012, 0.0012, 0.0024, 0.0062, 0.0008, 0.0009, 0.001, 0.0023, 0.0018, 0.0003, 0.0009, 0.0016},
    {0.0019, 0.002, 0.0141, 0.0037, 0.0004, 0.0015, 0.0022, 0.0029, 0.0014, 0.001, 0.0014, 0.0024, 0.0005, 0.0008, 0.0009, 0.0031, 0.0022, 0.0002, 0.0007, 0.0012},
    {0.0022, 0.0016, 0.0037, 0.0213, 0.0004, 0.0016, 0.0049, 0.0025, 0.001, 0.0012, 0.0015, 0.0024, 0.0005, 0.0008, 0.0012, 0.0028, 0.0019, 0.0002, 0.0006, 0.0013},
    {0.0016, 0.0004, 0.0004, 0.0004, 0.0119, 0.0003, 0.0004, 0.0008, 0.0002, 0.0011, 0.0016, 0.0005, 0.0004, 0.0005, 0.0004, 0.001, 0.0009, 0.0001, 0.0003, 0.0014},
    {0.0019, 0.0025, 0.0015, 0.0016, 0.0003, 0.0073, 0.0035, 0.0014, 0.001, 0.0009, 0.0016, 0.0031, 0.0007, 0.0005, 0.0008, 0.0019, 0.0014, 0.0002, 0.0007, 0.0012},
    {0.003, 0.0027, 0.0022, 0.0049, 0.0004, 0.0035, 0.0161, 0.0019, 0.0014, 0.0012, 0.002, 0.0041, 0.0007, 0.0009, 0.0014, 0.003, 0.002, 0.0003, 0.0009, 0.0017},
    {0.0058, 0.0017, 0.0029, 0.0025, 0.0008, 0.0014, 0.0019, 0.0378, 0.001, 0.0014, 0.0021, 0.0025, 0.0007, 0.0012, 0.0014, 0.0038, 0.0022, 0.0004, 0.0008, 0.0018},
    {0.0011, 0.0012, 0.0014, 0.001, 0.0002, 0.001, 0.0014, 0.001, 0.0093, 0.0006, 0.001, 0.0012, 0.0004, 0.0008, 0.0005, 0.0011, 0.0007, 0.0002, 0.0015, 0.0006},
    {0.0032, 0.0012, 0.001, 0.0012, 0.0011, 0.0009, 0.0012, 0.0014, 0.0006, 0.0184, 0.0114, 0.0016, 0.0025, 0.003, 0.001, 0.0017, 0.0027, 0.0004, 0.0014, 0.012},
    {0.0044, 0.0024, 0.0014, 0.0015, 0.0016, 0.0016, 0.002, 0.0021, 0.001, 0.0114, 0.0371, 0.0025, 0.0049, 0.0054, 0.0014, 0.0024, 0.0033, 0.0007, 0.0022, 0.0095},
    {0.0033, 0.0062, 0.0024, 0.0024, 0.0005, 0.0031, 0.0041, 0.0025, 0.0012, 0.0016, 0.0025, 0.0161, 0.0009, 0.0009, 0.0016, 0.0031, 0.0023, 0.0003, 0.001, 0.0019},
    {0.0013, 0.0008, 0.0005, 0.0005, 0.0004, 0.0007, 0.0007, 0.0007, 0.0004, 0.0025, 0.0049, 0.0009, 0.004, 0.0012, 0.0004, 0.0009, 0.001, 0.0002, 0.0006, 0.0023},
    {0.0016, 0.0009, 0.0008, 0.0008, 0.0005, 0.0005, 0.0009, 0.0012, 0.0008, 0.003, 0.0054, 0.0009, 0.0012, 0.0183, 0.0005, 0.0012, 0.0012, 0.0008, 0.0042, 0.0026},
    {0.0022, 0.001, 0.0009, 0.0012, 0.0004, 0.0008, 0.0014, 0.0014, 0.0005, 0.001, 0.0014, 0.0016, 0.0004, 0.0005, 0.0191, 0.0017, 0.0014, 0.0001, 0.0005, 0.0012},
    {0.0063, 0.0023, 0.0031, 0.0028, 0.001, 0.0019, 0.003, 0.0038, 0.0011, 0.0017, 0.0024, 0.0031, 0.0009, 0.0012, 0.0017, 0.0126, 0.0047, 0.0003, 0.001, 0.0024},
    {0.0037, 0.0018, 0.0022, 0.0019, 0.0009, 0.0014, 0.002, 0.0022, 0.0007, 0.0027, 0.0033, 0.0023, 0.001, 0.0012, 0.0014, 0.0047, 0.0125, 0.0003, 0.0009, 0.0036},
    {0.0004, 0.0003, 0.0002, 0.0002, 0.0001, 0.0002, 0.0003, 0.0004, 0.0002, 0.0004, 0.0007, 0.0003, 0.0002, 0.0008, 0.0001, 0.0003, 0.0003, 0.0065, 0.0009, 0.0004},
    {0.0013, 0.0009, 0.0007, 0.0006, 0.0003, 0.0007, 0.0009, 0.0008, 0.0015, 0.0014, 0.0022, 0.001, 0.0006, 0.0042, 0.0005, 0.001, 0.0009, 0.0009, 0.0102, 0.0015},
    {0.0051, 0.0016, 0.0012, 0.0013, 0.0014, 0.0012, 0.0017, 0.0018, 0.0006, 0.012, 0.0095, 0.0019, 0.0023, 0.0026, 0.0012, 0.0024, 0.0036, 0.0004, 0.0015, 0.0196}};

float **fmatrix_k1()
{
    int k, j;
    float marg[20];
    float sum;

    //initialize margin array
    for (int i = 0; i < 20; i++)
    {
        marg[i] = 0.0;
    }
    //initialize k1
    for (int i = 0; i < 20; i++)
    {
        for (int j = 0; j < 20; j++)
        {
            k1[i][j] = 0.0;
        }
    }
    //normalize matrix by marginal frequencies
    for (j = 0; j < 20; j++)
    {
        sum = 0;
        for (k = 0; k < 20; k++)
            sum += blm_qij[j][k];
        marg[j] = sum;
    }
    //calculate K1
    for (j = 0; j < 20; j++)
    {
        for (k = 0; k < 20; k++)
        {
            k1[j][k] = blm_qij[j][k] / (marg[j] * marg[k]);
            k1[j][k] = pow(k1[j][k], p_beta);
        }
    }

    return (k1);
}

float k2_prod(int *i1, int *i2, int start1, int start2, int k)
{
    float k2;
    int x, j1, j2;
    k2 = 1;
    for (x = 0; x < k; x++)
    {

        j1 = i1[x + start1];
        j2 = i2[x + start2];
        k2 *= k1[j1][j2];
    }
    return (k2);
}

float k3_sum(int *s1, int *s2, int l1, int l2)
{
    float k3, prod;
    int k;
    int start1, start2;

    k3 = 0.0;
    for (k = p_kmin; k <= p_kmax; k++)
    {
        for (start1 = 0; start1 <= l1 - k; start1++)
        {
            for (start2 = 0; start2 <= l2 - k; start2++)
            {
                prod = k2_prod(s1, s2, start1, start2, k);
                k3 += prod;
            }
        }
    }

    return (k3);
}

PyObject* calculate_k3_list(PyObject *self, PyObject *args)
{
    char* seq1;
    char* seq2;

    if (!PyArg_ParseTuple(args, "ss", &seq1, &seq2))
        return NULL;
    
    char *alphabet = "ARNDCQEGHILKMFPSTWYV";
    //declare peptide structs for input strings
    PEPTIDE pep1;
    PEPTIDE pep2;

    //Calculate K1 matrix
    fmatrix_k1();
    
    int length_s1;
    length_s1 = strlen(seq1);
    pep1.seq = seq1;
    pep1.len = length_s1;
    int *pep1_i = (int *)malloc(sizeof(int)*length_s1);

    //Simple loop to get position of current char in alphabet
    for (int j = 0; j < length_s1; j++)
    {
        char *e;
        e = strchr(alphabet, pep1.seq[j]);
        pep1_i[j] = (int)(e - alphabet);
    }
    pep1.i = pep1_i;
    //Calculate normalization matrix
    pep1.aff = k3_sum(pep1.i, pep1.i, pep1.len, pep1.len);


    int length_s2;
    length_s2 = strlen(seq2);
    pep2.seq = seq2;
    pep2.len = length_s2;
    int *pep2_i = (int *)malloc(sizeof(int)*length_s2);

    //Simple loop to get position of current char in alphabet
    for (int j = 0; j < length_s2; j++)
    {
        char *e;
        e = strchr(alphabet, pep2.seq[j]);
        pep2_i[j] = (int)(e - alphabet);
    }
    pep2.i = pep2_i;
    //Calculate normalization matrix
    pep2.aff = k3_sum(pep2.i, pep2.i, pep2.len, pep2.len);

    
    PyObject *res;
    float sco;
    sco = k3_sum(pep1.i, pep2.i, pep1.len, pep2.len) / sqrt(pep1.aff * pep2.aff);
    res = PyTuple_Pack(3, PyUnicode_FromString(pep1.seq), PyUnicode_FromString(pep2.seq), PyFloat_FromDouble(sco));

    return res;
}

//Python C API requires these to be defined
static PyMethodDef MatchMethods[] = {
    {"tcrmatch", calculate_k3_list, METH_VARARGS, "TCRMatch Python interface"},
    {NULL, NULL, 0, NULL}};

static struct PyModuleDef matchmodule = {
    PyModuleDef_HEAD_INIT,
    "tcrmatch",
    "TCRMatch Python interface",
    -1,
    MatchMethods};

PyMODINIT_FUNC PyInit_tcrmatch_c(void)
{
    return PyModule_Create(&matchmodule);
}