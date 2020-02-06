# NEEP
## Purpose
Null Empirically Estimated P-values (NEEP) is a non-parametric, high-throughput survival analysis meant for molecular expression data. The algorithm empirically finds the optimal cutoff (within a range) to separate expression values into *low* and *high*. The survival curves for this optimal cutoff are constructed. We calculate the MANTEL log-rank test for these optimally separated curves. The subsequent p-value distribution across all molecular objects will be skewed. To estimate the correct p-values, we sample 1M (or user specified) 'genes', which are a uniform permutation of the patient survival values from the provided clinical data. Then, p-values are estimated empirically and adjusted using the Benjamini-Hochberg approach.

Additional processing on NEEP output into protein interaction visualizations can be done using SINBAD: Survival-significant Isoform Networks By Altered Domain-inclusion (https://github.com/scwest/SINBAD).

Adding Makefile instructions shortly.

## Installation

```console
cd ~/neep
g++ -O3 -o neep neep.C util.C -std=c++11 -fopenmp
```
The ```-fopenmp``` option allows the computation of the null distribution to run in parallel on a multicore machine using the OpenMP library. It is optional but highly recommended as it speed things up.

## Usage

```console
cd ~/neep
./a.out -c <clinical filename> -e <expression filename> -o <output filename> -n <size of null distribution to sample> -t <minimum threshold to check> -u
```

The minimum threshold to check refers to the lowest acceptable split between low and high expression when conducting KM tests. The maximum threshold will automatically be 1-t. The ```-u``` (optional) option is for to force a uniform null distribution. We were seeing some weird activity from c++ method random_shuffle() and we found that manually assigning a uniform distribution generator we obtained nulls which matched those from the Python and Julia programming languages (previous versions of NEEP).

**Clinical file**

Each row is a patient. There is no header. *Days to event* is a combination of the normal *days to death* and *days to last followup* used in clinical survival analysis. A 0 in *event* means that there are no death has been recorded (censored values). 

```
<patient id>,<days to event>,<event>
...
```

**Expression file**
This file is the normal expression matrix, where the rows are the molecular objects of interest and the columns are the patients. The patients must be the same set from the clinical file and they must be in identical order. 

```
,<patient id 1>,<patient id 2>,...,<final patient id>
<molecular object 1>,<expression value>,<expression value>,...,<expression value>
...
```

**Output file**
The output will be in tab-separated format with a header.
Column 1: molecular object id
Column 2: best log-rank statistic
Column 3: null empirically estimated p-value (neep)
Column 4: adjusted neep
Column 5: direction (i.e. which survival curve is higher, *low* or *high*)

