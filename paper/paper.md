---
title: 'NEEP: null empirically estimated p-values for high-throughput molecular survival analysis'
tags:
  - Python
  - bioinformatics
  - survival analysis
  - RNA-seq
authors:
  - name: Sean West
    orcid: 0000-0002-1394-1999
    affiliation: 1
  - name: Hesham Ali
    affiliation: 1
  - name: Dario Ghersi
    orcid: 0000-0002-0630-0843
    affiliation: 1
affiliations:
 - name: School of Interdisciplinary Informatics, University of Nebraska at Omaha
   index: 1
date: 17 December 2019
bibliography: paper.bib
---

# Summary

When conducting survival analysis for molecular expression, a researcher has
two main options: Cox-Proportional Hazards (CPH) or Kaplan-Meier (KM). CPH uses regression to
calculate the survival-significance of each expression vector. It has two 
assumptions that are frequently violated in cancer research. First, CPH 
assumes that the censoring mechanism is not associated with patient survival.
When using patients from cancer research, those who survive to the end of 
the study are likely to have a higher rate of survival than those who did
not [@west2019]. Second, the core of CPH is the proportional hazards assumption. 
For molecular data, this assumes that the effect that a molecule has on survival
is constant over time. This cannot be the case as the pathology of cancer
changes across stages and time. 

KM is non-parametric and uses a log-rank test to check if patients with
low expression of the molecule of interest have altered survival rates 
from those with high expression. KM has the same two assumptions mentioned
for CPH. In addition, high-throughput KM survival analysis
using a single threshold has been shown to be sensitive to patient group
re-sampling [@sehgal2015]. However, we can get around the first assumption by conducting
the log-rank test along a range of splits, where the threshold that 
splits the expression-ordered list of patients into low- and high-expression
is tested at multiple points.  The second assumption is weaker in KM than
in CPH, but it is important that survival-significant molecules do not 
have KM curves that cross. 
 
Conducting KM across a range of thresholds for each molecules produces a
range of p-values, of which we could sample the minimum. However, 
taking the lowest p-value from a range will produce a non-uniform, 
right-skewed distribution of p-values. Since p-values should be uniform
under the null distribution, the skewed distribution cannot be used for
valid statistical analysis.
An equation was developed that could predict the correct p-values was
developed [@lausen1992]; however, it is not precise enough for p-value correction
procedures that are sensitive to very small p-value changes. Thus, we 
developed NEEP, which overcomes these issues by re-sampling permutations
of the patients to construct a null distribution in parallel. Using this
null distribution, NEEP transforms the p-values so they are statistically
valid. Finally, NEEP conducts False Discovery Rate (FDR) p-value correction and calculates 
effect sizes, the hazard ratio and the 1, 2, and 5 year mortality
ratios.


# Statement of Need 

**Research purpose**: `NEEP` offers non-parametric, high-throughput, and statistically valid survival analysis of molecular expression vectors.

**Problem solved**: `NEEP` permutes different orders of patients and calculates their p-values across a range in order to produce a `Null Empirically Estimated P-value' for each molecular expression vector, thus overcoming the assumptions of CPH survival analyiss and the issues of correct p-value estimation. 

**Target audience**: The target audience is anyone conducting molecular, high-throughput survival analysis that does not have confounding clinical variables and whose expression vectors may violate CPH assumptions. Full documentation is available in the [project repository](https://github.com/thecodingdoc/neep).


# Acknowledgements 

This project was partly funded by a
University of Nebraska Collaboration Initiative/
System Science Seed Grant to Sushil Kumar, Hesham Ali, and Dario Ghersi and
by the NIH AA026428 R21 grant to Sushil Kumar. The funder
website is https://nebraska.edu/collaboration-initiative. The funders had no role in project design,
data collection and analysis, decision to publish, or
preparation of the manuscript.

# References