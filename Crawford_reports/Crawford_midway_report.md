# Reproducing Davies et al. (_Nature Medicine_ 2017)

__Author__: David R. Crawford (drcbio@terpmail.umd.edu)  __Date__: October 22, 2018

## Plan of Investigation
------

### Determine Relevant Genomic Features
The goal is to reproduce findings presented in Supplementary tables 12 & 14.
This includes the stability analysis of genomic features (i.e. how often does is a feature
  assigned a non-zero coefficient across 100 cross-validation runs?) and the determination
  of model weights to use for predictions.

#### Experiment
text

#### Data, resources, and implementation
For this part I used the __R__ package referenced in Davies et al., __glmnet__.

This part is performed on data for the set of 77 _BRCA1_/_BRCA2_ carriers and 234 quiescent tumors.

Following Davies et al. 2017, the values are log normalized and then normalized to mean = 0 and SD = 1.
For these two normalization steps I used:
```
log_col = function(the_vector){
    the_vector1 = log(the_vector+1)
    return(the_vector1)
}
```
and
```
norm_col = function(the_vector){
    the_mean = mean(the_vector)
    the_sd = sd(the_vector)
    the_vector1 = the_vector-the_mean
    the_vector2 = the_vector1/the_sd
    return(the_vector2)
}
```

The outer folds divide the data into 10% for holdout validation and 90% for cross-validation.
For each outer fold, cross-validation is performed with 10 inner folds. I used the
following call for cross-validation:

```
cv.glmnet(x = train_mat, y = train_response, nfolds=10,
          lower.limits = 0, type.measure = "auc",
          family = "binomial", intercept = 0)
```
To force coefficients to be non-negative, __lower.limits__ is set to 0;
to optimize for AUC, __type.measure__ is set to "auc"; [explain!] __family__ is set to
"binomial"; and [explain!] __intercept__ is set to 0.


#### Validation
text


### Build AUC Plot for HRDetect
text

#### Experiment
text

#### Data, resources, and implementation
text

#### Validation
text

## Results, Conclusions, and Caveats
------
### Determine Relevant Genomic Features
For 100 runs the validation AUC had mean = 0.997 and sd = 0.00897. The results are as follows:

| Genomic Feature Name | Rounds with non-zero <br> coefficient (out of 100) | Coefficient <br> Mean | Coefficient <br> SD |
| --- | :---: | :---: | :---: |
| Deletions with microhomology | 100 | 0.698 | 0.564 |
| Substitutions signature 3 | 100 | 1.205 | 0.359 |
| Rearrangement signature RS5 | 91 | 0.149 | 0.102 |
| Rearrangement signature RS3 | 55 | 0.384 | 0.563 |
| HRD Index | 8 | 0.00764 | 0.0400 |
| Deletions at repeats | 1 | 0.00164 | 0.0164

Davies et al. (2017) use a different set of features (Supplementary table 14), adding "Substitutions signature 8"
and removing "Deletions at repeats", with respect to my table. The table below includes the coefficients they used in
their final training set (Supplementary table 14), and the non-zero counts for coefficients from their
stability analysis of genomic features trained on half of data from 77 carriers & 234 quiescent (Supplementary table 12).

| Genomic Feature Name | Rounds with non-zero <br> coefficient (out of 100) <br> [Davies] | Rounds with non-zero <br> coefficient (out of 100) <br> [Crawford] | Coefficient <br> Mean | Coefficient <br> SD |
| --- | :---: | :---: | :---: | :---: |
| Deletions with microhomology | 100 | 100 | 0.218 | 0.090 |
| Substitutions signature 3 | 99 | 100 | 2.096 | 3.555 |
| Rearrangement signature RS5 | 67 | 91 | 1.935 | 1.483 |
| Rearrangement signature RS3 | 61 | 55 | 1.260 | 1.657 |
| HRD Index | 92 | 8 | 2.195 | 0.750 |
| Substitutions signature 8 | 81 | 0 | 4.390 | 3.179 |
| Deletions (other) | 15 | 0 | - | - |
| Substitutions signature 5 | 13 | 0 | - | - |
| Substitutions signature 13 | 6 | 0 | - | - |
| Rearrangement signature RS1 | 2 | 0 | - | - |
| Deletions at repeats | 0 | 1 | - | - |

#### Midway Status
Although my analysis found a similar set of relevant genomic features, I'd like to look more closely and see
why I get a different set and why even for some shared items the non-zero round counts are so different.

I would also like to look into the differences in coefficients. It might be that I'm using different
arguments in my __cv.glmnet__ call.

One obvious difference is their inclusion of __Substitutions signature 8__, based on 81 rounds with non-zero coefficients,
 with a final coefficient of 4.390 (the largest for their features)!

### Build AUC Plot for HRDetect
#### Midway Status
This part is relatively straightforward once I've settled on a set of genomic features & their coefficients.
I just need to run __predict__ on the full data set using the consensus model (__glmnet.fit__) from the
training.



## References
------
1. Davies, H., Glodzik, D., Morganella, S., Yates, L. R., Staaf, J., Zou, X.,
Ramakrishna, M., Martin, S., Boyault, S., Sieuwerts, A. M., Simpson, P. T., King, T. A.,
Raine, K., Eyfjord, J. E., Kong, G., Borg, A., Birney, E., Stunnenberg, H. G.,
van de Vijver, M. J., Borresen-Dale, A.-L., Martens, J. W. M., Span, P. N., Lakhani, S. R.,
Vincent-Salomon, A., Sotiriou, C., Tutt, A., Thompson, A. M., Van Laere, S.,
Richardson, A. L., Viari, A., Campbell, P. J., Stratton, M. R., Nik-Zainal, S.
(2017) "HRDetect is a predictor of _BRCA1_ and _BRCA2_ deficiency based on mutational
signatures." _Nature Medicine_ 23(4): 517-525. [doi: 10.1038/nm.4292](10.1038/nm.4292)
