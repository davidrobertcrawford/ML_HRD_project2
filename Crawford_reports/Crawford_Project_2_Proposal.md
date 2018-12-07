# Extending Davies et al. (_Nature Medicine_ 2017)

__Author__: David R. Crawford (drcbio@terpmail.umd.edu)  __Date__: November 12, 2018

------
## Planned Experiment
The goal for the experiment is to train machine learning models on the breast cancer data used by Davies et al. 2017 using accepted train/test/validate practices. Davies et al. validated their HRDetect model using the original training set plus 60 additional control cases. This is contrary to the basic principle of validating a model with held-out data. A more appropriate approach would be to separate the data into a cross-validation set (say, 80%) and a final validation set (say, 20%). Under this approach, validation consists in testing the model's performance (e.g. the ROC AUC) when applied to data _that was not used for model training or cross-validation_. In addition, this experiment will include steps to balance classes to avoid overestimating accuracy due to imbalanced classes (e.g. if 90% of cases are positive, achieving high accuracy is as easy as predicting positive for every case).

The experiment will involve training multiple regression models, as done by Davies et al. and in Project 1, and support vector machines.

## Data & Resources
Data will be that used in Davies et al. I previously worked with this data for the Project 1 reproduction of Davies et al. Figure 4a & Supplemental Figures 12 & 14. Machine learning packages for __R__ will include __glmnet__ for multiple regression and __caret__ for support vector machines. I have performed machine learning with these packages before, so familiarity with their implementation methods will not be a hurdle. ROC analysis will be performed using __pROC__ and __plotROC__ packages for __R__.

## Validation Plan
Validation of the models will consist in measuring ROC AUC values for model predictions when run on a held-out validation set. As noted above, the hold-out set will be determined at the outset so that model validation is carried out on data not used for training. While an AUC of 98% (as found in Davies et al.) is an unrealistic expectation with standard validation practices, the goal for BRCA-ness prediction will be a model AUC of >80%, with >90% considered a great success.


-----
## References

1. Davies, H., Glodzik, D., Morganella, S., Yates, L. R., Staaf, J., Zou, X.,
Ramakrishna, M., Martin, S., Boyault, S., Sieuwerts, A. M., Simpson, P. T., King, T. A.,
Raine, K., Eyfjord, J. E., Kong, G., Borg, A., Birney, E., Stunnenberg, H. G.,
van de Vijver, M. J., Borresen-Dale, A.-L., Martens, J. W. M., Span, P. N., Lakhani, S. R.,
Vincent-Salomon, A., Sotiriou, C., Tutt, A., Thompson, A. M., Van Laere, S.,
Richardson, A. L., Viari, A., Campbell, P. J., Stratton, M. R., Nik-Zainal, S.
(2017) "HRDetect is a predictor of _BRCA1_ and _BRCA2_ deficiency based on mutational
signatures." _Nature Medicine_ 23(4): 517-525. [doi: 10.1038/nm.4292](10.1038/nm.4292)
