# Predicting BRCAness with Multiple Machine Learning Approaches (Revisiting Davies et al. 2017)

**Author**: David R. Crawford (drcbio@terpmail.umd.edu)
**Date**: December 7, 2018

## Plan of investigation

The goal for the experiment is to train machine learning models on the breast cancer data used by Davies et al. 2017 using accepted train/test/validate practices. Davies et al. validated their HRDetect model using the original training set plus 60 additional cases. This is contrary to the basic principle of validating a model with held-out data. A more appropriate approach would be to separate the data into a cross-validation set (say, 90%) and a final validation set (say, 10%). Under this approach, final validation consists in testing the model's performance (e.g. the AUROC) when applied to data _that was not used for model training or cross-validation_. In addition, this experiment will include steps to balance classes to avoid overestimating accuracy due to imbalanced classes (e.g. if 90% of cases are positive, achieving high accuracy is as easy as predicting positive for every case).

The experiment will involve training multiple regression models, as done by Davies et al. and in my Project 1, and two additional model types, support vector machines and random forests. I will also run on raw and on processed data [what difference does the processing make for us?]


### Machine Learning with Multiple Feature Sets and Multiple Models

#### Experiment
I propose to run machine learning models on data from Davies et al. [1]. I depart from Davies et al. in three ways:
1. I will use a standard train/test/validate structure so that cross-validated models are tested with a holdout set.
2. I will train models on the raw data, on the processed data, and on the combined raw and processed data.
3. In addition to generalized linear models (GLM), I will train random forest (RF) and support vector machine (SVM) models.


#### Data, resources, and implementation

Data will be that used in Davies et al. [1]. The data will be obtained from the course UMIACS server folder for Davies 2017: 'https://obj.umiacs.umd.edu/mutation-signature-explorer/publications/Davies2017/processed/'. Individual file paths are listed in the __Snakefile__. Data files include:
* _add-feat_: mutation signature counts (11 signatures); chromosomal rearrangement counts (6 rearrangements); three indel signatures; homologous recombination deficiency index; two Boolean variables indicating control status and inclusion in Davies et al.'s evaluation.
* _counts-brca_: raw counts for mutations in each of 96 trinucleotide categories
* _sample-key_: key for linking different patient and sample ids
* _raw-counts_: raw counts for indel and rearrangement events

All analysis will be implemented in __R__. I will write my own code for data processing and analysis.
* Random forest analysis will be performed using the __randomForest__ package [4].
* Support vector machine analysis will be performed using the __caret__ package [5].
* Generalized linear model analysis will be performed using the __glmnet__ package [2], as done by Davies et al. [1] and in my Project 1.
* ROC analysis will be performed using the __pROC__ package [3].
* General data manipulation tools will include __tidyverse__ packages [7] and the __data.table__ package [6].

We run 20 rounds to get a good number of rounds without taking too long... User can change number of rounds by setting value for __num_runs__ at the top of "main_code.R".

#### Machine Learning Parameters

##### Generalized Linear Model (GLM)
__cv.glmnet__

Family (binomial b/c case is positive or negative); measure is AUC b/c that's the ultimate project.

| Parameter                        | Setting           | Source            |
| -------------------------------- | ----------------- | ----------------- |
| nfolds                           |    10             |                   |
| lower.limits                     |    0              |                   |
| type.measure                     | "auc"             |                   |
| family                           | "binomial"        |                   |
| standardize                      | FALSE             |     source        |

__predict.glmnet__
Response gives us prediction of classes.

| Parameter                        | Setting           | Source            |
| -------------------------------- | ----------------- | ----------------- |
| type                             | "response"        |   source          |


##### Support Vector Machine (SVM)
__trainControl__

Repeated CV for repeated cross-validation. Class probs so we get classification probabilities we can use for AUROC calcs.

| Parameter                        | Setting           | Source            |
| -------------------------------- | ----------------- | ----------------- |
|  method                          |  "repeatedcv"     |                   |
|  number                          |   5               |                   |
| repeats                          | 5                 |                   |
| classProbs                       | TRUE              |                   |
| summaryFunction                  | twoClassSummary   |                   |
| savePredictions                  | TRUE              | source            |

__train__

Metric ROC since we want AUROC as measure.

| Parameter                        | Setting           | Source            |
| -------------------------------- | ----------------- | ----------------- |
|  method                          |  "svmRadial"      |                   |
| tuneLength                       | 5                 |   |
| metric                           | "ROC"             |   |
| verbose                          | FALSE             |  source           |

__predict__

Return probs so we can calc AUROC.

| Parameter                        | Setting           | Source            |
| -------------------------------- | ----------------- | ----------------- |
|  type                            |  "prob"           |   source          |


##### Random Forests (RF)
__predict__

Prob so we can calc AUROC. Norm.votes means we get the fraction of votes for the particular case. So, it's a frequentist probability in a sense.

| Parameter                        | Setting           | Source            |
| -------------------------------- | ----------------- | ----------------- |
|  type                            |  "prob"           |                   |
|  norm.votes                      |  TRUE             |  source           |




#### Validation

Validation of the models will consist in measuring AUROC values for model predictions when run on a held-out validation set. As noted above, the hold-out set will be determined at the outset so that model validation is carried out on data not used for training. While an AUC of 98% (as found in Davies et al.) is an unrealistic expectation with standard validation practices, the goal for BRCA-ness prediction will be a model AUC of >80%, with >90% considered a great success.





## Results, conclusions, and caveats

### Mutation signature discovery

We ran the experiment described above, and identified four mutation signatures...

We also attempted to...

![Reproduced Figure 1a from Kim, et al. [1]](fig1a-plot.jpg)
![Reproduced Figure 1b from Kim, et al. [1]](fig1b-plot.jpg)

In conclusion, according to the validation criteria defined above, we find that the results of the mutation signature discovery experiments from Kim, et al. [1] are...
The only two caveats to this conclusion are that:
1. text...
2. text...

### Future Directions

Mention validation on outside datasets (as Davies et al. did?)

## References

1. Davies, H., Glodzik, D., Morganella, S., Yates, L. R., Staaf, J., Zou, X.,
Ramakrishna, M., Martin, S., Boyault, S., Sieuwerts, A. M., Simpson, P. T., King, T. A.,
Raine, K., Eyfjord, J. E., Kong, G., Borg, A., Birney, E., Stunnenberg, H. G.,
van de Vijver, M. J., Borresen-Dale, A.-L., Martens, J. W. M., Span, P. N., Lakhani, S. R., Vincent-Salomon, A., Sotiriou, C., Tutt, A., Thompson, A. M., Van Laere, S.,
Richardson, A. L., Viari, A., Campbell, P. J., Stratton, M. R., Nik-Zainal, S.
(2017) "HRDetect is a predictor of _BRCA1_ and _BRCA2_ deficiency based on mutational
signatures." _Nature Medicine_ 23(4): 517-525. [doi: 10.1038/nm.4292](10.1038/nm.4292)
2. J. Friedman, T. Hastie, R. Tibshirani, N. Simon, B. Barasimhan, J. Qian (2018) "Package 'glmnet'. Lasso and Elastic-Net Regularized Generalized Linear Models." Version 2.0-16, 2018-03-12. CRAN:  https://CRAN.R-project.org/package=glmnet.
3. X. Robin, N. Turk, A. Hainard, N. Tiberti, F. Lisacek, J.-C. Sanchez, M. Muller, S. Siegbert (2018) "Package 'pROC'. Display and Analyze ROC Curves." Version 1.13.0, 2018-09-23. CRAN: https://CRAN.R-project.org/package=pROC.
4. L. Breiman, A. Cutler, A. Liaw, M. Wiener (2018) "Package 'randomForest'. Breiman and Cutler's Random Forests for Classification and Regression." Version 4.6-14, 2018-03-22. CRAN: https://CRAN.R-project.org/package=randomForest.
5. M. Kuhn (2018) "Package 'caret'. Classification and Regression Training." Version 6.0-81, 2018-11-20. CRAN:  https://CRAN.R-project.org/package=caret.
6. M. Dowle, A. Srinivasan, J. Gorecki, M. Chirico, P. Stetsenko, T. Short, S. Lianoglou, E. Antonyan, M. Bonsch, H. Parsonage (2018) "Package 'data.table'. Extension of 'data.frame'." Version 1.11.8, 2018-09-30. CRAN:  https://CRAN.R-project.org/package=data.table.
7. H. Wickham, RStudio (2017) "Package 'tidyverse'. Easily Install and Load the 'Tidyverse'." Version 1.2.1, 2017-11-14. CRAN:  https://CRAN.R-project.org/package=tidyverse.
