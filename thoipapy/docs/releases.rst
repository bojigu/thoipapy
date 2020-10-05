=================
THOIPApy releases
=================

1.0.0
-----
* added new prediction features including conservation calculated with rate4site
* removed blind test TMDs from the training dataset
* added a feature selection pipeline, including the removal of duplicate features, and selection of best predictive features
* added automatic tuning of machine-learning predictor
* excluded putative distant homologues from training in each iteration of leave-one-out validation
* added scripts for conducting bootstrapped t-test, comparing interface and non-interface residues
* extended the output for a single file/experiment to include the mean EC50 values for replicates with identical sample names (issue #8)

0.0.7
-----
* initial public release