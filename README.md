# euroformix_cli

## Overview

`euroformix_cli` is an commandline interface of the R package `euroformix` version 3.0.0 and above.


## Installation

``` r
# get the 3.0.0 version of euroformix from GitHub:
install.packages("devtools")
devtools::install_github("saschawilluweit/euroformix", ref = "v3.0.0")
```

If you find a bug, please file a minimal reproducible example in the
[issues](https://github.com/saschawilluweit/euroformix_clie/issues).


## Usage

Call `euroformix_cli.R` or `Rscript euroformix_cli.R` with the following commandline options
```
  -h --help                       Show this screen.
  -V --version                    Show version.
  -v --verbose                    Verbose output.
  --db_file <file>                Allele frequencies table.
  --samples_path <dir>            Path of sample files (one CSV file per sample with or without replicates).
  --sample_file <file>            Single sample file with or without replicates.
  --queried_file <file>           Genotypes of reference profiles in question (for use under Hp).
  --known_file <file>             Genotypes of known reference profiles (for use under Hp and Hd).
  --output_path <dir>             Path to place output files [default: results].
  --kit_database_file <kit_file>  Kit database file with all known kits (if not provided an internal one is used).
  --kit <name>                    Kit used for all given samples/replicates (in auto mode a separate column need to indicate the appropriate kit) [default: auto].
  --unknowns <n>                  Nubmber of unknows under prosecution hypotheses [default: 1].
  --fst <float>                   F_ST subpopulation/coancestry correction [default: 0.03].
  --dropin                        Consider Drop-In.
  --dropin_probability <float>    Probability of Drop-In [default: 0.05].
  --dropin_lambda <float>         Drop-In hyperparameter [default: auto].
  --stutter                       Consider stutter.
  --stutter_function <code>       Density function of stutter proportion parameter x [default: 1].
  --degradation                   Consider degradation.
  --analysis_threshold <n>        Analysis threshold in RFU (auto for autodetection = minimal observed peak height) [default: auto].
  --seed <n>                      Seed to reproduce randomness.
  --mle_steps <n>                 Maximum number of random evaluations nlm-optimizing routing [default: 4].
  --mle_iterations <n>            Maximum number of iterations for the MLE optimization [default: 100].
  --mcmc                          Calculate conservative LR estimate based on MCMC.
  --mcmc_iterations <n>           Number of MCMC iterations [default: 1000].
  --mcmc_delta <n>                A numerical parameter to scale with the covariance function Sigma [default: 10].
  --validation                    Do validation of assumption under Hp/Hd.
  --validation_alpha <float>      Rejection significance of peak-height-model under Hp/Hd [default: 0.01].
  --deconvolution                 Perform deconvolution under Hp/Hd.
  --parallel <n>                  Perfom calculatations in parallel at n cores (auto for autodetection) [default: auto].
  -Q --quiet                      Don't write any file or output anything except MLE LR (disables verbosity, MCMC, validation and deconvolution).
```


### Test

`euroformix_tests.R` is currently under development. It should test against `euroformix` and assure correct results.


### Contributing

We welcome contributions of all types\!

If you have never made a pull request to an R package before, `euroformix` is
an excellent place to start. Find an
[issue](https://github.com/saschawilluweit/euroformix_cli/issues/) with the **Beginner
Friendly** tag and comment that you’d like to take it on and we’ll help
you get started.

We encourage typo corrections, bug reports, bug fixes and feature
requests. Feedback on the clarity of the documentation is especially
valuable.
