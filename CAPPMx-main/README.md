Common Atoms Product Partition Model on Covariates for Generating
Synthetic Controls in Single Arm Clinical Trials
================
Noirrit Kiran Chandra
2023-06-12

An implementation of the CA-PPMx model proposed by (Chandra et al. 2023;
Müller, Chandra, and Sarkar 2023). Please cite our works.

### Installation

The package has been tried and tested in `Ubuntu` and `macOS`.

#### Required UNIX Packages

- `libgsl`
- `openmp`
- `R (>= 4.2.2)`

#### Installation from Github

``` r
devtools::install_github("noirritchandra/CAPPMx", build_vignettes = TRUE)
```

Must set `build_vignettes = TRUE` to install the `vignette` files
containing illustrations.

``` r
library(CAPPMx)
```

Check the `vignette` files of the package for the reference manual:

``` r
vignette("Intro_to_CAPPMx","CAPPMx")
```

### Glioblastoma Data

A resampled version of the Glioblastoma dataset used by Chandra et al.
(2023) is included in the `CAPPMx` package which can be accessed using
the following:

``` r
data("MDACC_reproduced")
head(MDACC_reproduced)
```

    ##   Surgery Reason Histologic Grade EOR Gender ATRX MGMT CT SOC RT Dose KPS Age
    ## 1              1                1   1      1   NA    0  0   1      NA   0   0
    ## 2              1                1   0      0    1    0  0   1       1   0   1
    ## 3              1                1   0      1   NA    0  0   1       0   0   1
    ## 4              1                1   1      1    1    0  1   1       1   0   0
    ## 5              1                1   1      1    1    0  1   1       1   0   0
    ## 6              1                1   0      1   NA    0  1   1       1   0   1
    ##     endpts surv_inds
    ## 1 30.71429     FALSE
    ## 2 78.00000     FALSE
    ## 3 10.56548      TRUE
    ## 4 64.70833     FALSE
    ## 5 44.85119     FALSE
    ## 6 47.42857     FALSE

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-chandra_GBM22" class="csl-entry">

Chandra, N. K., A. Sarkar, J de Groot, Y. Yuan, and P. and Müller. 2023.
“Bayesian Nonparametric Common Atoms Regression for Generating Synthetic
Controls in Clinical Trials.” *arXiv Preprint arXiv:2201.00068*.

</div>

<div id="ref-mueller_cam23" class="csl-entry">

Müller, P., N. K. Chandra, and A. Sarkar. 2023. “Bayesian Approaches to
Include Real-World Data in Clinical Studies.” *Philosophical
Transactions of the Royal Society A: Mathematical, Physical and
Engineering Sciences* 381: 20220158.

</div>

</div>
