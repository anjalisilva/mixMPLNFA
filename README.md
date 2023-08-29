
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mixMPLNFA

Mixtures of Multivariate Poisson-Log Normal Factor Analyzers for
Clustering Count Data

<!-- badges: start -->

![GitHub commit activity
(branch)](https://img.shields.io/github/commit-activity/y/anjalisilva/MPLNClust/master)
<!-- badges: end -->

## Description

`mixMPLNFA` is an R package for performing clustering using parsimonious
mixtures of multivariate Poisson-log normal factor analyzers family
(PMPLNFA) via variational Gaussian approximations. It was developed for
count data, with clustering of RNA sequencing data as a motivation.
However, the clustering method may be applied to other types of count
data. The package provides functions for parameter estimation via a
variational Gaussian approximation with Expectation-Maximization (EM)
algorithm. Information criteria (AIC, BIC, AIC3 and ICL) utilized for
model selection.

## Installation

To install the latest version of the package:

``` r
require("devtools")
devtools::install_github("anjalisilva/mixMPLNFA", build_vignettes = TRUE)
library("mixMPLNFA")
```

To run the Shiny app (under construction):

``` r
mixMPLNFA::runMixMPLNFA()
```

## Overview

To list all functions available in the package:

``` r
ls("package:mixMPLNFA")
```

`MPLNClust` contains 7 functions.

1.  ***PMPLNFAClust*** for carrying out clustering of count data using
    mixtures of MPLN via variational expectation-maximization
2.  ***mplnVisualizeLine*** for visualizing clustering results as line
    plots
3.  ***runMixMPLNFA*** is the shiny implementation of *PMPLNFAClust*
    (under construction)

## Details

Mixture model-based clustering methods can be over-parameterized in
high-dimensional spaces, especially as the number of clusters increases.
Subspace clustering allows to cluster data in low-dimensional subspaces,
while keeping all the dimensions and by introducing restrictions to
mixture parameters (Bouveyron and Brunet, 2014). Restrictions are
introduced to the model parameters with the aim of obtaining
parsimonious models, which are sufficiently flexible for clustering
purposes. Since the largest contribution of free parameters is through
the covariance matrices, it is a natural focus for the introduction of
parsimony.

In 2008, a family of eight parsimonious Gaussian mixture models (PGMMs;
[McNicholas and Murphy,
2008](https://link.springer.com/article/10.1007/s11222-008-9056-0)) were
introduced with parsimonious covariance structures. In 2019, a
model-based clustering methodology using mixtures of multivariate
Poisson-log normal distribution (MPLN; [Aitchison and Ho,
1989](mixMPLNFA)) was developed to analyze multivariate count
measurements by [Silva et al.,
2019](https://pubmed.ncbi.nlm.nih.gov/31311497/). In this work, a family
of mixtures of MPLN factor analyzers that is analogous to the PGMM
family is developed, by considering the constraints. This family is
referred to as the parsimonious mixtures of MPLN factor analyzers family
(PMPLNFA).

### Variational-EM Framework for Parameter Estimation

[Subedi and Browne, 2020](https://doi.org/10.1002/sta4.310) proposed a
variational Gaussian approximation that alleviates challenges of MCMC-EM
algorithm. Here the posterior distribution is approximated by minimizing
the Kullback-Leibler (KL) divergence between the true and the
approximating densities. A variational-EM based framework is used for
parameter estimation.

## Model Selection and Other Details

Four model selection criteria are offered, which include the Akaike
information criterion (AIC; Akaike, 1973), the Bayesian information
criterion (BIC; Schwarz, 1978), a variation of the AIC used by Bozdogan
(1994) called AIC3, and the integrated completed likelihood (ICL;
Biernacki et al., 2000).

Starting values (argument: *initMethod*) play an important role to the
successful operation of this algorithm. There maybe issues with
singularity, in which case altering initialization method may help.

## Shiny App

The Shiny app employing ***PMPLNFAClust*** could be run and results
could be visualized:

``` r
mixMPLNFA::runMixMPLNFA()
```

## Tutorials

For tutorials and plot interpretation, refer to the vignette (under
construction):

``` r
browseVignettes("mixMPLNFA")
```

## Citation for Package

``` r
citation("mixMPLNFA")
```

Finite Mixtures of Multivariate Poisson-Log Normal Factor Analyzers for
Clustering Count Data. (2023) *Unpublished*.

``` r
A BibTeX entry for LaTeX users is

  @Article{,
    title = {Finite Mixtures of Multivariate Poisson-Log Normal
Factor Analyzers for Clustering Count Data},
    year = {2023}
  }
```

## References

- [Aitchison, J. and C. H. Ho (1989). The multivariate Poisson-log
  normal distribution.
  *Biometrika.*](https://www.jstor.org/stable/2336624?seq=1)

- [Silva, A., S. J. Rothstein, P. D. McNicholas, and S. Subedi (2019). A
  multivariate Poisson-log normal mixture model for clustering
  transcriptome sequencing data. *BMC
  Bioinformatics.*](https://pubmed.ncbi.nlm.nih.gov/31311497/)

- [Subedi, S., R.P. Browne (2020). A family of parsimonious mixtures of
  multivariate Poisson-lognormal distributions for clustering
  multivariate count data. *Stat.*
  9:e310.](https://doi.org/10.1002/sta4.310)

- [McNicholas, P. D., and T. B. Murphy (2008). Parsimonious Gaussian
  mixture models. *Statistics and Computing.* 18,
  285–296.](https://link.springer.com/article/10.1007/s11222-008-9056-0)

## Authors

- Anjali Silva (<anjali@alumni.uoguelph.ca>).
- Andrea Payne (<andreapayne@cmail.carleton.ca>).
- Sanjeena Dang (<sanjeena.dang@carleton.ca>).

## Maintainer

- Anjali Silva (<anjali@alumni.uoguelph.ca>).

## Contributions

`mixMPLNFA` repository welcomes issues, enhancement requests, and other
contributions. To submit an issue, use the [GitHub
issues](https://github.com/anjalisilva/mixMPLNFA/issues).

## Acknowledgments

- Dr. Marcelo Ponce, SciNet HPC Consortium, University of Toronto, ON,
  Canada for all the computational support.
- Early work was funded by Natural Sciences and Engineering Research
  Council of Canada (Subedi) and Queen Elizabeth II Graduate Scholarship
  (Silva).
- Later work was supported by the Postdoctoral Fellowship award from the
  Canadian Institutes of Health Research (Silva) and the Canada Natural
  Sciences and Engineering Research Council grant 400920-2013 (Subedi).
