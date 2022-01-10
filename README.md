
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PMPLNFA

Parsimonious Finite Mixtures of Multivariate Poisson-Log Normal Factor
Analyzers for Clustering Count Data

<!-- badges: start -->

![GitHub commit activity
(branch)](https://img.shields.io/github/commit-activity/y/anjalisilva/MPLNClust/master)
<!-- badges: end -->

## Description

Mixture model-based clustering methods can be over-parameterized in
high-dimensional spaces, especially as the number of clusters increases.
Subspace clustering allows to clus- ter data in low-dimensional
subspaces, while keeping all the dimensions and by introducing
restrictions to mixture parameters (Bouveyron and Brunet, 2014).
Restrictions are introduced to the model parameters with the aim of
obtaining parsimonious models, which are sufficiently flexible for
clustering purposes. Since the largest contribution of free parameters
is through the covariance matrices, it is a natural focus for the
introduction of parsimony.

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

## References

-   [Aitchison, J. and C. H. Ho (1989). The multivariate Poisson-log
    normal distribution.
    *Biometrika.*](https://www.jstor.org/stable/2336624?seq=1)

-   [Silva, A., S. J. Rothstein, P. D. McNicholas, and S. Subedi (2019).
    A multivariate Poisson-log normal mixture model for clustering
    transcriptome sequencing data. *BMC
    Bioinformatics.*](https://pubmed.ncbi.nlm.nih.gov/31311497/)

-   [Subedi, S., R.P. Browne (2020). A family of parsimonious mixtures
    of multivariate Poisson-lognormal distributions for clustering
    multivariate count data. *Stat.*
    9:e310.](https://doi.org/10.1002/sta4.310)

-   [McNicholas, P. D., and T. B. Murphy (2008). Parsimonious Gaussian
    mixture models. *Statistics and Computing.* 18,
    285–296.](https://link.springer.com/article/10.1007/s11222-008-9056-0)

## Authors

-   Anjali Silva (<anjali@alumni.uoguelph.ca>).
-   Sanjeena Dang (<sanjeena.dang@carleton.ca>).

## Maintainer

-   Anjali Silva (<anjali@alumni.uoguelph.ca>).

## Contributions

`PMPLNFA` repository welcomes issues, enhancement requests, and other
contributions. To submit an issue, use the [GitHub
issues](https://github.com/anjalisilva/mixMPLNFA/issues).

## Acknowledgments

-   Dr. Marcelo Ponce, SciNet HPC Consortium, University of Toronto, ON,
    Canada for all the computational support.
-   This work was funded by Natural Sciences and Engineering Research
    Council of Canada, Queen Elizabeth II Graduate Scholarship, and
    Arthur Richmond Memorial Scholarship.
