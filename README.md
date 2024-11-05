
<!-- badges: start -->

[![R-CMD-check](https://github.com/EpiForeSITE/branching_process/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/EpiForeSITE/branching_process/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/EpiForeSITE/branching_process/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/EpiForeSITE/branching_process/actions/workflows/pkgdown.yaml)
<!-- badges: end -->

## Infectious disease outbreak quantification using a branching process model with negative binomial offspring distribution

This package provides functions that quantify infectious disease outbreaks
using a branching process, a stochastic process in which each individual in
generation n produces a random number of individuals in generation n+1,
continuing for some number of generations or until there are no individuals
remaining.

The random number of next-generation individuals produced by each individual
is drawn from the offspring distribution, a discrete probability distribution
with non-negative range. To model infectious disease outbreaks, it is common
to use a negative binomial offspring distribution, parameterized by the mean
`R` and dispersion parameter `k`. This parameterization is equivalent to
using mu = R and size = k in R's "NegBinomial", e.g. dnbinom(x, mu=R, size=k) would give the density,  i.e. the probability of exactly x transmissions from one individual.

The functions in the package can be used to quantify the risk posed by individual importers of a novel transmissible
pathogen to a generic population, including intervention effects. They can also be used for transmission parameter estimation, e.g. via maximum likelihood, for observed outbreak clusters, such as the basic reproduction number, a dispersion parameter quantifying variance in transmission, and a post-control reproduction number.

Many functions in the package were used in the following publications.

- Toth D, Gundlapalli A, Khader K, Pettey W, Rubin M, Adler F, Samore M
  (2015). Estimates of outbreak risk from new introductions of Ebola
  with immediate and delayed transmission control. Emerg Infect Dis,
  21(8), 1402-1408. <https://doi.org/10.3201/eid2108.150170>.

- Toth D, Tanner W, Khader K, Gundlapalli A (2016). Estimates of the
  risk of large or long-lasting outbreaks of Middle East respiratory
  syndrome after importations outside the Arabian Peninsula. Epidemics,
  16, 27-32. <https://doi.org/10.1016/j.epidem.2016.04.002>

## Installing the package

To install the package, you can use the following code:

``` r
devtools::install_github("EpiForeSITE/branching_process")
```

## To cite the package in publications

To cite the package in publications, please use:

    ## To cite package 'branchingprocess' in publications use:
    ## 
    ##   Toth D (????). _Branching Process Outbreak Simulator_. R package
    ##   version 0.0-9, <https://epiforesite.github.io/branching_process/>.
    ## 
    ## A BibTeX entry for LaTeX users is
    ## 
    ##   @Manual{,
    ##     heather = {To cite branchingprocess in publications use:},
    ##     title = {Branching Process Outbreak Simulator},
    ##     author = {Damon Toth},
    ##     note = {R package version 0.0-9},
    ##     url = {https://epiforesite.github.io/branching_process/},
    ##   }
