
<!-- README.md is generated from README.Rmd. Please edit that file -->

# elu2

<!-- badges: start -->
<!-- badges: end -->

The goal of elu2 is to …

## Installation

## Installation

You can install the development version of `elu2` from
[GitHub](https://github.com/Elusan/elu2) with:

``` r
# install.packages("devtools")
devtools::install_github("Elusan/elu2")
#> Using github PAT from envvar GITHUB_PAT. Use `gitcreds::gitcreds_set()` and unset GITHUB_PAT in .Renviron (or elsewhere) if you want to use the more secure git credential store instead.
#> Downloading GitHub repo Elusan/elu2@HEAD
#> ── R CMD build ─────────────────────────────────────────────────────────────────
#> * checking for file ‘/private/var/folders/h2/gtrfj58d751fknk0r1n1sct00000gn/T/Rtmp8uqgex/remotes5731648d45d0/Elusan-elu2-acb5a6a/DESCRIPTION’ ... OK
#> * preparing ‘elu2’:
#> * checking DESCRIPTION meta-information ... OK
#> * checking for LF line-endings in source and make files and shell scripts
#> * checking for empty or unneeded directories
#> Omitted ‘LazyData’ from DESCRIPTION
#> * building ‘elu2_0.1.0.tar.gz’
#> Installing package into '/private/var/folders/h2/gtrfj58d751fknk0r1n1sct00000gn/T/RtmpTqTuj4/temp_libpath54fb61cba57f'
#> (as 'lib' is unspecified)
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(elu2)
## basic example code
hello()
#> [1] "Hello from elu2!"
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
