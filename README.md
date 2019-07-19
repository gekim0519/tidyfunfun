
## `tidyfunfun`

### FDA Functions for tf-Class

This package extends functions for FDA in **`tidyfun`**. Specfically,
four functions were made for conducting **fpca** and
**function-on-scalar regression** on **`tfd`** data types.

Functional principal components analysis on tf data types is implemented
in **`fpca.tfd`**.

Fitting function for Function-on-Scalar Regression for cross-sectional
data is implemented in **`ols_cs_tfd`**, **`gibbs_cs_fpca_tfd`**,
**`gibbs_cs_wish_tfd`**.

The functions are modified from the **`refund`** package
([link](https://github.com/refunders/refundable)) to take in a dataframe
containing `tfd` columns directly as an input.

## Installation

Installation to tidyfunfun

``` r
# install.packages("devtools")
devtools::install_github('gekim0519/tidyfunfun')
```

Installation to tidyfun

``` r
devtools::install_github("fabian-s/tidyfun")
```
