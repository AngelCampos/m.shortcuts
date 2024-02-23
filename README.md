
# m.shortcuts

<!-- badges: start -->

[![](https://img.shields.io/badge/devel%20version-0.1.1-blue.svg)](https://github.com/AngelCampos/m.shortcuts)
<!-- badges: end -->

## Description

R package with sourced, original, and wrappers to recurrent functions I
like and find recurrently useful.

## Installation

You can install the development version from
[GitHub](https://github.com/AngelCampos/m.shortcuts) typing in the
following commands in the R console:

``` r
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
devtools::install_github("AngelCampos/m.shortcuts")
```

## Functions list

- check_figs(): Make a table of figures in ‘figs’ dir, generated in .R
  scripts and plot in .Rmd scripts
- ggBoxplot(): Shortcut for a ggplot2 boxplot
- intervalHeatmap(): Heatmap over 3 dimensions. It divides a vector over
  a set number of intervals over two dimensions. Useful to understand
  how values distribute across 2 dimensions.
- lsos(): Improved list of objects. List the ten biggest objects in the
  global environment by default.
- multiplot(): An alternative to plot several plots, but it does not
  work to save to file.
- rmTmp(): Remove all objects in the global environment that start with
  ‘tmp’.
