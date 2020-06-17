---
title: "rcsurplus1d"
output: html_document
---

## Purpose

Package `rcsurplus1d` lets you explore a novel surplus
production model for fisheries stock assessment. The model employs a production
function that differs from the canonical logistic (Schaefer) and Gompertz (Fox)
functions, but can still be related to the Pella-Tomlinson formulation. 
In this package you can use your own data to compare the four approaches, or you can
explore the well-studied Namibian hake fishery dataset.

## Getting started

If you don't have the software OpenBUGS, first download and install it (http://www.openbugs.net/w/Downloads).

Then, install three packages in the following order:

```{r}
devtools::install_github("rtlemos/rcrandom")
devtools::install_github("rtlemos/rcvirtual")
devtools::install_github("rtlemos/rcsurplus1d")
```

Load the third one:
```{r}
library(rcsurplus1d)
```

Create your instance of the package class, and fire up the graphical user interface (GUI) on your web browser:
```{r}
myinstance <- rcsurplus1d()
myinstance$gui()
```

To load your own data, navigate to the tab "Input"; on the left, under "Data input", you will see "Browse...". Click on it to load your CSV file, which *must* have the columns "year", "catch", "effort", and "cpue".

And that is it! Please read the `About` tabs, to understand the theoretical concepts and how to navigate the GUI. Additional information about the model can be found in [Rankin and Lemos (2015)](https://www.researchgate.net/publication/279962333_An_alternative_surplus_production_model?ev=prf_high).

## Further analyses

Inside the temporary folder of your `R` session (`tempdir()`), you will find several files that may be useful for troubleshooting and running further analyses. The following short list is generated when we only run a single chain for the alternative model:

* `data.txt` -> your fisheries dataset in a format that OpenBUGS understands
* `alternative_model.txt` -> OpenBUGS model specification
* `inits1.txt` -> initial values for OpenBUGS to start the model fit
* `log.txt` -> log of model fit, including possible error messages; also includes summary statistics for model parameters and deviance information.
* `CODAindex.txt` -> key to interpret `CODAchain1.txt`: indicates the first and last index of each model parameter in the chain
* `CODAchain1.txt` -> vectorized MCMC array associated with the model fit
