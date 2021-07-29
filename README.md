## Stata DCA

Functions to perform Decision Curve Analysis (DCA) using Stata.

The functions within the dca package include detailed help files (access the help files with help dca, or help stdca after installation).
For detailed vignettes and examples of DCA in action using Stata, R and SAS, visit [decisioncurveanalysis.org](decisioncurveanalysis.org).
Use the following command to install the collection of DCA functions in Stata. 

```stata
net install dca, from("https://raw.github.com/ddsjoberg/dca.stata/master/") replace
```

## Unit Testing

1. Confirm the working directory set in `unit-testing/tests.do` is correct.
1. Run the file `unit-testing/tests.do`. This will...
    - Source the current versions of `dca.do` and `stdca.do`.
    - Run each of the unit testing files located in the `unit-testing/tests` folder.
1. Confirm the results in the Stata console window.

## Release History

#### (development version)

* Added prevalence()` option to the `dca` functions. Users working with case-control data can now specify the population prevalence.

* Addded a unit testing infrastructure for the package.

#### v1.0.0 (2017-12-06)

* Initial release.