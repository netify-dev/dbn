# DBN validation suite

These files run a small simulation study to spot-check parameter recovery
and posterior-predictive calibration. They are not part of the CRAN
test suite because each fit takes 1-10 seconds and the checks are
statistical (not deterministic).

To run:

```r
library(dbn)
library(testthat)
Sys.setenv(DBN_VALIDATION = "true")
testthat::test_dir(system.file("validation", package = "dbn"))
```
