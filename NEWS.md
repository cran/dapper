# dapper 1.1.1

* Added a randomized-response vignette demonstrating a complete analysis.
* Clarified the model callback requirements, including the conditional i.i.d.
  assumption for latent records, and corrected the documentation of draw counts
  and acceptance statistics.
* Corrected the discrete Gaussian probability mass formula and DOI links, and
  added plotting to the shared sampler example.
* `ddnorm()` and `rdnorm()` now require an integer `mu`. `rdlaplace()` now
  requires a positive integer `scale`; `ddlaplace()` continues to accept any
  positive `scale`. Added tests for these input restrictions.

# dapper 1.1.0

* `dapper_sample()` now additionally returns the mean acceptance rate for each latent record.
* updated naming scheme for `new_privacy`.
* fixed minor typos and updates links in documentation.

# dapper 1.0.1

* `st_f()` assertion bug fix. Incorrect argument order in one of the checks.
* `st_f()` description typo corrected.
* `ddnorm()` now returns the log unnormalized density when the log option is set to `TRUE`.

# dapper 1.0.0

* Initial CRAN submission.
