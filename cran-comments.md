## Submission

This is an update of blockmodels (1.1.5, December 2021) to version 1.2.1.

* **Change of maintainer**: Julien Chiquet <julien.chiquet@inrae.fr> takes over the
  maintenance from Jean-Benoist Leger <jbleger@hds.utc.fr>, who agrees and remains an author.
  Jean-Benoist Leger will confirm by email.
* New features: covariates on the nodes, missing values in the network.
* Bug fixes, some of which change the results (see NEWS.md).

## Test environments

* local: Ubuntu 24.04, R 4.6.2
* GitHub Actions: ubuntu-latest (R devel, release, oldrel-1), macos-latest (release),
  windows-latest (release)
* win-builder: R devel (2026-11-05 r90840) and R release (4.6.2)
* mac-builder was unavailable (HTTP 502) at the time of submission: macOS is covered by
  GitHub Actions (macos-latest, R release)

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE: new maintainer.

## Reverse dependencies

We checked the 8 reverse dependencies on CRAN (expSBM, graphclust, GREMLINS, mimiSBM,
missSBM, robber, sand, sbm), comparing R CMD check results with the CRAN version and
with this version of blockmodels.

* 7 packages: no new problem.
* sbm: the tests checking that the predictions of sbm match, to 1e-12, those of
  blockmodels for the Bernoulli models with covariates fail (12 failures). blockmodels 1.2.1 fixes these
  predictions (they were exact only for hard memberships), so the test now detects the
  same bug in sbm. We maintain sbm too: sbm 1.0.0, which requires blockmodels (>= 1.2.1)
  and fixes its predictions accordingly, will be submitted as soon as this version is
  accepted.
