# Contributing

First of all, thank you for your interest in the spNetwork package! Contributions and suggestions are welcome.

## How to report a Bug

if you have found a bug, please, create an associated issue and illustrate the bug with a minimal 
reproductible example. Do not forget to specify your R version, the session informations and your OS.

## How to Propose a Change to spNetwork

### Prerequisites

Before making any pull request, you must file an issue and make sure that @JeremyGelb agrees with 
what you plan to develop. This is just a way to ensure that what you propose is doable considering
the actual architecture of the package and to respect its spirit.

### Pull request process

Guidelines for contributors:

* Create a named branch (not master) which merges cleanly with master.
* Limit the scope of the pull request to a clearly defined improvement or fix.
* For naming functions, use name_separated_by_underscore.
* spNetwork relies on the spatial class from the package sf and NOT sp anymore.
* Avoid to include new dependencies.
* Document (more is better than less) your code.
* Use [roxygen2](https://cran.r-project.org/package=roxygen2) for documentation.
* Use [testthat](https://cran.r-project.org/package=testthat). Contributions with test cases included are easier to accept.
* After you merge a pull request, delete the associated branch if you can. You will see an option for this on the GitHub pull request page.

Every pull request will be reviewed by @JeremyGelb. He will either merge it or provide feedback on issues to be addressed.

### Code of Conduct

Please note that the spNetwork is released with a [Contributor Code of Conduct](CONDUCT.md). By contributing to this project, you agree to abide by its terms.
