# contrastanalysis

Contrast Analysis with Equivalence Testing for Residual Contrasts.

Implements contrast analysis integrated with Campbell (2024) TOST equivalence testing. Uses the Intersection-Union Test (IUT; Berger, 1982) principle.

## Installation

### From a local folder

After downloading or cloning this package:

```r
# Install devtools if you don't have it
install.packages("devtools")

# Generate documentation and install
devtools::document("path/to/contrastanalysis")
devtools::install("path/to/contrastanalysis")
```

### From GitHub

```r
devtools::install_github("Y45HV1N/contrastanalysis")
```

## Quick Start

```r
library(contrastanalysis)

# Between-subjects: linear dose-response across 4 groups
dat <- data.frame(
  group = rep(c("A","B","C","D"), each = 50),
  score = c(rnorm(50,2,3), rnorm(50,4,3), rnorm(50,6,3), rnorm(50,8,3))
)

result <- contrast_analysis(
  data       = dat,
  dv         = "score",
  conditions = "group",
  hypothesis = c(A = 0, B = 1, C = 2, D = 3),
  delta      = 0.3
)
```

## Features

- **Between-subjects and within-subjects** designs (k ≥ 3 groups).
- **NHST** on the contrast of interest + **TOST equivalence tests** on residuals.
- Two equivalence bound types: standardized *d* (`delta_type = "dz"`) or percentage of model variance (`delta_type = "share_signal"`).
- **Forest plots** and **variance decomposition plots**.
- **APA-formatted output** ready for copy-paste.

## Documentation

- **Online tutorial:** [https://Y45HV1N.github.io/contrastanalysis/](https://Y45HV1N.github.io/contrastanalysis/)
- **Vignette (after installing):** `vignette("tutorial", package = "contrastanalysis")`
- **Help page:** `?contrast_analysis`

## Citation

If you use this package in your research, please cite it as:

> Seetahul, Y. (2026). *contrastanalysis: Contrast Analysis with Equivalence Testing for Residual Contrasts*. R package version 1.0.0. https://github.com/Y45HV1N/contrastanalysis

Or in BibTeX:

```bibtex
@Manual{seetahul2026contrastanalysis,
  title  = {contrastanalysis: Contrast Analysis with Equivalence Testing for Residual Contrasts},
  author = {Yashvin Seetahul},
  year   = {2026},
  note   = {R package version 1.0.0},
  url    = {https://github.com/Y45HV1N/contrastanalysis}
}
```

In R, after installing the package, you can also retrieve the citation with:

```r
citation("contrastanalysis")
```

## References

- Berger, R. L. (1982). Multiparameter hypothesis testing and acceptance sampling. *Technometrics*, *24*(4), 295-300.
- Campbell, H. (2024). Equivalence testing for linear regression. *Psychological Methods*, *29*(1), 88–98. https://doi.org/10.1037/met0000596
- Lakens, D. (2017). Equivalence tests: A practical primer. *Social Psychological and Personality Science*, *8*(4), 355-362.
