# XRank
# Afterburner for limma fit objects, providing posterior LFC and ranking.
Good for
- Ranking genes that change the most.
- Estimating true LFC, for example for GSEA.

To install from R (using devtools), copy-paste the following:

```R
devtoolsInstalled = library("devtools", logical.return = T)
if ( !devtoolsInstalled ) {
  install.packages("devtools")
  library("devtools")
}
install_github('ChristofferFlensburg/XRank')
library(XRank)
?XRank
```

The algorithm provides two useful statistics in addition to the standard limma output: XRank and best.guess. These statistics are added as new entries in the fit data.frame in the same format as p-values or log fold changes.

XRank is the posterior expectation value of the rank by absolute LFC. This statistic is meant to replace p-values as the ranking statistic when repporting p-values. P-values rank by how strongly the null hypothesis is refuted, and is thus biased towards large and expressed genes that have more statistical power. XRank takes the magnitude of the log fold change into account, and ranks the large true log fold changes first, while avoiding false positives.

best.guess is a posterior median of the log fold change. It correlates better with true log fold change than does fit$coefficient (measures LFC) in both SEQC data and simulations. This makes the statistic particularly useful for downstream gene set analysis, such as GSEA. This statistic is unsuitable for ranking in top tables, as it is ranks null genes high much more often than a ranking by p-values or XRank.
