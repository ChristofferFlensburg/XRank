# XRank
# Afterburner for limma fit objects, providing posterior LFC and ranking.
Good for
- Ranking genes that change the most.
- Estimating true LFC, for example for GSEA.

To install from R (using devtools), copy-paste the following:

devtoolsInstalled = library("devtools", logical.return = T)
if ( !devtoolsInstalled ) {
  install.packages("devtools")
  library("devtools")
}
install_github('ChristofferFlensburg/XRank')
library(XRank)
?XRank
