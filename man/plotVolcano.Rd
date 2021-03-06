% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/XRank.R
\name{plotVolcano}
\alias{plotVolcano}
\title{Plots the LFC and p-value of a fit object.}
\usage{
plotVolcano(fit, coef = 1, print = 20, main = NA, shrink = T,
  line = "fdr", xlim = "default", ylim = "default",
  xlab = "corrected LFC", ylab = "-log10(p-value)", specialGenes = c(),
  col = "red", ...)
}
\arguments{
\item{fit}{A fit object that has gone through XRank.}

\item{coef}{The columns of the fit object to be plotted. Default 1.}

\item{print}{Non-negative integer. The number of top ranked genes to have
their names printed on the plot.}

\item{shrink}{Logical. Whether to print the higher ranked genes in a smaller font.}

\item{line}{character. Where to plot a horisontal support line. 'nGenes' plots a
line at p = 1/nrow(fit), above which there will be one false positive
on average. 'fdr' draws the line at the p corresponding to fdr = 5\%,
above which 5\% of the genes are false positives on average.
'fdr' reverts to 'nGenes' if no gene has fdr <= 5\%. Default 'fdr'.}

\item{specialGenes}{vector of rownames of fit. These genes will be highlighted no matter
where on the plot they are.}

\item{col}{The colour of the best.guess dots. Default 'red'.}

\item{...}{Remaining arguments are passed on to plot(...).}
}
\description{
Plots the LFC and p-value of a fit object.
}
\details{
This function plots a standard volcano plot: The -log10(p-value) as function
            of the LFC. It also superimpose the posterior estimates of the LFC by XRank:
            best.guess in red. The best.guess will in general be strongly shrunk towards
            LFC = 0 for large p-values, and leave small p-values essentially unchanged.
            See XRank for more information about these statistics.
}
\examples{
\donttest{
#Set up a (highly artificial) count matrix.
counts = matrix(rpois(1000*6, 100), ncol=6,
         dimnames=list(paste0('gene', 1:1000), c('a', 'b', 'c', 'D', 'E', 'F')))

#set up experimental design
group = c(rep('lower', 3), rep('upper', 3))
design = model.matrix(~0+group)
colnames(design) = gsub('^group', '', colnames(design))
contrasts = limma::makeContrasts('upper-lower', levels=colnames(design))

#run voom and limma
fit = limma::voom(counts, design, plot=T)
fit = limma::lmFit(fit, design=design)
fit = limma::contrasts.fit(fit, contrasts)
fit = limma::eBayes(fit)
fit = XRank(fit, plot=F)
plotVolcano(fit)
}
}

