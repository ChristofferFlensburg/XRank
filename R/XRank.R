

#This function calculates the prior from the data supplied in the fit object.
#The bins are denser around LFC = 0, so default number of bins should be enough.
getPriors = function(fit, coefs = 0, nbins = 50, plot=T, mode='empirical') {
  #if no coefficients specified, get priors for all.
  if ( length(coefs) == 1 && coefs == 0 ) coefs = colnames(fit$coefficients)
  N = length(coefs)

  #find maximum LFC, to know the range to take data over.
  range = lapply(coefs, function(coef) 1.1*max(abs(fit$coefficients[,coef])))
  names(range) = coefs

  #Set the breaks for the bins in the histogram.
  #Smaller bins close to 0, where accuracy is more important,
  #and density is higher
  breaks = lapply(coefs, function(coef)
    (((-nbins):(nbins))/nbins)^2*sign((-nbins):(nbins))*range[[coef]] )
  names(breaks) = coefs

  #Now get the histograms, which will be the priors
  priors = lapply(coefs, function(coef) {
    if ( mode == 'posterior' & 'best.guess' %in% names(fit) ) hist(fit$best.guess[,coef], breaks = breaks[[coef]], plot=plot)
    else hist(fit$coefficients[,coef], breaks = breaks[[coef]], plot=plot)
  })
  names(priors) = coefs

  #Return the priors
  return(priors)
}

#A function generating some data evenly distributed in LFC up to
#the supplied maximum LFC, and then transforms it into a flat prior.
getFlatPrior = function(max = 15, nbins = 100, plot=T) {
  #Set the breaks for the bins in the histogram.
  #Smaller bins close to 0, where accuracy is more important,
  #and density is higher
  breaks = (((-nbins):(nbins))/nbins)^2*
    sign((-nbins):(nbins))*max

  #fake a flat dataset
  flat = runif(100000, -max, max)

  #Now get the histogram, which will be the prior
  prior = hist(flat, breaks = breaks, plot=plot)

  #Return the priors
  return(prior)
}

#' Plots the LFC and p-value of a fit object.
#'
#' @param fit   A fit object that has gone through XRank.
#' @param coef  The columns of the fit object to be plotted. Default 1.
#' @param print Non-negative integer. The number of top ranked genes to have
#'              their names printed on the plot.
#' @param shrink Logical. Whether to print the higher ranked genes in a smaller font.
#' @param line character. Where to plot a horisontal support line. 'nGenes' plots a
#'             line at p = 1/nrow(fit), above which there will be one false positive
#'             on average. 'fdr' draws the line at the p corresponding to fdr = 5\%,
#'             above which 5\% of the genes are false positives on average.
#'             'fdr' reverts to 'nGenes' if no gene has fdr <= 5\%. Default 'fdr'.
#' @param specialGenes  vector of rownames of fit. These genes will be highlighted no matter
#'             where on the plot they are.
#' @param col  The colour of the best.guess dots. Default 'red'.
#' @param ...  Remaining arguments are passed on to plot(...).
#'
#' @details This function plots a standard volcano plot: The -log10(p-value) as function
#'             of the LFC. It also superimpose the posterior estimates of the LFC by XRank:
#'             best.guess in red. The best.guess will in general be strongly shrunk towards
#'             LFC = 0 for large p-values, and leave small p-values essentially unchanged.
#'             See XRank for more information about these statistics.
#'
#' @examples
#' \donttest{
#' #Set up a (highly artificial) count matrix.
#' counts = matrix(rpois(1000*6, 100), ncol=6,
#'          dimnames=list(paste0('gene', 1:1000), c('a', 'b', 'c', 'D', 'E', 'F')))
#'
#' #set up experimental design
#' group = c(rep('lower', 3), rep('upper', 3))
#' design = model.matrix(~0+group)
#' colnames(design) = gsub('^group', '', colnames(design))
#' contrasts = limma::makeContrasts('upper-lower', levels=colnames(design))
#'
#' #run voom and limma
#' fit = limma::voom(counts, design, plot=T)
#' fit = limma::lmFit(fit, design=design)
#' fit = limma::contrasts.fit(fit, contrasts)
#' fit = limma::eBayes(fit)
#' fit = XRank(fit, plot=F)
#' plotVolcano(fit)
#' }
#'
#' @import limma
#' @export
plotVolcano = function(fit, coef=1, print = 20, main = NA, shrink=T, line='fdr',
  xlim = 'default', ylim = 'default',
  xlab='corrected LFC', ylab='-log10(p-value)', specialGenes=c(), col='red', ...) {
  #get the name of the column if only index provided
  if ( is.numeric(coef) ) coef = colnames(fit)[coef]
  print=min(print, nrow(fit))
  #get the FCs and p values
  bgs = fit$best.guess[,coef]
  cos = fit$coefficients[,coef]
  ps = fit$p.value[, coef]

  #decide for a range on the x-axis
  xmax = max(1, 1.5*max(abs(bgs)))
  xlim = c(-xmax, xmax)

  if ( is.na(main) ) main = coef
  
  #plot the points and segments
  plot(1, type='n', xlim = xlim, ylim=c(0, max(-log10(ps))), xlab=xlab, ylab=ylab, main = main, ...)
  if ( line == 'nGenes' ) segments( 2*xlim[1], log10(nrow(fit)), 2*xlim[2], log10(nrow(fit)), cex=1, col='orange', lty=2)
  if ( line == 'fdr' ) {
    fdr = p.adjust(ps, method='fdr')
    cut = -log10(max(ps[fdr <= 0.05], 1/nrow(fit)))
    segments(2*xlim[1], cut, 2*xlim[2], cut, cex=1, col='orange', lty=2)
  }
  segments(bgs, -log10(ps), cos, -log10(ps), cex=0.4, col=rgb(0,0,0,0.15))
  points(cos, -log10(ps), pch=16, cex=0.5)
  points(bgs, -log10(ps), pch=16, cex=0.7, col=col)

  legend('bottomright', c('measured', 'best guess'), pch=16, pt.cex=c(0.5, 0.8), col=c('black', 'red'))

  #add labels for top DE tags.
  if ( print > 0 ) {
    names = rownames(fit)
    tops = which(cos > 0)[order(fit$XRank[cos > 0,coef])[1:print]]
    bots = which(cos < 0)[order(fit$XRank[cos < 0,coef])[1:print]]

    ymax = -log10(min(ps))
    if ( shrink ) cex = 0.6 + 0.4*(print:1)/print
    else cex = 0.8
    shift = (0.5+1.5*cex)*ymax/100
    text(bgs[tops], -log10(ps[tops])+shift, names[tops], col='red', cex=cex)
    text(bgs[bots], -log10(ps[bots])+shift, names[bots], col='red', cex=cex)
    special = which(rownames(fit) %in% specialGenes)
    if ( length(special) > 0 ) {
      points(bgs[special], -log10(ps[special]), col=rgb(0, 0.8, 1), cex=1, pch=16)
      text(bgs[special], -log10(ps[special])+ymax/50, names[special], col=rgb(0, 0.8, 1), cex=1)
    }
  }
}

findGLFCandPV = function(dist, breaks, p) {
  dist = as.numeric(dist)
  breaks = as.numeric(breaks)
  #integrated probability distribution
  sums = cumsum(as.numeric(dist))

  #the break that bring the integral over p
  i = max(which(sums < p))

  #interpolate to find the x between the bins best fitting the limit.
  x0 = breaks[i+1] + (p - sums[i])*(breaks[i+2] - breaks[i+1])/dist[i+1]

  #integrated probability to be below 0
  p0 = interpolate(c(-Inf, breaks), sums, 0)
  #find smallest side (above or below 0).
  #multiply by 2 for double side test.
  p0 = 2*min(p0, 1-p0)

  return(c(x0, p0))
}

getPosterior = function(r, d=0.1, prior = 'flat', mode='normal',
  df=4, t=1, verbose=F, zoom = F, binsize = 0) {
  #just a shorthand notation.
  v = verbose

  #get the flat prior if needed.
  if ( is.character(prior) &&  prior == 'flat' ) {
    max = 2*(abs(log2(r))+d+0.5)
    prior = getFlatPrior(max = max, plot=F)
  }
  
  #set the grid for the prior.
  x = prior$mids
  #The density, not depending on binwidth.
  priorD = prior$density

  #calcualte binsize on the grid
  if ( !is.list(binsize) ) {
    nbins = length(x)
    binsize = prior$breaks[2:(nbins+1)] - prior$breaks[1:nbins]
  }

  #set the prior probability distribution of the LFC on the grid.
  #Both the density (useful for multiplying probabilities) and
  #"counts" (useful for integrating) are kept track of.
  if ( mode == 'normal' ) {
    flatD = dnorm(x, log2(r), d)
  }
  else if ( mode == 't' ) {
    flatD = dt((x-log2(r))*t/abs(log2(r)), df, 0)
  }

  #find the normalised (over all bins) posterior prob dist
  #Multiply the densities to get the posterior probability density.
  postD = flatD*priorD
  #Transformt to counts for easier integration.
  postC = postD*binsize
  #Normalise to give probabilities.
  postCNorm = postC/sum(postC)

  return(postCNorm)
}

posteriors = function(fit, mode='t', priors='empirical', FDR = F,
                     coefs = 0, plot=T, quiet = F, cpus=1) {
  if ( coefs[1] == 0 ) coefs = 1:ncol(fit)
  if ( is.numeric(coefs) ) coefs = colnames(fit)[coefs]

  #get the priors
  if ( is.character(priors) && priors == 'empirical' ) {
    #call the emprical prior function.
    if ( !quiet ) cat('Preparing empirical priors... ')
    priors = getPriors(fit, coefs = coefs, plot=plot, mode='empirical')
    if ( !quiet ) cat('done.\n')
  }
  else if ( is.character(priors) && priors == 'posterior' ) {
    #call the emprical prior function.
    if ( !quiet ) cat('Preparing empirical priors... ')
    priors = getPriors(fit, coefs = coefs, plot=plot, mode='posterior')
    if ( !quiet ) cat('done.\n')
  }
  else if ( is.character(priors) && priors == 'flat' ) {
    #call the flat prior function.
    max = sapply(1:length(coefs), function(i) 1.1*max(abs(fit$coefficients[,i])))
    priors = lapply( 1:length(coefs), function(i) getFlatPrior(max = max[i], plot=plot))
  }
  else if ( is.list(priors) ) {
    #if user supplied prior, do a sanity check.
    if ( length(priors) != length(coefs) ) {
      cat(length(priors),' priors, and ', length(coefs),
          ' coefficients. Need matching numbers.\n',sep='')
      return(-1)
    }
    else if ( !quiet ) cat('Using provided priors.\n')
  }
  else {
    cat('Could not make sense of prior input. Please use \'empirical\', \'flat\' or provide a list of prepared priors.\n')
    return(-1)
  }
  
  if ( !quiet ) cat('Calculating posteriors:\n')

  #calculate binsizes once for all
  getBinsize = function(prior) {
    nbins = length(prior$mids)
    prior$breaks[2:(nbins+1)] - prior$breaks[1:nbins]
  }
  binsize = lapply(priors, getBinsize) 

  #loop through the genes and coefficients, calculating the safe
  #ratio for each by calling the GR function.
  #The uncertainty in LFC is taken from the t-statistics or the
  #s2.post and stdev.unscaled of the fit object, depending on the input option.
  getPosteriors = function(rep, gene) {
    getPosterior(2^fit$coefficients[gene, rep], prior = priors[[rep]],
                 mode=mode, df = fit$df.total[gene], t=abs(fit$t[gene,rep]),
                 d = sqrt(fit$s2.post[gene])*fit$stdev.unscaled[gene, rep],
                 binsize = binsize[[rep]])
  }
  
  getMatrix = function(rep) {
    if ( !quiet ) cat(rep, '...', sep='')
    do.call(rbind,
            mclapply(1:length(rownames(fit)), function(gene) getPosteriors(rep, gene), mc.cores=cpus))
  }

  data = lapply(coefs, getMatrix)
  
  if ( !quiet ) cat('done.\n')

  #name the columns, for prettier output.
  names(data) = coefs
  for ( coef in coefs ) {
    rownames(data[[coef]]) = rownames(fit)
    colnames(data[[coef]]) = priors[[coef]]$mids
  }

  return(list(data = data, prior=priors, fit = fit))
}

interpolates = function(x, y, x0s) {
  y0s = sapply(x0s, function(x0) interpolate(x, y, x0))
  return(y0s)
}

interpolate = function(x, y, x0) {
  if ( x0 < min(x) ) return(y[1])
  if ( x0 > max(x) ) return(y[length(y)])
  if (x0 %in% x) return(y[which(x == x0)])
  
  lowI = max(which(x < x0))
  highI = min(which(x > x0))
  lowX = x[lowI]
  highX = x[highI]
  lowY = y[lowI]
  highY = y[highI]
  ret = lowY + (x0-lowX)/(highX-lowX)*(highY-lowY)
  return(ret)
}


postRank = function(posts, quiet=F, cpus=1) {
  if ( !quiet ) cat('Calculating expected ranks...')

  #set up return matrix
  retRank = matrix(0, ncol=length(posts$data), nrow=nrow(posts$data[[1]]))
  retMeanRank = matrix(0, ncol=length(posts$data), nrow=nrow(posts$data[[1]]))
  colnames(retRank) = colnames(retMeanRank) = names(posts$data)
  rownames(retRank) = rownames(retMeanRank) = rownames(posts$data[[1]])

  for ( coef in names(posts$data) ) {
    if ( !quiet ) cat(coef, '..', sep='')
    mx = posts$data[[coef]]
    mx = mx[,1:50] + mx[,100:51]
    xs = posts$prior[[coef]]$mids
    xs = xs[1:50]
    F = colSums(mx)

    ranks = cumsum(F) - F/2
    binsize = F
    cumPost = do.call(rbind, mclapply(1:nrow(mx), function(gene) cumsum(mx[gene,]), mc.cores=cpus))
    
    currentRank = 1
    ranking = rep(1, nrow(mx))
    for (i in 1:length(ranks)) {
      if ( floor(ranks[i]) < currentRank ) next
      n = floor(ranks[i]) - currentRank + 1
      tops = order(-cumPost[,i])[1:n]
      ranking[currentRank:(currentRank+n-1)] = tops
      cumPost[tops,] = cumPost[tops,]*0
      currentRank = currentRank + n
    }
    
    meanRank = sapply(1:nrow(mx), function(gene) sum(ranks*mx[gene,]))
    if ( T | coef == names(posts$data)[1] ) {
      retRank[,coef] = ranking
      retMeanRank[,coef] = meanRank
      #retRank = ranking
      #retMeanRank = meanRank
    }
    else {
      retRank = cbind(retRank, ranking)
      retMeanRank = cbind(retMeanRank, meanRank)
    }
  }

  if ( !quiet ) cat('done.\n')
  return(list(ranking = retRank, XRank = retMeanRank))
}


bestGuess = function(fit, coefs = 0, quiet=F, cpus=1) {
  if ( !quiet ) cat('Calculating best guess...')
  fit$best.guess = do.call(cbind,
    lapply(coefs, function(coef) {
      if ( !quiet ) cat(coef, '..', sep='')
      unlist(mclapply(1:nrow(fit), function(gene)
                      findGLFCandPV(fit$posterior[[coef]][gene,],
                                    fit$prior[[coef]]$breaks, 0.5)[1], mc.cores=cpus))
    }
           )
    )
  colnames(fit$best.guess) = names(fit$posterior)
  rownames(fit$best.guess) = rownames(fit)
  if ( !quiet ) cat('done.\n')
  return(fit)
}


#' Calculates posterior LFC and ranks on expectation value of rank of true LFC.
#'
#' @param fit A fit object generated from limma and voom.
#' @param coefs The columns of the fit object to be analysed.
#'   Defaults to 0 which means all columns.
#' @param keepPosterior Logical. Whether to store the numerical posterior distributions
#'   for each gene. By default a 100 column matrix, so can take some space.
#'   If stored, it allows for plots of the posterior distributions by gene.
#' @param verbose Logical. Whether a few rows of output are generated to track progress.
#' @param plot Logical. Whether a plot is generated at the end of analysis.
#' @param cpus Integer. The number of parallel cpus to use.
#' @return The \code{fit} object, with new columns:
#'   'best.guess': The posterior LFC estimate.
#'   'XRank': The expectation value of the posterior on rank on true absolute LFC.
#'            Ranking by this statistic generally outperforms other stats.
#'
#' @details Calculates a posterior PDF on the LFC based on an empirical prior. The median
#'          of this PDF is stored as the 'best.guess' and has been shown to correlate better
#'          with the true LFC than 'coefficients'. The posterior PDF on LFC is transformed into a PDF on
#'          rank by true absolute LFC. The expectation value of this PDF is stored as
#'          XRank, and is intended to be used for ranking DE genes in way to put the largest
#'          true absolute LFCs on the top of the list, as opposed to 'p.values' or 'lods'
#'          that rank the most strongly refused null hypothesis.
#'
#' @examples
#' \donttest{
#' #Set up a (highly artificial) count matrix.
#' counts = matrix(rpois(1000*6, 100), ncol=6,
#'          dimnames=list(paste0('gene', 1:1000), c('a', 'b', 'c', 'D', 'E', 'F')))
#'
#' #set up experimental design
#' group = c(rep('lower', 3), rep('upper', 3))
#' design = model.matrix(~0+group)
#' colnames(design) = gsub('^group', '', colnames(design))
#' contrasts = limma::makeContrasts('upper-lower', levels=colnames(design))
#'
#' #run voom and limma
#' fit = limma::voom(counts, design, plot=T)
#' fit = limma::lmFit(fit, design=design)
#' fit = limma::contrasts.fit(fit, contrasts)
#' fit = limma::eBayes(fit)
#' fit = XRank(fit, plot=F)
#' plotVolcano(fit)
#' }
#'
#' @import limma
#' @export
XRank = function(fit, coefs = 0, keepPosterior = T, verbose=F, plot=F, cpus=1) {
  if ( coefs == 0 ) coefs = colnames(fit)
  if ( is.numeric(coefs) ) coefs = colnames(fit)[coefs]
  require(parallel)
  posts = posteriors(fit, coefs = coefs, quiet= !verbose, plot=plot, cpus=cpus, prior='empirical')
  ranks = postRank(posts, quiet = !verbose, cpus=cpus)

  if ( keepPosterior ) fit$posterior = posts$data
  fit$prior = posts$prior
  fit$XRank = as.matrix(ranks$XRank)
  fit = bestGuess(fit, coefs = coefs, quiet = !verbose, cpus=cpus)
  
  if ( plot ) plotVolcano(ret)

  return(fit)
}
