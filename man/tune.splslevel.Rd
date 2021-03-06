\name{tune.splslevel}
\encoding{UTF-8}
\alias{tune.splslevel}


\title{Tuning functions for multilevel sPLS method}

\description{
For a multilevel spls analysis, the tuning criterion is based on the maximisation of the correlation between the components from both data sets
}

\usage{
tune.splslevel(X, Y, multilevel, ncomp = NULL,
mode = "regression",
test.keepX = rep(ncol(X), ncomp),
test.keepY = rep(ncol(Y), ncomp),
already.tested.X = NULL,
already.tested.Y = NULL)

}

\arguments{
\item{X}{numeric matrix of predictors. \code{NA}s are allowed.}
\item{Y}{\code{if(method = 'spls')} numeric vector or matrix of continuous responses (for multi-response models) \code{NA}s are allowed.}
\item{multilevel}{Design matrix for multilevel analysis (for repeated measurements) that indicates the repeated measures on each individual, i.e. the individuals ID. See Details.}
\item{ncomp}{the number of components to include in the model.}
\item{mode}{character string. What type of algorithm to use, (partially) matching
one of \code{"regression"}, \code{"canonical"}, \code{"invariant"} or \code{"classic"}.}
\item{test.keepX}{numeric vector for the different number of variables to test from the \eqn{X} data set}
\item{test.keepY}{numeric vector for the different number of variables to test from the \eqn{Y} data set}
\item{already.tested.X}{Optional, if \code{ncomp > 1} A numeric vector indicating the number of variables to select from the \eqn{X} data set on the firsts components.}
\item{already.tested.Y}{Optional, if \code{ncomp > 1} A numeric vector indicating the number of variables to select from the \eqn{Y} data set on the firsts components.}


}

\details{
For a multilevel spls analysis, the tuning criterion is based on the maximisation of the correlation between the components from both data sets
}

\value{

\item{cor.value}{correlation between latent variables}

}

\references{
mixOmics article:

Rohart F, Gautier B, Singh A, Lê Cao K-A. mixOmics: an R package for 'omics feature selection and multiple data integration. PLoS Comput Biol 13(11): e1005752


}



\author{Kim-Anh Lê Cao, Benoit Gautier, Francois Bartolo, Florian Rohart.}

\seealso{\code{\link{splsda}}, \code{\link{predict.splsda}} and http://www.mixOmics.org for more details.}

\examples{
data(liver.toxicity)
# note: we made up those data, pretending they are repeated measurements
repeat.indiv <- c(1, 2, 1, 2, 1, 2, 1, 2, 3, 3, 4, 3, 4, 3, 4, 4, 5, 6, 5, 5,
6, 5, 6, 7, 7, 8, 6, 7, 8, 7, 8, 8, 9, 10, 9, 10, 11, 9, 9,
10, 11, 12, 12, 10, 11, 12, 11, 12, 13, 14, 13, 14, 13, 14,
13, 14, 15, 16, 15, 16, 15, 16, 15, 16)
summary(as.factor(repeat.indiv)) # 16 rats, 4 measurements each

# this is a spls (unsupervised analysis) so no need to mention any factor in design
# we only perform a one level variation split
design <- data.frame(sample = repeat.indiv)

tune.splslevel(X = liver.toxicity$gene,
Y=liver.toxicity$clinic,
multilevel = design,
test.keepX = c(5,10,15),
test.keepY = c(1,2,5),
ncomp = 1)

}

\keyword{regression}
\keyword{multivariate}
