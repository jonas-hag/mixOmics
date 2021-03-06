\name{image}
\encoding{UTF-8}
\alias{image.tune.rcc}

\title{Plot the cross-validation score.}

\description{
This function provide a image map (checkerboard plot) of the cross-validation score obtained
by the \code{tune.rcc} function.
}

\usage{
\method{image}{tune.rcc}(x, col = heat.colors, \ldots)
}

\arguments{
\item{x}{object returned by \code{estim.regul}.}
\item{col}{a character string specifying the colors function to use:
\code{\link{terrain.colors}}, \code{\link{topo.colors}},
\code{\link{rainbow}} or similar functions. Defaults to \code{\link{heat.colors}}.}
\item{\ldots}{not used currently.}
}

\details{
\code{image.estim.regul} creates an image map of the matrix
\code{object$mat} containing the cross-validation
score obtained by the \code{estim.regul} function.
Also a color scales strip is plotted.
}

\author{Sébastien Déjean and Ignacio González.}

\seealso{\code{\link{tune.rcc}}, \code{\link{image}}.}

\examples{
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene

## this can take some seconds
cv.score <- tune.rcc(X, Y, validation = "Mfold", plot = FALSE)
image(cv.score)
}

\keyword{dplot}
\keyword{hplot}
\value{none}
