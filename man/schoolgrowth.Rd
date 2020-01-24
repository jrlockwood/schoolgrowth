\name{schoolgrowth}
\alias{schoolgrowth}

\title{Smoothing School-Level Growth Measures Using Empirical Best Linear Prediction}

\description{Uses student-level growth scores for particular years,
  grades and subjects to compute direct estimates, and empirical
  best linear predictors (EBLPs), of specified school-level aggregate
  growth targets.}

\usage{
schoolgrowth(d, target = NULL, target_contrast = NULL, control = list())
}

\arguments{
  \item{d}{
    TODO
  }

  \item{target}{
    TODO
  }

  \item{target_contrast}{
    TODO
  }

  \item{control}{
    TODO
  }
}

\details{
  TODO describe student-level model \code{Y = mu + alpha + theta + error}.
}

\value{A list with elements

  \itemize{
    \item{item1: }{TODO}

    \item{item2: }{TODO}
  }

  The fields of \code{TODO} are:

  \itemize{
    \item{field1: }{TODO}

    \item{field2: }{TODO}
  }
}

\references{
  Lockwood J.R, Castellano K.E. and McCaffrey D.F. (2020).
  \dQuote{Improving accuracy and stability of aggregate student
    growth measures using best linear prediction},
  \emph{Unpublished}
}

\author{
  J.R. Lockwood \email{jrlockwood@ets.org}
}

\examples{
data(growthscores)
}

\keyword{models}
