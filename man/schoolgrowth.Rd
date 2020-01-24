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
print(head(growthscores))
print(c(nrec = nrow(growthscores), nstudent = length(unique(growthscores$stuid)), nschool = length(unique(growthscores$school))))

## data have B=6 blocks:
print(unique(growthscores[,c("year","grade","subject")]))

## fit model
##
## NOTE: jackknife step skipped with 'jackknife=FALSE' control option
## only to minimize computation time for the example
m <- schoolgrowth(growthscores, target = c(years="final", subjects="math", grades="all", weights="n"), control=list(quietly=TRUE, mse_blp_chk=FALSE, jackknife=FALSE))
print(names(m))

## summary information for blocks and block pairs
print(m$dblock)
print(nrow(m$dblockpairs)) ## B*(B+1)/2
print(head(m$dblockpairs))

## summary information for each school
print(length(m$dsch))
print(names(m$dsch[[1]]))

## miscellaneous meta-data
print(m$modstats)

## summary of student patterns
print(m$tab_patterns)

## estimated variance components
print(m$G)
print(m$R)

## direct, EBLP and hybrid estimates of target for each school
print(m$aggregated_growth)
}

\keyword{models}
