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
    A dataframe consisting of student-level growth scores.  Each record
    of \code{d} provides the growth score for a given student from a
    given \dQuote{block}, where a block is defined as the combination of
    a given year, a given grade-level, and a given subject.  The fields
    of \code{d} must consist of at least the following: \code{stuid},
    \code{school}, \code{grade}, \code{year}, \code{subject} and
    \code{Y}.  None of these fields can have missing values.  The field
    \code{stuid} is a student identifier that links growth scores from
    the same student in different blocks.  The field \code{school} is a
    school identifier that, for each record, indicates which school is
    associated with the growth score for the given student and block.
    The fields \code{grade}, \code{year}, and \code{subject} are the
    grade-level, year and subject associated with the given growth
    score.  The field \code{Y} is the growth score, which must be
    numeric.  The field \code{year} must be either numeric, or in a
    format that can be coerced to numeric.
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
  TODO describe student-level model \code{Y = mu + alpha + theta +
    error}. TODO note that code forces the end-user to think about and confront
   strata directly, rather than having an automatic solution implemented that
   might not do what someone wants.  the code will not run without
   connectedness.  TODO note that jackknife relies on random partitioning of schools into batches,
   and for replicability, user will have to set.seed() outside the function
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
print(c(nrecords = nrow(growthscores),
        nstudent = length(unique(growthscores$stuid)),
        nschool  = length(unique(growthscores$school))))

## data have B=12 blocks:
print(unique(growthscores[,c("year","grade","subject")]))

## fit model
##
## NOTE: jackknife step skipped with 'jackknife=FALSE' control option
## only to minimize computation time for the example
m <- schoolgrowth(growthscores,
     target = c(years="final", subjects="math", grades="all", weights="n"),
     control=list(quietly=TRUE, mse_blp_chk=FALSE, jackknife=FALSE))

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

## direct, EBLP and hybrid estimates of target for each school,
## along with estimated MSEs
print(m$aggregated_growth)

## changing the target
\dontrun{
## math, final year, grade 4:
tmp <- schoolgrowth(growthscores,
          target = c(years="final", subjects="math", grades="4", weights="n"))

## math and ela, both years, grade 4:
tmp <- schoolgrowth(growthscores,
          target = c(years="1,2", subjects="ela,math", grades="4", weights="n"))

## ela, both years, grades 5,6 equally weighted rather than student-weighted:
tmp <- schoolgrowth(growthscores,
          target = c(years="1,2", subjects="ela", grades="5,6", weights="equal"))
print(tmp$dsch[[1]]$weights)
}

}

\keyword{models}
