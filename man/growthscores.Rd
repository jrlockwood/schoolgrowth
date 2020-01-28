\name{growthscores}
\alias{growthscores}
\docType{data}
\title{Example growth score data}
\description{
  Example growth score data in long format.  Growth scores in two
  subjects (mathematics and English language arts) are provided for
  students in grades 4-6 for each of two consecutive school years.  Each
  record provides the growth score for a particular student in a
  particular year and subject.  Repeated measurements from the same
  student across years and/or subjects are linked with a unique student
  identifier.
}
\usage{data(growthscores)}
\format{
  A data frame with 5797 observations and 6 fields:
  \describe{
    \item{\code{stuid}}{Unique identifier for each student}
    \item{\code{school}}{Unique identifier for each school}
    \item{\code{grade}}{Grade level of the student in the given year (4,5 or 6)}
    \item{\code{year}}{Year (1 or 2)}
    \item{\code{subject}}{Subject ('math' or 'ela')}
    \item{\code{Y}}{Growth score for the given student in the given year
      and subject}
  }
}

\source{
  Anonymous
}

\keyword{datasets}
