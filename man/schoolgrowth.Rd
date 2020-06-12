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
    given 'block', where a block is defined as the combination of a
    given year, a given grade-level, and a given subject. The dataframe
    \code{d} must have at least the following named fields:
    \code{stuid}, \code{school}, \code{grade}, \code{year},
    \code{subject}, and \code{Y}.  None of these fields can have missing
    values.  The field \code{stuid} is a student identifier that links
    growth scores from the same student in different blocks.  The field
    \code{school} is a school identifier that, for each record,
    indicates which school is associated with the growth score for the
    given student and block.  The fields \code{grade}, \code{year}, and
    \code{subject} are the grade-level, year and subject associated with
    the given growth score.  The field \code{Y} is the growth score,
    which must be numeric.  The field \code{year} must be either
    numeric, or in a format that can be coerced to numeric.
  }

  \item{target}{
    A named character vector of length at most 4 that is
    used to define the target estimand for each school. The function
    computes the direct estimate and EBLP of this target for each school
    (see Details), provided that the school has at least some data in
    the required blocks.
    
    The names of \code{target} must be among the following: 'years',
    'subjects', 'grades', and 'weights'.  Any elements that are missing
    are set to defaults as follows: the 'years' element is set to
    \code{"final"}, the 'subjects' element is set to \code{"all"}, the
    'grades' element is set to \code{"all"}, and the weights element is
    set to \code{"n"}.  Thus if \code{target} is left \code{NULL} by the
    user, the default target estimand for each school will be the
    average growth in the final year of the modeled data, across all
    grades and subjects, where each grade and subject will be weighted
    in proportion to the number of students with growth scores in this
    grade and subject for each school during the final year of data.

    The valid values for the 'years' element of \code{target} are
    \code{"final"}, corresponding to the final year of growth data in
    \code{d}, \code{"all"}, referring to all years, or a character
    string of the form \code{"y1,y2,..."} where 'y1','y2', etc are
    specific values of years that occur in \code{d}.  For example, if
    the growth data have years '1' through '6' and it is desired to
    define the target based on the final two years of data, the value of
    \code{target["years"]} would be \code{"5,6"}.

    The valid values for the 'subjects' element of \code{target} are
    \code{"all"}, referring to all subjects, or a character string of
    the form \code{"s1,s2,..."} where 's1','s2', etc are specific values
    of subjects that occur in \code{d}.

    The valid values for the 'grades' element of \code{target} are
    \code{"all"}, referring to all grades, or a character string of the
    form \code{"g1,g2,..."} where 'g1','g2', etc are specific values of
    grades that occur in \code{d}.

    The combination of the 'years', 'subjects' and 'grades' elements of
    \code{target} define which blocks are part of the target estimand
    for each school.  The 'weights' element of \code{target} specifies
    how these blocks should be weighted for each school to define the
    target estimand. The valid values for the 'weights' element of
    \code{target} are \code{"equal"} or \code{"n"}.  The value of
    \code{"equal"} implies that for a given school, each block that is
    part of the target and has at least one observed growth score will
    receive equal weight.  The value of \code{"n"} implies that for a
    given school, each blocks that is part of the target and has at
    least one observed growth score will receive weight proportional to
    the number of observed growth scores in that block.  For example, if
    the target estimand is defined as math growth in the final year
    across grades 4-6, this target consists of three blocks.  Suppose a
    school has 5,10, and 15 student-level growth scores in these three
    blocks, respectively.  Then setting \code{target["weights"] =
    "equal"} would define the target estimand for this school by putting
    weight \code{1/3} on each of the three blocks, whereas setting
    \code{target["weights"] = "n"} would define the target estimand for
    this school by putting weights \code{1/6}, \code{1/3} and \code{1/2}
    on the three blocks.
    
    See the Examples section for examples for specifying \code{target}.
  }

  \item{target_contrast}{
    An optional named character vector, analogous to \code{target} with
    identical syntax, that is used to specify blocks that will be given
    negative weight in the target estimand.  Thus the combination of
    \code{target} and \code{target_contrast} can be used to define
    target estimands that represent contrasts among blocks.

    See the Examples section for an example for specifying
    \code{target_contrast}.
   }
   
  \item{control}{
    An optional named list of control arguments. Named elements that are
    not among the following list will be ignored.

    \itemize{
      
      \item{quietly: }{
	Logical. If TRUE, will print progress messages for each
	step of the estimation process.  Default is FALSE.
      }	

      \item{school_nmin: }{
	Integer. The minimum number of growth scores required for a
	school to be included in the analysis.  Current behavior is to
	halt with an error message if any school has fewer than
	\code{school_nmin} associated growth scores. Default is 1.
      }

      \item{pattern_nmin: }{
	Integer. The minimum number of students required to define a
	pattern for the pattern-mixture specification of the model (see
	Details).  Patterns with fewer than pattern_nmin students will
	be collapsed with other patterns and modeled with a shared set
	of mean parameters. Default is 100.  If \code{pattern_nmin} is
	the character string \code{"min"}, the function will
	find the smallest value that does not result in
	collinearity between the block-by-pattern dummy variables and
	the school dummy variables, and collapse any patterns with
	counts fewer than this value.
      }
      
      \item{blockpair_student_nmin: }{
	Integer. The minimum number of students required to estimate a
	residual covariance (i.e., an off-diagonal element of the matrix
	\eqn{R}{R}; see Details) for a pair of blocks.  Pairs of blocks
	with fewer than \code{blockpair_student_nmin} students who have
	a growth score in each of the two blocks will have the
	corresponding element of \eqn{R}{R} set to 0.  If any individual
	block has fewer than \code{blockpair_student_nmin} student
	growth measures, the function will halt with an error message.
	Default is 100.
      }

      \item{blockpair_school_nmin: }{
	Integer. The minimum number of schools required to estimate a
	block random effect covariance (i.e., an off-diagonal element of
	the matrix \eqn{G}{G}; see Details) for a pair of blocks.  Pairs
	of blocks with fewer than \code{blockpair_school_nmin} schools
	with at least one growth score in each of the two blocks will
	have the corresponding element of \eqn{G}{G} set to 0.  If any
	individual block has fewer than \code{blockpair_school_nmin}
	schools with growth measures, the function will halt with an
	error message.  Default is 10.
      }
      
      \item{eig.tol: }{
	Numeric. Used to decide when eigenvalues of estimated values of
	\eqn{R}{R} and \eqn{G}{G} are sufficiently small to be treated as
	zero.  Eigenvalues lambda_k are treated as zero when the ratio
	of lambda_k to the largest eigenvalue of the matrix is less than
	or equal to eig.tol.  Default is 1e-06.
      }

      \item{keepDiag_first: }{
	Logical. Used in the coercion of estimated values of \eqn{R}{R}
	and \eqn{G}{G} to be positive semi-definite (PSD).  If
	TRUE, the initial attempt to coerce the matrix to be PSD
	will do so treating the diagonal elements as fixed. If
	FALSE, diagonal elements may be altered in the coercion
	to PSD.  Default is FALSE.
      }

      \item{return_d: }{
	Logical.  If TRUE, the dataset of growth measures
	\code{d} is returned, along with additional fields added by the
	\code{schoolgrowth} function.  Default is FALSE.
      }

      \item{return_schjack: }{
	Logical.  If TRUE, and if jackknifing is used, returns
	extra jackknife information for each school.  This can result
	in large objects if the number of schools and jackknife batches
	are both large.  Default is TRUE.
      }
      
      \item{mse_blp_chk: }{
	Logical.  If TRUE, conducts an auxiliary calculation to
	double-check the estimated mean squared error (MSE) of the
	EBLPs.  Default is FALSE.
      }

      \item{jackknife: }{
	Logical. If TRUE, uses jackknifing of schools to estimate
	a second-order correction term to the estimated MSE of the
	EBLPs.  This term is intended to account for the extra error in
	EBLPs due to the estimation of the matrix \eqn{G}{G}.  The
	computation of this term will increase both computation time and
	RAM footprint, but will result in estimated MSEs that tend to be
	closer to true MSEs.  Default is TRUE.
      }

      \item{jackknife_J: }{
	Can be either the character string \code{"max"}, or a positive
	integer ranging from 1 to the total number of schools with
	growth measures in at least two blocks.  It is the number of
	batches of schools used in the jackknife procedure for EBLP MSE
	estimation if \code{control$jackknife = TRUE}. If
	\code{jackknife_J} is \code{"max"}, standard delete-1
	jackknifing of schools is used.  In some settings this can
	result in long computation times and/or large RAM requirements,
	and a warning will be provided.  If \code{jackknife_J} is
	smaller than the number of schools with growth measures in at
	least two blocks, the \code{schoolgrowth} function creates
	\code{jackknife_J} batches of schools to conduct jackknife
	estimation of the MSEs.  Batches are approximately
	equally-sized, and are determined at random, so that the user
	will need to call \code{set.seed} outside of the function in
	order to reproduce the estimated MSEs.  If \code{jackknife_J}
	exceeds 100, a warning will be issued about possible RAM
	requirements.  The default value is the minimum of 50 and the
	number of schools with growth measures in at least two blocks.
      }	
      
      \item{regularize_Ghat: }{
	Logical. An experimental option for using shrinkage to reduce
	error in the elements of the estimated value of \eqn{G}{G}.
	Default is FALSE.
      }
      
      \item{R: }{
	A member of the class \code{dsCMatrix} of sparse, symmetric
        matrices from the \code{\link{Matrix}} package.  If
        \code{control$R} is supplied, it is treated as the fixed value
        of the matrix \eqn{R}{R} for the entire estimation procedure,
        rather than \eqn{R}{R} being estimated from the growth measures
        in \code{d}.  The value of \code{control$R} must have row names
        and column names equal to the block labels created by the
        \code{schoolgrowth} function. A typical usage may be for
        sensitivity analysis or simulation studies in which the
        \code{schoolgrowth} function is called once, and then subsequent
        calls to the function uses the value of \eqn{R}{R}, or a
        perturbation of it, from the initial fit.  Default is
        \code{NULL} so that \eqn{R}{R} is estimated from the supplied
        data.
      }

      \item{G: }{
	Analogous to \code{control$R}, but for the variance-covariance
	matrix \eqn{G}{G} of the block random effects for each school.  If
	\code{control$G} is supplied, it must have the appropriate row
	names and column names, and must satisfy the sum-to-zero
	constraints implied by the model (see Details and References).
	Default is \code{NULL} so that \eqn{G}{G} is estimated from the
	supplied data.
      }
    }
  }
}

\details{
  Details on the statistical model and estimation procedures are
  provided in Lockwood, Castellano and McCaffrey (2020).  Each record of
  \code{d} provides the growth score for a given student from a given
  'block', where a block is defined as the combination of a given year,
  a given grade-level, and a given subject.  These data are used to
  estimate the parameters of a linear mixed-effects model.

  Let \eqn{B}{B} be the total number of blocks represented in the data
  and \eqn{S}{S} be the total number of schools represented in the data.
  Students are first partitioned into "patterns" defined by the set of
  blocks for which they have observed growth scores. The linear
  mixed-effects model then assumes that the growth score
  \eqn{Y_{i,b}}{Y(i,b)} for a given student \eqn{i}{i} in a given block
  \eqn{b}{b} depends on a pattern-specific, block-specific mean
  parameter, a fixed effect for the school \eqn{s}{s} to which student
  \eqn{i}{i} is associated for block \eqn{b}{b}, a random effect for the
  combination of school \eqn{s}{s} and block \eqn{b}{b}, and a residual
  error. Residual errors across blocks for the same student are allowed
  to be correlated, but are assumed to be independent across students.
  School-by-block random effects are allowed to be correlated within
  schools, are constrained to satisfy sum-to-zero constraints within
  schools for model identifiability, and are assumed to be independent
  across schools. The model parameters are the block-by-pattern means,
  the \eqn{s=1,\ldots,S}{s=1,...,S} school fixed effects, the \eqn{(B
  \times B)}{(B x B)} variance-covariance matrix \eqn{G}{G} of the
  school-by-block random effects, and the \eqn{(B \times B)}{(B x B)}
  variance-covariance matrix \eqn{R}{R} specifying the variances and
  covariances of the residual errors.  The matrix \eqn{G}{G} has rank at
  most \eqn{B-1}{B-1} due to the sum-to-zero constraints on the
  school-by-block random effects.

  The \code{schoolgrowth} function computes estimates of these
  parameters using the moment-based estimation procedure detailed in
  Lockwood, Castellano and McCaffrey (2020).  It then uses these
  parameter estimates, in conjunction with the observed growth scores,
  to construct the EBLP of the target estimand for each school.  The
  target estimand for each school is a linear combination of fixed
  effects and random effects.  It represents a hypothetical average
  growth for a school that would be observed if there were infinitely
  many growth scores observed for the school for each of the blocks that
  are part of the target.  The MSE of the EBLP is estimated for each
  school using a first-order plug-in approximation, and if
  \code{control$jackknife = TRUE}, a second-order term.

  The \code{schoolgrowth} function checks the relationship between
  block-by-pattern indicator (dummy) variables and school indicator
  variables to ensure that the design matrix based on these two sets of
  variables has full column rank.  It halts with an error message if
  this condition is not met.  The condition will not be met if there
  exists a proper subset of block-by-pattern combinations, and a proper
  subset of schools, with the property that a growth score is linked to
  a block-by-pattern combination in the first subset if and only if it
  is also linked to a school in the second subset. For example, such a
  situation could arise in the analysis of data from grades 4-8 where
  each school in the sample served either grades 4-5, or grades 6-8,
  with no schools serving more extensive grade ranges.

  The EBLP operations requires \eqn{R}{R} and \eqn{G}{G} to be PSD.
  The estimated values of these matrices may not be PSD without
  modifications.  The \code{schoolgrowth} function coerces these
  matrices to be PSD using an adaptation of the \code{\link{nearPD}}
  function in the \code{\link{Matrix}} package to accommodate fixed
  zeros in the estimated matrices.  The default behavior for the
  coercion of each of \eqn{R}{R} and \eqn{G}{G} is to first use
  \code{\link{nearPD}} with no restrictions.  If this does not converge,
  a second attempt is made using \code{\link{nearPD}} with the diagonal
  values of the matrix held constant.  The order of these two steps can
  be switched by setting \code{control$keepDiag_first=TRUE}.  If neither
  step converges, the a PSD approximation is constructed using a
  spectral decomposition and setting all eigenvalues that are less than
  (\code{control$eig.tol} times the maximum eigenvalue) to zero.
}

\value{
  A list, with elements specified below.  Most of these elements contain
  meta-data provided for diagnostic purposes and will not be needed by
  many users.  The most important element is
  \code{aggregated_growth}. We describe this element first, and then the
  remaining elements.
  
  \itemize{
    \item{aggregated_growth: }{
      A dataframe with one record per school providing the EBLP estimate
      of the target estimand for each school, along with other
      information.  Fields include: \code{school}, the school ID
      variable from \code{d}; \code{gconfig}, a character string
      summarizing the "grade configuration" of the school, defined here
      as the unique set of grades with growth scores linked to the
      school over all records in \code{d}; \code{ntotal}, the total
      number of growth scores linked to the school over all records in
      \code{d}; \code{ntarget}, the number of growth scores linked to
      the school in blocks that are part of the target;
      \code{ncontrast}, the number of growth scores linked to the school
      in blocks that are part of the contrast (if \code{target_contrast}
      is not \code{NULL}); \code{est.direct}, the "direct estimate" of
      the target estimand, equal to the appropriately-weighted average
      of the growth scores in the target blocks; \code{mse.direct}, the
      estimated MSE of the direct estimate; \code{est.blp}, the EBLP of
      the target estimand; \code{mse.blp}, the estimated MSE of the
      EBLP; \code{est.hybrid}, equal to \code{est.blp} if \code{mse.blp}
      is less than or equal to \code{mse.direct} and otherwise equal to
      \code{est.direct}; \code{mse.hybrid}, equal to the minimum of
      \code{mse.blp} and \code{mse.direct}; and \code{prmse.direct},
      equal to 1 - (\code{mse.blp} / \code{mse.direct}), the estimated
      proportional reduction in MSE for the EBLP relative to the direct
      estimate.

      For schools that do not have growth scores in any of the blocks
      that are part of the target estimand, the fields
      \code{est.direct}, \code{mse.direct}, \code{est.blp},
      \code{mse.blp}, \code{est.hybrid}, \code{mse.hybrid}, and
      \code{prmse.direct} are \code{NA}.
    }
    
    \item{control: }{
      The value of \code{control} used during estimation, including
      defaults for elements that are not specified by the user in the
      function call.
    }

    \item{target: }{
      The value of \code{target} used during estimation.
    }

    \item{target_contrast: }{
      The value of \code{target_contrast} used during estimation
      (\code{NULL} if \code{target_contrast} is not used).
    }
    
    \item{dblock: }{
      A dataframe with \eqn{B}{B} rows providing block-level meta-data.
      There is one row for each block.  Fields include the year, the
      grade level, the subject, the block label, a numeric block ID
      variable ranging from \eqn{1}{1} to \eqn{B}{B}, the number of
      growth scores in the data for each block, and a logical variable
      indicating whether each block is part of the target.  If
      \code{target_contrast} is specified, \code{dblock} also includes a
      logical variable indicating whether each block is included in the
      contrast.
    }

    \item{dblockpairs: }{
      A dataframe with \eqn{B(B+1)/2}{B(B+1)/2} rows providing meta-data
      for blocks and all pairwise combinations of blocks.  Rows are in
      column-major order for a symmetric matrix. Fields include the
      block labels for each member of the pair, the numeric block ID for
      each member of the pair, the number of schools with growth scores
      in both members of the pair, the number of students with growth
      score in both members of the pair, and information about the
      estimated matrices \eqn{G}{G} and \eqn{R}{R}.  The field
      \code{G_est} is a logical variable indicating whether each pair of
      blocks has an estimated random-effect covariance parameter (vs the
      parameter being set to 0).  The field \code{Graw} contains the
      estimated covariance parameters prior to PSD adjustment, whereas
      the field \code{G} contains the estimated covariance parameters
      after PSD adjustment.  The field \code{R_est}, \code{Rraw} and
      \code{R} are analogous for the error variance-covariance matrix.
    }

    \item{dsch: }{
      A list with one element per school containing meta-data for each
      school used during the estimation procedure.  Fields that are most
      likely to be of interest to users include: \code{school}, the
      school ID; \code{oblocks}, a vector of block IDs containing the
      set of blocks for which the school has associated growth measures;
      \code{nblock}, the length of \code{oblocks}; \code{N}, a \eqn{(B
      \times B)}{(B x B)} sparse symmetric matrix indicating the counts
      of growth scores in each block for the school (diagonal elements)
      and the counts of students with growth scores in both elements of
      each pair of blocks (off-diagonal elements); \code{mu}, the
      estimated school fixed effect; \code{var_muhat}, the estimated
      variance of \code{mu}; \code{tab}, a dataframe containing
      block-level average growth scores and block-level EBLPs for the
      school, along with other information; \code{R_sb}, a sparse
      symmetric matrix estimating the variance-covariance matrix of the
      student-level sampling errors of the block-level average growth
      scores for the school; \code{mse_blp}, the estimated MSE of the
      block-level EBLPs for the school; and \code{weights}, a matrix
      providing the weights that each block received in the computation
      of the direct estimate and EBLP for the school.  The
      \code{weights} element is \code{NULL} for schools with no growth
      measures in any of the target blocks.
    }

    \item{bhat_ols: }{
      The estimated block-by-pattern means and school fixed effects
      obtained during first-stage OLS estimation.  Note that the
      estimated school fixed effects in \code{bhat_ols} will not
      generally correspond to the final estimated school fixed effects
      (the \code{mu} elements of \code{dsch}), which are obtained using
      generalized least squares in a later estimation step.
    }

    \item{modstats: }{
      A list containing miscellaneous summary statistics from the data
      and estimation procedure.  Elements include: \code{ntot}, the
      total number of growth scores used in the estimation procedure
      (i.e., \code{nrow(d)}); \code{nstu}, the number of unique students
      in \code{d}; \code{nsch}, the number of unique schools in
      \code{d}; \code{varY}, the sample variance of the growth scores in
      \code{d}; \code{estimated_variance_among_schools}, the estimate of
      the variance among the school fixed effects; and
      \code{estimated_percvar_among_schools}, the estimate of the
      percentage of variance in the growth scores explained by the
      school fixed effects.
    }

    \item{tab_patterns: }{
      A dataframe summarizing the student patterns.  Fields include
      \code{pattern}, a character string consisting of 0s and 1s
      indicating the blocks in which the pattern does or does not have
      observed growth scores; \code{pcount}, the number of students
      represented by the pattern; and \code{cpattern}, equal to
      \code{pattern} if the pattern is maintained during the estimation
      procedure, or \code{"collapsed"} if the pattern was collapsed
      during the estimation procedure.
    }

    \item{G: }{
      The estimated variance-covariance matrix \eqn{G} of the
      school-by-block random effects, after any PSD adjustment. It is a
      member of the class \code{dsCMatrix} of sparse, symmetric matrices
      from the \code{\link{Matrix}} package.
    }

    \item{R: }{
      The estimated variance-covariance matrix \eqn{R} of the residual
      errors, after any PSD adjustment. It is a member of the class
      \code{dsCMatrix} of sparse, symmetric matrices from the
      \code{\link{Matrix}} package.
    }

    \item{d: }{
      If \code{control$return_d = TRUE}, a dataframe similar to \code{d}
      but with additional fields created by the function.
    }
    
    \item{G_jack: }{
      If \code{control$jackknife = TRUE}, a list of length equal to the
      number of jackknife batches used in the EBLP MSE estimation
      procedure.  Each element is a member of the class \code{dsCMatrix}
      of sparse, symmetric matrices from the \code{\link{Matrix}}
      package, equal to the estimated value of \eqn{G} for each
      jackknife batch.
    }

    \item{Gstar_jack: }{
      If \code{control$jackknife = TRUE}, a list of length equal to the
      number of jackknife batches used in the EBLP MSE estimation
      procedure.  Each element is a member of the class \code{dsCMatrix}
      of sparse, symmetric matrices from the \code{\link{Matrix}}
      package, equal to the estimated value of \code{G*} for each
      jackknife batch.  Each is a \eqn{((B-1) \times (B-1))}{((B-1) x
      (B-1))} matrix corresponding to the estimable parameters of
      \eqn{G} in light of the sum-to-zero constraints.  See Lockwood,
      Castellano and McCaffrey (2020) for details.
    }

    \item{varhat_G: }{
      If \code{control$jackknife = TRUE}, a jackknife variance estimate
      of each element of the estimated value of \eqn{G}.  The estimated
      variances are organized as a \eqn{(B \times B)}{(B x B)} matrix,
      so that the \eqn{(i,j)}{(i,j)} element of \code{varhat_G} is the
      jackknife variance estimate of the \eqn{(i,j)}{(i,j)} element of
      \eqn{G}.
    }
    
  }
}

\references{
  Lockwood J.R, Castellano K.E. and McCaffrey D.F. (2020).
  \dQuote{Improving accuracy and stability of aggregate student growth
  measures using best linear prediction}, \emph{Unpublished}.
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

\dontrun{
## changing the target
##
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

## defining a target that is the difference in math growth between years
## 2 and 1, where growth in each year is student-weighted across grades:
tmp <- schoolgrowth(growthscores,
  target =          c(years="2", subjects="math", grades="all", weights="n"),
  target_contrast = c(years="1", subjects="math", grades="all", weights="n"))
print(tmp$dsch[[1]]$weights)
}

}

\keyword{models}
