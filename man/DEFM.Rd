% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R, R/defm-package.R
\name{DEFM}
\alias{DEFM}
\alias{new_defm}
\alias{defm}
\alias{init_defm}
\alias{print_stats}
\alias{nterms_defm}
\alias{nrow_defm}
\alias{ncol_defm_y}
\alias{ncol_defm_x}
\alias{nobs_defm}
\alias{morder_defm}
\title{Discrete Exponential Family Model (DEFM)}
\usage{
init_defm(m)

print_stats(m, i = 0L)

nterms_defm(m)

nrow_defm(m)

ncol_defm_y(m)

ncol_defm_x(m)

nobs_defm(m)

morder_defm(m)

new_defm(id, Y, X, order = 1)
}
\arguments{
\item{m}{An object of class \code{DEFM}.}

\item{i}{An integer scalar indicating which set of statistics to print (see details.)}

\item{id}{Integer vector of length \code{n}. Observation ids, for example,
person id.}

\item{Y}{0/1 matrix of responses of \code{n_y} columns and \code{n} rows.}

\item{X}{Numeric matrix of covariates of size \code{n_x} by \code{n}.}

\item{order}{Integer. Order of the markov process, by default, 1.}
}
\value{
An external pointer of class \code{DEFM.}

\itemize{
\item \code{nterms_defm} returns the number of terms in the model.
}

\itemize{
\item \code{nrow_defm} returns the number of rows in the model.
}

\itemize{
\item \code{ncol_defm_y} returns the number of output variables in
the model.
}

\itemize{
\item \code{ncol_defm_x} returns the number of covariates in the model.
}

\itemize{
\item \code{nobs_defm} returns the number of observations (events) in the
model.
}

\itemize{
\item \code{morder_defm} returns the order of the Markov process.
}

An external pointer of class \code{DEFM.}
}
\description{
Discrete Exponential Family Models (DEFMs) are models from the exponential
family that deal with discrete data. Here, we deal with binary arrays which
can be used to represent, among other things, networks and multinomial binary
Markov processes.

Discrete Exponential Family Models (DEFMs) are models from the exponential
family that deal with discrete data. Here, we deal with binary arrays which
can be used to represent, among other things, networks and multinomial
binary Markov processes.
}
\details{
The \code{print_stats} function prints the supportset of the ith type
of array in the model.
}
\examples{
# Loading Valente's SNS data
data(valentesnsList)

mymodel <- new_defm(
  id = valentesnsList$id,
  Y = valentesnsList$Y,
  X = valentesnsList$X,
  order = 1
)

# Adding the intercept terms and a motif from tobacco to mj
term_defm_logit_intercept(mymodel)
term_defm_transition_formula(mymodel, "{y1, 0y2} > {y1, y2}")

# Initialize the model
init_defm(mymodel)

# Fitting the MLE
}
\references{
Vega Yon, G. G., Pugh, M. J., & Valente, T. W. (2022). Discrete Exponential-Family Models for Multivariate Binary Outcomes (arXiv:2211.00627). arXiv. \url{https://arxiv.org/abs/2211.00627}
}
\seealso{
\code{\link[=defm_mle]{defm_mle()}} for maximum likelihood estimation and \code{\link[=loglike_defm]{loglike_defm()}}
for the log-likelihood function.
}
