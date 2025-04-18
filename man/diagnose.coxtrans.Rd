% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coxtrans.R
\name{diagnose.coxtrans}
\alias{diagnose.coxtrans}
\title{Diagnose Cox Transfer Model's Optimization Process}
\usage{
\method{diagnose}{coxtrans}(object, ...)
}
\arguments{
\item{object}{An object of class \code{coxstream}.}

\item{...}{Additional arguments (currently unused).}
}
\description{
Diagnose Cox Transfer Model's Optimization Process
}
\details{
This function produces two plots:
\itemize{
\item Residuals Convergence: Plots the evolution of primal and dual residuals
along with their tolerance levels.
\item Loss Decomposition: Plots the negative log-likelihood, total loss, and
penalty term.
}
}
