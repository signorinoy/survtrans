% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{preprocess}
\alias{preprocess}
\title{Preprocess Survival Data}
\usage{
preprocess(formula, data, group, offset)
}
\arguments{
\item{formula}{A formula specifying the survival model
(e.g., Surv(time, status) ~ covariates).}

\item{data}{A data frame containing the variables referenced in the formula.}

\item{group}{Optional grouping vector.}

\item{offset}{Optional offset vector.}
}
\value{
A list with the following components:
\describe{
\item{x}{A matrix of standardized covariates.}
\item{time}{A vector of survival times, sorted in descending order.}
\item{status}{A vector of censoring indicators, corresponding to the sorted
survival times.}
\item{group}{A factor vector representing the group classification for each
sample, sorted by time.}
\item{offset}{A numeric vector representing offsets for each sample, sorted
by time.}
}
}
\description{
This function preprocesses survival data for a Cox transfer analysis. It
performs several steps including extracting the response and covariates from
a model frame, standardizing the covariates, ensuring the validity of offsets
and grouping, and sorting the data by survival time in descending order.
}
\details{
The function first creates a model frame from the provided formula and data,
then extracts the response variable, which should contain the survival time
and censoring status. The covariates are extracted using a model matrix and
standardized using the scale function. The function also checks the
\code{group} and \code{offset} arguments, assigning default values if they
are not provided. Finally, all components (time, status, covariates, group,
and offset) are sorted in descending order based on the survival time.
}
