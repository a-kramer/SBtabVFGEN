% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sbtab_utilities.R
\name{time.series}
\alias{time.series}
\title{Create a Time Series Simulation Experiment}
\usage{
\method{time}{series}(
  outputValues,
  outputTimes = as.double(1:dim(outputValues)[2]),
  errorValues = 0.05 * outputValues + 0.05 * max(outputValues),
  inputParameters = c(),
  initialTime = as.double(min(outputTimes)),
  initialState,
  events = NA
)
}
\arguments{
\item{outputValues}{(data.frame) measured data, to be replicated by the
simulation}

\item{outputTimes}{(vector) time values at which the outputValues were
measured}

\item{errorValues}{(data.frame) an estimate of the measurement noise, if
available}

\item{inputParameters}{(vector) a parameter vector that the model needs to
operate (can be some default value to be changed later)}

\item{initialTime}{t0 for the time series experiment: y(t0) = y0}

\item{initialState}{(vector) initial values of the state variables, at t0}
}
\value{
list with these quantities as list items
}
\description{
Given the constituents of a time series simulation experiment,
return a list that rgsl will understand as a simulation.
}
\details{
The output of a real experiment can be a function of the state
variables, not necessarily the state variable trajectories
themselves. We assume that parameter estimation of some sort will
happen.

If an estimate of the measurement error is not available, then
--strictly speaking-- the data is useless: the noise must be
understood as unbounded. But, instead, we assume that the case is
very complex and the user knows what to do about it. We
automatically set the noise to 5% (relative to each value) and
another 5% of the largest value as an estimate of scale, this
represents the absolute error:

tl;dr      default_error = 5% REL + 5% MAXVAL.

The user should replace these values with something, or
simply not use them, if the application goes beyond testing.
}
