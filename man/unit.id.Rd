% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/units.R
\name{unit.id}
\alias{unit.id}
\title{Converts a unit to a string that works as an identifier}
\usage{
unit.id(unit.str, prnt = FALSE)
}
\arguments{
\item{unit.str}{the original string representastion of that unit}

\item{prnt}{logical switch: if TRUE, the name will be printed.}
}
\value{
unit.id string
}
\description{
Some formats require a name for a unit definition. This functions
creates a name from a unit, converting math/symbols to text. The
returned value should work as an SBML unit id.
}
\examples{
> unit.id("s^9")
[1] "s_to_the_power_of_9"

> unit.id("cm^2")
[1] "cm_square"

> unit.id("1/s")
}
