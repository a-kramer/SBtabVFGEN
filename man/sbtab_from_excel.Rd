% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sbtab_utilities.R
\name{sbtab_from_excel}
\alias{sbtab_from_excel}
\title{A parser for Excel files with SBtab document structure}
\usage{
sbtab_from_excel(excel.file = dir(pattern = "[.]xlsx?$")[1])
}
\arguments{
\item{excel.file}{a string (file's name)}
}
\value{
SBtab a list of tables (data.frames), one per ods sheet
    SBtab[['TableName']] is a data.frame, comment(SBtab) is
    the name of the document
}
\description{
This function uses the readxl package to read the file.  SBtab
sheets are themselves tables dedicated to a specific type of model
property: Reaction, Compound, Parameter, etc.  The first column
will be used as the row.names. The second row will be used as
column names, the first row is used to determine the document-name
and table-name.
}
\examples{
model.sbtab<-sbtab_from_excel('model.xlsx')
}
\keyword{import}
