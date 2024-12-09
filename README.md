# SBtabVFGEN

Convert a Model written in [SBtab](https://www.sbtab.net/), saved as a
series of `tsv` files or alternatively an [Open Document
Spreadsheet](https://www.documentfoundation.org/) `ods` to a
[VFGEN](https://warrenweckesser.github.io/vfgen/) vector field file
`vf`.

This project is supported by [EBRAINS](https://ebrains.eu)
infrastructure and the [Human Brain
Project](https://www.humanbrainproject.eu), with more detail in [acknowledgements](./ACKNOWLEDGMENTS.md).

## Install

Using `remotes`:
```R
remotes::install_github("a-kramer/SBtabVFGEN")
```
You may have to check `.libPaths()` to verify that it includes a path
that you have permission to write to (this is just generally the case,
not just for this package).

Currently, this will work on platforms that have
[R](https://www.r-project.org/). But, any user who needs SBML output
must install the libSBML package for R. 

This is entirely optional, but SBML files may be useful to share
models with others. This is not easy, as `libSBML` is not a
[cran](https://cran.r-project.org/) package, nor is it on
github/lab. But browsing the [SBML
repository](https://sourceforge.net/p/sbml/libsbml/) on sourceforge
makes it possible to find the right version of the [R
interface](https://sourceforge.net/projects/sbml/files/libsbml/5.18.0/stable/R%20interface/).

Download the `tar.gz` file `libSBML_*.tar.gz` in the appropriate
version, and the install it as a package using

```sh
$ R CMD INSTALL liblibSBML_*.tar.gz
```

## Purpose

This model conversion tool can be used by scientists working in the
field of _systems biology_ and all adjacent fields that work with
_ordinary differential equation_ (ODE) models. 

It can be helpful when collaborating with other researchers as it keeps
the model separate from any programming language choice. The user writes the model
in SBtab form, a simple, human readable format; afterwards this SBtab
model can be converted to an ODE and further processed via `vfgen`.

The final result is code for the ODE _right hand side_ function and
analytical _jacobian function_ (among other things) in the chosen
programming language.

This tool prepares a model _M_ for use in numerical analysis application
such as parameter estimation:

```
  User written:                           generated
  +-----------+      +-----------+      +------------+      +----------+
  |           |      |           |      |  (CVODE)   |      |          |
  | SBtab (M) +--+-->+   VFGEN   +----->+  ODE code  +----->+   MCMC   |
  |           |  |   |           |      | +jacobian  |      |  (e.g.)  |
  +-----------+  |   +-----------+      +------------+      +----------+
                 |
            +----+--------------+
            |                   |
            | sbtab_to_vfgen()  |
            |                   |
            +-------------------+
```

The above sketch is an illustration of this tools location within a
larger workflow (context).

As an alternative to VFGEN, we have written a less powerful tool that
just creates C and R code (while VFGEN covers many langauges):
[RPN-derivative](icpm-kth/RPN-derivative)
([forked](andrei-k/RPN-derivative)), which includes the shell script
`sh/ode.sh`. The shell script uses the `derivative` program (same
repository, written in C) to calculate model-jacobians analytically.

## Usage Example

Within an interactive R session called from the folder that contains
the tsv files:
```R
library(SBtabVFGEN)
model.tsv <- dir(pattern="[.]tsv$");
model.sbtab <- sbtab_from_tsv(model.tsv)
sbtab_to_vfgen(model.sbtab)
```

## Shell Interface via `Rscript`

An additional R [script](./R/sbtab_to_vfgen) can be called from the commandline; it uses `sbtab_to_vfgen.R` internally:
```bash
$ alias sbtab_to_vfgen='.../path/to/R/sbtab_to_vfgen'
$ sbtab_to_vfgen *.tsv
```
This will work if the `.tsv` files have acceptable SBtab content.

## Other Output Formats 

As a by-product, [sbtab_to_vfgen()](./R/sbtab_to_vfgen.R) also produces
a `mod` file intended as a starting point for use in
[neuron](https://neuron.yale.edu/neuron/). This is not the primary
purpose of this function.

If _libSBML_ is installed and R bindings available, an attempt will be
made to produce an SBML file (see further below).

The type of _systems biology_ models that we have in mind are
(plain) ordinary differential equations (ODEs). 

_NEURON_ comprehensively models biochemistry and electrophysiology
(e.g. membrane action potentials). These two different simulations
have to be coupled (neuron does that). This aspect is missing from the
SBtab files we typically use, so the produced `mod` file is really
only an initial point; the user has to change and adapt the file to
the intended purpose and make it work inside of neuron. The user must
also be aware of NEURONs units and make the necessary unit
conversions at some point.

There is no automatic conversion of units (yet).


## VFGEN Output

[VFGEN](https://github.com/WarrenWeckesser/vfgen) is a very useful
tool that reformats an ODE (given in vfgen's xml format) and
convert it into various programming languages, including right hand
side functions for two ODE solvers in C:
[gsl](https://www.gnu.org/software/gsl/doc/html/ode-initval.html) and
[cvode](https://computing.llnl.gov/projects/sundials/cvode). While
doing so, it uses [GiNaC](https://ginac.de/) to calculate the model's
Jacobian analytically (among other things). 

The R script [sbtab_to_vfgen.R](./sbtab_to_vfgen.R) converts an SBtab
model to vfgen's xml format (with normal math, as string attributes),
among others (file will end in `vf`).

Biological models don't necessarily map uniquely onto ODE models, a
compound can be a state variable or an algebraic assignment or a
constant, this has to be inferred a bit from the SBtab files.

## SBtab

[SBtab](https://www.sbtab.net/) is a tabular format for biochemical
models (as in Systems Biology). It is ~perhaps~ easier to understand
than `sbml` and can be parsed/worked via shell scripts (e.g. with line
oriented tools such as `sed` and `awk`) due to its tabular nature
(stored as e.g. tab separated value text files).

In addition to the official upstream documentation, we have summarised
the SBtab entries that this script can use in
[sbtab.md](./docs/sbtab.md). The SBtab specification does not go into
detail about many use-cases, so some interpretation on our part was
needed. The SBtab files we create may not adhere perfectly to the
official specs and similarly SBtab files created by official SBtab
software may not work here. This will probably improve over time, but
currently we use no code/software from the SBtab authors (upstream).

## Open Document Format, Gnumeric and Spreadsheets in General

Even though `.tsv` files are more fundamental types (have least amount
of prerequisites), `.ods` files keep all the sheets in one file and
are slightly more convenient. [Gnumeric](http://www.gnumeric.org/) is
a spreadsheet software that handles both `ods` and `tsv` files fairly
well (it also has its own format `.gnumeric`)

Conversion between spreadsheet formats like `.ods`, `.gnumeric` and
`.tsv` files is very convenient using `ssconvert`, a part of
[gnumeric](http://www.gnumeric.org/). The shell scripts
[ods_to_tsv.sh](./ods_to_tsv.sh) and [tsv_to_ods.sh](./tsv_to_ods.sh)
in this repository are an example of `ssconvert` usage.

An SBtab document can be imported from an open document spreadsheet
(`.ods`) directly using the
[readODS](https://cran.r-project.org/web/packages/readODS/index.html)
package:

```R
library("SBtabVFGEN")
model.ods <- "examplemodel.ods" 
if (file.exists(model.ods){
 model.sbtab <- sbtab_from_ods(model.ods)
 sbtab_to_vfgen(model.sbtab)
}
```

The result is written to several files (`.vf`,`.mod`, and
`.xml`). Some other results with additional information are also
created.

In either case, whether TSV or ODS was used, `model.sbtab` will be a
list of `data.frame` objects.

Other spreadsheet programs such as _google spreadsheets_ and _libre
office_ export to `tsv` _one sheet at a time_ (with no easy
workarounds) and lack an option to export _N_ sheets into _N_
files. Gnumeric's `ssconvert` command does.

# Systems Biology Markup Language (SBML level 2 version 4)

The program `sbtab_to_vfgen.R` also produces an `.xml` file in the _Systems Biology Markup Language_ (SBML).
This is only done, if _libsbml_ is installed with `R` bindings, like this:

```bash
$ R CMD INSTALL libSBML_5.18.0.tar.gz
```

If this check: `if (requireNamespace(libSBML))` succeeds, then the
scripts attempts to make an sbml file. SBML is a format that has
units, and the units defined in SBtab are forwarded to SBML. The
formats are very different with regard to unit handling and math
generally. The method we use to parse human readble text units is
described in [units.md](./docs/units.md).

[Here](./docs/libsbml.md) is a small (incomplete) list of libsbml
functions in R (that we used). An auto-generated [full
list](./libSBMLused.md), without comments, is also present.

There is a more modern level of SBML, level 3, but this package lacks
the ability to generate this format.

