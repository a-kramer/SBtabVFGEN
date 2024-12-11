# SBtab

We use _Systems Biology tables_ (SBtab) because this format can be
extended to contain additional information more easily (such as
experimental data, and conditions under which data was measured).

To make a conversion feasible, we decided on a set of columns and
tables (some specified by `TableName` some also by `TableType`
according to the official specification) which have to be present for
the conversion to work.

In contrast to the [official
documentation](https://www.sbtab.net/sbtab/default/documentation.html),
we need all tables to be kept in their own respective `.tsv` file (not
all tables in one huge tsv file), or different sheets of the same
`ods` file (the official SBtab project uses `.xlsx`).

We don't use any of the code from the original [SBtab
authors](https://www.sbtab.net/sbtab/default/team.html).

The following Sections have more information on specific tables and
their required columns (in addition to the obvious `!ID` and `!Name`
columns).

## Text files, Unicode and ASCII

Sometimes, spreadsheet software introduces non-ascii Unicode characters such as
`−` ('MINUS SIGN' U+2212, in html this is «\&minus;») into the
document. They should be replaced by the ascii character `-`
('HYPHEN-MINUS' U+002D) after export to text-based formats. And
similarly for other Unicode characters, unless they appear in comments
(or generally unparsed content). _Minus_ and _Hyphen_ can look quite
similar: `−-` (depending on the chosedn font), but the ascii _hyphen_ is the
character that programming languages understand as _subtraction_.

You can check your files for non-ascii characters like this:
```bash
grep -P '[^[:ascii:]]' *.tsv
#OR
grep -n '[^a-z_A-Z[:digit:][:punct:][:space:]]' *.tsv
# automatic minus to hyphen replacement:
sed -i 's/−/-/g' *.tsv
```
Neither `grep` nor `egrep` define the `[:ascii:]` character class
without the `-P` option for _perl regular expressions_. If such
characters appear outside of Formulas and IDs, they may be harmless.

Of course you can use [perl](https://www.perl.org/) directly, or
anything else that has regular expressions.

## ID and Name of Quantites

All tables require a unique `!ID` column (the ID can be seen as a key
for associative arrays; _dictionaries_, _hash tables_). The
`!Name` column must be unique as well and the entries should work as
variable names in the language that you plan to convert the model to
(in some formats the rules for `Name` entries are more lenient than
`ID` strings). We see no good reason to have both IDs and Names, so
using the same string for both is fine.

The script in this repository uses the `make.names()` function on this
column. This will make them unique, but almost certainly break the
model, the reason is that human error will lead to infromative error
messages later on instead of wrong simulation results. Using the same
ID multiple times may otherwise go unnoticed and merge two distinct
entities (species/parameters/erc.). This will also replace `.` with
`_` in all names (dots are allowed in R variable names, but often
illegal in other languages).

## Scale

Many numbers can be given in a specified scale (like `log`), these
numbers will be converted to linear scale when a model file is written
to file.

Let a quantity `y` be measured in unit `M` (`y` is a number
followed by a unit, `y/M` is just a number), and `!Scale` be set to
`log10`, then the number you write in the `![Default]Value` column is
`z=log10(y/M)`. The script will do the inverse to generate the model
and pass the unit on to `.mod` files. Here are some examples:

|     z | y         | scale indicator     |
|------:|:---------:|:--------------------|
|  `-5` | $10^{-5}$ | `log10`             |
|       | 1.0E-5    | `base-10 logarithm` |
|------:|:---------:|:--------------------|
| `1.2` | 3.32      | `log`               |
|       | 3.32      | `ln`                |
|       | 3.32      | `natural logarithm` |
|------:|:---------:|:--------------------|
| `1.6` | 1.6       | `lin`               |
|       | 1.6       | `linear`            |

## Compound

This table defines the compounds that are supposed to be modeled by
state variables and are subject to change by the reactions in the
systems.

|        Column | Values                       | Comment                                                                                                       |
|--------------:|:----------------------------:|:--------------------------------------------------------------------------------------------------------------|
|        !Scale | log, log10, linear           | and some variants of these                                                                                    |
| !InitialValue | a number                     | (per unit) in the above scale                                                                                 |
|         !Unit | the unit of the above number | as it would be in linear scale                                                                                |
|  !SteadyState | `TRUE`                       | this compound should reach a steady state in at least one scenario and you want to know whether this happened |
|               | `FALSE`                      | it is not important whether or not this compound reaches steady state                                         |
|   !Assignment | `Name` or `ID`               | this field will assign a pre-defined algebraic expression to the compound                                     |

The conversion script will make a file called
`[…]SuggestedOutput.tsv`, it will have lines that can be used to check
whether a compound has reached steady state (or not), this is done for
each compound that has `!SteadyState` marked as `TRUE`. If that output
is close to `0`, then steady state was reached (it's the sum of all
fluxes for the compound in question).

Others columns are unused but may be informative to the user, or others.

### Compound Assignments

In some cases, a compound's amount or concentration is not supposed to
be governed by reactions (kinetic laws, stoichiometry) but rather by a
fixed (time-dependent) value. In SBML this is called a boundary
condition. If the species is supposed to be constant, the field
`!IsConstant` can be set to `TRUE`; otherwise, you can assign the
value of an expression, listed in the `Expression` table to this
compound. The `!Assignment` field can contain the name of an
`Expression`. In SBML, a rule will be created in the `listOfRules`,
that rule will target the boundary condition species:

```xml
      <species    id="PKC_active" 
                name="PKC_active_value" 
         compartment="Comp1" 
initialConcentration="0" 
      substanceUnits="substance" 
   boundaryCondition="true"/>
<!-- ......... -->
<!-- and later -->
<!-- ......... -->
      <assignmentRule variable="PKC_active">
        <math xmlns="http://www.w3.org/1998/Math/MathML">
          <apply>
            <plus/>
            <ci> PKC_DAG_AA_p </ci>
            <ci> PKC_Ca_memb_p </ci>
            <ci> PKC_Ca_AA_p </ci>
            <ci> PKC_DAG_memb_p </ci>
            <ci> PKC_basal_p </ci>
            <ci> PKC_AA_p </ci>
          </apply>
        </math>
      </assignmentRule>

```

In the other formats, `vf` and `mod`, the relationship is much simpler:
```xml
<Expression Name="PKC_active_value" Description="defined expression Ex0" Formula="PKC_DAG_AA_p+PKC_Ca_memb_p+PKC_Ca_AA_p+PKC_DAG_memb_p+PKC_basal_p+PKC_AA_p"/>
<Expression Name="PKC_active" Description="defined expression S11" Formula="PKC_active_value"/>
```
which will lead to code such as:
```matlab
PKC_active_value = PKC_DAG_AA_p+PKC_Ca_memb_p+PKC_Ca_AA_p+PKC_DAG_memb_p+PKC_basal_p+PKC_AA_p;
PKC_active = PKC_active_value;
```

A blank cell, `NONE`, `FALSE`, or `NO` means that there is no
Assignment for this species/compound. Empty cells can be tricky if export or
import functions merge multiple delimiters (so `\t\t` is not
recognized as an empty cell).

A `TRUE` value in `!IsConstant` and meaningful assignments are
mutually exclusive and may lead to weird results.

## Parameters

| Column | Values  | Comment |
| -----: | :-----: | :------ |
| `!Scale` | `log`, `log10`, `linear` | some aliases of these are possible (such as `base-10 logarithm`)|
| `!DefaultValue` | a number | in above scale, normalized to the unit of measurement, possibly subject to fitting/sampling |
| `!Std` | a number | standard deviation / uncertainty of this parameter |
| `!Min` and `!Max` | numbers | respectively, used if `!Std` is not present |

The columns `!Std` and `!Min`/`!Max` are only used in
sampling/optimization, the model conversion is unaffected by them, the
DefaultValue is passed on to the model files (if there is a place to
put them).

## Reactions

The column `!ReactionFormula` determines the stoichiometry of the
model, the `!KineticLaw` column determines the flux of the given
reaction. Both er required and are standard columns in SBtab.

| Column | Values  | Comment |
| -----: | :-----: | :------ |
| !KineticLaw | e.g. `kf*A*B-kr*AB` | the flux, as a math expression |
| !ReactionFormula | e.g. `A+B<=>AB` | so, `AB` will increase and both `A` and `B` will decrease by this reaction whenever the flux is positive |

Since the kinetic law determines the reversibility of the reaction, the
column `!IsReversible` is not necessary, but if you determine the
kinetics based on the law of mass action it may be important for you
to have that column as a reminder (for when you are auto generating
the `!KineticLaw` column, which this script doesn't do).

## Input

The input parameters to the model that distinguish different
experiments. These quantities are known and can be influenced by the
people who are performing an experiment (or rather the real
counterparts of these quantities can be influenced). These play the
roles of (additional) parameters, but a different kind of parameter
than in the Parameter table. Experiments are supposed to have the same
parameters of the normal kind and different parameters of the input
kind.

The `!DefaultValue` column serves the same role as with normal
parameters, but these are set to known values during an _experiment_
(these have to be known values).

## Output

The outputs are _observable quantities_ of this system; what is and
isn't an output depends on what you can measure (or have knowledge
about). Outputs are usually converted to functions in the target
language. _Experimental Data_ and _Outputs_ are intimately related as
the outputs are the model's equivalent of the data and in some way
those can be compared to one another.

It is possible to include the measured data in other sheets, that data should be
stored together with an estimate of the measurement noise
levels. Regardless of the nature of the noise and underlying
distributions we use `!ErrorName` to indicate which column (the one
that has this name) is storing information about this measurement
error.

| Column | Values  | Comment |
| -----: | :-----: | :------ |
| `!ErrorName` | a string | indicates the column in data sheets that hold the measurement error of an observable |
| `!ErrorType` | not used | this is for the user |
| `!ProbDist`  | a string | the probability distribution of the noise model (currently unused); this is for humans to read |
| `!Formula`   | a math expression | the right hand side of the assignment |

Many data columns may share the same Error column. This is useful if
you have only a very rough estimate of the noise anyway, and the
outputs are in the same number range so using the same standard
deviation (etc.) for all data points seems good enough.

The data sheets are not used by `sbtab_to_vfgen()`, but they are used for Parameter Estimation, e.g. via
[MCMC](https://github.com/a-kramer/mcmc_clib).

## Expression

These will be local variables of the model. These variables store a
one line mathematical expression for reusability of the resulting
value.

Expressions are assignments that are calculated repeatedly each
time the ODEs right hand side is called (before fluxes are calculated).

The `!Formula` column stores a string math expression that will serve
as the right hand side of the assignment.

## Experiments

This table holds the mapping between input parameters and data
sheets. It determines the conditions under which a data set should be
replicated using the model. The _conditions_ als include _events_ that
need to happen to replicate an experiment.

|   Column | Values        | Comment |
| -------: | :-----------: | :------ |
|  `!Type` | `Time␣Series` |  indicates that the data is a `t->output` mapping |
|          | `Dose␣Response`| data sheet is an input/output curve, i.e. `input->output` mapping |
|`>some_id`| a number      | sets the input parameters for this experiment |
| `!Event` | a table name  | the name of an event table that holds time instantaneous model state changes |

This is not used by this converter, but useful for parameter fitting
and interpretation of _the input_ and _output_ concepts.

## Events

Events are sudden changes in state variables or inputs at
speciefied points in time.

The table `Transformation` defines what sudden transformation the
model is capable of. For each experiment an additional _event
schedule_ table assigns time points to transformation events.

```
!!SBtab TableName='Transformation' Document='MyModel'
 !ID  >A        >B     >C
 TF1  A-B*dose  B      C
 TF2  A         B/2    C
 TF3  A*0.1     B*0.1  C*0.1
```

Each row describes a named transformation, where the column names make
a reference to the affected variable. Row `TF1` means `A <- A-B*dose`
(other variables are unaffected).

Events will interrupt the ODE solver in the rgsl package (not other
solvers), change the values of input parameters, or state variables
according to the available named tranformations (see above). After the
transformation the solver is re-initialized to continue from there.

```
!!SBtab TableName='EventScheduleAlpha' Document='myModel'
  !TimePoint  !Time  !Tranformation !Dose
 Event0Time0  10.0   TF1            0.1
 Event0Time0  11.0   TF2            0.2
         ...  ...    ...
```

The event schedule creates a time sequence for transformation events.

```
!!SBtab TableName='Experiments' Document='myModel'
               !ID  !Type          !Time  !T0    !Event              !Citation
    Smith2019Fig2E  Time Series    600.0  -100.0  EventScheduleAlpha  https//doi.org/10.[...]
               ...  ...            ...    ...     ...                ...
```

And finally the table of experiments will assign each exeriment (rows)
the appropriate sequence of transformation events. Event formulas may
also refer to the `dose` of a transformation, and a `!Dose` column in
the schedule can assign each time-point an event dose that is local to
that specific time and event.

### Old Event System

An older event type assumed linear transformations (and affine
transformations). This is discouraged now.

Eventtables have a `!Time` column and an effect column, where the
header combines a target and operation: `>OPERATION:TARGET`. An example:

| !TimePoint | !Time | >ADD:Ca |
|-----------:|:-----:|:--------|
|        TP0 | 100   | 1e3     |
|        TP0 | 200   | 1e3     |
|        TP0 | 300   | 1e3     |


In this example, we add 1000 units of Ca, every 100 time-units. The
possible operations are: `SET`,`ADD`,`SUB`,`MUL`,`DIV`.
