.create.sbml.unit <- function(sbml,unit,id){
	n <- nrow(unit)
	unitdef <- libSBML::Model_createUnitDefinition(sbml);
	libSBML::UnitDefinition_setId(unitdef,id);
	for (i in 1:n){
		u <- libSBML::UnitDefinition_createUnit(unitdef);
		Kind <- switch(unit$kind[i],litre="UNIT_KIND_LITRE",metre="UNIT_KIND_METRE",second="UNIT_KIND_SECOND",mole="UNIT_KIND_MOLE",kilogram="UNIT_KIND_KILOGRAM",gram="UNIT_KIND_GRAM","UNIT_KIND_DIMENSIONLESS")
		libSBML::Unit_setKind(u, Kind);
		libSBML::Unit_setExponent(u,unit$exponent[i]);
		libSBML::Unit_setMultiplier(u,unit$multiplier[i]);
		libSBML::Unit_setScale(u,unit$scale[i]);
	}
}

.unique.units <- function(sbml,Unit){
	U <- unique(Unit)
	uid <- unit.id(U,prnt=TRUE)
	unit <- list()
	for (i in 1:length(uid)){
		unit[[i]] <- unit.from.string(U[i])
		.create.sbml.unit(sbml,unit[[i]],uid[i])
	}
	names(unit) <- uid
	return(unit)
}

.sbtab.parameter.to.sbml <- function(sbml,Parameter){
	all.uid <- unit.id(Parameter$Unit)
	num.parameters <- nrow(Parameter)
	Name <- row.names(Parameter)

	if ("Value" %in% names(Parameter)){
		Value <- as.numeric(Parameter$Value)
	} else if ("DefaultValue" %in% names(Parameter)) {
		Value <- as.numeric(Parameter$DefaultValue)
	} else {
		stop("Parameters have no Value column")
	}
	for (i in 1:num.parameters){
		message(sprintf("adding SBtab parameter %i: «%s»",i,Name[i]))
		para <- libSBML::Model_createParameter(sbml)
		libSBML::Parameter_setId(para,Name[i])
		libSBML::Parameter_setName(para,Name[i])

		this.unit.id <- all.uid[i]
		message(sprintf("unit of parameter %i (%s): «%s»",i,Name[i],this.unit.id))
		libSBML::Parameter_setUnits(para,this.unit.id)
		libSBML::Parameter_setValue(para,Value[i])
	}
}

.sbtab.compound.to.sbml <- function(sbml,Compound,CompName,SubstanceUnitID){
	all.uid <- unit.id(Compound$Unit)
	num.species <- nrow(Compound)
	Name <- row.names(Compound)
	for (i in 1:num.species){
		print(Name[i])
		sp <- libSBML::Model_createSpecies(sbml)
		this.unit.id <- all.uid[i]
		message(sprintf("unit of species %i (%s): «%s»",i,Name[i],this.unit.id))
		libSBML::Species_setId(sp, Name[i])

		libSBML::Species_setCompartment(sp, CompName)
		libSBML::Species_setUnits(sp,SubstanceUnitID)
		iv <- as.double(Compound$InitialValue[i])
		if (any(is.na(iv))){
			iv[is.na(iv)] <- 0.0
			message("some initial values were set to 0 because they were not numbers.")
		}
		libSBML::Species_setInitialConcentration(sp, iv)
		libSBML::Species_setHasOnlySubstanceUnits(sp, FALSE)
		Ai <- Compound$Assignment[i]
		if (Compound$IsConstant[i]){
			libSBML::Species_setBoundaryCondition(sp,"true")
			libSBML::Species_setName(sp, Name[i])
		} else if (!grepl("^(0|FALSE|NONE|NULL)?$",Ai)){
			libSBML::Species_setBoundaryCondition(sp,"true")
			libSBML::Species_setName(sp, Ai)
		} else {
			libSBML::Species_setName(sp, Name[i])
		}
	}
}

.sbtab.constant.to.sbml <- function(sbml,Constant){
	all.uid <- unit.id(Constant$Unit)
	num.constants <- nrow(Constant)
	print(num.constants)
	Name <- row.names(Constant)
	for (i in 1:num.constants){
		message(sprintf("adding SBtab Constant %i: «%s»",i,Name[i]))
		print(Name[i])
		const <- libSBML::Model_createParameter(sbml)
		libSBML::Parameter_setId(const,Name[i])
		libSBML::Parameter_setName(const,Name[i])
		this.unit.id <- all.uid[i]
		message(sprintf("unit of constant %i (%s): «%s»",i,Name[i],this.unit.id))
		libSBML::Parameter_setUnits(const, this.unit.id)
		libSBML::Parameter_setValue(const, as.numeric(Constant$Value[i]))
	}
}

.sbtab.reaction.to.sbml <- function(sbml,Reaction,Compartment){
	num.reactions <- nrow(Reaction)
	Name <- row.names(Reaction)
	OneCompartmentModel=FALSE;
	Compartment.id <- row.names(Compartment)
	if (all(is.na(Compartment)) || nrow(Compartment)==1){
		OneCompartmentModel=TRUE;
	}
	AB <- strsplit(Reaction$Formula,split="<=>")
	for (i in 1:num.reactions){
		message(sprintf("adding reaction %i: «%s»",i,Name[i]))
		reaction <- libSBML::Model_createReaction(sbml);

		libSBML::Reaction_setId(reaction,Name[i]);
		libSBML::Reaction_setName(reaction,Name[i]);
		libSBML::Reaction_setReversible(reaction, Reaction$IsReversible[i]);

		message(sprintf("converting flux to mathml: «%s»",Reaction$Flux[i]))
		kl <- libSBML::Reaction_createKineticLaw(reaction);
		Flux <- Reaction$Flux[i]
		if (OneCompartmentModel){
			Flux <- sprintf("(%s)*%s",Flux,Compartment.id[1]);
			message("reaction rates in SBML need to have the unit «substance/time» rather than «concentration/time».\nSo conventional kinetics need to be multiplied by compartment volumes.");
			message("Since this model has only one Compartment, all kinetics will be multiplied by its volume for the SBML model.");
		}
		astMath <- libSBML::parseL3Formula(Flux);
		libSBML::KineticLaw_setMath(kl, astMath);

		A <- ftsplit(AB[[i]][1],"+")
		l <- nzchar(A) & !grepl("^(NULL|NIL|Ø)$",A)
		A <- A[l]
		B <- ftsplit(AB[[i]][2],"+")
		l <- nzchar(B) & !grepl("^(NULL|NIL|Ø)$",B)
		B <- B[l]

		message("[SBML] Reactants: ")
		print(A)
		message("[SBML] Products: ")
		print(B)
		message("---")
		for (a in A){
			r <- regexec(pattern="([0-9]*)\\s*[*]?\\s*(\\w+)",text=a)
			m <- unlist(regmatches(a,r))
			print(m)
			RefName <- m[3]
			spr  <-  libSBML::Reaction_createReactant(reaction)
			libSBML::SimpleSpeciesReference_setSpecies(spr, RefName)
			if (nzchar(m[2])){
				libSBML::SpeciesReference_setStoichiometry(spr,as.numeric(m[2]))
			}
		}
		for (b in B){
			r <- regexec(pattern="([0-9]*)\\s*[*]?\\s*(\\w+)",text=b)
			m <- unlist(regmatches(b,r))
			print(m)
			RefName <- m[3]
			spr  <-  libSBML::Reaction_createProduct(reaction)
			libSBML::SimpleSpeciesReference_setSpecies(spr, RefName)
			if (nzchar(m[2])){
				libSBML::SpeciesReference_setStoichiometry(spr,as.numeric(m[2]))
			}
		}
	}
}

.sbtab.expression.to.sbml <- function(sbml,Expression,CompartmentName,Compound){
	all.uid <- unit.id(Expression$Unit)
	num.species <- nrow(Expression)
	ExpressionName <- row.names(Expression)
	for (i in 1:num.species){
		print(ExpressionName[i])
		this.unit.id <- all.uid[i]
		message(sprintf("unit of expression %i (%s): «%s»",i,ExpressionName[i],this.unit.id))
		F <- Expression$Formula[i]
		if (grepl("^\\s*[+-]?[0-9]+\\s*$",F)){
			message("this expression is a constant, mapping it to a parameter")
			par <- libSBML::Model_createParameter(sbml)
			libSBML::Parameter_setId(par,ExpressionName[i]);
			libSBML::Parameter_setName(par,ExpressionName[i]);
			libSBML::Parameter_setValue(par, as.numeric(F));
			libSBML::Parameter_setUnits(par, this.unit.id);
		} else {
			pat <- paste0("^",ExpressionName[i],"$")
			if (any(grepl(pat,Compound$Assignment))){
				j <- grep(pat,Compound$Assignment)
				stopifnot(length(j)==1)
				VarName <- row.names(Compound)[j]
				if (all(is.na(VarName))) {
					print(ExpressionName[i])
					print(pat)
					print(j)
					print(row.names(Compound)[j])
					print(row.names(Compound))
					stop("Assignment not understood correctly.")
				}
			} else {
				sp <- libSBML::Model_createSpecies(sbml);
				libSBML::Species_setId(sp, ExpressionName[i]);
				libSBML::Species_setName(sp, ExpressionName[i]);
				libSBML::Species_setCompartment(sp, CompartmentName);
				libSBML::Species_setBoundaryCondition(sp, "true")
				VarName <- ExpressionName[i]
			}
			astMath <- libSBML::parseL3Formula(F);
			rule <- libSBML::Model_createAssignmentRule(sbml)
			libSBML::Rule_setVariable(rule,VarName)
			libSBML::Rule_setMath(rule,astMath)
		}
	}
}

.sbtab.output.to.sbml <- function(sbml,Output,CompName,OutputUnitID){
	all.uid <- unit.id(Output$Unit)
	num.species <- nrow(Output)
	Name <- row.names(Output)
	for (i in 1:num.species){
		message(sprintf("output %i:",i))
		print(Name[i])
		sp <- libSBML::Model_createSpecies(sbml)
		this.unit.id <- all.uid[i]
		message(sprintf("unit of output %i (%s): «%s»",i,Name[i],this.unit.id))
		libSBML::Species_setId(sp, Name[i])

		libSBML::Species_setCompartment(sp, CompName)
		libSBML::Species_setUnits(sp,OutputUnitID)
		libSBML::Species_setInitialConcentration(sp, 0.0)
		libSBML::Species_setBoundaryCondition(sp,"true")
		libSBML::Species_setName(sp, Name[i])
		F <- Output$Formula[i]
		astMath <- libSBML::parseL3Formula(F)
		message(sprintf("output %i has formula:",i))
		print(F)
		message(libSBML::formulaToString(astMath))
		rule <- libSBML::Model_createAssignmentRule(sbml)
		libSBML::Rule_setVariable(rule,Name[i])
		libSBML::Rule_setMath(rule,astMath)
	}
}

default_substance_unit <- function(Compound){
	comp_units  <- table(Compound$Unit)
	u_name <- names(sort(comp_units,decreasing=TRUE))[1]
	unit <- unit.from.string(u_name)
	w <- which(unit$kind=='litre')
	if (length(w)==1 && nrow(unit)>1){
		return(unit[-w,])
	} else {
		return(unit)
	}
}

default_time_unit <- function(Parameter){
	par_units  <- Parameter$Unit
	for (u in par_units) {
		unit <- unit.from.string(u)
		w <- which(unit$kind=='second')
		if (length(w)==1){
			u_time<-unit[w,]
			u_time$exponent <- 1.0
			return(u_time)
		}
	}
	warning("Unable to determine default time unit from parameter units.")
	print(table(par_units))
	return(unit.from.string("s"))
}

default_volume_unit <- function(Compartment){
	C_units  <- Compartment$Unit
	for (u in C_units) {
		unit <- unit.from.string(u)
		## try liters
		w <- which(unit$kind=='litre')
		if (length(w)==1){
			return(unit[w,])
		}
		## try cubic meters
		w <- which(unit$kind=='metre' && unit$exponent==3)
		if (length(w)==1){
			return(unit[w,])
		}
	}
	warning("Unable to determine default volume unit from compartment size units.")
	print(table(C_units))
	return(unit.from.string("s"))
}

#' Create an SBML model from Tables
#'
#' This function creates a character vector that contains the
#' biological model in SBML form. The return value can be written into
#' a file to create an SBML document.
#'
#' The input is an SBtab model as data.frames
#'
#' @param H model header, a Model identifying string
#' @param Defaults A data.frame with default units
#' @param Constant Table of Constants for this model [optional]
#' @param Parameter Table of Parameters
#' @param Input Table of input parameters, treated the same as normal Parameters
#' @param Expression will be interpreted and transformed into somthing that SBML understands (optional)
#' @param Reaction Table of reactions
#' @param Compound Table of species (reacting compounds)
#' @param Output Table of output functions (a model of measurement/observables)
#' @param Comp Table of Compartments
#' @return SBML model as a character vector
.make.sbml <- function(H,Defaults,Constant=NULL,Parameter,Input=NULL,Expression=NULL,Reaction,Compound,Output,Comp){
	Doc <- libSBML::SBMLDocument(2,4); # (LEVEL,VERSION)
	sbml  <-  libSBML::SBMLDocument_createModel(Doc)
	libSBML::Model_setId(sbml,H)
	default.id <- row.names(Defaults)
	## default units
	.create.sbml.unit(sbml,default_substance_unit(Compound),"substance")
	.create.sbml.unit(sbml,default_time_unit(Parameter),"time")
	.create.sbml.unit(sbml,default_volume_unit(Comp),"volume")
	.create.sbml.unit(sbml,unit.from.string("micrometer^2"),"area")
	.create.sbml.unit(sbml,unit.from.string("micrometer"),"length")
	##
	u.units <- .unique.units(sbml,c(Expression$Unit,Compound$Unit,Parameter$Unit,Constant$Unit,Input$Unit,Output$Unit,Comp$Unit))
	CompName <- row.names(Comp)
	print(CompName)
	Comp.id <- row.names(Comp)
	for (i in 1:nrow(Comp)){
		this.unit.id <- unit.id(Comp$Unit[i])
		comp <- libSBML::Model_createCompartment(sbml)
		libSBML::Compartment_setId(comp,Comp.id[i])
		libSBML::Compartment_setName(comp,CompName[i])
		libSBML::Compartment_setUnits(comp,this.unit.id)
		libSBML::Compartment_setVolume(comp,as.numeric(Comp$Size[i]))
	}
	if(!is.null(Constant)){
		.sbtab.constant.to.sbml(sbml,Constant)
	}
	.sbtab.compound.to.sbml(sbml,Compound,Comp.id[1],"substance")
	.sbtab.parameter.to.sbml(sbml,Parameter)
	if (!is.null(Input)){
		.sbtab.parameter.to.sbml(sbml,Input)
	}
	.sbtab.reaction.to.sbml(sbml,Reaction,Comp)
	if (!is.null(Expression)){
		.sbtab.expression.to.sbml(sbml,Expression,Comp.id[1],Compound)
	}
	if (!is.null(Output)){
		.sbtab.output.to.sbml(sbml,Output,Comp.id[1],"substance")
	}
	return(Doc)
}
