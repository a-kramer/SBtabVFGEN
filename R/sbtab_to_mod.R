#' Writes MOD variable definitions
#'
#' If the number of definitions is low, they all appear in one line
#' Otherwise, this prints one definition per line
#'
#' @param Prefix a prefix to print for every line
#' @param Table a data.frame with variable definitions; we use only the row.names
#' @param Suffix a suffix to print for every line
#' @return a character vector with strings in MOD format
OneOrMoreLines <- function(Prefix,Table,Suffix){
	if (!is.null(Table) && nrow(Table)>0){
		Names <- sprintf("%s %s %s",Prefix,paste0(row.names(Table),collapse=", "),Suffix)
	}else{
		return(character())
	}
	if (nchar(Names)<76){
		return(Names)
	}else{
		return(sprintf("%s %s %s",Prefix,row.names(Table),Suffix))
	}
}

#' A unit string compatible with NEURON
#'
#' NEURON has additional rules for units, this converts a unit-string
#' to comply with the rules.
#' @param unit a string
#' @return string with MOD unit (NEURON)
NeuronUnit<-function(unit){
	unit <- gsub("^[ ]*1/","/",unit)
	unit <- gsub("[()]","",unit)
	unit <- gsub("[ ]*[*][ ]*","-",unit)
	return(unit)
}

.make.mod <- function(H,Constant,Parameter,Input,Expression,Reaction,Compound,Output,ODE,ConLaw=NULL){
	Mod <- list()
	fmt <- list(const="\t%s = %s (%s) : a constant",
			    par="\t%s = %g (%s): %s",
			    input="\t%s  = %g (%s) : an input",
			    total="\t%s = %g : the total amount of a conserved sub-set of states",
			    ConservationLaw="\t%s = %s : conservation law",
			    expression="\t%s : a pre-defined algebraic expression",
			    flux="\t%s : a flux, for use in DERIVATIVE mechanism",
			    comment="\t: Compound %s with initial condition %s had derivative %s, but is calculated by conservation law.",
			    state="\t%s (%s) : a state variable",
			    ode="\t%s' = %s : affects compound %s",
			    reactigon="\t %s <-> %s (%s, %s)",
			    output="\t%s = %s : Output ID %s",
			    assignment="\t%s = %s : assignment for expression %s")
	##    Mod[["header"]] <- "TITLE Mod file for componen"
	Mod[["TITLE"]] <- sprintf("TITLE %s",H)
	Mod[["COMMENT"]] <- sprintf("COMMENT\n\tautomatically generated from an SBtab file\n\tdate: %s\nENDCOMMENT",date())

	Range <- character()
	Range <- c(Range,OneOrMoreLines("\tRANGE",Input,": input"))
	Range <- c(Range,OneOrMoreLines("\tRANGE",Output,": output"))
	Range <- c(Range,OneOrMoreLines("\tRANGE",Expression,": assigned"))
	Range <- c(Range,OneOrMoreLines("\tRANGE",Compound,": compound"))

	Mod[["NEURON"]] <- c("NEURON {",
			             sprintf("\tSUFFIX %s : OR perhaps POINT_PROCESS ?",H),
			             Range,
			             ": USEION ca READ cai VALENCE 2 : sth. like this may be needed for ions you have in your model",
			             "}")
	##
	l <- grepl("0$",row.names(Compound))
	if (any(l)) {
		message(sprintf("possibly problematic names: %s.",paste0(row.names(Compound)[l],collapse=", ")))
		stop("Compound names should not end in '0'. \nNEURON will create initial values by appending a '0' to the state variable names. \nShould your variables contain var1 and var10, NEURON will also create var10 and var100 from them, which will be in conflict with your names. Please choose different names in the source files (SBtab).")
	}
	## Conservation Laws
	if (is.null(ConLaw) || all(is.na(ConLaw))){
		ConservationLaw <- NULL
		ConservationInput <- NULL
		CName <- NULL
		nLaws <- 0
	} else {
		k <- ConLaw$Eliminates
		CName <- row.names(Compound)[k]
		ConservationInput <- sprintf(fmt$total,ConLaw$ConstantName,ConLaw$Constant)
		F <- sprintf("%s - (%s)",ConLaw$ConstantName,ConLaw$Formula)
		nLaws <- length(F)
		ConservationLaw <- sprintf(fmt$ConservationLaw,CName,F)
	}
	Mod[["CONSTANT"]] <- c("CONSTANT {",
			               sprintf(fmt$const,row.names(Constant),Constant$Value, NeuronUnit(Constant$Unit)),
			               "}")
	l=grepl("/.*\\<second\\>",Parameter$Unit)
	if (any(l)){
		Parameter$Unit[l] <- gsub("\\<second\\>","millisecond",Parameter$Unit[l])
		Parameter$Value[l] <- Parameter$Value[l]/1e3;
	}
	r_of_k <- sapply(row.names(Parameter),
		\(p) {
			ptt <- paste0("\\b",p,"\\b")
			return(
				ifelse(any(grepl(ptt,Reaction$Flux)),paste(Reaction$Formula[grep(ptt,Reaction$Flux)],collapse=", "),"does not appear in reactions.")
			)
		}
	)
	Mod[["PARAMETER"]] <- c("PARAMETER {",
			                sprintf(fmt$par,row.names(Parameter),Parameter$Value, NeuronUnit(Parameter$Unit), r_of_k),
			                sprintf(fmt$input,row.names(Input),Input$DefaultValue, NeuronUnit(Input$Unit)),
			                ConservationInput,
			                "}")

	## Expressions and Reaction Fluxes
	Mod[["ASSIGNED"]] <- c("ASSIGNED {",
			               "\ttime (millisecond) : alias for t",
			               sprintf(fmt$expression,row.names(Expression)),
			               sprintf(fmt$flux,row.names(Reaction)),
			               sprintf("\t%s : computed from conservation law",CName),
			               sprintf("\t%s : an observable",row.names(Output)),
			               "}")
	Assignment <- sprintf(fmt$assignment,row.names(Expression),Expression$Formula,row.names(Expression))
	nC <- dim.data.frame(Compound)
	CName <- row.names(Compound)
	STATE=vector(mode="character",length=nC[1])
	DERIVATIVE=vector(mode="character",length=nC[1])
	IVP=vector(mode="character",length=nC[1])
	for (i in 1:nC[1]){
		if (nLaws>0 && i %in% ConLaw$Eliminates){
			message(sprintf("MOD: StateVariable %s will be commented out as it was already defined as a Mass Conservation Law Expression.",CName[i]))
			STATE[i] <- sprintf("\t: %s is calculated via Conservation Law",CName[i])
			DERIVATIVE[i] <- sprintf(fmt$comment, CName[i], Compound$InitialValue[i],ODE[i])
			IVP[i] <- sprintf("\t: %s cannot have initial values as it is determined by conservation law",CName[i])
		}else{
			STATE[i] <- sprintf(fmt$state,CName[i],NeuronUnit(Compound$Unit[i]))
			Right.Hand.Side <- sub("^[[:blank:]]*[+]","",ODE[i]) # remove leading plus signs, if present
			DERIVATIVE[i] <- sprintf(fmt$ode,CName[i], Right.Hand.Side, CName[i])
			IVP[i] <- sprintf("\t %s = %s : initial condition",CName[i],Compound$InitialValue[i])
		}
	}
	Mod[["EXPRESSION"]] <- c("PROCEDURE assign_calculated_values() {",
			                 "\ttime = t : an alias for the time variable, if needed.",
			                 ConservationLaw,
			                 Assignment,
			                 sprintf("\t%s = %s : flux expression %s",row.names(Reaction),Reaction$Flux,row.names(Reaction)),
			                 "}")

	Mod[["STATE"]] <- c("STATE {",STATE,"}")
	Mod[["INITIAL"]] <- c("INITIAL {",IVP,"}")
	Mod[["BREAKPOINT"]] <- c("BREAKPOINT {","\tSOLVE ode METHOD cnexp",
			                 "\tassign_calculated_values() : procedure",
			                 "}")
	Mod[["DERIVATIVE"]] <- c("DERIVATIVE ode {",DERIVATIVE,"}")
	## Output Functions
	Mod[["FUNCTION"]] <- c("PROCEDURE observables_func() {",sprintf(fmt$output,row.names(Output),Output$Formula,row.names(Output)),"}")
	return(Mod)
}
