.GetConservationLaws <- function(N){
	if (requireNamespace("pracma",quietly=TRUE)){
		M <- pracma::null(t(N))
		if (is.null(M)){
			return(NULL)
		} else if (all(dim(M)>1)){
			M <- t(pracma::rref(t(M)))
		} else {
			M <- M/max(M)
		}
		nr=M
		count=0
		f <- c(2,3,5,7)
		while (norm(nr-round(nr),type="F") > 1e-6 && count<length(f)){
			count <- count+1
			message(sprintf("nullspace is not represented by integers. \nTo make the mass conservation more readable, we multiply them by %i and round.",f[count]))
			nr <- nr*f[count]
		}
		Laws <- round(nr)
		n <- dim(Laws)[2]
		if (n>1){
			Laws <- Laws[,n:1] # reverse order of laws to avoid dependency issues between the laws.
		}
	} else {
		Laws=NULL
	}
	print(Laws)
	return(Laws)
}

AppendAmounts <- function(S,Quantity,QuantityName,Separator){
	n <- length(QuantityName)
	aq <- abs(as.numeric(Quantity))
	s <- paste(aq,QuantityName,sep="*",collapse=Separator)
	S <- paste0(S,s)
	return(S)
}

.GetLawText <- function(Laws,CompoundName,InitialValue){
	print(Laws)
	nLaws <- dim(Laws)[2]
	print(nLaws)
	nC <- length(CompoundName)
	I <- 1:nC
	ConLaw <- list()
	ConLaw[["Constant"]] <- vector(length=nLaws)
	ConLaw[["Eliminates"]] <- vector(length=nLaws)
	ConLaw[["Formula"]] <- vector(length=nLaws,mode="character")
	ConLaw[["ConstantName"]] <- vector(length=nLaws,mode="character")
	ConLaw[["Text"]] <- vector(length=nLaws,mode="character")
	a <- vector(length=nC,mode="logical")
	a[] <- TRUE # allowed replacement index...everything that is TRUE may be replaced by a conservation law expression
	for (j in 1:nLaws){
		l <- as.integer(round(Laws[,j]))
		p <- l>0
		n <- l<0
		LawTextP <- AppendAmounts("",l[p],CompoundName[p],"+")
		LawTextN <- AppendAmounts("",l[n],CompoundName[n],"-")
		##message(sprintf("length of LawTextN: %i («%s»)\n",nchar(LawTextN),LawTextN))
		if (nzchar(LawTextN)){
			LawText <- paste(LawTextP,LawTextN,sep="-")
		} else {
			LawText <- LawTextP
		}
		## MaximumInitialValue <- max(InitialValue[(p|n)&a]) # the thing to be replaced has to appear in the law, but not replaced by previous laws.

		u <- which.max(InitialValue[(p|n)&a])
		k <- I[(p|n)&a][u]
		##print(k)
		##message(sprintf("Replacing compound %i («%s») by Conservation Law Expression.\n",k,CompoundName[k]))
		a[k] <- FALSE # this compound may not be replaced in the future
		Const <- as.numeric(l %*% InitialValue)
		ConLaw$ConstantName[j] <- sprintf("%s_ConservedConst",CompoundName[k])
		ConLaw$Text[j] <- paste(ConLaw$ConstantName[j],sub("1[*]","",LawText),sep=" = ")
		m <- (1:nC != k)

		ConLaw$Constant[j] <- Const
		ConLaw$Eliminates[j] <- k

		FormulaP <- AppendAmounts("", l[p&m], CompoundName[p & m],"+")
		FormulaN <- AppendAmounts("", l[n&m], CompoundName[n & m],"-")
		if (nzchar(FormulaN)){
			Formula <- paste(FormulaP,FormulaN,sep="-")
		}else{
			Formula <- FormulaP
		}
		ConLaw$Formula[j] <- gsub("1[*]","",Formula)
		message(LawText)
		message(sprintf("This will comment out compound %i («%s», initial value: %g), Conserved Constant = %f\n",k,CompoundName[k],InitialValue[k],Const))
	}
	return(ConLaw)
}

PrintSteadyStateOutputs <- function(Compound,ODE,Reaction,document.name){
	ss <- Compound$SteadyState
	RN <- row.names(Reaction)
	N <- length(RN)
	if (any(ss)){
		CName <- row.names(Compound)[ss]
		ODE <- ODE[ss]
		## for working SBML, replace all flux names with the actual flux expressions:
		for (j in 1:N){
			ODE <- gsub(sprintf("\\<%s\\>",RN[j]),sprintf("(%s)",Reaction$Flux[j]),ODE)
		}
		## ODE <- gsub("^[+]","",ODE)
		header <- character()
		header[1] <- sprintf("!!SBtabSBtabVersion='1.0'\tTableName='Output' TableTitle='These Outputs describe how well the SteadyState has been achieved' TableType='Quantity' Document='%s'",document.name)
		header[2] <- sprintf("!ID\t!Name\t!Comment\t!ErrorType\t!Unit\t!ProbDist\t!Formula")
		Name <- sprintf("%s_NetFlux",CName)
		SuggestedMeasureOfEquilibrium <- c(header,sprintf("%s\t%s\tmeasures deviation from steady state\t%s\tWeight\tnM\tNormal\t%s",Name,Name,ODE))
		ssfname <- paste0(document.name,"_SteadyStateMetrics.tsv")
		cat(SuggestedMeasureOfEquilibrium,sep="\n",file=ssfname)
	}
}

.GetLogical <- function(Column){
	n <- length(Column)
	LC <- vector(mode = "logical", length = n)
	l10 <- grepl("^1$|^0$",Column)
	lTF <- grepl("^T(RUE)?$|^F(ALSE)?$",toupper(Column))
	LC[lTF] <- as.logical(Column[lTF])
	LC[l10] <- as.logical(as.numeric(Column[l10]))
	return(LC)
}

.GetReactions <- function(SBtab){
	Formula <- SBtab[["Reaction"]][["!ReactionFormula"]]
	Flux <- SBtab[["Reaction"]][["!KineticLaw"]]
	IsReversible <- SBtab[["Reaction"]][["!IsReversible"]]
	Reaction <- data.frame(Formula,Flux,IsReversible,row.names=row.names(SBtab$Reaction))
	return(Reaction)
}

.GetConstants <- function(SBtab){
	if ("Constant" %in% names(SBtab)){
		n <- nrow(SBtab$Constant)
		Value <- sbtab_quantity(SBtab$Constant)
		if ("!Unit" %in% names(SBtab$Constant)){
			Unit <- SBtab$Constant[["!Unit"]]
		} else {
			Unit <- rep("1",n)
		}
		Constant <- data.frame(Value,Unit,row.names=row.names(SBtab$Constant))
		return(Constant)
	} else {
		message("There is no «Constant» Table in this model.")
		return(NULL)
	}
}

## this will at first be for logical vectors, not general yet
.OptionalColumn <- function(SBtab,Name,mode="logical"){
	n <- nrow(SBtab)
	if (Name %in% names(SBtab)){
		Column <- switch(mode,
			             logical=.GetLogical(SBtab[[Name]]),
			             numeric=as.numeric(SBtab[[Name]]),
			             as.character(SBtab[[Name]])
			             )
	} else {
		Column <- vector(mode,length=n)
	}
	return(Column)
}

.GetCompounds <- function(SBtab){
	nComp <- nrow(SBtab[["Compound"]])
	Name <- row.names(SBtab[["Compound"]])
	message("compound names:")
	print(Name)
	## replace possible non-ascii "-"
	CleanIV <- gsub("−","-", SBtab[["Compound"]][["!InitialValue"]])
	InitialValue <- as.numeric(CleanIV);
	SteadyState <- .OptionalColumn(SBtab[["Compound"]],"!SteadyState","logical")
	Unit <- SBtab[["Compound"]][["!Unit"]]
	message("Units: ")
	print(Unit)
	message("---")
	Assignment <- .OptionalColumn(SBtab[["Compound"]],"!Assignment","character")
	IsConstant <- .OptionalColumn(SBtab[["Compound"]],"!IsConstant","logical")
	Interpolation <- .OptionalColumn(SBtab[["Compound"]],"!Interpolation","logical")
	Compound <- data.frame(InitialValue,SteadyState,Unit,IsConstant,Assignment,Interpolation,row.names=Name)
	return(Compound)
}

.GetExpressions <- function(SBtab){
	if ("Expression" %in% names(SBtab)){
		Formula <- SBtab[["Expression"]][["!Formula"]]
		Unit <- SBtab[["Expression"]][["!Unit"]]
		Expression <- data.frame(Formula,Unit,row.names=row.names(SBtab$Expression))
		return(Expression)
	} else {
		message("There is no «Expression» Table in this model.")
		return(NULL)
	}
}

.GetParameters <- function(SBtab){
	nPar <- nrow(SBtab[["Parameter"]]);
	if ("!Scale" %in% names(SBtab[["Parameter"]])){
		Scale <- SBtab[["Parameter"]][["!Scale"]]
	}else{
		Scale <- vector(mode="character",len=nPar)
		Scale[] <- "linear"
	}
	Name <- row.names(SBtab[["Parameter"]])
	Value <- sbtab_quantity(SBtab[["Parameter"]])
	if (any(grepl("log",Scale))) {
		l <- grepl("^(((natural|base-e)?[ ]*log(arithm)?)|ln)$",Scale)
		Value[l] <- exp(Value[l])
		l <- grepl("^(((decadic|base-10)[ ]*logarithm)|log10)$",Scale)
		Value[l] <- 10^Value[l]
	}
	Unit <- SBtab[["Parameter"]][["!Unit"]]
	Parameter <- data.frame(Value=Value,Unit=Unit,row.names=Name)
	return(Parameter)
}

.GetOutputs <- function(SBtab){
	Name <- row.names(SBtab[["Output"]])
	Formula <-  SBtab[["Output"]][["!Formula"]]
	Unit <- SBtab[["Output"]][["!Unit"]]
	Output <- data.frame(Formula,Unit,row.names=Name)
	return(Output)
}

.GetCompartments <- function(SBtab){
	if ("Compartment" %in% names(SBtab)){
		Name <- row.names(SBtab[["Compartment"]])
		Size <-  SBtab[["Compartment"]][["!Size"]]
		Unit <- SBtab[["Compartment"]][["!Unit"]]
	} else {
		Name <- 'Compartment'
		Size <- 1.0
		Unit <- 'liter'
	}
	Comp <- data.frame(Size,Unit,row.names=Name)
	return(Comp)
}


.GetInputs <- function(SBtab){
	if ("Input" %in% names(SBtab)){
		if ("!ConservationLaw" %in% names(SBtab[["Input"]])){
			Disregard <- .GetLogical(SBtab[["Input"]][["!ConservationLaw"]])
			message("Some input parameters may be earlier detected Conservation Law constants: ")
			print(Disregard)
			message("---")
		} else {
			n <- nrow(SBtab[["Input"]])
			Disregard <- vector(mode="logical",length=n)
		}
		Name <- row.names(SBtab$Input)[!Disregard]
		DefaultValue <- sbtab_quantity(SBtab$Input)[!Disregard]
		Unit <- SBtab$Input[["!Unit"]][!Disregard]
		##
		Input <- data.frame(DefaultValue,Unit,row.names=Name)
		return(Input)
	} else {
		message("There is no «Input» Table in this model.")
		return(NULL)
	}
}

.GetDefaults <- function(SBtab){
	if ("Defaults" %in% names(SBtab)){
		Name <- rownames(SBtab$Defaults)
		Unit <-  SBtab$Defaults[["!Unit"]]
	} else {
		Name <- c("time","substance","volume","area","length")
		Unit <- c("millisecond","millimole","liter","nanometer^2","nanometer")
	}
	Defaults <- data.frame(Unit,row.names=Name)
	return(Defaults)
}


UpdateODEandStoichiometry <- function(Term,Compound,FluxName,Expression,Input){
	l <- length(Term)
	J <- vector(mode="integer",len=l)
	C <- vector(mode="integer",len=l)
	##
	## TODO: j and n must be vectors, otherwise only the last matching Compound will be returned
	for (i in 1:l){
		## find possible factors within string
		xb <- unlist(strsplit(trimws(Term[i]),"[* ]"))
		##
		if (length(xb)>1){
			n <- round(as.numeric(xb[1]))
			compound <- make.cnames(xb[2])
		} else {
			compound <- make.cnames(xb[1])
			n <- 1
		}
		cat(sprintf("%i × %s",n,compound))
		if (compound %in% row.names(Compound)){
			j <- as.numeric(match(compound,row.names(Compound)))
			cat(sprintf("\t\t\t(%s is compound %i)\n",compound,j))
		} else if (compound %in% row.names(Expression)){
			j <- (-1)
			cat(sprintf("\t\t\t«%s» is a fixed expression, it has no influx. ODE will be unaffected, but the expression may be used in ReactionFlux calculations\n",compound))
		} else if (compound %in% row.names(Input)){
			j <- (-1)
			cat(sprintf("\t\t\t«%s» is an input parameter (a parameter that represents a constant concentration of a substance outside of the model's scope), it has no influx. ODE will be unaffected, but the expression may be used in ReactionFlux calculations\n",compound))
		} else if (compound %in% c("null","NULL","NIL","NONE","NA","Ø","[]","{}")) {
			cat(sprintf("\t\t\t«%s» (Ø) is a placeholder to formulate degradation in reaction formulae.\n",compound))
			j <- (-2)
		} else {
			stop(sprintf("\t\t\t«%s» is neither in the list of registered compounds nor is it an expression.\n",compound))
		}
		J[i] <- j
		C[i] <- n
	}
	return(list(n=C,compound=compound,J=J))
}

NFlux <- function(n,RName){
	if (n>1){
		NF <- paste(as.character(n),RName,sep="*")
	} else if (n==1) {
		NF <- RName
	} else {
		print(n)
		stop("weird n.")
	}
	return(NF)
}

ParseReactionFormulae <- function(Compound,Reaction,Expression,Input){
	message(class(Reaction$Formula))
	lhs_rhs <- strsplit(as.vector(Reaction$Formula),"<=>")

	nC <- dim.data.frame(Compound)
	nR <- dim.data.frame(Reaction)
	## stoichiometry matrix:
	N <- matrix(0,nrow=nC[1],ncol=nR[1])
	##print(N)
	ODE<-vector(mode="character",length=nC[1])
	## Names
	RName <- row.names(Reaction)
	CName <- row.names(Compound)
	lhs <- vector(mode="character",length=nR[1])
	rhs <- vector(mode="character",length=nR[1])
	##
	for (i in 1:nR[1]){
		line=lhs_rhs[[i]]
		lhs[i]=trimws(line[1])
		rhs[i]=trimws(line[2])
		cat(sprintf("Reaction %i:",i))
		cat(sprintf("line (a->b): «%s» ←→ «%s»\n",line[1],line[2]))
		a <- unlist(strsplit(line[1],"[+]"))
		b <- unlist(strsplit(line[2],"[+]"))
		cat(" where a: ")
		print(a)
		cat("   and b: ")
		print(b)

		## the following two «for» blocks (1,2) operate by adding
		## things to N and ODE. I don't see how to make them into a
		## function without copying ODE and N a lot into that function
		## and back into the caller; Weirdly the <<- operator did not
		## work at all. But perhaps <<- is also bad practice.

		## 1
		cat("Products:\n")
		L <- length(b);
		Term <- UpdateODEandStoichiometry(b,Compound,RName[i],Expression,Input)
		for (k in 1:L){
			j <- Term$J[k]
			if (j>0){
			    ODE[j] <- paste(ODE[j],NFlux(Term$n[k],RName[i]),sep="+")
			    N[j,i] <- N[j,i] + Term$n[k]
			}
		}
		## 2
		cat("Reactants:\n")
		L <- length(a);
		Term <- UpdateODEandStoichiometry(a,Compound,RName[i],Expression,Input)
		for (k in 1:L){
			j <- Term$J[k]
			if (j>0){
			    ODE[j] <- paste(ODE[j],NFlux(Term$n[k],RName[i]),sep="-")
			    N[j,i] <- N[j,i] - Term$n[k]
			}
		}
	}
	message(sprintf("Number of compounds:\t%i\nNumber of Reactions:\t%i",nC[1],nR[1]))
	ModelStructure <- list(ODE=ODE,Stoichiometry=N,LHS=lhs,RHS=rhs)
	return(ModelStructure)
}

PrintConLawInfo <- function(ConLaw,CompoundName,document.name){
	nLaws <- length(ConLaw$Text)
	header<-character(length=2)
	header[1] <- sprintf("!!SBtab\tDocument='%s' TableName='Suggested Input' TableTitle='The model has conservation laws that were automatically determined, these are the conserved constants' TableType='Quantity'",document.name)
	header[2] <- sprintf("!ID\t!Name\t!DefaultValue\t!Unit\t!ConservationLaw\t!Comment")
	SuggestedParameters <- c(header,sprintf("CLU%i\t%s\t%g\tnM\tTRUE\t%s",1:nLaws,ConLaw$ConstantName,ConLaw$Constant,ConLaw$Text))
	infname <- paste0(document.name,"_SuggestedInput.tsv")
	cat(SuggestedParameters,sep="\n",file=infname)

	header[1] <- sprintf("!!SBtab\tDocument='%s' TableName='Suggested Output' TableTitle='Automatically determined conservation laws remove state variables, these outputs make them observable' TableType='Quantity'",document.name)
	header[2] <- sprintf("!ID\t!Name\t!Comment\t!ErrorName\t!ErrorType\t!Unit\t!ProbDist\t!Formula")
	k <- ConLaw[["Eliminates"]]
	SuggestedOutput=c(header,sprintf("YCL%i\t%s_mon\tmonitors implicit state\tSD_YCL%i\tnot applicable\tnM\tnone\t%s",1:nLaws,CompoundName[k],1:nLaws,CompoundName[k]))
	outfname <- paste0(document.name,"_SuggestedOutput.tsv")
	cat(SuggestedOutput,sep="\n",file=outfname)
	message(sprintf("If you'd like to monitor omitted compounds, add this to the Output table: %s\n",outfname))
}

.make.vfgen <- function(H,Constant,Parameter,Input,Expression,Reaction,Compound,Output,ODE,ConLaw){
	vfgen <- list()
	fmt <- list(const=" <Constant Name=\"%s\" Description=\"constant\" Value=\"%s\"/>",
			    par=" <Parameter Name=\"%s\" Description=\"independent parameter\" DefaultValue=\"%g\"/>",
			    input=" <Parameter Name=\"%s\" Description=\"input parameter\" DefaultValue=\"%s\"/>",
			    total=" <Parameter Name=\"%s\" Description=\"conserved quantity eliminates %s as a state variable\" DefaultValue=\"%f\"/>",
			    ConservationLaw=" <Expression Name=\"%s\" Description=\"derived from conservation law %i\" Formula=\"%s\"/>",
			    expression=" <Expression Name=\"%s\" Description=\"defined expression\" Formula=\"%s\"/>",
			    flux=" <Expression Name=\"%s\" Description=\"flux\" Formula=\"%s\"/>",
			    comment="<!-- <StateVariable Name=\"%s\" Description=\"removed compound\" DefaultInitialCondition=\"%s\" Formula=\"%s\"/> -->",
			    ode=" <StateVariable Name=\"%s\" Description=\"compound\" DefaultInitialCondition=\"%s\" Formula=\"%s\"/>",
			    output=" <Function Name=\"%s\" Description=\"output\" Formula=\"%s\"/>")
	vfgen[["header"]] <- "<?xml version=\"1.0\" ?>"
	vfgen[["model"]] <- sprintf("<VectorField Name=\"%s\" Description=\"model created by an R script «sbtab_to_vfgen.R» (https://github.com/a-kramer/SBtabVFGEN)\">",H)
	## Constants
	vfgen[["const"]] <- sprintf(fmt$const,row.names(Constant),Constant$Value)
	## Parameters
	vfgen[["par"]] <- sprintf(fmt$par,row.names(Parameter),Parameter$Value)
	## Inputs
	vfgen[["input"]] <- sprintf(fmt$input,row.names(Input),Input$DefaultValue)
	## Conservation Laws
	if (is.null(ConLaw) || any(is.na(ConLaw))){
		vfgen[["ConservationLaw"]] <- NULL
		vfgen[["ConservationInput"]] <- NULL
		nLaws <- 0
	} else {
		k <- ConLaw$Eliminates
		CName <- row.names(Compound)[k]
		vfgen[["ConservationInput"]] <- sprintf(fmt$total,ConLaw$ConstantName,CName,ConLaw$Constant)
		F <- sprintf("%s - (%s)",ConLaw$ConstantName,ConLaw$Formula)
		nLaws <- length(F)
		vfgen[["ConservationLaw"]] <- sprintf(fmt$ConservationLaw,CName,c(1:nLaws),F)
	}
			                            # Expressions and Reaction Fluxes
	vfgen[["expression"]] <- sprintf(fmt$expression,row.names(Expression),Expression$Formula)
	vfgen[["flux"]] <- sprintf(fmt$flux,row.names(Reaction),Reaction$Flux)
			                            # ODE right-hand-sides
	nC <- dim.data.frame(Compound)
	CName <- row.names(Compound)
	for (i in 1:nC[1]){
		if (nLaws>0 && i %in% ConLaw$Eliminates){
			cat(sprintf("StateVariable %s will be commented out as it was already defined as a Mass Conservation Law Expression.",CName[i]))
			vfgen[["ode"]][i] <- sprintf(fmt$comment,CName[i], Compound$InitialValue[i], ODE[i])
		}else{
			vfgen[["ode"]][i] <- sprintf(fmt$ode,CName[i], Compound$InitialValue[i], ODE[i])
		}
	}
	## Output Functions
	vfgen[["function"]] <- sprintf(fmt$output,row.names(Output),Output$Formula)
	vfgen[["endmodel"]] <- "</VectorField>"
	return(vfgen)
}

.write.txt <- function(H,Constant,Parameter,Input,Expression,Reaction,Compound,Output,ODE,ConLaw){
	files<-c("StateVariables.txt","Parameters.txt","ODE.txt")
	if (!is.null(Constant)) {
		write.table(Constant$Value,row.names=row.names(Constant),col.names=FALSE,sep='\t',file="Constants.txt",quote=FALSE)
		files<-c(files,"Constants.txt")
	}
	write.table(Parameter$Value,row.names=row.names(Parameter),col.names=FALSE,sep='\t',file="Parameters.txt",quote=FALSE)
	if (!is.null(Input)){
		write.table(Input$DefaultValue,row.names=row.names(Input),col.names=FALSE,sep='\t',append=TRUE,file="Parameters.txt",quote=FALSE)
	}
	if (!is.null(Expression)){
		write.table(Expression$Formula,row.names=row.names(Expression),col.names=FALSE,sep='\t',file="Expressions.txt",quote=FALSE)
		files<-c(files,"Expressions.txt")
	}
	##
	N <- dim(Compound)[1]
	if (!is.null(ConLaw) && !any(is.na(ConLaw)) && is.list(ConLaw)) {
		ConLaw <- as.data.frame(ConLaw)
		k <- ConLaw$Eliminates
		CName <- row.names(Compound)[k]
		write.table(ConLaw[,c('Constant')],row.names=row.names(ConLaw),col.names=FALSE,sep='\t',append=TRUE,file="Parameters.txt",quote=FALSE)
		F <- sprintf("(%s - (%s))",ConLaw$ConstantName,ConLaw$Formula)
		names(F) <- CName
		write.table(F,row.names=TRUE,col.names=FALSE,sep='\t',append=TRUE,file="Expressions.txt",quote=FALSE)
		if (!("Expressions.txt" %in% files)){
			files<-c(files,"Expressions.txt")
		}
		i <- (-k)
	} else {
		i <- seq(N)
	}
	CNames<-row.names(Compound)
	RNames<-row.names(Reaction)
	write.table(Compound[i,c('InitialValue','Unit')],row.names=CNames[i],col.names=FALSE,sep='\t',file="StateVariables.txt",quote=FALSE)
	write.table(Reaction[i,'Flux'],row.names=RNames[i],col.names=FALSE,sep='\t',append=TRUE,file="Expressions.txt",quote=FALSE)
	if (!is.null(Output)) {
		write.table(Output[,'Formula'],row.names=row.names(Output),col.names=FALSE,sep='\t',file="OutputFunctions.txt",quote=FALSE)
		files<-c(files,"OutputFunctions.txt")
	}
	ODE<-data.frame(rhs=ODE,row.names=row.names(Compound))
	write.table(ODE[i,],row.names=FALSE,col.names=FALSE,sep='\t',file="ODE.txt",quote=FALSE)
	tar(paste0(H,".tar.gz"),files=files,compression="gzip")
	zip(paste0(H,".zip"),files=files)
}


#' SBtab content exporter
#'
#' This function interprets the SBtab content and exports it to three
#' possible formats: (1) vfgen's `vf` format (xml); (2) NEURON's mod
#' file; (3) SBML if the libSBML package can be loaded.
#'
#' @param SBtab a list as returned from any of the import functions
#'     [sbtab_from_tsv()] or [sbtab_from_ods()]
#' @param cla a logical value that determines whether conservation law
#'     analysis should be performed. If this is TRUE (default) the
#'     resulting ODE model will be reduced to a state variable sub-set
#'     that is independent (in the sense of linear algebraic
#'     relationships).
#' @return if conservation law analysis is turned on, this function
#'     returns the results, otherwise NA; if the model cannot be reduced, this
#'     also returns NULL.
#' @export
sbtab_to_vfgen <- function(SBtab,cla=TRUE){
	options(stringsAsFactors = FALSE)
	## message("The names of the SBtab list:")
	## message(cat(names(SBtab),sep=", "))
	document.name <- comment(SBtab)
	cat(sprintf("Document Name: %s.\n",document.name))
	cat(sprintf("SBtab has %i tables.\n",length(SBtab)))
	cat("The names of SBtab[[1]]:\n")
	cat(colnames(SBtab[[1]]),sep=", ")
	cat("\n")
	Reaction <- .GetReactions(SBtab)
	Constant <- .GetConstants(SBtab)
	Expression <- .GetExpressions(SBtab)
	Compound <- .GetCompounds(SBtab)
	Parameter <- .GetParameters(SBtab)
	Output <- .GetOutputs(SBtab)
	Input <- .GetInputs(SBtab)
	Defaults <- .GetDefaults(SBtab)
	Comp  <- .GetCompartments(SBtab)
	if (requireNamespace("libSBML",quietly=TRUE)){
		## definitely without conservation laws
		SBML <- .make.sbml(document.name,Defaults,Constant,Parameter,Input,Expression,Reaction,Compound,Output,Comp)
		file.xml <- sprintf("%s.xml",document.name)
		libSBML::writeSBML(SBML,file.xml);
		stopifnot(file.exists(file.xml))
		## fix sbml's problem with the time variable
		sbml.text <- readLines(file.xml)
		sbml.text <- gsub('<ci>\\s*t(ime)?\\s*</ci>','<csymbol encoding="text" definitionURL="http://www.sbml.org/sbml/symbols/time"> time </csymbol>',sbml.text)
		cat(sbml.text,sep="\n",file=file.xml)
		message(sprintf("The sbml file has been named «%s».",file.xml))
	}
	## some biological compounds are better represented as expressions/assignments or constants
	## most will be state variables
	if ("IsConstant" %in% names(Compound)){
		IsConstant <- Compound$IsConstant
		message(sprintf("class(IsConstant): %s.\n",class(IsConstant)))
		CC <- Compound[IsConstant,]
		NewExpression <- data.frame(Formula=CC$InitialValue,Unit=CC$Unit,row.names=row.names(CC))
		print(NewExpression)
		Expression <- rbind(Expression,NewExpression)
		print(Expression)
		Compound <- Compound[!IsConstant,]
		print(row.names(Expression))
	}
	message("---")
	if ("Assignment" %in% names(Compound)){
		A <- Compound$Assignment
		U <- Compound$Unit
		l <- vector(mode="logical",len=length(A))
		F <- vector(mode="character",len=length(A))
		for (i in 1:length(A)){
			a <- A[i]
			if (a==""){
				ex<-NA
			}else{
				ex <- charmatch(a,rownames(Expression))
			}
			if (!any(is.na(ex))){
				l[i] <- TRUE
				F[i] <- row.names(Expression[ex,])
				message(sprintf("Compound «%s» is mapped to expression %i «%s» (matched by ID).\n",a,ex,Expression[ex,"Name"]))
			} else if (!any(is.na(Expression[a,"Formula"]))){
				l[i] <- TRUE
				F[i] <- a
				message(sprintf("Compound «%s» is mapped to expression «%s» (matched by Name).\n",a,Expression[a,"Name"]))
			}
		}
		if (any(l)){
			NewExpression <- data.frame(ID=Compound[l,"ID"],Formula=F[l],Unit=U[l],row.names=row.names(Compound[l,]))
			print(NewExpression)
			Expression <- rbind(Expression,NewExpression)
			Compound <- Compound[!l,]
		}
	}
	##
	ModelStructure <- ParseReactionFormulae(Compound,Reaction,Expression,Input)
	ODE <- ModelStructure$ODE
	Reaction[["lhs"]] <- ModelStructure$lhs
	Reaction[["rhs"]] <- ModelStructure$rhs

	if (cla){
		Laws <- .GetConservationLaws(ModelStructure$Stoichiometry)
		if (!is.null(Laws)){
			nLaws <- dim(Laws)[2]
			N <- ModelStructure$Stoichiometry
			message("Stoichiometric Matrix:")
			print(N)
			message("---")
			message(sprintf("Conservation Law dimensions:\t%i × %i\n",dim(Laws)[1],dim(Laws)[2]))
			message(sprintf("To check that the conservation laws apply: norm(t(StoichiometryMatrix) * ConservationLaw == %6.5f)",norm(t(N) %*% Laws,type="F")))
			ConLaw <- .GetLawText(Laws,row.names(Compound),Compound$InitialValue)
			PrintConLawInfo(ConLaw,row.names(Compound),document.name)
			if (require("hdf5r")){
				f5 <- h5file("ConservationLaws.h5",mode="w")
				f5[["ConservationLaws"]] <- t(Laws); # hdf5 will transpose this matrix again
				f5[["/Stoichiometry"]] <- N
				f5[["/Description"]]<-ConLaw$Text
				f5[["/Document"]]<-document.name
				f5[["/Constant"]]<-ConLaw$Constant
				f5[["/ConstantName"]]<-ConLaw$ConstantName
				f5[["/EliminatedCompounds"]]<-(ConLaw$Eliminates-1); # this is intended for C
				h5close(f5)
			} else {
				rownames(Laws)<-row.names(Compound)
				write.table(t(Laws),file="ConservationLaws.tsv",sep="\t",col.names=FALSE,row.names=FALSE)
				colnames(N)<-row.names(Reaction)
				rownames(N)<-row.names(Compound)
				write.table(N,file="Stoichiometry.tsv",sep="\t",col.names=FALSE,row.names=FALSE)
			}
		} else {
			ConLaw <- NULL
		}
	} else {
		Laws <- NULL
		ConLaw <- NA
	}
	##
	PrintSteadyStateOutputs(Compound,ODE,Reaction,document.name)
	H <- make.cnames(document.name)
	.write.txt(H,Constant,Parameter,Input,Expression,Reaction,Compound,Output,ODE,ConLaw)
	vfgen <- .make.vfgen(H,Constant,Parameter,Input,Expression,Reaction,Compound,Output,ODE,ConLaw)
	fname<-sprintf("%s.vf",H)
	cat(unlist(vfgen),sep="\n",file=fname)
	message(sprintf("The vf content was written to: %s\n",fname))

	Mod <- .make.mod(H,Constant,Parameter,Input,Expression,Reaction,Compound,Output,ODE,ConLaw)
	fname<-sprintf("%s.mod",H)
	cat(unlist(Mod),sep="\n",file=fname)
	message(sprintf("The mod content was written to: %s\n",fname))

	return(ConLaw)
}
