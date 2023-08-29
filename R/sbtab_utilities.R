#' Make C compatible names
#'
#' Uses make.names internally, but replaces dots with underscores.
#'
#' @param Labels a character vector with words that need to be turned
#'     into names.
#' @export
#' @return a vector with unique names, can be used as variable names
#'     in C
make.cnames <- function(Labels){
	Names <- gsub("'([^']*)'","lsquo\\1rsquo",trimws(Labels))
	Names <- gsub('"([^"]*)"',"ldquo\\1rdquo",Names)
	Names <- gsub("&","and",Names,fixed=TRUE)
	Names <- gsub("|","or",Names,fixed=TRUE)
	Names <- gsub("'","prime",Names)
	Unique.Names <- gsub("[.]","_",make.names(Names, unique = TRUE, allow_ = TRUE))
	return(Unique.Names)
}

#' Fixed-token-split a scalar string into parts
#'
#' This function splits a string like strsplit, but it removes
#' whitespace from the result (on both ends) and the split token is
#' taken literally
#'
#' @param str a string (character vector of length 1)
#' @param s a split token
#' @param re defaults to FALSE, if TRUE s is treated as a regular
#'     expression
#' @return a character vector of the components without leading or
#'     trailing whitespace
#' @export
#' @examples ftsplit(" A + 2*B ","+")
#' [1] "A"   "2*B"
#' @examples x<-c('a+b','c+d','1 + 2'); lapply(x,ftsplit,'+')
#' [[1]]
#' [1] "a" "b"
#' [[2]]
#' [1] "c" "d"
#' [[3]]
#' [1] "1" "2"
#' @examples this also works, but mixes up the components:
#' x<-c('a+b','c+d+1','1 / 2'); ftsplit(x,'+')
#' [1] "a"	 "b"	 "c"	 "d"	 "1"	 "1 / 2"
ftsplit <- function(str,s=" ",re=FALSE){
	s<-trimws(unlist(strsplit(str,s,fixed=!re)))
	l<-nzchar(s)
	return(s[l])
}

#' Update a named vector with values from a table
#'
#' This function can be used in the circumstance that you already have
#' a default vector and want to update some of its entries from a
#' data.frame with more specific values. The table describes the
#' circumstamces of these more specific values and has columns named
#' like the vector elements. The function returns a matrix with one
#' column per scenario, updated using the table. This is used while
#' parsing SBtab tables.
#'
#' The returned value M has several copies of vector v
#' with some values changed accoring to a table (Qantity Matrix) If
#' the table has columns that name members of v (also named) they will
#' be used. The M will have as many columns as the table has rows.
#'
#' @param v the vector to update
#' @param Table a table with column names partially matching those in
#'     v
#' @return a matrix with various versions of v (columns) one per
#'     setting described in data.frame Table. The names can have a ">"
#'     prefix in the names (see SBtab rules)
#' @export
#' @examples
#' > v<-c(1,2,3)
#' > names(v)<-c('a','b','c')
#'
#' > data<-data.frame(row.names=c('low','med','high'),
#'   b=c(0.5,2.5,5.5),
#'   comment=c('b < 1','close to default','b > 2×default'))
#'
#'        b   comment
#' low  0.5   b < 1
#' med  2.5   close to default
#' high 5.5   b > 2×default
#'
#' > update_from_table(v,data)
#'   low med high
#' a   1   1   1
#' b   2   2   2
#' c   3   3   3
update_from_table <- function(v,Table,prefix=">"){
	if (is.null(v) || is.null(Table)) return(NULL)
	N <- names(v)
	stopifnot(is.data.frame(Table))
	n <- nrow(Table)
	M <- matrix(v,nrow=length(v),ncol=n)
	rownames(M)<-names(v)
	l <- startsWith(names(Table),prefix)
	T<-Table[l]
	names(T)<-sub(prefix,"",names(Table[l]))
	NT <- N[N %in% names(T)]
	M[NT,] <- t(T[NT])
	colnames(M)<-rownames(Table)
	return(M)
}

#' read vector from sbtab Quantity table
#'
#' This function extracts a vector from an sbtab table, and names the
#' elements according to the row names of the table. The value column can be named any of these:
#' Value, DefaultValue, InitialValue, Mean, Median.
#'
#' @param Table the sbtab table imported via sbtab_from_{ods,tsv}
#'     (either)
#' @return a vector with names corresponding to !ID and values taken
#'     from !DefaultValue, !Value, or !InitialValue
#' @keywords sbtab quantity
#' @export
sbtab_quantity <- function(Table){
	colNames <- names(Table)
	l <- grepl("^!((Default|Initial)?Value|Mean|Median)$",colNames)
	if (any(l)){
		v <- as.double(Table %1% l)
		names(v) <- rownames(Table)
	} else {
		stop("Table has no Value column.")
	}
	return(v)
}

#' Read a property from the SBtab header
#'
#' This function retrieves a property from an sbtab header:
#' PropertyName='PropertyValue' (properties are these key=value
#' pairs).
#'
#' @param sbtab.header the first line of an SBtab file
#' @param key the left side of each key=value pair
#' @return value of the key=value pair
sbtab.header.value <- function(sbtab.header,key='Document'){
	m <- unlist(sbtab.header %~% sprintf("%s='([^']+)'",key))
	if (length(m)>0){
		property <- m[2] # so the first experssion in parentheses
	} else {
		warning(sprintf("property «%s» not set in SBtab header: «%s».",key,header))
		property <- NULL
	}
	return(property)
}

#' A parser for .ods files with SBtab document structure
#'
#' This uses the readODS package to read the file.
#' SBtab sheets are themselves tables dedicated to a specific type of
#' model property: Reaction, Compound, Parameter, etc.
#'
#' The SBtab content is not interpreted in any way.
#' @param ods.file a string (file's name)
#' @param verbose if FALSE, nothing is printed on screen (with cat)
#' @return SBtab a list of tables (data.frames), one per ods sheet
#'     SBtab[['TableName']] retrives a data.frame comment(SBtab) is
#'     the name of the document
#'
#' @keywords import
#' @examples
#' model.sbtab<-sbtab_from_ods('model.ods')
#' @export
sbtab_from_ods <- function(ods.file,verbose=TRUE){
	M <- readODS::read.ods(ods.file)
	lM <- length(M)
	SBtab <- list(length=lM)
	header <- paste0(M[[1]][1,],collapse='')
	document.name <- sbtab.header.value(header,'Document')
	table.name <- vector(length=lM)
	for (i in 1:lM){
		header <- paste0(M[[i]][1,],collapse='')
		table.name[i] <- sbtab.header.value(header,'TableName')
		SBtab[[i]] <- M[[i]][-c(1,2),-1]
		names(SBtab[[i]]) <- M[[i]][2,]
		row.names(SBtab[[i]]) <- as.character(M[[i]][-c(1,2),1])
	}
	names(SBtab) <- table.name
	if (verbose) cat("Tables found: ",paste(table.name,collapse=', '),"\n")
	comment(SBtab) <- document.name
	return(SBtab)
}

#' A parser for a bunch of .tsv files with SBtab document content
#'
#' This function reads the files (not using any package).
#' SBtab sheets are themselves tables dedicated to a specific type of
#' model property: Reaction, Compound, Parameter, etc.
#'
#' The SBtab content is not interpreted in any way.
#' @param tsv.file a character vector (file names, one per sheet),
#'     defaults to all tsv files in the current directory.
#' @param verbose if FALSE, nothing is printed with cat()
#' @return SBtab a list of tables, one per file in tsv.file list
#'     SBtab[['Reaction']] retrieves the table of reactions, a
#'     data.frame comment(SBtab) retrieves the SBtab document name
#' @keywords import
#' @examples
#' model.sbtab<-sbtab_from_tsv(dir(pattern='.*[.]tsv$'))
#' @export
sbtab_from_tsv <- function(tsv.file=dir(pattern='[.]tsv$'),verbose=TRUE){
	SBtab <- list()
	header <- readLines(tsv.file[1],n=1)
	document.name <- sbtab.header.value(header,"Document")
	if (verbose) printf("[tsv] file[1] «%s» belongs to Document «%s»\n\tI'll take this as the Model Name.\n",tsv.file[1],document.name)
	for (f in tsv.file){
		header <- readLines(f,n=1)
		TableName <- sbtab.header.value(header,'TableName')
		SBtab[[TableName]] <- read.delim(f,as.is=TRUE,skip=1,check.names=FALSE,comment.char="%",blank.lines.skip=TRUE,row.names=1)
		attr(SBtab[[TableName]],"TableName") <- TableName
	}
	comment(SBtab) <- document.name
	return(SBtab)
}

#' collect all variable names
#'
#' this function collects all varibale names in SBtab, here we use the
#' !ID column (because the !Name column can be something that isn't a
#' programming language kind of variable-name).
#'
#' If the react flag is FALSE then really all names are collected,
#' otherwise only those that can appear in reaction formulae and
#' kinetic laws of reactions.
#'
#' @param tab a list of lists with the table content
#' @param react if TRUE, this function only collects the IDs of
#'     Compouds, Parameters and Expressions
#' @return a vector of all names
all.vars <- function(tab,react=TRUE){
	allvars <- c(row.names(tab$Parameter),row.names(tab$Compound))
	tNames <- names(tab)
	if (react){
		l <- tNames %in% c("Input","Constant","Expression")
		tNames <- tNames[l]
	}
	for (T in tNames){
		allvars <- c(allvars,row.names(tab[[T]]))
	}
	return(trimws(unlist(allvars)))
}

#' Checks that all SBtab references are valid
#'
#' A reference is a column name that starts with ">": >Calcium is a
#' reference to a compound (probably) called Calcium (or a Parameter,
#' etc.)
#'
#' @param tab a list of tables
#' @param allvars a list of all variable IDs
all.refs.valid <- function(tab,allvars=all.vars(tab,reac=FALSE)){
	r <- TRUE
	for (T in tab){
		cNames <- colnames(T)
		l <- grepl("^[~>][a-zA-Z]\\w*$",cNames)
		refs <- cNames[l] %s% "[~>]"
		valid.refs <- refs %in% allvars
		if (!all(valid.refs)){
			message(sprintf("Table %s contains invalid references:\n",attr(T,"TableName")))
			message(sprintf("%30s\n",refs[!valid.refs]))
			r <- FALSE
			warning("Not all refs are valid.")
		}
	}
	return(r)
}

#' Are the IDs the same as the Names
#'
#' This is not mandatory, but there could be errors if somewhere in
#' the document a Name is used inside of mathematical expressions
#' while the ID was needed. This can help find those mistakes.
#'
#' @param T one table (data.frame)
id.eq.name <- function(T){
	if ("!Name" %in% names(T)){
		ID <- row.names(T)
		Name <- T[["!Name"]]
		l <- ID == Name
		if (any(!l)){
			message("these IDs and Names are different:")
			message(sprintf("%30s    %s\n","ID","Name"))
			message(sprintf("%30s != %s\n",ID[!l],Name[!l]))
		}
	} else {
		l <- NA
	}
	return(l)
}

#' validate spelling in SBtab document
#'
#' This function checks the SBtab document for consistency.
#' It finds misspelled varibale names in the reaction kinetics.
#'
#' @param tab a list of lists as returned by `sbtab_from_tsv()`
#' @export
sbtab.valid <- function(tab){
	stopifnot("Reaction" %in% names(tab))
	vars <- ftsplit(gsub("[-()+*/]+"," ",tab$Reaction[["!KineticLaw"]]),"[ ]+",re=TRUE)
	av <- all.vars(tab)
	l <- vars %in% av
	if (any(!l)) {
		warning("These variables appear in the reaction table, but are not defined in the rest of the document")
		message(sprintf("%30s\n",vars[!l]))
		return(FALSE)
	}
	## point out where ID and Name are different, maybe that is a problem
	for (T in tab){
		l <- id.eq.name(T)
		if (!any(is.na(l)) && !all(l)){
			warning(sprintf("Not all IDs are equal to the Name attribute in «%s», but maybe this is on purpose.\n",attr(T,"TableName")))
		}
	}
	if (all.refs.valid(tab)){
		message("All internal references are valid (~REF and >REF)\n")
	} else {
		return(FALSE)
	}
	return(TRUE)
}

#' Create a Time Series Simulation Experiment
#'
#' Given the constituents of a time series simulation experiment,
#' return a list that rgsl will understand as a simulation.
#'
#' The output of a real experiment can be a function of the state
#' variables, not necessarily the state variable trajectories
#' themselves. We assume that parameter estimation of some sort will
#' happen.
#'
#' If an estimate of the measurement error is not available, then
#' --strictly speaking-- the data is useless: the noise must be
#' understood as unbounded. But, instead, we assume that the case is
#' very complex and the user knows what to do about it. We
#' automatically set the noise to 5% (relative to each value) and
#' another 5% of the largest value as an estimate of scale, this
#' represents the absolute error:
#'
#' tl;dr      default_error = 5% REL + 5% MAXVAL.
#'
#' The user should replace these values with something, or
#' simply not use them, if the application goes beyond testing.
#'
#' @export
#' @param outputValues (data.frame) measured data, to be replicated by the
#'     simulation
#' @param outputTimes (vector) time values at which the outputValues were
#'     measured
#' @param errorValues (data.frame) an estimate of the measurement noise, if
#'     available
#' @param inputParameters (vector) a parameter vector that the model needs to
#'     operate (can be some default value to be changed later)
#' @param initialTime t0 for the time series experiment: y(t0) = y0
#' @param initialState (vector) initial values of the state variables, at t0
#' @return list with these quantities as list items
time.series <- function(outputValues,outputTimes=as.double(1:dim(outputValues)[2]),errorValues=0.05*outputValues+0.05*max(outputValues),inputParameters,initialTime=as.double(min(outputTimes)),initialState,events=NA){
	outNames <- names(outputValues)
	names(outputValues) <- outNames %s% ">"
	if (is.null(inputParameters))
		inputParameters <- c()
	if (is.null(events) || all(is.na(events))){
		experiment <- list(
			outputValues=outputValues,
			errorValues=errorValues,
			input=inputParameters,
			initialTime=as.double(initialTime),
			initialState=initialState,
			outputTimes=outputTimes
		)
	} else {
		experiment <- list(
			outputValues=outputValues,
			errorValues=errorValues,
			input=inputParameters,
			initialTime=as.double(initialTime),
			initialState=initialState,
			outputTimes=as.double(outputTimes),
			events=events
		)
	}
	return(experiment)
}

sbtab.events <- function(ename,tab){
	if (is.na(ename)) {
		return(NULL)
	}
	if ("Input" %in% names(tab)){
		ni <- nrow(tab$Input)
	} else {
		ni <- 0
	}
	if (nzchar(ename) && ename %in% names(tab)){
		n.sv <- nrow(tab$Compound)
		n.par <- nrow(tab$Parameter) + ni
		n.t <- nrow(tab[[ename]])

		tf <- list(
			state=list(
				A=array(1.0,dim=c(n.sv,1,n.t)),
				b=array(0.0,dim=c(n.sv,1,n.t))
			),
			param=list(
				A=array(1.0,dim=c(n.par,1,n.t)),
				b=array(0.0,dim=c(n.par,1,n.t))
			)
		)
		rownames(tf$state$b) <- row.names(tab$Compound)
		rownames(tf$state$A) <- row.names(tab$Compound)
		par.names <- c(row.names(tab$Parameter),row.names(tab$Input))
		rownames(tf$param$b) <- par.names
		rownames(tf$param$A) <- par.names

		colNames <- names(tab[[ename]])
		n <- nrow(tab[[ename]])

		event.time <- as.double(tab[[ename]][["!Time"]])

		for (i in grep("^>[A-Z]{3}:.+$",colNames)){
			m <- colNames[i] %~% ">([A-Za-z]{3}):(.+)$"
			op <- tolower(m[[1]][2])
			quantity <- m[[1]][3]
			kind <- ifelse(quantity %in% row.names(tab$Compound),"state","param")
			if (op == "add") {
				tf[[kind]]$b[quantity,1,] <- as.numeric(tab[[ename]][[i]])
			} else if (op == "sub") {
				tf[[kind]]$b[quantity,1,] <- (-1.0)*as.numeric(tab[[ename]][[i]])
			} else if (op == "mul") {
				tf[[kind]]$A[quantity,1,] <- as.numeric(tab[[ename]][[i]])
			} else if (op == "div") {
				tf[[kind]]$A[quantity,1,] <- 1.0/as.numeric(tab[[ename]][[i]])
			} else if (op == "set") {
				tf[[kind]]$A[quantity,1,] <- 0.0
				tf[[kind]]$b[quantity,1,] <- as.numeric(tab[[ename]][[i]])
			} else {
				stop(sprintf("unknown operation in event table %s: %s\n",ename,op))
			}
		}
	} else {
		return(NULL)
	}
	return(list(time=event.time,tf=tf))
}

#' Read data from SBtab
#'
#' This function reads all datasets from SBtab Documents that are
#' referenced in the Experiment table and returns them as a list of
#' data.frames
#'
#' The returned list always describes simulation experiments:
#' time series experiments, possibly with scheduled events.
#' @param tab a list of data.frames with SBtab content
#' @export
sbtab.data <- function(tab){
	l <- grepl("([tT]able([oO]f)?)?[eE]xperiments?|[mM]easurements?|[sS]imulations?",names(tab))
	if (any(l) && sum(l)==1) {
		E <- tab %1% l
	} else {
		warning("There must be exactly one Table of Experiments.")
		print(names(tab)[l])
		return(NA)
	}
	out.id <- row.names(tab$Output)
	if ("Input" %in% names(tab)){
		input.id <- row.names(tab$Input)
	} else {
		input.id <- NULL
	}
	n <- dim(E)[1]
	experiments <- list()
	l <- grepl("!([eE]xperiment)?[tT]ype",names(E))
	if (any(l)){
		type <- E %1% l
		time.series <- grepl("[Tt]ime[_ ][Ss]eries",type,ignore.case=TRUE)
		dose.response <- grepl("[Dd]ose[_ ][Rr]esponse",type,ignore.case=TRUE)
		dose.sequence <- grepl("[Dd]ose[ _][Ss]equence",type,ignore.case=TRUE)
	} else {
		time.series <- rep(TRUE,n) # by default
		dose.response <- !time.series
		dose.sequence <- !time.series
	}

	t0 <- E %1% grepl("^!([tT]0|[Ii]nitialTime)$",names(E))
	v <- sbtab_quantity(tab$Compound)
	initVal <- update_from_table(v,E)
	if ("Input" %in% names(tab)){
		v <- sbtab_quantity(tab$Input)
		input <- update_from_table(v,E)
	} else {
		input <- NULL
	}
	id <- row.names(E)

	if ("!Event" %in% names(E)){
		event.names <- E[["!Event"]]
	} else {
		event.names <- rep(NA,n)
	}

	for (i in 1:n){
		stopifnot(id[i] %in% names(tab))
		tNames <- names(tab[[id[i]]])
		l <- grepl(">([a-zA-Z][^ ]*)",tNames)
		m <- dim(tab[[id[i]]])[1]

		v.out <- rep(NA,length(out.id))
		names(v.out) <- out.id

		if (is.character(event.names[i])){
			events <- sbtab.events(event.names[i],tab)
		} else {
			events <- NA
		}

		if (time.series[i]){
			oTime <- as.double(tab[[id[i]]][["!Time"]])
			initialTime <- ifelse(is.null(t0),min(oTime),t0[i])
			OUT <- as.data.frame(t(update_from_table(v.out,tab[[id[i]]],prefix=">")))
			ERR <- as.data.frame(t(update_from_table(v.out,tab[[id[i]]],prefix="~")))
			experiments[[id[i]]] <- time.series(
				outputValues=OUT,
				errorValues=ERR,
				inputParameters=input[,i],
				initialTime=as.double(initialTime),
				initialState=initVal[,i],
				outputTimes=oTime,
				events=events
			)
		} else if (dose.response[i]){
			if (is.null(input)){
				u <- NULL
			} else {
				u <- update_from_table(input[,i],tab[[id[i]]])
			}
			OUT <- as.data.frame(t(update_from_table(v.out,tab[[id[i]]],prefix=">")))
			ERR <- as.data.frame(t(update_from_table(v.out,tab[[id[i]]],prefix="~")))
			outTime <- E[["!Time"]]
			initialTime <- ifelse(is.null(t0),0.0,t0[i])
			dose.id <- row.names(tab[[id[i]]])
			for (j in 1:m){
				experiments[[dose.id[j]]] <- time.series(
					outputValues=OUT[j,],
					errorValues=ERR[j,],
					inputParameters=u[,j],
					initialTime=initialTime,
					initialState=initVal[,i],
					outputTimes=outTime[i],
					events=events
				)
			}
		}
	}
	return(experiments)
}
