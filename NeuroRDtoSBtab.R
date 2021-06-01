make.cnames <- function(Labels){
    ## this makes names that work in C, where dots have special meanings.
    Unique.Names <- gsub("[.]","_",make.names(trimws(Labels), unique = TRUE, allow_ = TRUE),useBytes=TRUE)
    return(Unique.Names)
}


attr <- function(name='name',text){
    r <- sprintf('%s\\s*=\\s*"([^"]+)"',name)
    M <- regmatches(text,regexec(pattern=r,text))
    a <- unlist(lapply(M,function(m) m[2]))
    return(a)
}

tag.attr <- function(tag,name,text){
    r <- sprintf("<%s[ ]",tag)
    l <- grepl(pattern=r,text)
    return(attr(name,text[l]))
}


tagged.value <- function(tag='.*Rate',text){
    r <- sprintf("<%s>\\s*([0-9e.+-]+)\\s*</%s>",tag,tag)
    M <- regmatches(text,regexec(pattern=r,text))
    a <- as.numeric(unlist(lapply(M,function(m) m[2])))
    return(a)
}


NeuroRD <- function(file){
    xml <- readLines(file);
    ## Reactions:
    i <- grep(pattern='<Reaction[ >]',xml)
    j <- grep(pattern='</Reaction>',xml)
    #print(i)
    #print(j)
    L <- length(i)
    stopifnot(L==length(j))
    ReactionName <- vector(mode='character',length=L)
    kf <- numeric(L)
    kr <- numeric(L)
    Q10 <- numeric(L)
    ReactantID <- list(length=L)
    ProductID <- list(length=L)
    for (n in 1:L){
        r <- seq(i[n],j[n])
        R <- xml[r]
        R.flat <- paste0(R,collapse=' ')
        message(R.flat)
        ReactionName[n]  <- attr('id',R.flat)
        message(sprintf("id: «%s»",ReactionName[n]))
        ReactantID[[n]] <- tag.attr('Reactant','specieID',R)
        ProductID[[n]] <- tag.attr('Product','specieID',R)
        kf[n] <- tagged.value('forwardRate',R.flat) 
        kr[n] <- tagged.value('reverseRate',R.flat)
        Q10[n] <- tagged.value('Q10',R.flat)        
    }
    ReactionName <- make.cnames(ReactionName)
    names(kf) <- ReactionName
    names(kr) <- ReactionName
    names(Q10) <- ReactionName
    Reaction <- list(Name=ReactionName,Reactant=ReactantID,Product=ProductID,kf=kf,kr=kr,Q10=Q10)
    
    ##names(Reaction) <- ReactionName
    ## Species:
    i <- grep(pattern='<Spec',xml)
    L <- length(i)
    SpeciesName <- make.cnames(attr('id',xml[i]))
    kdiff <- as.numeric(attr('kdiff',xml[i]))
    kdiffunit <- attr('kdiffunit',xml[i])
    ##species <- data.frame(row.names=SpeciesName,kdiff=kdiff,kdiffunit=kdiffunit)

    ## Initial Values:
    i <- grep(pattern='<InitialConditions',xml)
    j <- grep(pattern='</InitialConditions',xml)
    r <- seq(i,j)
    IC <- xml[r]
    id <- make.cnames(tag.attr('NanoMolarity','specieID',IC))
    val <- as.numeric(tag.attr('NanoMolarity','value',IC))
    names(val) <- id
    species <- data.frame(row.names=SpeciesName,kdiff=kdiff,kdiffunit=kdiffunit,value=val[SpeciesName])
    return(list(Reaction=Reaction,Species=species))
}

.write.compound.table <- function(NeuroRD,DocumentName='MODEL'){
    species <- NRD$Species
    Header <- sprintf("!!SBtab Document='%s' TableName='Compound' TableType='Compound' TableTitle='Species in this Model'",DocumentName)
    header <- sprintf("!ID\t!Name\t!InitialValue\t!Unit")
    s <- paste(row.names(species),row.names(species),species$value,'nanomole/liter',sep='\t')
    cat(c(Header,header,s),file='Compound.tsv',sep="\n")
}

.write.reaction.table <- function(NeuroRD,K,DocumentName='MODEL'){
    reactions <- NRD$Reaction
    Header <- sprintf("!!SBtab Document='%s' TableName='Reaction' TableType='Reaction' TableTitle='Reactions in this Model'",DocumentName)
    header <- sprintf("!ID\t!Name\t!KineticLaw\t!ReactionFormula\t!IsReversible")
    n <- length(reactions$Name)
    Kinetic <- vector(mode='character',length=n)
    Formula <- vector(mode='character',length=n)
    IsReversible <- vector(mode='logical',length=n)
    for (i in 1:n){
        kf <- reactions$kf[[i]]
        kr <- reactions$kr[[i]]
        kfid <- sprintf('kf_%s',reactions$Name[[i]])
        krid <- sprintf('kr_%s',reactions$Name[[i]])
        L <- paste(reactions$Reactant[[i]],collapse='+')
        R <- paste(reactions$Product[[i]],collapse='+')
        F <- paste(L,R,sep='<=>')
        L <- paste(kfid,paste(reactions$Reactant[[i]],collapse='*'),sep='*')
        R <- paste(krid,paste(reactions$Product[[i]],collapse='*'),sep='*')
        IsReversible[i] <- (!is.na(kr) && kr > 0.0)
        if (IsReversible[i]){
            K <- paste(L,R,sep='-')
            IsReversible[i] <- TRUE
            print(K)
        } else {
            K <- L            
        }
        Formula[i] <- F
        Kinetic[i] <- K
    }
    r <- paste(reactions$Name,reactions$Name,Kinetic,Formula,IsReversible,sep='\t')
    cat(c(Header,header,r),file='Reaction.tsv',sep="\n")
    
}

.write.parameter.table <- function(NeuroRD,DocumentName="MODEL"){

    kf <- NRD$Reaction[['kf']]
    nf <- unlist(lapply(NRD$Reaction[['Reactant']],length))
    names(nf) <- names(kf)
    unitf <- trimws(gsub(pattern='nanomole\\^0 liter\\^0',
                         replacement='',
                         sprintf('nanomole^%i liter^%i second^%i',1-nf,nf-1,-1)
                         )
                    ) # i+nf = 1
    kr <- NRD$Reaction[['kr']]
    nr <- unlist(lapply(NRD$Reaction[['Product']],length))
    names(nr) <- names(kr)
    unitr <- trimws(gsub(pattern='nanomole\\^0 liter\\^0',
                         replacement='',
                         sprintf('nanomole^%i liter^%i second^%i',1-nr,nr-1,-1)
                         )
                    ) # i+nr = 1
    Header <- sprintf("!!SBtab Document='%s' TableName='Parameter' TableType='Quantity' TableTitle='reaction rate coefficients'",DocumentName)
    header <- sprintf("!ID\t!Name\t!Value\t!Unit")
    kfid <- sprintf('kf_%s',names(kf))
    krid <- sprintf('kr_%s',names(kr))
    sf <- paste(kfid,kfid,kf,unitf,sep='\t')
    sr <- paste(krid,krid,kr,unitr,sep='\t')
    cat(c(Header,header,sf,sr),file='Parameter.tsv',sep="\n",append=FALSE)
    return(data.frame(kf=kf,kr=kr,nf=nf,nr=nr,unitf=unitf,unitr=unitr,row.names=names(kf)))
}

.write.defaults <- function(NRD,DocumentName='MODEL'){
    Header <- sprintf("!!SBtab Document='%s' TableName='Defaults' TableType='Quantity' TableTitle='default units'",DocumentName)
    header <- sprintf("!ID\t!Name\t!Unit")
    u<-c(paste('time','time','millisecond',sep='\t'),
         paste('volume','volume','liter',sep='\t'),
         paste('area','area','um^2',sep='\t'),
         paste('length','length','um',sep='\t'),
         paste('substance','substance','nanomole',sep='\t')
         )
    cat(c(Header,header,u),file='Defaults.tsv',sep="\n",append=FALSE)
    
}

NeuroRD_to_SBtab <- function(file){
    NRD <- NeuroRD(file)
    .write.compound.table(NRD)
    K <- .write.parameter.table(NRD)
    .write.reaction.table(NRD,K)
    .write.defaults(NRD)
}
