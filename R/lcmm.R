lcmm <- function(fixed,mixture,random,subject,classmb,ng=1,idiag=FALSE,nwg=FALSE,link="linear",intnodes=NULL,epsY=0.5,data,B,convB=0.0001,convL=0.0001,convG=0.0001,maxiter=500,nsim=100,prior)
{

m <- match.call()
if(missing(fixed)) stop("The argument Fixed must be specified in any model")
if(missing(data)){ stop("The argument data should be specified and defined as a data.frame")}


attr.fixed <- attributes(terms(fixed))
depvar <- as.character(attr.fixed$variables[2])
Y0 <- data[,depvar]
minY0 <- min(Y0)
maxY0 <- max(Y0) 

if(all.equal((maxY0-minY0),0) == T){
	stop("All the values of the dependent variable are the same. No estimation can be performed in that case.")
}
if((any(is.na(Y0))==TRUE)){
	stop("The dependent variable should not contain any missing value")
}

if(length(grep("-",unlist(strsplit(link,split="")))) > 2){
	stop("Please check and revise the 'link' argument according to the format given in the help.")
}



################################# pas de separateur "-", uniquement pour les splines


if(all.equal(length(grep("-",unlist(strsplit(link,split="")))),0)==T){

	if (!(link %in% c("linear","beta","thresholds","splines"))){
		stop("The only available link functions in lcmm are 'linear', 'beta', 'splines' and 'thresholds' functions.")
	}else{
		nbzitr0 <- switch(link,"linear"=2,"beta"=2,"splines"=5,"thresholds"=2)
		idlink0 <- switch(link,"linear"=0,"beta"=1,"splines"=2,"thresholds"=3)
		ntrtot0 <- switch(link,"linear"=2,"beta"=4, "splines"= (nbzitr0 + 2), "thresholds"= as.integer(maxY0-minY0)) 
		if(all.equal(link,"splines")==T){
			link <- "splines"
			type <- "equi"	
		}
		
	} 	
	if ((link %in% c("thresholds"))){


		############################## PARTIE A COMPLETER POUR ORDINAL #####################
		if(!(all.equal(minY0,as.integer(minY0))==T) | !(all.equal(maxY0,as.integer(maxY0))==T)|!all(Y0 %in% minY0:maxY0)){
			stop("With the threshold link function, the longitudinal outcome must be discrete")
		}

		IND <- sort(unique(Y0))
		ide0 <- rep(0,as.integer(maxY0-minY0))
		ide0[(IND-minY0+1)] <- 1	
	
	#######################################################################
	}
	zitr <- rep(0,nbzitr0)
	zitr[1] <- minY0
	zitr[nbzitr0] <- maxY0
	
}


################################# Avec un seul separateur "-", uniquement pour les splines
if(all.equal(length(grep("-",unlist(strsplit(link,split="")))),1)==T){
	stop("The number and location of the nodes only apply for the splines link function. For 'thresholds', 'linear' and 'beta' links, no nodes are required. For the splines link function, both the number and the type of location for the nodes should be specified (ex: 5-manual-splines for 5 manual nodes)")
}


if(all.equal(length(grep("-",unlist(strsplit(link,split="")))),2)==T){

	if(any(unlist(strsplit(link,"-")) %in% c("linear","beta","thresholds"))){
		stop("The number and location of the nodes only apply for the 'splines' link function. For 'thresholds', 'linear' and 'beta' links, no nodes are required.")
	}
	
### Verification de l'ordre de replissage du link
	if(!(unlist(strsplit(link,"-"))[3] %in% c("splines"))){
		stop("When defining a link function using splines, the third part of the 'link' argument must include only 'splines' (ex: 5-equi-splines for 5 equidistant nodes)")
	}
	
	if(!(unlist(strsplit(link,"-"))[2] %in% c("equi","manual","quant"))){
		stop("When defining a link function using splines, the second part of the 'link' argument must include only 'equi' 'manual' 'quant' for equidistant manual or quantile nodes")
	}
	
	nbzitr0 <- as.integer(unlist(strsplit(link,"-"))[1])
	ntrtot0 <- nbzitr0 + 2
	idlink0 <- 2 
	type <- unlist(strsplit(link,"-"))[2]	   
	link <- "splines" 
	if((nbzitr0-2) < 0) stop("At least 2 nodes should be specified for the splines link function.")
}


if (all.equal(idlink0,2)==T){

	zitr <- rep(0,nbzitr0)
	zitr[1] <- minY0
	zitr[nbzitr0] <- maxY0
	
	if(all.equal("manual",type)==T){
		if (is.null(intnodes)){
		stop("If 'manual' option is specified for the splines link function, intnodes argument should include the list of interior nodes")
		}else{            
		if(!(all.equal(length(intnodes),(nbzitr0-2))==T)==T) stop("Intnodes does not include the correct number of interior nodes")     
		intnodes <- sort(intnodes)
		if(intnodes[1]<=zitr[1]|intnodes[nbzitr0-2]>=zitr[nbzitr0])stop("Intnodes are not inside the boundaries of the marker")     
		zitr[2:(nbzitr0-1)] <- intnodes[1:(nbzitr0-2)]
		}
	}    
	if(all.equal("quant",type)==T){
		pas <-c(1:(nbzitr0-2))/(nbzitr0-1) 
		zitr[2:(nbzitr0-1)] <- quantile(sort(Y0),probs=pas)
	}        	       
	if(all.equal("equi",type)==T){
		pas=as.double(maxY0-minY0)/as.double(nbzitr0-1)
		for(i in 2:(nbzitr0-1)){
		zitr[i] <- zitr[i-1]+pas
		}
	} 
}      


if (idlink0==1) {
if (epsY<=0) {
epsY <- 0.5
cat("Argument 'epsY' should be a definite positive real. It is changed to the default value of 0.5. \n") 
}
}



link <- as.character(link)
### appel des differents modeles selon la valeur de l'argument link
result <- switch(link
,"linear"=.Contlcmm(fixed=fixed,mixture=mixture,random=random,subject=subject,classmb=classmb,ng=ng,idiag=idiag,nwg=nwg,data=data,B=B,convB=convB,convL=convL,convG=convG,prior=prior,maxiter=maxiter,epsY=epsY,idlink0=idlink0,ntrtot0=ntrtot0,nbzitr0=nbzitr0,zitr=zitr,nsim=nsim,call=m)

,"beta"=.Contlcmm(fixed=fixed,mixture=mixture,random=random,subject=subject,classmb=classmb,ng=ng,idiag=idiag,nwg=nwg,data=data,B=B,convB=convB,convL=convL,convG=convG,prior=prior,maxiter=maxiter,epsY=epsY,idlink0=idlink0,ntrtot0=ntrtot0,nbzitr0=nbzitr0,zitr=zitr,nsim=nsim,call=m)

,"splines"=.Contlcmm(fixed=fixed,mixture=mixture,random=random,subject=subject,classmb=classmb,ng=ng,idiag=idiag,nwg=nwg,data=data,B=B,convB=convB,convL=convL,convG=convG,prior=prior,maxiter=maxiter,epsY=epsY,idlink0=idlink0,ntrtot0=ntrtot0,nbzitr0=nbzitr0,zitr=zitr,nsim=nsim,call=m)
,"thresholds"=.Ordlcmm(fixed,mixture,random,subject,classmb,ng,idiag,nwg,data,B,convB,convL,convG,prior,maxiter,call=m)
)
##### a faire pour ordinal ###########
}




