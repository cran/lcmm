
predictL.multlcmm <- function(x,newdata,na.action=1,...){

if(missing(newdata)) stop("The argument newdata should be specified")
if(missing(x)) stop("The argument x should be specified")
if (!inherits(x, "multlcmm")) stop("use only with \"multlcmm\" objects")
# ad 2/04/2012 Xnames2
if (!all(x$Xnames2 %in% colnames(newdata))) {
stop(paste(c("newdata should at least include the following covariates: ","\n",x$Xnames2),collapse=" "))}
if (!inherits(newdata, "data.frame")) stop("newdata should be a data.frame object")


if(x$conv==1|x$conv==2) {

if(!(na.action%in%c(1,2)))stop("only 1 for 'na.omit' or 2 for 'na.fail' are required in na.action argument")


### Traitement des donnees manquantes
newdata <- newdata[,x$Xnames2]
newdata <- as.data.frame(newdata)
colnames(newdata) <- x$Xnames2

linesNA <- apply(newdata,2,function(v) which(is.na(v)))
linesNA <- unique(unlist(linesNA))

if(length(linesNA) & na.action==2) stop("newdata contains missing values")
if(length(linesNA) & na.action==1) 
{
 newdata <- as.data.frame(newdata[-linesNA,])
 colnames(newdata) <- x$Xnames2[-1]
}

##pour les facteurs
olddata <- eval(x$call$data)
termes <- x$Xnames[-1]

 #cas où une variable dans le dataset du modèle est un facteur
 for(v in x$Xnames2[-1])
 {                                                          
  if (is.factor(olddata[,v]))
  {
   mod <- levels(olddata[,v])
   if (!(levels(as.factor(newdata[,v])) %in% mod)) stop(paste("invalid level in factor", v))
   newdata[,v] <- factor(newdata[,v], levels=mod)

   for(m in mod)
   {
    termes <- gsub(paste(v,m,sep=""),v,termes)
   }
   
  }
 }

 #cas où on a factor() dans l'appel de la fonction
 dans_appel <- c(all.names(x$call$fixed),all.names(x$call$random),all.names(x$call$mixture),all.names(x$call$classmb))
 ind_factor <- which(dans_appel=="factor")
 if(length(ind_factor))
 {
  nom.factor <- dans_appel[ind_factor+1]
  for (v in nom.factor)
  {
   mod <- levels(as.factor(olddata[,v]))
   if (!all(levels(as.factor(newdata[,v])) %in% mod)) stop(paste("invalid level in factor", v))
   newdata[,v] <- factor(newdata[,v], levels=mod)
   
   factorv <- paste("factor(",v,")",sep="")
   for(m in mod)
   {
    termes <- gsub(paste(factorv,m,sep=""),v,termes,fixed=TRUE)
   }
  }
 }


###matrice avec toutes les var et toutes les interactions
newdata1 <- model.matrix(as.formula(paste("~",paste(termes,collapse="+"))),data=newdata)
#remettre les termes dans le bon ordre
 Xnames <- x$Xnames[-1]
 z <- grep("factor\\(",Xnames) 
 if (length(z))
 {
  Xnames <- gsub("factor\\(","",Xnames)
  Xnames[z] <- gsub("\\)","",Xnames[z])
 }
 newdata1 <- newdata1[,c("(Intercept)",Xnames)]

###calcul des predictions
   
X1 <- NULL
X2 <- NULL
b1 <- NULL
b2 <- NULL

kk<-0
for(k in 1:length(x$idg0)){
if(x$idg0[k]==1){
X1 <- cbind(X1,newdata1[,k])
if (k==1) b1 <- c(b1,0)
if (k>1) {
place <- x$N[1]+kk
b1 <- c(b1,x$best[place+1])
kk <- kk+1
}
}

if(x$idg0[k]==2){
X2 <- cbind(X2,newdata1[,k])
if (k==1){
place1 <- x$N[1]+kk+1
place2 <- x$N[1]+kk+x$ng-1
b2 <- rbind(b2,c(0,x$best[place1:place2]))
kk <- kk+x$ng-1
}
if (k>1){
place1 <- x$N[1]+kk+1
place2 <- x$N[1]+kk+x$ng
b2 <- rbind(b2,x$best[place1:place2])
kk <- kk+x$ng}
}
}

Ypred<-matrix(0,length(newdata1[,1]),x$ng)
for(g in 1:x$ng){
if(length(b1) != 0){
Ypred[,g]<- X1 %*% b1
}
if(length(b2) != 0){
Ypred[,g]<- Ypred[,g] + X2 %*% b2[,g]
}
}

colnames(Ypred) <- paste("pred_class",1:x$ng,sep="")
if (x$ng==1) colnames(Ypred) <- "pred"

res <- Ypred
}
else{
cat("Output can not be produced since the program stopped abnormally.")
res <- NA
}

 res
}


predictL <- function(x,...) UseMethod("predictL")
