postprob.Jointlcmm <-
function(x,...){
if (!inherits(x, "Jointlcmm")) stop("use only with \"Jointlcmm\" objects")

if(x$specif[[3]]==1) { 
cat("Postprob function can only be used when ng > 1 \n")
} else {
if(x$conv==1|x$conv==2) {

classif<-NULL
classifY <- NULL
cl.table<-NULL
pprob <- x$pprob[,-1]
pprobY <- x$pprobY[,-1]
for (g in 1:x$specif[[3]]) {
temp<- subset(pprob,pprob[,1]==g)
tempY<- subset(pprobY,pprobY[,1]==g)
temp1<-apply(temp[,2:(x$specif[[3]]+1)],2,mean)
temp1Y<-apply(tempY[,2:(x$specif[[3]]+1)],2,mean)
cl.table<-rbind(cl.table,temp1)
classif<-cbind(classif,length(temp[,1]))
classifY<-cbind(classifY,length(tempY[,1]))
}

classif <- rbind(classif,100*classif/x$specif[[2]])
classifY <- rbind(classifY,100*classifY/x$specif[[2]])

rownames(cl.table)<-paste("class",1:x$specif[[3]],sep="")
colnames(cl.table)<-paste("prob",1:x$specif[[3]],sep="")
colnames(classif)<-paste("class",1:x$specif[[3]],sep="")
rownames(classif)<-c("N","%")
colnames(classifY)<-paste("class",1:x$specif[[3]],sep="")
rownames(classifY)<-c("N","%")

if(sum(is.na(pprob[,2:(x$specif[[3]]+1)]))==0){
cat(" \n")
cat("Posterior classification based on longitudinal and time-to-event data:", "\n")
print(round(classif,2))
cat(" \n")

cat("Posterior classification table:", "\n")
cat("     --> mean of posterior probabilities in each class", "\n")
print(round(cl.table,4))
#print(cl.table)
cat(" \n")
}else{
cat("Error in the computation of posterior class-membership probabilities given all the information")
}


if(sum(is.na(pprobY[,2:(x$specif[[3]]+1)]))==0){
cat(" \n")
cat("Posterior classification based only on longitudinal data:", "\n")
print(round(classifY,2))
cat(" \n")
}else{
cat("Error in the computation of posterior class-membership probabilities given the repeated measures of the marker")
}
}else{
cat("Output can not be produced since the program stopped abnormally.")
}
}
}







postprob <- function(x,...) UseMethod("postprob")
