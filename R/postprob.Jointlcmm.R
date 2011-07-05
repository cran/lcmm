postprob.Jointlcmm <-
function(x,...){
if (!inherits(x, "Jointlcmm")) stop("use only with \"Jointlcmm\" objects")

if(x$ng==1) { 
cat("Postprob function can only be used when ng > 1 \n")
} else {
if(x$conv==1|x$conv==2) {

classif<-NULL
classifY <- NULL
cl.table<-NULL
for (g in 1:x$ng) {
temp<- subset(x$pprob,x$pprob[,2]==g)
tempY<- subset(x$pprobY,x$pprobY[,2]==g)
temp1<-apply(temp[,3:(x$ng+2)],2,mean)
temp1Y<-apply(tempY[,3:(x$ng+2)],2,mean)
cl.table<-rbind(cl.table,temp1)
classif<-cbind(classif,length(temp[,2]))
classifY<-cbind(classifY,length(tempY[,2]))
}

classif <- rbind(classif,classif/x$ns)
classifY <- rbind(classifY,classifY/x$ns)

rownames(cl.table)<-paste("class",1:x$ng,sep="")
colnames(cl.table)<-paste("prob",1:x$ng,sep="")
colnames(classif)<-paste("class",1:x$ng,sep="")
rownames(classif)<-c("N","%")
colnames(classifY)<-paste("class",1:x$ng,sep="")
rownames(classifY)<-c("N","%")

if(sum(is.na(x$pprob[,3:(x$ng+2)]))==0){
cat(" \n")
cat("Posterior classification based on longitudinal and time-to-event data:", "\n")
print(classif)
cat(" \n")

cat("Posterior classification table:", "\n")
cat("     --> mean of posterior probabilities in each class", "\n")
print(signif(cl.table,3))
cat(" \n")
}else{
cat("Error in the computation of posterior class-membership probabilities given all the information")
}


if(sum(is.na(x$pprobY[,3:(x$ng+2)]))==0){
cat(" \n")
cat("Posterior classification based only on longitudinal data:", "\n")
print(classifY)
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
