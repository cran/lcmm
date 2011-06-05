postprob.lcmm <-
function(x,...){
if (!inherits(x, "lcmm")) stop("use only with \"lcmm\" objects")
if(x$ng==1) stop("use only with  ng > 1")

classif<-NULL
cl.table<-NULL
for (g in 1:x$ng) {
temp<- subset(x$pprob,x$pprob[,2]==g)
temp1<-apply(temp[,3:(x$ng+2)],2,mean)
cl.table<-rbind(cl.table,temp1)
classif<-cbind(classif,c(as.integer(length(temp[,2])),as.double((length(temp[,2])/x$ns*100))))
}

rownames(cl.table)<-paste("class",1:x$ng,sep="")
colnames(cl.table)<-paste("prob",1:x$ng,sep="")
colnames(classif)<-paste("class",1:x$ng,sep="")
rownames(classif)<-c("N","%")

cat(" \n")
cat("Posterior classification:", "\n")
print(classif)
cat(" \n")

cat("Posterior classification table:", "\n")
cat("     --> mean of posterior probabilities in each class", "\n")
print(signif(cl.table,3))
cat(" \n")
}



postprob <- function(x,...) UseMethod("postprob")
