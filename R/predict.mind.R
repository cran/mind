library(tm)

predict.mind<-function(object,data,...)

{
  if (class(object)!="mind") stop(paste("no applicable method for 'predict.mind' applied to an object of class '",class(object),"'",sep=""))

    y<-unique(object$beta[,1])
    
    x<-unique(object$beta[,2])
    x <- removeNumbers(x)
    x<-unique(x[-1])
    
    z<-do.call(rbind,lapply(object$r_effect[[1]],colnames))[,1]
    
    data[x] <- lapply(data[x] , factor)
    
    X_p<-model.matrix(~.,as.data.frame(data[,x]))
    
    pred<-list()
    
    for (i in y)
    {
      b<-object$beta[object$beta[,1]==i,]                
      
      bb<-colnames(object$sigma_e)[ncol(object$sigma_e)]

      p_fixed1<-list()
      
      if (bb!="broadarea")
      {
      for (v in names(b)[which(!(names(b)%in%c("y","fixed")))])
      {
      k<-which(names(b)[which(!(names(b)%in%c("y","fixed")))]==v)
      
      p_fixed1[[k]]<-as.matrix(X_p[data[,bb]==unique(data[,bb])[k],])%*%as.matrix(as.numeric(b[,v]))
      }
      p_fixed1<-do.call(rbind,p_fixed1)
      }else{
        p_fixed1<-as.matrix(X_p)%*%as.matrix(as.numeric(b$beta_ba_1))
             }
      
      r<-as.data.frame(data[,z])
      colnames(r)<-z
      ncolr<-ncol(r)
      X_r<-list()
      
      for (j in 1:length(object$r_effect[[1]]))
      {
      X_r[[j]]<-as.data.frame(merge(r,as.data.frame(object$r_effect[[1]][j]),by=z[j],all.x=T))
      }
      
      if (j>1){
        X_r<-Reduce('+', X_r)}else{
          X_r<-as.data.frame(X_r[[1]])
        }
      #X_r<-X_r[,(ncol(X_r)-length(y)+1):ncol(X_r)]
      
      pred[[i]]<-data.frame(X_r[,i]+p_fixed1)
    }
    pred<-do.call(cbind,pred)
    colnames(pred)<-y
    
    return(as.data.frame(pred))
    NextMethod("mind")
}

#