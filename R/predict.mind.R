
predict.mind<-function(object,data,type ="proj",dir_s=NULL,dir_cov=NULL,...)

{
  if (class(object)!="mind") stop(paste("no applicable method for 'predict.mind' applied to an object of class '",class(object),"'",sep=""))
  if (!(type%in%c("eblup","proj","synth"))) stop("'arg'should be one of 'eblup', 'proj', 'synth'")
  if (type=="eblup" & is.null(dir)) stop("when 'eblup' is selected 'dir' must be provided. See 'Details' for more info")
      
  if (type=="synth"){
        y<-unique(object$beta[,1])
        
        x<-unique(object$beta[,2])
        x <- removeNumbers(x)
        x<-unique(x[-1])
        
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
                pred[[i]]<-data.frame(p_fixed1)
        }
        pred<-do.call(cbind,pred)
        colnames(pred)<-y
      }
  if (type=="proj"){
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
      XX_r<-list()
      for (k in 1:length(object$r_effect))
           {
      for (j in 1:length(object$r_effect[[1]]))
      {
      X_r[[j]]<-as.data.frame(merge(r,as.data.frame(object$r_effect[[k]][j]),by=z[j],all.x=T))
      }
      
      if (j>1){
        X_r<-Reduce('+', X_r)}else{
          X_r<-as.data.frame(X_r[[1]])
        }
      X_r[is.na(X_r)]<-0
      XX_r[[k]]<-X_r
      }
      XX_r<- Reduce('+', XX_r)
      
      pred[[i]]<-data.frame(XX_r[,i]+p_fixed1)
    }
    pred<-do.call(cbind,pred)
    colnames(pred)<-y
  }
  if (type=="eblup"){
        dom<-colnames(object$EBLUP)[1]
        
        y<-unique(object$beta[,1])
        
        x<-unique(object$beta[,2])
        x <- removeNumbers(x)
        x<-unique(x[-1])
        
        z<-do.call(rbind,lapply(object$r_effect[[1]],colnames))[,1]
        
        data[x] <- lapply(data[x] , factor)
        
        bb<-colnames(object$sigma_e)[ncol(object$sigma_e)]
        if (bb!="broadarea") ba<-bb else ba<-NULL
        
        tt_names<-names(data[!(names(data)%in%c(dom,x,z,ba))])
        appo <- aggregate(as.formula(paste(tt_names,"~",dom,"+",ba,"+",paste(x,collapse = "+"),"+",paste(z,collapse ="+"),sep="")),data=data,sum)
        appo1<-merge(appo,dir_cov,all.x=T)
        appo1[is.na(appo1)]<-0
        appo1$PD<-appo1[,tt_names]-appo1[,names(dir_cov[!(names(dir_cov)%in%c(dom,x,z))])]
        appo2 <- aggregate(as.formula(paste("PD","~",ba,"+",paste(x,collapse = "+"),sep="")),data=appo1,sum)
        
        X_p<-model.matrix(~.,as.data.frame(appo2[,x]))*appo2$PD
        pred<-list()
        

                for (i in y)
        {
          b<-object$beta[object$beta[,1]==i,]                
          
          p_fixed1<-list()
          appo2b<-list()
          
          if (bb!="broadarea")
          {
            for (v in names(b)[which(!(names(b)%in%c("y","fixed")))])
            {
              k<-which(names(b)[which(!(names(b)%in%c("y","fixed")))]==v)
              
              p_fixed1[[k]]<-as.matrix(X_p[appo2[,bb]==unique(appo2[,bb])[k],])%*%as.matrix(as.numeric(b[,v]))
              appo2b[[k]]<-appo2[appo2[,bb]==unique(appo2[,bb])[k],]
              }
            p_fixed1<-do.call(rbind,p_fixed1)
            appo2b<-do.call(rbind,appo2b)
          }else{
            appo2b<-NULL
            appo2b<-appo2
            p_fixed1<-as.matrix(X_p)%*%as.matrix(as.numeric(b$beta_ba_1))
          }

          appo3<-cbind(appo2b,p_fixed1)
          colnames(appo3)[ncol(appo3)-1]<-"PD_tot"
          appo4<-merge(appo1[,c(dom,x,z,"PD",ba)],appo3,
                       by = c(x,ba),
                       all.x=T)
          appo4$p_fixed_fin<-appo4$p_fixed1*appo4$PD/appo4$PD_tot
          appo5<-appo4
          

          for (j in 1:length(object$r_effect[[1]]))
          {
            xxx<-appo4
            for (k in 1:length(object$r_effect))
            {
            xxx<-as.data.frame(merge(xxx,as.data.frame(as.data.frame(object$r_effect[[k]][j])[,c(z[j],i)]),by=z[j],all.x=T))
            }
            xxx[is.na(xxx)]<-0
            if (k>1)
            {
              xxx[,i]<-rowSums(xxx[,(ncol(xxx)-k+1):ncol(xxx)])
            }
              appo5<-merge(appo5,xxx[,c(dom,x,z[j],i)],all.x=T)
            rm(xxx)
            appo5[,paste(i,j,sep="_")]<-appo5[,i]*appo5$PD
            appo5[,i]<-NULL
            }
          appo5b<-dir_s[,c(dom,x,z,i)]
          colnames(appo5b)[ncol(appo5b)]<-paste(i,"_d",sep="")
          appo6<-merge(appo5,appo5b,by=c(dom,x,z),all.x=T)
          appo6[is.na(appo6)]<-0
          cn<-colnames(appo5)[(ncol(appo5)-j+1):ncol(appo5)]
          
          appo7<-merge(appo4[,c(dom,x,z,"p_fixed_fin")],appo5[,c(dom,x,z,cn)])
          appo7<-merge(appo7,appo6[,c(dom,x,z,paste(i,"_d",sep=""))],all.x=T)
          appo7[is.na(appo7)]<-0          
          pred[[i]]<-data.frame(rowSums(appo7[,!(names(appo7)%in%c(dom,x,z))]))
          
          rm("appo6","appo5","appo5b","appo4","appo3","appo2b")
          gc()
          
                }
        pred<-do.call(cbind,pred)
        colnames(pred)<-y
        pred<-cbind(appo7[,c(dom,z,x)],pred)
  }
  
    return(as.data.frame(pred))
    NextMethod("mind")
}

#