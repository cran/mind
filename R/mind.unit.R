mind.unit<-function(formula,dom,data,universe,weights=NA,broadarea=NA,
               max_iter=200,max_diff=1e-05,phi_u0=0.05,
               REML=TRUE)
  {

  rho_u0=0.05

z_i_list<-list()
stima_omega_XZ_eblup_dom<-list()
stima_omega_XZ_proj_dom<-list()
stima_omega_X_proj_dom<-list()
stima_omega_Z_proj_dom<-list()
mse_BLUP<-list()
Z_piuu<-list()
cv_BLUP<-list()
u_omegaa<-list()
n_d<-list()
r_effect<-list()
beta_omega_broad<-list()
sigma_e<-list()
sigma_u_fin<-list()
ICC<-list()
mod_perf<-list()

id<-"supp_id"
universe[,id]<-1:nrow(universe)
data[,id]<-1:nrow(data)

if (is.na(weights))
{
  weights<-"weights"
  data[,weights]<-1
}

data_all<-data
if (is.na(broadarea))
{
  broadarea<-"broadarea"
  universe[,broadarea]<-1
}

allvars<-all.vars(formula)

# Trovo le y

y_y<-trimws(as.vector(unlist(genXtract(as.character(formula),"cbind(",")"))))
if(length(y_y)==0)
{
  y_y<-allvars[1]
  ncols<-ncol(data_all)
  data_all<-cbind(data_all,model.matrix(as.formula(paste("~-1+factor(",y_y,")",sep="")),data_all))
  y_y<-paste(allvars[1],1:(ncol(data_all)-ncols),sep="_")
  colnames(data_all)[(ncols+1):ncol(data_all)]<-y_y
  allvars<-c(y_y,allvars[-1])
}else{
  y_yy<-vector()
  for (i in 1:length(allvars)){
    y_yy[i]<-as.vector(regmatches(y_y, gregexpr(allvars[i], y_y)))
  }
  y_y<-unique(unlist(y_yy));rm(y_yy)
}
n_y<-length(y_y)

# Trovo le z

z_z<-trimws(as.vector(unlist(genXtract(as.character(formula),"|",")"))))
n_z<-length(z_z)

# Trovo le x

x_x<-allvars[which(!(allvars%in%c(y_y,z_z)))]
n_x<-length(x_x)

myfun<-function(x){length(unique(x))}
x_i<-apply(data_all[x_x],2,myfun)

int1<-grepl( "- 1",as.character(formula)[3],fixed = T)

if (!any(x_i==1) & int1==F)
{
  data_all$intercept<-1
  universe$intercept<-1
  formula<-update(formula,~.+factor(intercept))
  allvars<-c(allvars,"intercept")
  x_x<-allvars[which(!(allvars%in%c(y_y,z_z)))]
  n_x<-length(x_x)
  x_i<-apply(data_all[x_x],2,myfun)
  }

pos_id<-match(id,names(data))
pos_dom<-match(dom,names(data))
pos_w<-match(weights,names(data))

#

univ_all<-universe
bba<-sort(unique(univ_all[,broadarea]))
macro<- aggregate(as.formula(paste(broadarea,"~",dom,sep="")),univ_all,unique)

  for (ba in 1:length(unique(macro[,broadarea])))

  {

data<-data_all[,match(c(id,dom,y_y,z_z,names(sort(x_i)),weights),names(data_all))]
data<-data[data[,dom]%in%macro[macro[,broadarea]==bba[ba],dom],]

data_xz<-data[,c(3:ncol(data))]
for (i in 1 : n_y){data_xz[,c(i)]<-data_xz[,c(i)]*data_xz[,ncol(data_xz)]}

f<-".~-1"
  app<-paste(names(sort(x_i)), sep="", collapse="+")
    f<-paste(f,app,sep="+")

f1<-""
  app<-paste(z_z, sep="", collapse="+")
  f1<-paste(f1,app,sep="+")

f<-paste(f,f1,sep="")
f<-as.formula(f)

######## QUI!
# app<-paste(names(sort(x_i)), sep="", collapse="+")
# ff<-as.formula(paste("~",app,sep=""))
# if (any(x_i==1)) {beta_name<-model.matrix(ff,as.data.frame(data),contrasts.arg = lapply(as.data.frame(data), contrasts, contrasts=TRUE))
# } else {X_p<-model.matrix( update(ff,~.-1),as.data.frame(P_x_appo),contrasts.arg = lapply(as.data.frame(P_x_appo), contrasts, contrasts=TRUE))}
########
data_yq<-as.data.frame((data[,c(3:(2+n_y))]^2)*data[,c(ncol(data))])

data_yq<-setnames(data_yq,old=colnames(data_yq), new = c(paste("yq",1:n_y,sep="")))
data_yq$w_s<-1
data_xz<-cbind(data_xz,data_yq)

A_xz<-aggregate(data=data_xz,f,sum)
A_xz_tot<-A_xz[,c(1:(n_x+n_z+n_y),((ncol(A_xz)-n_y-1):ncol(A_xz)))]
A_xz<-A_xz[,c(1:(n_x+n_z+n_y),((ncol(A_xz)-n_y-1)))]
names(A_xz)[1:(n_x+n_z)]<-paste("delta",1:(n_x+n_z),sep="")

A_xz_s<-A_xz_tot[,c(1:(n_x+n_z+n_y))]
A_xz_s[,weights]<-A_xz_tot[,weights]

names(A_xz_s)[1:(n_x+n_z)]<-paste("delta",1:(n_x+n_z),sep="")

A_xz_q<-A_xz_tot[,c(1:(n_x+n_z), (ncol(A_xz_tot)-3):(ncol(A_xz_tot)-1),(n_x+n_z+n_y+1))]
names(A_xz_q)[1:(n_x+n_z)]<-paste("delta",1:(n_x+n_z),sep="")

app0<-".~-1"
r=1
M_xz_riga<-NULL
app1<-paste(paste("+factor(delta",r,sep=""),")",sep="")
f1<-paste(app0,app1,sep="")

app1<-paste(paste("+factor(delta",r,sep=""),")",sep="")
fr<-paste(app0,app1,sep="")
fr<-as.formula(fr)
Abk<-aggregate(data=A_xz,fr,sum)

Abk<-Abk[,c(1,(ncol(Abk)-n_y):ncol(Abk))]
names(Abk)[1]<-c("riga")

if(int1==F){
for (j in 1 : n_x){
  app2<-paste(paste("+factor(delta",j,sep=""),")",sep="")
  frc<-paste(f1,app2,sep="")
  frc<-as.formula(frc)
  app2<-""

  Abk<-aggregate(data=A_xz,frc,sum)

  if (r==j) (Abk<-cbind(Abk[,c(1)],Abk[,c(1,ncol(Abk))]))
  Abk<-Abk[,c(1,2,ncol(Abk))]
  names(Abk)[1:2]<-c("riga","colonna")

  m_riga<-as.data.frame(with(data=Abk,table(riga)))[1]
  m_colonna<-as.data.frame(with(data=Abk,table(colonna)))[1]
  m_riga<-as.numeric(as.vector(unlist(m_riga)))
  m_colonna<-as.numeric(as.vector(unlist(m_colonna)))
  R<-length(m_riga)
  C<-length(m_colonna)
  delta_rc<-NULL
  for (i  in 1 : length(m_colonna)){
    delta_rc<-rbind(delta_rc,cbind(m_riga,rep(m_colonna[i],length(m_riga))))}
  delta_rc<-as.data.frame(delta_rc)
  names(delta_rc)<-c("riga","colonna")

    Abk<-merge(delta_rc,Abk,by=c("riga","colonna"),all.x=T)
  Abk[,weights]<-ifelse(is.na(Abk[,weights]),0,Abk[,weights])
  if (C<=2){
    Abk<-Abk[order(Abk$colonna,Abk$riga),]
    M_xzAbk<-matrix(Abk[,weights], ncol=C, nrow=R)
  }
  if (C>2) {
    Abk<-Abk[order(Abk$riga,Abk$colonna),]
    M_xzAbk<-matrix(Abk[,weights], ncol=C, nrow=R, byrow=T)}
  if (j!=1) (M_xzAbk<-M_xzAbk[,-1])

  dim_s<-dim(as.data.frame(M_xzAbk))
  if (dim_s[2]==1) M_xz_riga<-cbind(M_xz_riga,t(as.data.frame(M_xzAbk)))
  if (dim_s[2]>1) M_xz_riga<-cbind(M_xz_riga,M_xzAbk)
  if (r==1 & j==n_x) (M_x<-M_xz_riga)
  if (r>1 & j==n_x) M_x<-rbind(M_x,M_xz_riga[-1,])
}
M_x_r1<-M_x
M_x=NULL

app0<-".~-1"
for (r in 2 : n_x){
  M_xz_riga<-NULL
  app1<-paste(paste("+factor(delta",r,sep=""),")",sep="")
  f1<-paste(app0,app1,sep="")

  app1<-paste(paste("+factor(delta",r,sep=""),")",sep="")
  fr<-paste(app0,app1,sep="")
  fr<-as.formula(fr)
  Abk<-aggregate(data=A_xz,fr,sum)

  Abk<-Abk[,c(1,(ncol(Abk)-n_y):ncol(Abk))]
  names(Abk)[1]<-c("riga")

  for (j in 2 : n_x){
    app2<-paste(paste("+factor(delta",j,sep=""),")",sep="")
    frc<-paste(f1,app2,sep="")
    frc<-as.formula(frc)
    app2<-""

    Abk<-aggregate(data=A_xz,frc,sum)

    if (r==j) (Abk<-cbind(Abk[,c(1)],Abk[,c(1,ncol(Abk))]))
    Abk<-Abk[,c(1,2,ncol(Abk))]
    names(Abk)[1:2]<-c("riga","colonna")

    m_riga<-as.data.frame(with(data=Abk,table(riga)))[1]
    m_colonna<-as.data.frame(with(data=Abk,table(colonna)))[1]
    m_riga<-as.numeric(as.vector(unlist(m_riga)))
    m_colonna<-as.numeric(as.vector(unlist(m_colonna)))
    R<-length(m_riga)
    C<-length(m_colonna)
    delta_rc<-NULL
    for (i  in 1 : length(m_colonna)){
      delta_rc<-rbind(delta_rc,cbind(m_riga,rep(m_colonna[i],length(m_riga))))}
    delta_rc<-as.data.frame(delta_rc)
    names(delta_rc)<-c("riga","colonna")

        Abk<-merge(delta_rc,Abk,by=c("riga","colonna"),all.x=T)
    Abk[,weights]<-ifelse(is.na(Abk[,weights]),0,Abk[,weights])
    if (C<=2){
      Abk<-Abk[order(Abk$colonna,Abk$riga),]
      M_xzAbk<-matrix(Abk[,weights], ncol=C, nrow=R)
    }
    if (C>2) {
      Abk<-Abk[order(Abk$riga,Abk$colonna),]
      M_xzAbk<-matrix(Abk[,weights], ncol=C, nrow=R, byrow=T)}
    dim_s2<-ncol(M_xzAbk)
    if (j!=1 & dim_s2==2) (M_xzAbk<-as.data.frame(M_xzAbk[,2]))
    if (j!=1 & dim_s2>2) (M_xzAbk<-M_xzAbk[,-1])

        if (j==2) (M_xz_riga<-M_xzAbk)
    if (j>2)  (M_xz_riga<-cbind(M_xz_riga,M_xzAbk))

    if (r==1 & j==n_x) (M_x<-M_xz_riga)
    if (r>1 & j==n_x) M_x<-rbind(M_x,M_xz_riga[-1,])
  }
}

  M_x_fin<-as.data.frame(cbind(t(M_x_r1)[-1],M_x))

   M_x_r1<-as.data.frame(M_x_r1)
  colnames(M_x_fin)<-names(M_x_r1)
  M_x_fin<-rbind(M_x_r1,M_x_fin)
  M_x<-as.matrix(M_x_fin)
}else{
  for (j in 1 : n_x){
    app2<-paste(paste("+factor(delta",j,sep=""),")",sep="")
    frc<-paste(f1,app2,sep="")
    frc<-as.formula(frc)
    app2<-""

    Abk<-aggregate(data=A_xz,frc,sum)

    if (r==j) (Abk<-cbind(Abk[,c(1)],Abk[,c(1,ncol(Abk))]))
    Abk<-Abk[,c(1,2,ncol(Abk))]
    names(Abk)[1:2]<-c("riga","colonna")

    m_riga<-as.data.frame(with(data=Abk,table(riga)))[1]
    m_colonna<-as.data.frame(with(data=Abk,table(colonna)))[1]
    m_riga<-as.numeric(as.vector(unlist(m_riga)))
    m_colonna<-as.numeric(as.vector(unlist(m_colonna)))
    R<-length(m_riga)
    C<-length(m_colonna)
    delta_rc<-NULL
    for (i  in 1 : length(m_colonna)){
      delta_rc<-rbind(delta_rc,cbind(m_riga,rep(m_colonna[i],length(m_riga))))}
    delta_rc<-as.data.frame(delta_rc)
    names(delta_rc)<-c("riga","colonna")

    Abk<-merge(delta_rc,Abk,by=c("riga","colonna"),all.x=T)
    Abk[,weights]<-ifelse(is.na(Abk[,weights]),0,Abk[,weights])
    if (C<=2){
      Abk<-Abk[order(Abk$colonna,Abk$riga),]
      M_xzAbk<-matrix(Abk[,weights], ncol=C, nrow=R)
    }
    if (C>2) {
      Abk<-Abk[order(Abk$riga,Abk$colonna),]
      M_xzAbk<-matrix(Abk[,weights], ncol=C, nrow=R, byrow=T)}
    if (j!=1) (M_xzAbk<-M_xzAbk[,-1])

    dim_s<-dim(as.data.frame(M_xzAbk))
    if (dim_s[2]==1) M_xz_riga<-cbind(M_xz_riga,t(as.data.frame(M_xzAbk)))
    if (dim_s[2]>1) M_xz_riga<-cbind(M_xz_riga,M_xzAbk)
    if (r==1 & j==n_x) (M_x<-M_xz_riga)
    if (r>1 & j==n_x) M_x<-rbind(M_x,M_xz_riga[-1,])
  }
if (n_x==1){
  M_x<-as.matrix(M_x)
}else{
  M_x_r1<-M_x
  M_x=NULL
  for (r in 2 : n_x){
    M_xz_riga<-NULL
    app1<-paste(paste("+factor(delta",r,sep=""),")",sep="")
    f1<-paste(app0,app1,sep="")

    app1<-paste(paste("+factor(delta",r,sep=""),")",sep="")
    fr<-paste(app0,app1,sep="")
    fr<-as.formula(fr)
    Abk<-aggregate(data=A_xz,fr,sum)

    Abk<-Abk[,c(1,(ncol(Abk)-n_y):ncol(Abk))]
    names(Abk)[1]<-c("riga")

    for (j in 2 : n_x){
      app2<-paste(paste("+factor(delta",j,sep=""),")",sep="")
      frc<-paste(f1,app2,sep="")
      frc<-as.formula(frc)
      app2<-""

      Abk<-aggregate(data=A_xz,frc,sum)

      if (r==j) (Abk<-cbind(Abk[,c(1)],Abk[,c(1,ncol(Abk))]))
      Abk<-Abk[,c(1,2,ncol(Abk))]
      names(Abk)[1:2]<-c("riga","colonna")

      m_riga<-as.data.frame(with(data=Abk,table(riga)))[1]
      m_colonna<-as.data.frame(with(data=Abk,table(colonna)))[1]
      m_riga<-as.numeric(as.vector(unlist(m_riga)))
      m_colonna<-as.numeric(as.vector(unlist(m_colonna)))
      R<-length(m_riga)
      C<-length(m_colonna)
      delta_rc<-NULL
      for (i  in 1 : length(m_colonna)){
        delta_rc<-rbind(delta_rc,cbind(m_riga,rep(m_colonna[i],length(m_riga))))}
      delta_rc<-as.data.frame(delta_rc)
      names(delta_rc)<-c("riga","colonna")

      Abk<-merge(delta_rc,Abk,by=c("riga","colonna"),all.x=T)
      Abk[,weights]<-ifelse(is.na(Abk[,weights]),0,Abk[,weights])
      if (C<=2){
        Abk<-Abk[order(Abk$colonna,Abk$riga),]
        M_xzAbk<-matrix(Abk[,weights], ncol=C, nrow=R)
      }
      if (C>2) {
        Abk<-Abk[order(Abk$riga,Abk$colonna),]
        M_xzAbk<-matrix(Abk[,weights], ncol=C, nrow=R, byrow=T)}
      dim_s2<-ncol(M_xzAbk)
      if (j!=1 & dim_s2==2) (M_xzAbk<-as.data.frame(M_xzAbk[,2]))
      if (j!=1 & dim_s2>2) (M_xzAbk<-M_xzAbk[,-1])

      if (j==2) (M_xz_riga<-M_xzAbk)
      if (j>2)  (M_xz_riga<-cbind(M_xz_riga,M_xzAbk))

      if (r==1 & j==n_x) (M_x<-M_xz_riga)
      if (r>1 & j==n_x) M_x<-rbind(M_x,M_xz_riga[-1,])
    }
  }

  M_x_r1<-as.data.frame(M_x_r1)
  M_x_fin<-as.data.frame(cbind(t(M_x_r1[,-c(1:sort(x_i)[1])]),M_x))
  #colnames(M_x_fin)<-names(M_x_r1)
  M_x_fin<-rbind(M_x_r1,M_x_fin)
  M_x<-as.matrix(M_x_fin)
}
    }

rm("Abk","app","app1","app2","C","delta_rc","m_colonna",
   "m_riga","M_xz_riga","M_xzAbk","R")

for (r in (n_x+1) : (n_x+n_z)){
  M_xz_riga<-NULL
  app1<-paste(paste("+factor(delta",r,sep=""),")",sep="")
  f1<-paste(app0,app1,sep="")

  app1<-paste(paste("+factor(delta",r,sep=""),")",sep="")
  fr<-paste(app0,app1,sep="")
  fr<-as.formula(fr)
  Abk<-aggregate(data=A_xz,fr,sum)

  Abk<-Abk[,c(1,(ncol(Abk)-n_y):ncol(Abk))]
  names(Abk)[1]<-c("riga")

  for (j in (n_x+1) : (n_x+n_z)){
    app2<-paste(paste("+factor(delta",j,sep=""),")",sep="")
    frc<-paste(f1,app2,sep="")
    frc<-as.formula(frc)
    app2<-""

    Abk<-aggregate(data=A_xz,frc,sum)

    if (r==j) (Abk<-cbind(Abk[,c(1)],Abk[,c(1,ncol(Abk))]))
    Abk<-Abk[,c(1,2,ncol(Abk))]
    names(Abk)[1:2]<-c("riga","colonna")

    m_riga<-as.data.frame(with(data=Abk,table(riga)))[1]
    m_colonna<-as.data.frame(with(data=Abk,table(colonna)))[1]
    m_riga<-as.numeric(as.vector(unlist(m_riga)))
    m_colonna<-as.numeric(as.vector(unlist(m_colonna)))
    R<-length(m_riga)
    C<-length(m_colonna)
    delta_rc<-NULL
    for (i  in 1 : length(m_colonna)){
      delta_rc<-rbind(delta_rc,cbind(m_riga,rep(m_colonna[i],length(m_riga))))}
    delta_rc<-as.data.frame(delta_rc)
    names(delta_rc)<-c("riga","colonna")

        Abk<-merge(delta_rc,Abk,by=c("riga","colonna"),all.x=T)
    Abk[,weights]<-ifelse(is.na(Abk[,weights]),0,Abk[,weights])
    if (C==2){
      Abk<-Abk[order(Abk$colonna,Abk$riga),]
      M_xzAbk<-matrix(Abk[,weights], ncol=C, nrow=R)
    }
    if (C>2) {
      Abk<-Abk[order(Abk$riga,Abk$colonna),]
      M_xzAbk<-matrix(Abk[,weights], ncol=C, nrow=R, byrow=T)}

        M_xz_riga<-cbind(M_xz_riga,M_xzAbk)
    if (r==(n_x+1) & j==(n_x+n_z)) (M_z<-M_xz_riga)
    if (r>(n_x+1) & j==(n_x+n_z)) M_z<-rbind(M_z,M_xz_riga)
  }}

rm("Abk","app1","app2","C","delta_rc","m_colonna",
   "m_riga","M_xz_riga","M_xzAbk","R")

for (r in 1 : n_x){
  M_xz_riga<-NULL
  app1<-paste(paste("+factor(delta",r,sep=""),")",sep="")
  f1<-paste(app0,app1,sep="")


  for (j in (n_x+1) : (n_x+n_z)){
    app2<-paste(paste("+factor(delta",j,sep=""),")",sep="")
    frc<-paste(f1,app2,sep="")
    frc<-as.formula(frc)
    app2<-""


    Abk<-aggregate(data=A_xz,frc,sum)

    if (r==j) (Abk<-cbind(Abk[,c(1)],Abk[,c(1,ncol(Abk))]))
    Abk<-Abk[,c(1,2,ncol(Abk))]
    names(Abk)[1:2]<-c("riga","colonna")

    m_riga<-as.data.frame(with(data=Abk,table(riga)))[1]
    m_riga<-as.numeric(as.vector(unlist(m_riga)))
    m_colonna<-as.data.frame(with(data=Abk,table(colonna)))[1]
    m_colonna<-as.numeric(as.vector(unlist(m_colonna)))

    R<-length(m_riga)
    C<-length(m_colonna)
    delta_rc<-NULL
    for (i  in 1 : length(m_colonna)){
      delta_rc<-rbind(delta_rc,cbind(m_riga,rep(m_colonna[i],length(m_riga))))}
    delta_rc<-as.data.frame(delta_rc)
    names(delta_rc)<-c("riga","colonna")

        Abk<-merge(delta_rc,Abk,by=c("riga","colonna"),all.x=T)
    Abk[,weights]<-ifelse(is.na(Abk[,weights]),0,Abk[,weights])
    if (C==2){
      Abk<-Abk[order(Abk$colonna,Abk$riga),]
      M_xzAbk<-matrix(Abk[,weights], ncol=C, nrow=R)
    }
    if (C>2) {
      Abk<-Abk[order(Abk$riga,Abk$colonna),]
      M_xzAbk<-matrix(Abk[,weights], ncol=C, nrow=R, byrow=T)}

        M_xz_riga<-cbind(M_xz_riga,M_xzAbk)
    if (r==1 & j==(n_x+n_z)) (M_xz<-M_xz_riga)
    if (r>1 & j==(n_x+n_z)) M_xz<-rbind(M_xz,M_xz_riga[-1,])
  }}

rm("Abk","app1","app2","C","delta_rc","m_colonna",
   "m_riga","M_xz_riga","M_xzAbk","R")

for (r in 1 : n_x){
  M_xz_riga<-NULL
  app1<-paste(paste("+factor(delta",r,sep=""),")",sep="")
  f1<-paste(app0,app1,sep="")

  for (j in (n_x+1) : (n_x+n_z)){
    app2<-paste(paste("+factor(delta",j,sep=""),")",sep="")
    frc<-paste(f1,app2,sep="")
    frc<-as.formula(frc)
    app2<-""

    Abk<-aggregate(data=A_xz_s,frc,sum)

    if (r==j) (Abk<-cbind(Abk[,c(1)],Abk[,c(1,ncol(Abk))]))
    Abk<-Abk[,c(1,2,ncol(Abk))]
    names(Abk)[1:2]<-c("riga","colonna")

    as.data.frame(with(data=Abk,table(riga)))[1]
    m_riga<-as.data.frame(with(data=Abk,table(riga)))[1]
    m_colonna<-as.data.frame(with(data=Abk,table(colonna)))[1]
    m_riga<-as.numeric(as.vector(unlist(m_riga)))
    m_colonna<-as.numeric(as.vector(unlist(m_colonna)))
    R<-length(m_riga)
    C<-length(m_colonna)
    delta_rc<-NULL
    for (i  in 1 : length(m_colonna)){
      delta_rc<-rbind(delta_rc,cbind(m_riga,rep(m_colonna[i],length(m_riga))))}
    delta_rc<-as.data.frame(delta_rc)
    names(delta_rc)<-c("riga","colonna")

       Abk<-merge(delta_rc,Abk,by=c("riga","colonna"),all.x=T)
    Abk[,weights]<-ifelse(is.na(Abk[,weights]),0,Abk[,weights])
    if (C==2){
      Abk<-Abk[order(Abk$colonna,Abk$riga),]
      M_xzAbk<-matrix(Abk[,weights], ncol=C, nrow=R)
    }
    if (C>2) {
      Abk<-Abk[order(Abk$riga,Abk$colonna),]
      M_xzAbk<-matrix(Abk[,weights], ncol=C, nrow=R, byrow=T)}

        M_xz_riga<-cbind(M_xz_riga,M_xzAbk)
    if (r==1 & j==(n_x+n_z)) (M_xz_s<-M_xz_riga)
    if (r>1 & j==(n_x+n_z)) M_xz_s<-rbind(M_xz_s
                                          ,M_xz_riga[-1,])
  }
  }
rm("Abk","app1","app2","C","delta_rc","m_colonna",
   "m_riga","M_xz_riga","M_xzAbk","R")

app0<-".~-1"

m_xy<-NULL

for (r in 1 : n_x){
  app1<-paste(paste("+factor(delta",r,sep=""),")",sep="")
  f1<-paste(app0,app1,sep="")

  app1<-paste(paste("+factor(delta",r,sep=""),")",sep="")
  fr<-paste(app0,app1,sep="")
  fr<-as.formula(fr)
  Abk<-aggregate(data=A_xz,fr,sum)

  Abk<-Abk[,c(1,(ncol(Abk)-n_y):ncol(Abk))]
  names(Abk)[1]<-c("riga")
  #blocco m_y1
  Y=Abk[,c(2:(ncol(Abk)-1))]
  if (r>1) (Y<-Y[-1,])
  y=as.numeric(t(Y))
  m_y=as.matrix(y)
  m_xy<-rbind(m_xy,m_y)
}

m_xy_q<-NULL

for (r in 1 : n_x){
  app1<-paste(paste("+factor(delta",r,sep=""),")",sep="")
  f1<-paste(app0,app1,sep="")

  app1<-paste(paste("+factor(delta",r,sep=""),")",sep="")
  fr<-paste(app0,app1,sep="")
  fr<-as.formula(fr)
  Abk<-aggregate(data=A_xz_q,fr,sum)

  Abk<-Abk[,c(1,(ncol(Abk)-n_y):ncol(Abk))]
  names(Abk)[1]<-c("riga")
  #blocco m_y1
  Y=Abk[,c(2:(ncol(Abk)-1))]
  if (r>1) (Y<-Y[-1,])
  y=as.numeric(t(Y))
  m_y=as.matrix(y)
  m_xy_q<-rbind(m_xy_q,m_y)
}

m_zy<-NULL
for (r in (n_x+1) : (n_x+n_z)){
  M_xz_riga<-NULL
  app1<-paste(paste("+factor(delta",r,sep=""),")",sep="")
  f1<-paste(app0,app1,sep="")

  app1<-paste(paste("+factor(delta",r,sep=""),")",sep="")
  fr<-paste(app0,app1,sep="")
  fr<-as.formula(fr)
  Abk<-aggregate(data=A_xz,fr,sum)

  Abk<-Abk[,c(1,(ncol(Abk)-n_y):ncol(Abk))]
  names(Abk)[1]<-c("riga")
  Y=Abk[,c(2:(ncol(Abk)-1))]
  y=as.numeric(t(Y))
  m_y=as.matrix(y)
  m_zy<-rbind(m_zy,m_y)
}

data_xzd<-data[,c(2:ncol(data))]
data_xzd$n<-1


app0<-".~-1"
appd<-paste("+factor(",dom,")",sep="")
app1<-paste(app0,appd,sep="")

fx=""
fx<-paste(paste("+factor(",names(sort(x_i)),")",sep=""),collapse="",sep="")

fz<-""
fz<-paste(paste("+factor(",z_z,")",sep=""),collapse="",sep="")

fd<-paste(app1,fx,fz,sep="")
fd<-as.formula(fd)

f<-paste(app0,fx,fz,sep="")
f<-as.formula(f)

A_xzd<-aggregate(data=data_xzd,fd,sum)
A_xzd<-A_xzd[,c((1:(n_x+n_z+1)),(ncol(A_xzd)-1),ncol(A_xzd))]
names(A_xzd)[1]<-dom
names(A_xzd)[2:(n_x+n_z+1)]<-paste("delta",1:(n_x+n_z),sep="")

A_xzd_new<-A_xzd[,-(ncol(A_xzd)-1)]


universe<-univ_all[univ_all[,broadarea]==bba[ba],]
ww<-names(universe)[which(!(names(universe)%in%c(id,dom,x_x,z_z,broadarea)))]
names(universe)[names(universe)==ww]<-weights
universe<-universe[,c(id,dom,names(sort(x_i)),z_z,weights)]
dom_aux<-aggregate(as.formula(paste(weights,"~",
                          paste(dom,paste(z_z,collapse = "+"),sep="+"),sep=""))
                                ,universe,sum)
z_i<-apply(universe[z_z],2,myfun)
z_i_list[[ba]]<-z_i

univ_xzd<-universe[,c(2:ncol(universe))]

app0<-".~-1"
appd<-paste("+factor(",dom,")",sep="")
app1<-paste(app0,appd,sep="")

fx=""
fx<-paste(paste("+factor(",names(sort(x_i)),")",sep=""),collapse="",sep="")

fz<-""
fz<-paste(paste("+factor(",z_z,")",sep=""),collapse="",sep="")

fd<-paste(app1,fx,fz,sep="")
fd<-as.formula(fd)

f<-paste(app0,fx,fz,sep="")
f<-as.formula(f)

U_xzd<-aggregate(data=univ_xzd,fd,sum)
U_xzd<-U_xzd[,c((1:(n_x+n_z+1)),ncol(U_xzd))]
names(U_xzd)[1]<-dom
names(U_xzd)[2:(n_x+n_z+1)]<-paste("delta",1:(n_x+n_z),sep="")
#length(with(data=U_xzd,table(dom)))

names(U_xzd)[ncol(U_xzd)]<-c("N")
XZ_dom<-(aggregate(U_xzd[,ncol(U_xzd)], list(dom=U_xzd[,dom]),sum))[1]
names(XZ_dom)[1]<-dom
xz_dom<-(aggregate(A_xzd[,ncol(A_xzd)], list(dom=A_xzd[,dom]),sum))[1]
names(xz_dom)[1]<-dom

XZ_dom1<-XZ_dom
xz_dom1<-xz_dom


P_xzd<-merge(U_xzd,A_xzd,all.x=T)
P_xzd[,weights]<-ifelse(is.na(P_xzd[,weights])==TRUE,0,P_xzd[,weights])
P_xzd$n<-with(data=P_xzd,ifelse(is.na(n)==TRUE,0,n))


f3<-""
f0<-".~-1"
for (i in 1 : (n_x+n_z)){
  app<-paste(paste("+factor(delta",i,sep=""),")",sep="")
  f3<-paste(f3,app,sep="")
}
fxz<-paste(f0,f3,sep="")
fxz<-as.formula(fxz)

P_xz<-aggregate(data=P_xzd,fxz,sum)

P_xz<-P_xz[,c((1:(n_x+n_z)),(ncol(P_xz)-2),ncol(P_xz))]
names(P_xz)[1:(n_x+n_z)]<-paste("delta",1:(n_x+n_z),sep="")

f1<-""
f0<-".~-1"
f00<-"~-1"

for (i in 1 : (n_x)){
  app<-paste(paste("+factor(delta",i,sep=""),")",sep="")
  f1<-paste(f1,app,sep="")}

fx<-paste(f0,f1,sep="")
fx<-as.formula(fx)

fxx<-paste(f00,f1,sep="")
fxx<-as.formula(fxx)

P_x<-aggregate(data=P_xzd,fx,sum)
P_x<-P_x[,c((1:n_x),(ncol(P_x)-2),ncol(P_x))]
names(P_x)[1:n_x]<-paste("delta",1:n_x,sep="")
P_x_appo<-P_x[,1:n_x]

P_x<-cbind(p_x=c(1:nrow(P_x)),P_x)

if (any(x_i==1)) {X_p<-model.matrix(~.,as.data.frame(P_x_appo[,-1]),contrasts.arg = lapply(as.data.frame(P_x_appo[,-1]), contrasts, contrasts=TRUE))
} else {X_p<-model.matrix(~.-1,as.data.frame(P_x_appo),contrasts.arg = lapply(as.data.frame(P_x_appo), contrasts, contrasts=TRUE))}

X_p=as.data.frame(X_p)
X_p=as.matrix(X_p)
X_p=(kronecker(X_p, diag(n_y), FUN = "*"))

f2<-""
f0<-".~-1"
f00<-"~-1"
for (i in (n_x+1) : (n_x+n_z)){
  app<-paste(paste("+factor(delta",i,sep=""),")",sep="")
  f2<-paste(f2,app,sep="")
}
fz<-paste(f0,f2,sep="")
fz<-as.formula(fz)

fzz<-paste(f00,f2,sep="")
fzz<-as.formula(fzz)

P_z<-aggregate(data=P_xz,fz,sum)

P_z_appo<-data.frame(P_z[,1:n_z])
P_z<-P_z[,c((1:n_z),(ncol(P_z)-1),ncol(P_z))]

names(P_z)[1:n_z]<-paste("delta",(n_x+1):(n_x+n_z),sep="")
P_z<-cbind(p_z=c(1:nrow(P_z)),P_z)

P_z_appo<-as.data.frame(P_z_appo)
Z_p<-model.matrix(~.-1,P_z_appo,contrasts.arg = lapply(P_z_appo, contrasts, contrasts=FALSE))

Z_p=as.data.frame(Z_p)
Z_p=as.matrix(Z_p)
Z_p=(kronecker(Z_p, diag(n_y), FUN = "*"))

# if (loop==2) {

P_xd_appo<-as.data.frame(P_xzd[,c(2:(n_x+1))])
P_zd_appo<-as.data.frame(P_xzd[,c((n_x+2):(n_x+n_z+1))])


if (any(x_i==1)) {X_pd<-model.matrix(~.,as.data.frame(P_xd_appo[,-1]),contrasts.arg = lapply(as.data.frame(P_xd_appo[,-1]), contrasts, contrasts=TRUE))
} else {X_pd<-model.matrix(~.-1,as.data.frame(P_xd_appo),contrasts.arg = lapply(as.data.frame(P_xd_appo), contrasts, contrasts=TRUE))}

P_xzd$r<-P_xzd$N-P_xzd$n

X_pd=as.data.frame(X_pd)

N_pd=X_pd*P_xzd$N
n_pd=X_pd*P_xzd$n
r_pd=X_pd*P_xzd$r

N_pd=cbind(P_xzd[,1],N_pd)
n_pd=cbind(P_xzd[,1],n_pd)
r_pd=cbind(P_xzd[,1],r_pd)

colnames(N_pd)[1]<-c("dom")
X_piu<-aggregate(N_pd[,-1], by=list(dom=N_pd$dom),sum)

colnames(n_pd)[1]<-c("dom")
x_piu<-aggregate(n_pd[,-1], by=list(dom=n_pd$dom),sum)

colnames(r_pd)[1]<-c("dom")
Xr_piu<-aggregate(r_pd[,-1], by=list(dom=r_pd$dom),sum)

X_piu<-X_piu[,-1]
x_piu<-x_piu[,-1]
Xr_piu<-Xr_piu[,-1]

X_piu<-as.matrix(X_piu)
x_piu<-as.matrix(x_piu)
Xr_piu<-as.matrix(Xr_piu)

P_zd_appo<-as.data.frame(P_zd_appo)
if (any(z_i==1)) {
  Z_pd<-model.matrix(~.,P_zd_appo[,z_i!=1],contrasts.arg = lapply(P_zd_appo[,z_i!=1], contrasts, contrasts=FALSE))
} else {
  Z_pd<-model.matrix(~.-1,P_zd_appo,contrasts.arg = lapply(P_zd_appo, contrasts, contrasts=FALSE))
}


Z_pd=as.data.frame(Z_pd)

N_pd=Z_pd*P_xzd$N
n_pd=Z_pd*P_xzd$n
r_pd=Z_pd*P_xzd$r

N_pd=cbind(P_xzd[,1],N_pd)
n_pd=cbind(P_xzd[,1],n_pd)
r_pd=cbind(P_xzd[,1],r_pd)

colnames(N_pd)[1]<-c("dom")
Z_piu<-aggregate(N_pd[,-1], by=list(dom=N_pd$dom),sum)

colnames(n_pd)[1]<-c("dom")
z_piu<-aggregate(n_pd[,-1], by=list(dom=n_pd$dom),sum)

colnames(r_pd)[1]<-c("dom")
Zr_piu<-aggregate(r_pd[,-1], by=list(dom=r_pd$dom),sum)

Z_piu<-Z_piu[,-1]
z_piu<-z_piu[,-1]
Zr_piu<-Zr_piu[,-1]

Z_piu<-as.matrix(Z_piu)
z_piu<-as.matrix(z_piu)
Zr_piu<-as.matrix(Zr_piu)



#}

xz_piu_dom<-aggregate(as.formula(paste(weights,"~",dom,sep="")),A_xzd,sum)

xz_piu_dom1<-xz_piu_dom[xz_piu_dom[,dom]%in%macro[macro[,broadarea]==bba[ba],dom],]

for (i in 1:length(z_z))
{
  if (!identical(dom_aux[,dom],dom_aux[,z_z[i]]))
  {
    tot_com1<-dom_aux[,c(z_z[i],weights)]
  }
}

M_x_kr<-kronecker(M_x, diag(n_y), FUN = "*")
M_z_kr<-kronecker(M_z, diag(n_y), FUN = "*")
M_xz_kr<-kronecker(M_xz, diag(n_y), FUN = "*")
M_xz_s_kr<-kronecker(t(M_xz_s), diag(n_y), FUN = "*")
M_x_inv=ginv(M_x)
M_x_inv_kr=kronecker(M_x_inv, diag(n_y), FUN = "*")
I_C<-diag(n_y)
Xr_piu_kr<-kronecker(Xr_piu, diag(n_y), FUN = "*")
Zr_piu_kr<-kronecker(Zr_piu, diag(n_y), FUN = "*")
X_piu_kr<-kronecker(X_piu, diag(n_y), FUN = "*")
Z_piu_kr<-kronecker(Z_piu, diag(n_y), FUN = "*")


scelta_z<-rep(1,n_z)
scelta_z<-as.data.frame(scelta_z)
scelta_z<-t(scelta_z)
colnames(scelta_z)<-paste("fatt",1:n_z,sep="")

#elenco_file_dist<-list(NA,"file2")


Q<-NULL
Q_cum<-NULL
Q_cum_app<-0

for (i in 1 : n_z){
  app<-length(table(A_xz[,c(n_x+i)]))
  Q_cum_app=Q_cum_app+app
  Q<-cbind(Q,app)
  Q_cum<-cbind(Q_cum,Q_cum_app)
}
Q<-as.data.frame(Q)
names(Q)<-paste("Q",1:n_z,sep="")
Q_cum<-as.data.frame(Q_cum)
names(Q_cum)<-paste("Q_cum",1:n_z,sep="")
Q_kr<-Q*n_y
Q_cum_kr<-Q_cum*n_y
base<-P_xz[,c((n_x+1):(n_x+n_z))]

n_iter<-max_iter

for (i in 1:n_z){
  app<-t(as.data.frame(rep(1,n_y)*phi_u0))
  colnames(app)<-paste("phi_u",1:n_y,sep="")
  rownames(app)<-"j"
  if (i==1) (phi_u<-app)
  if (i>1)  (phi_u<-rbind(phi_u,app))
}


for (i in 1:n_z){
  app<-t(as.data.frame(rep(1,n_y)*rho_u0))
  colnames(app)<-paste("rho_u",1:n_y,sep="")
  rownames(app)<-"j"
  if (i==1) (rho_u<-app)
  if (i>1)  (rho_u<-rbind(rho_u,app))
}


beta_ols_x<-M_x_inv_kr%*%m_xy


I_C<-diag(n_y)
for (c in 1 : n_y)
{
  i_cx<-as.data.frame(rep(I_C[,c],nrow(M_x)))
  i_cz<-as.data.frame(rep(I_C[,c],nrow(M_z)))
  m_xy_c<-m_xy*i_cx
  m_zy_c<-m_zy*i_cz
  m_xy_qc<-m_xy_q[seq(c,length(with(data=A_xz,table(delta1)))*n_y,n_y)]
  m_yc<-sum(m_xy_qc)
  sigma_0c<-(m_yc-t(m_xy_c)%*%beta_ols_x)*(sum(data$w)-nrow(M_x))^(-1)

  if (c==1) (sigma_0<-sigma_0c)
  if (c>1) (sigma_0<-rbind(sigma_0,sigma_0c))
}
sigma_0<-as.data.frame(t(sigma_0))
colnames(sigma_0)<-paste("sigma_0",1:n_y,sep="")
sigma_00<-sigma_0
sigma_00$iter<-0

  for (i in 1:n_z){
    app<-as.data.frame(rep(1,n_y)*phi_u0*sigma_0)
    colnames(app)<-paste("sigma_u",1:n_y,sep="")
    rownames(app)<-"j"
    if (i==1) (sigma_u<-app)
    if (i>1)  (sigma_u<-rbind(sigma_u,app))
  }

phi_u_start<-cbind(phi_u, fatt=c(1:n_z), iter=rep(0,n_z))
sigma_u_start<-cbind(sigma_u, fatt=c(1:n_z), iter=rep(0,n_z))
rho_u_start<-cbind(rho_u, fatt=c(1:n_z), iter=rep(0,n_z))
sigma_0_start<-sigma_00

mod_lmer<-list()

 for (i in 1:length(y_y))
{
 if (any(x_i==1)){
    formula_new<-update(formula,
                              as.formula(paste((y_y)[i],"~.- ",paste("factor(",names(x_i[x_i==1]),")",sep=""),sep="")))
  }else{formula_new<-update(formula,as.formula(paste((y_y)[i],"~.-1",sep="")))
    }
mod_lmer[[i]]<-lmer(formula_new,data=data)
}


I_C<-diag(n_y)

for (m in 1 : n_y)
{
app<-kronecker(mod_lmer[[m]]@beta, I_C[m,], FUN = "*")
if (m==1) (beta_omega_lmer<-app)
if (m>1)  (beta_omega_lmer<-beta_omega_lmer+app)
}

for (m in 1 : n_y)
{
app1<-as.data.frame(ranef(mod_lmer[[m]]))
app1<-app1[order(app1$grpvar),]
app1<-app1$condval
app<-kronecker(app1, I_C[m,], FUN = "*")
  if (m==1) (u_omega_lmer<-app)
  if (m>1)  (u_omega_lmer<-u_omega_lmer+app)
}


for (m in 1 : n_y)
{
vc_c<-as.data.frame(VarCorr(mod_lmer[[m]]))
  app<-vc_c[order(vc_c$grp),][-1,4]/vc_c[order(vc_c$grp),][1,4]
  if (m==1) (phi_u_lmer2<-app)
  if (m>1)  (phi_u_lmer2<-cbind(phi_u_lmer2,app))
}
colnames(phi_u_lmer2)<-paste("phi_u_lmer",1:n_y,sep="")
phi_u_lmer<-phi_u_lmer2

for (m in 1 : n_y)
{
vc_c<-as.data.frame(VarCorr(mod_lmer[[m]]))
mod_lmer[[m]]@theta
app<-vc_c[nrow(vc_c),4]
if (m==1) (sigma_0_lmer<-app)
if (m>1) (sigma_0_lmer<-cbind(sigma_0_lmer,app))
}
colnames(sigma_0_lmer)<-paste("sigma_0_lmer",1:n_y,sep="")

for (m in 1 : n_z)
{
app<-t(t(sigma_0_lmer)*t(phi_u_lmer)[,m])
if (m==1) (sigma_u_lmer<-app)
if (m>1) (sigma_u_lmer<-rbind(sigma_u_lmer,app))
}

colnames(sigma_u_lmer)<-paste("sigma_u_lmer",1:n_y,sep="")

# if (modelli==1) {n_iter<-1
#                  sigma_0<-sigma_0_lmer
#                  sigma_0<-cbind(sigma_0,sigma_0_start[,ncol(sigma_0_start)])
#                  colnames(sigma_0)<-c(paste("sigma_0",1:n_y,sep=""),"iter")
#
#                  colf<-ncol(phi_u_start)
#                  phi_u<-phi_u_lmer
#                  phi_u<-cbind(phi_u,phi_u_start[,(colf-1):colf])
#                  colnames(phi_u)<-c(paste("phi_u",1:n_y,sep=""),"fatt","iter")
#
#
#                  cols<-ncol(sigma_u_start)
#                  sigma_u<-sigma_u_lmer
#                  sigma_u<-cbind(sigma_u,sigma_u_start[,(cols-1):cols])
#                  colnames(sigma_u)<-c(paste("sigma_u",1:n_y,sep=""),"fatt","iter")
#
#                  }

# if (modelli==2) {n_iter<-max_iter
#                  sigma_0<-sigma_0_lmer
#                  sigma_0<-cbind(sigma_0,sigma_0_start[,ncol(sigma_0_start)])
#                  colnames(sigma_0)<-c(paste("sigma_0",1:n_y,sep=""),"iter")
#
#                  colf<-ncol(phi_u_start)
#                  phi_u<-phi_u_lmer
#                  phi_u<-cbind(phi_u,phi_u_start[,(colf-1):colf])
#                  colnames(phi_u)<-c(paste("phi_u",1:n_y,sep=""),"fatt","iter")
#
#
#                  cols<-ncol(sigma_u_start)
#                  sigma_u<-sigma_u_lmer
#                  sigma_u<-cbind(sigma_u,sigma_u_start[,(cols-1):cols])
#                  colnames(sigma_u)<-c(paste("sigma_u",1:n_y,sep=""),"fatt","iter")
#
#                  }

#if (modelli==3) {
n_iter<-max_iter
sigma_0<-t(rep(1,n_y)*0.15)
sigma_0<-cbind(sigma_0,sigma_0_start[,ncol(sigma_0_start)])
colnames(sigma_0)<-c(paste("sigma_0",1:n_y,sep=""),"iter")

colf<-ncol(phi_u_start)
phi_u<-matrix(0.05,nrow = n_z,ncol = n_y)
phi_u_start<-as.data.frame(phi_u_start)
phi_u<-cbind(phi_u,phi_u_start[,(colf-1):colf]) ##### PROBLEMA
colnames(phi_u)<-c(paste("phi_u",1:n_y,sep=""),"fatt","iter")

cols<-ncol(sigma_u_start)
sigma_u<-sigma_u_lmer
sigma_u<-cbind(sigma_u,sigma_u_start[,(cols-1):cols])
colnames(sigma_u)<-c(paste("sigma_u",1:n_y,sep=""),"fatt","iter")


#}

i<-1
dif<-0.9

while ( i<=n_iter && dif>=max_diff) {

  for (j in 1 : n_z){
    scelta<-as.numeric(scelta_z[j])
    sigma_ej<-as.numeric(sigma_0[1:n_y])
    phi_uj<-as.numeric(phi_u[j,1:n_y])
    rho_uj<-as.numeric(rho_u[j,1:n_y])

        sigma_uj<-as.numeric(sigma_u[j,1:n_y])

    if (scelta==1)
    {
      diag_phi_uj<-diag(phi_uj)
      diag_sigma_uj<-diag(sigma_uj)
      Omegaj<-kronecker(diag(Q[j]), diag(n_y), FUN = "*")
      phij_Omegaj<-kronecker(diag(Q[j]), diag_phi_uj, FUN = "*")
      sigmaj_Omegaj<-kronecker(diag(Q[j]), diag_sigma_uj, FUN = "*")
    }

    if (j==1) {Omegaj_tot<-Omegaj
    Omega<-phij_Omegaj
    Omegaa<-sigmaj_Omegaj
    }

    if (j>1) {Omegaj_tot<-bdiag(Omegaj_tot,Omegaj)
    Omega<-bdiag(Omega,phij_Omegaj)
    Omegaa<-bdiag(Omegaa,sigmaj_Omegaj)
    }

    Omega<-as.matrix(Omega)
    Omegaa<-as.matrix(Omegaa)
  }

    T_star<-ginv(M_z_kr+solve(Omega))

    M_x_omega<-M_x_kr-M_xz_kr%*%T_star%*%t(M_xz_kr)
  M_x_inv_omega<-ginv(M_x_omega)
  m_xy_omega<-m_xy-M_xz_kr%*%T_star%*%m_zy
  beta_omega<-M_x_inv_omega%*%m_xy_omega

    u_omega1<-m_zy-t(M_xz_kr)%*%beta_omega
  u_omega<-T_star%*%u_omega1

  for (c in 1 : n_y)
  {
    i_cx<-as.data.frame(rep(I_C[,c],nrow(M_x)))
    i_cz<-as.data.frame(rep(I_C[,c],nrow(M_z)))
    m_xy_c<-m_xy*i_cx
    m_zy_c<-m_zy*i_cz
    m_xy_qc<-m_xy_q[seq(c,length(with(data=A_xz,table(delta1)))*n_y,n_y)]

        m_yc<-sum(m_xy_qc)
    sigma_0c<-(m_yc-t(m_xy_c)%*%beta_omega-t(m_zy_c)%*%u_omega)*(nrow(data)-nrow(M_x))^(-1)
    if (REML==FALSE) (sigma_0c<-(m_yc-t(m_xy_c)%*%beta_omega-t(m_zy_c)%*%u_omega)*nrow(data)^(-1))
    if (c==1) (sigma_0<-sigma_0c)
    if (c>1) (sigma_0<-rbind(sigma_0,sigma_0c))
  }

  sigma_0<-t(sigma_0)
  sigma_0<-as.data.frame(sigma_0)
  colnames(sigma_0)<-paste("sigma",1:n_y,sep="")
  sigma_0$iter<-i

  T_s<-T_star+T_star%*%t(M_xz_kr)%*%M_x_inv_omega%*%M_xz_kr%*%T_star

if (REML==FALSE) (T_s=T_star)
T_s[,1]-T_star[,1]

  diag_sigma_0<-diag(sigma_0[1:n_y])
  diag_sigma_0_kr<-kronecker(diag(sum(Q)), diag_sigma_0, FUN = "*")
  diag_sigma_0_inv<-solve(diag_sigma_0)
  diag_sigma_0_inv_kr<-kronecker(diag(sum(Q)), diag_sigma_0_inv, FUN = "*")
  Omega_D<-diag_sigma_0_inv_kr%*%Omegaa
  Omega_D[,1]-Omega[,1]
  ZSZ<-M_z_kr-t(M_xz_kr)%*%M_x_inv_kr%*%M_xz_kr
  ZSZD<-ZSZ%*%Omega

  T_w_<-diag(sum(Q)*n_y)+ZSZD


  ZWZ<-M_z_kr
  W_=diag(sum(Q)*n_y)+ZWZ%*%Omega
  T_star_Omega<-Omega%*%ginv(T_star)
  T_star_Omega[,1]-W_[,1]
  sum(diag(solve(T_star_Omega)))
  sum(diag(solve(W_)))

   for (j in 1 : n_z){
     if (j==1) (a1<-as.numeric(Q_cum_kr[j]-Q_cum_kr[j]+1))
     if (j>1)  (a1<-as.numeric(Q_cum_kr[j-1]+1))
     a2<-as.numeric(Q_cum_kr[j])
     Omegaj_inv<-solve(Omegaj_tot[a1:a2,a1:a2])
     uj<-u_omega[a1:a2]
     T_sj<-T_s[a1:a2,a1:a2]

     for (c in 1 : n_y)
     {
       i_uj_c<-as.data.frame(rep(I_C[,c],Q[j]))
       I_C_c<-I_C*I_C[,c]
       F_jc<-kronecker(diag(Q[j]),I_C_c, FUN = "*")

       uj_c<-as.vector(uj*i_uj_c)
       uj_c<-unlist(uj_c)


        phi_jc_add1<-as.numeric(unlist(sum(diag(T_sj%*%Omegaj_inv%*%F_jc))))
        phi_jc_add2<- (1/sigma_0[c])*as.numeric((t(uj_c)%*%(Omegaj_inv%*%F_jc)%*%uj_c))
        phi_jc<-Q[j]^(-1)*(phi_jc_add1+phi_jc_add2)


       if (c==1) (phi_jc_iter<-phi_jc)
       if (c>1) (phi_jc_iter<-cbind(phi_jc_iter,phi_jc))
     }

     colnames(phi_jc_iter)<-paste("phi_j",1:n_y,sep="")
     phi_jc_iter$fatt<-j
     phi_jc_iter$iter<-i

     if (j==1) (phi_j_iter<-phi_jc_iter)

     if (j>1) (phi_j_iter<-rbind(phi_j_iter,phi_jc_iter))

   }


  if (i==1)
  {
  phi_j<-phi_j_iter
  sigma_j<-sigma_0

  diff_phi<-max(abs(phi_j[phi_j$iter==1,1:n_y]-phi_u[,1:n_y]))
  diff_sigma<-max(abs(sigma_0[sigma_0$iter==1,1:n_y]-sigma_00[sigma_00$iter==0,1:n_y]))}

  if (i>1)
  {
  phi_j<-rbind(phi_j,phi_j_iter)
  sigma_j<-rbind(sigma_j,sigma_0)
  diff_phi<-max(abs(phi_j[phi_j$iter==i,1:n_y]-phi_j[phi_j$iter==(i-1),1:n_y]))
  diff_sigma<-max(abs(sigma_j[sigma_j$iter==i,1:n_y]-sigma_j[sigma_j$iter==(i-1),1:n_y]))}

  dif<-max(diff_phi,diff_sigma)

  sigma_0<-sigma_j[sigma_j$iter==i,1:n_y]
  phi_u<-phi_j[phi_j$iter==i,1:n_y]

  i=i+1

}

# if (modelli==1) {
#   phi_j[phi_j$iter==max(phi_j$iter),c(1:n_y)]=phi_u_lmer
#   sigma_j[sigma_j$iter==max(sigma_j$iter),c(1:n_y)]=sigma_0_lmer
# }


stima_omega_X_proj<-X_piu_kr%*%beta_omega

u_omega_sample<-u_omega
r<-sum(apply(data[z_z],2,myfun))
r1<-apply(data[z_z[1]],2,myfun)


u_omega_mat<-matrix(u_omega_sample,nrow=r,ncol=n_y,byrow=TRUE)

domeff<-data.frame(dom=univ_xzd[,dom],z1=univ_xzd[,z_z[1]])
domeff<-aggregate(domeff$z1,by=list(dom=domeff$dom),unique)
colnames(domeff)[2]<-"z1"
domeff_c<-data.frame(dom=data_xzd[,dom],z1=data_xzd[,z_z[1]])
domeff_c<-aggregate(domeff_c$z1,by=list(dom=domeff_c$dom),unique)
colnames(domeff_c)[2]<-"z1"
u_omega_mat_z1<-cbind(z1=unique(domeff_c$z1),u_omega_mat[1:r1,])
u_omega_mat_z1<-merge(domeff,u_omega_mat_z1,by="z1",all.x=T)

u_omega_ext<-u_omega_mat_z1[,-1]
u_omega_ext1<-u_omega_ext[,-1]
u_omega_ext1[is.na(u_omega_ext1)]<-0
u_omega_ext1<-as.matrix(u_omega_ext1)
u_omega=u_omega_ext1
if (r>r1) (u_omega=rbind(u_omega_ext1,u_omega_mat[(r1+1):r,]))
u_omega<-as.numeric(t(u_omega))
u_omega<-as.matrix(u_omega)
if (length(unique(data_xzd[,z_z[1]]))<length(unique(data_xzd[,dom]))) (u_omega=u_omega_sample)

stima_omega_Z_proj<-Z_piu_kr%*%u_omega
stima_omega_XZ_proj<-stima_omega_X_proj+stima_omega_Z_proj

stima_omega_X_eblup0<-Xr_piu_kr%*%beta_omega
stima_omega_Z_eblup0<-Zr_piu_kr%*%u_omega
stima_omega_XZ_eblup0<-stima_omega_X_eblup0+stima_omega_Z_eblup0

u_omegaa[[ba]]<-u_omega

sample_y<-aggregate(data[3:(3+n_y-1)], list(dom=data[,dom]),sum)
yy=as.numeric(t(sample_y[,-1]))
yy=as.matrix(yy)

domains2<-as.data.frame(macro[macro[,broadarea]==bba[ba],dom])
colnames(domains2)<-dom

yy_mat<-matrix(yy,nrow=nrow(sample_y),ncol=n_y,byrow=TRUE)

yy_ext<-cbind(data.frame(dom=xz_piu_dom1[,1]),yy_mat)
names(yy_ext)[1]<-dom
yy_ext<-merge(domains2,yy_ext,by=dom,all.x=TRUE)
yy_ext[is.na(yy_ext)]<-0
yy_ext<-as.numeric(t(yy_ext[,-1]))
yy_ext<-as.matrix(yy_ext)

yy_sample<-yy
yy<-yy_ext

stima_omega_XZ_eblup<-yy+stima_omega_XZ_eblup0

stima_omega_XZ_eblup<-matrix(stima_omega_XZ_eblup,nrow=nrow(X_piu),ncol=n_y,byrow=TRUE)
stima_omega_XZ_proj<-matrix(stima_omega_XZ_proj,nrow=nrow(X_piu),ncol=n_y,byrow=TRUE)
stima_omega_X_proj<-matrix(stima_omega_X_proj,nrow=nrow(X_piu),ncol=n_y,byrow=TRUE)
stima_omega_Z_proj<-matrix(stima_omega_X_proj,nrow=nrow(X_piu),ncol=n_y,byrow=TRUE)


domains<-macro[macro[,broadarea]==bba[ba],dom]

beta_omega_broad[[ba]]<-beta_omega

stima_omega_XZ_eblup_dom[[ba]]<-cbind(data.frame(dom=domains),stima_omega_XZ_eblup)
stima_omega_XZ_proj_dom[[ba]]<-cbind(data.frame(dom=domains),stima_omega_XZ_proj)
stima_omega_X_proj_dom[[ba]]<-cbind(data.frame(dom=domains),stima_omega_X_proj)

#############################
############# MSE
#############################

univ_xzd$n<-1
univ_xzd$w_inv<-univ_xzd[,weights]^(-1)
data_xzd$w_inv<-data_xzd[,weights]^(-1)

N_d<-aggregate(univ_xzd[,c((ncol(univ_xzd)-2):ncol(univ_xzd))], by=list(dom=univ_xzd[,dom]),sum)
colnames(N_d)<-c("dom","W","N","W_inv")
n_dd<-aggregate(data_xzd[,c((ncol(data_xzd)-2):ncol(data_xzd))], by=list(dom=data_xzd[,dom]),sum)
r_d<-merge(N_d,n_dd, by=("dom"),all.x = T)
colnames(r_d)[5]<-c("w")
r_d$w[is.na(r_d$w)]<-0
r_d$w_inv[is.na(r_d$w_inv)]<-0
r_d$n[is.na(r_d$n)]<-0

r_d$W_r=r_d$W-r_d$w
r_d$W_inv_r=r_d$W_inv-r_d$w_inv
r_d$N_r=r_d$N-r_d$n

ndom<-ncol(Zr_piu)-(ncol(Zr_piu)-nrow(Zr_piu))
n_mod_fatt<-ncol(Zr_piu)

Zr_piu1<-Zr_piu[,c(1:r1)]
if (r>r1) {Zr_piu2=Zr_piu[,c((r1+1):r)]}
TT<-as.matrix(T_star)
Zr_piu_dom<-cbind(domains2,Zr_piu1)


if (r>r1) {
  Zr_piu2<-cbind(dom=domeff$dom, Zr_piu2)
  colnames(Zr_piu2)[1]<-dom
}
if (r>r1) {Zr_piu_dom<-merge(Zr_piu2, Zr_piu_dom, by=dom, all.x=TRUE)}

Zr_piu_dom<-Zr_piu_dom[,-1]
Zr_piu_dom<-as.matrix(Zr_piu_dom)
Zr_piu_dom_kr<-kronecker(Zr_piu_dom, diag(n_y), FUN = "*")
Zr_piu_dom_kr<-as.matrix(Zr_piu_dom_kr)

rho_j<-phi_j

M_z_omega_ML<-M_z_kr-M_z_kr%*%T_star%*%M_z_kr
M_xz_omega<-M_xz_s_kr%*%M_x_inv_omega%*%t(M_xz_s_kr)
I<-diag(sum(Q_kr))
L<-I-M_z_omega_ML%*%Omega
M_z_omega_REML<-M_z_omega_ML-L%*%M_xz_omega%*%t(L)
n=nrow(data)

k=ncol(M_x_inv_omega)/n_y
I_C<-diag(n_y)

  M_z_omega<-M_z_omega_REML

  for (j in 1 : n_z){
    if (j==1) (a1<-as.numeric(Q_cum_kr[j]-Q_cum_kr[j]+1))
    if (j>1)  (a1<-as.numeric(Q_cum_kr[j-1]+1))
    a2<-as.numeric(Q_cum_kr[j])
    Omegaj_inv<-solve(Omegaj_tot[a1:a2,a1:a2])
    T_sj<-T_s[a1:a2,a1:a2]

    for (c in 1 : n_y)
    {
      sigma<-sigma_j[sigma_j$iter==max(sigma_j$iter),c]
      phi<-phi_j[(phi_j$iter==max(phi_j$iter) & phi_j$fatt==j),c]
      I_C_c<-I_C*I_C[,c]
      F_jc<-kronecker(diag(Q[j]),I_C_c, FUN = "*")

      for (jj in 1 : n_z){
        a<-as.numeric(Q[jj])
        block<-matrix(0,nrow=a,ncol=a)
        block<-kronecker(block,I_C_c, FUN = "*")
        if (jj==j)  (block<-kronecker(diag(a),I_C_c, FUN = "*"))
        if (jj==1) F0_jc<-block
        if (jj>1)  F0_jc<-bdiag(F0_jc,block)
      }

      II_11_j<-(n-k)/sigma^2

mat_jc<-M_z_omega[a1:a2,a1:a2]%*%F_jc
      M_zj._omega<-M_z_omega[a1:a2,a1:a2]
      tr_jc=sum(diag(as.matrix(M_zj._omega)%*%F_jc))
      II_12_jc<-tr_jc/sigma



      if (c==1) {I_11_j<-II_11_j
      I_12_jc<-II_12_jc}

      if (c>1)  {I_11_j<-cbind(I_11_j,II_11_j)
      I_12_jc<-cbind(I_12_jc,II_12_jc)}

      if (c==n_y) {I_11_j=as.vector(I_11_j)
      I_11_j=diag(I_11_j)
      I_12_jc=I_12_jc=as.vector(I_12_jc)
      I_12_jc=diag(I_12_jc)}

    }

    colnames(I_12_jc)<-paste("I_12_j",1:n_y,sep="")
    if (j==1) {I_12_j<-I_12_jc}

    if (j>1)   {I_12_j<-cbind(I_12_j,I_12_jc)}


  }

    for (h in 1 : n_z){
    for (j in 1 : n_z){
      if (j==1) (a1<-as.numeric(Q_cum_kr[j]-Q_cum_kr[j]+1))
      if (j>1)  (a1<-as.numeric(Q_cum_kr[j-1]+1))
      a2<-as.numeric(Q_cum_kr[j])
      Omegaj_inv<-solve(Omegaj_tot[a1:a2,a1:a2])
      T_sj<-T_s[a1:a2,a1:a2]

      if (h==1) (a1_h<-as.numeric(Q_cum_kr[h]-Q_cum_kr[h]+1))
      if (h>1)  (a1_h<-as.numeric(Q_cum_kr[h-1]+1))
      a2_h<-as.numeric(Q_cum_kr[h])
      Omegah_inv<-solve(Omegaj_tot[a1_h:a2_h,a1_h:a2_h])
      T_sh<-T_s[a1_h:a2_h,a1_h:a2_h]
      T_shj<-T_s[a1_h:a2_h,a1:a2]
      T_sjh<-T_s[a1:a2,a1_h:a2_h]

      for (c in 1 : n_y)
      {
        sigma<-sigma_j[sigma_j$iter==max(sigma_j$iter),c]
        phij<-phi_j[(phi_j$iter==max(phi_j$iter) & phi_j$fatt==j),c]
        phih<-phi_j[(phi_j$iter==max(phi_j$iter) & phi_j$fatt==h),c]

        I_C_c<-I_C*I_C[,c]
        F_jc<-kronecker(diag(Q[j]),I_C_c, FUN = "*")
        F_hc<-kronecker(diag(Q[h]),I_C_c, FUN = "*")
        calcolo1=3

                if (calcolo1==3)
        {M_zj._omega<-M_z_omega[a1_h:a2_h,a1:a2]
        M_zh._omega<-M_z_omega[a1:a2,a1_h:a2_h]
        M_zjh._omega<-as.matrix(M_zj._omega)%*%F_jc%*%as.matrix(M_zh._omega)%*%F_hc
        II_22_jh_c=as.numeric(unlist(sum(diag(as.matrix(M_zjh._omega)))))
        }  else  {r_jc<-phi^(-1)*as.numeric(unlist(sum(diag(T_sj%*%Omegaj_inv%*%F_jc))))
        r_jhc<-as.numeric(unlist(sum(diag(T_shj%*%Omegaj_inv%*%F_jc%*%T_sjh%*%Omegah_inv%*%F_hc))))

        if (j==h) {II_22_jh_c<-phij^(-2)*((Q[j]-2*r_jc)*1+phij^(-2)*phih^(-2)*r_jhc)}
        else {II_22_jh_c<-phij^(-2)*((Q[j]-2*r_jc)*0+phij^(-2)*phih^(-2)*r_jhc)}
        }

        if (c==1) {I_22_jh_c<-II_22_jh_c}

        if (c>1)  {I_22_jh_c<-cbind(I_22_jh_c,II_22_jh_c)}

        if (c==n_y) {I_22_jh_c=as.vector(I_22_jh_c)
        I_22_jh_c=diag(I_22_jh_c)}

      }

            colnames(I_22_jh_c)<-paste("I_22_hj",1:n_y,sep="")
      if (j==1)  {I_22_jh<-I_22_jh_c
      } else {I_22_jh<-cbind(I_22_jh,I_22_jh_c)}

      if (h==1 & j==n_z)   {I_22<-I_22_jh}
      if (h>1 & j==n_z)    {I_22<-rbind(I_22,I_22_jh)}
      }
  }

 I_1<-t(cbind(I_11_j,I_12_j))
  I_2<-rbind(I_12_j,I_22)
  I_ANOVA<-cbind(I_1,I_2)
  I_ANOVA<-as.matrix(I_ANOVA)
  mat_Inf=0.5*I_ANOVA
  inv_Inf<-ginv(mat_Inf)
  inf_rho<-0*diag(n_y*n_z)
  inv_Inf1<-as.matrix(bdiag(inv_Inf,inf_rho))

  mat_Inf0<-as.matrix(bdiag(mat_Inf,inf_rho))
  inv_Inf0<-ginv(mat_Inf0)

  de_Omega_phi_tot<-Omegaj_tot

  for (j in 1 : n_z){
    if (j==1) (a1<-as.numeric(Q_cum_kr[j]-Q_cum_kr[j]+1))
    if (j>1)  (a1<-as.numeric(Q_cum_kr[j-1]+1))
    a2<-as.numeric(Q_cum_kr[j])

    if (j==1) {Omega_0<-0*Omega[a1:a2,a1:a2]}
    if (j>1)  {Omega_0<-bdiag(Omega_0,0*Omega[a1:a2,a1:a2])}
    Omega_0_inv_phi<-as.matrix(Omega_0)
    Omega_0_inv_rho<-as.matrix(Omega_0)
  }

  for (j in 1 : n_z){

    if (j==1) (a1<-as.numeric(Q_cum_kr[j]-Q_cum_kr[j]+1))
    if (j>1)  (a1<-as.numeric(Q_cum_kr[j-1]+1))
    a2<-as.numeric(Q_cum_kr[j])

    if (j==1) {Omega_0<-0*Omega[a1:a2,a1:a2]}
    if (j>1)  {Omega_0<-bdiag(Omega_0,0*Omega[a1:a2,a1:a2])}
    Omega_0_inv_phi<-as.matrix(Omega_0)
    Omega_0_inv_rho<-as.matrix(Omega_0)
  }
  for (j in 1 : n_z){
    scelta<-as.numeric(scelta_z[j])
    sigma_ej<-as.numeric(sigma_j[sigma_j$iter==max(sigma_j$iter),-ncol(sigma_j)])
    phi_uj<-as.numeric(phi_j[(phi_j$iter==max(phi_j$iter) & phi_j$fatt==j),1:(ncol(phi_j)-2)])
    rho_uj<-as.numeric(rho_j[(rho_j$iter==max(rho_j$iter) & rho_j$fatt==j),1:(ncol(rho_j)-2)])

    if (j==1) (a1<-as.numeric(Q_cum_kr[j]-Q_cum_kr[j]+1))
    if (j>1)  (a1<-as.numeric(Q_cum_kr[j-1]+1))
    a2<-as.numeric(Q_cum_kr[j])
    Omegaj<-Omega[a1:a2,a1:a2]
    Omegaj_inv<-ginv(Omegaj)
    de_Omegaj_phij<-Omegaj_tot[a1:a2,a1:a2]

    if (scelta==1)
    {
      diag_sigma_uj<-diag(sigma_uj)
      de_Omegaj_rho<-0*kronecker(diag(Q[j]), diag_sigma_uj, FUN = "*")
      de_Omega_rho<-de_Omegaj_rho
    }

    de_Omega_inv_phij<--Omegaj_inv%*%de_Omegaj_phij%*%Omegaj_inv
    de_Omega_inv_rhoj<--Omegaj_inv%*%de_Omegaj_rho%*%Omegaj_inv

    Omega_0_inv_phi[a1:a2,a1:a2]=as.matrix(de_Omega_inv_phij)
    Omega_0_inv_rho[a1:a2,a1:a2]=as.matrix(de_Omega_inv_rhoj)

    if (j==1) {de_Omega_inv_phi<-Omega_0_inv_phi
    de_Omega_inv_rho<-(-1)*Omega_0_inv_rho}
    if (j>1) {de_Omega_inv_phi<-rbind(de_Omega_inv_phi,Omega_0_inv_phi)
    de_Omega_inv_rho<-rbind(de_Omega_inv_rho,Omega_0_inv_rho)}

   }

  de_Omega_inv_gamma<-rbind(de_Omega_inv_phi,de_Omega_inv_rho)

  B=inv_Inf1[(n_y+1):(2*n_z*n_y+n_y),(n_y+1):(2*n_z*n_y+n_y)]

  A<-Zr_piu_dom_kr%*%TT
  A_kr<-kronecker(A, diag(n_z*2), FUN = "*")
  Delta_inv_alpha=-A_kr%*%de_Omega_inv_gamma%*%TT

  Sigma_s_star<-M_z_kr+M_z_kr%*%Omega%*%M_z_kr
  G3_a<-Delta_inv_alpha%*%Sigma_s_star%*%t(Delta_inv_alpha)

  D<-nrow(domains2)

  amp=2*n_z*n_y

for (rr1 in 1:D){
    initial.time_mse <- proc.time()
    ar<-amp*(rr1-1)+1
    br<-amp*rr1

    c1=rr1
    ac<-amp*(c1-1)+1
    bc<-amp*c1
    G3_rc<-G3_a[ar:br,ac:bc]
    for (c in 1:n_y){
      I_C_c<-I_C*I_C[,c]

      F_c<-kronecker(diag(2*n_z),I_C_c, FUN = "*")


      B_c<-B%*%F_c
      if (c==1) (G3=sigma_ej[c]*sum(diag(G3_rc%*%F_c%*%B_c)))
      if (c>1)  (G3<-cbind(G3,sigma_ej[c]*sum(diag(G3_rc%*%F_c%*%B_c))))
      if (c==n_y) {G3=as.data.frame(G3)}
    }

    if (rr1==1) (G3_fin=G3)
    if (rr1>1)  (G3_fin=cbind(G3_fin,G3))
    final.time_mse <- proc.time()
    time_mse<-final.time_mse-initial.time_mse
  }
    G3_fin=diag(G3_fin)

sigma_ej<-as.numeric(sigma_j[sigma_j$iter==max(sigma_j$iter),-ncol(sigma_j)])
sigma_C<-sigma_j[sigma_j$iter==max(sigma_j$iter),-ncol(sigma_j)]

sigma_e[[ba]]<-sigma_j[sigma_j$iter==max(sigma_j$iter),-ncol(sigma_j)]
sigma_u_fin[[ba]]<-rowSums(t(phi_u)*sigma_ej) #phi in riga ha effetti in colonna modalit? y
ICC[[ba]]<-sigma_u_fin[[ba]]/(sigma_u_fin[[ba]]+sigma_ej)

phi_C<-phi_j[(phi_j$iter==max(phi_j$iter) & phi_j$fatt==n_z),-c(ncol(phi_j),(ncol(phi_j)-1))]
Diag_sigma_C<-kronecker(diag(nrow(Z_piu)),diag(sigma_C), FUN = "*")

G1<-Diag_sigma_C%*%Zr_piu_dom_kr%*%T_star%*%t(Zr_piu_dom_kr)

G2_1<-Xr_piu_kr-Zr_piu_dom_kr%*%T_star%*%t(M_xz_kr)
G2<-Diag_sigma_C%*%G2_1%*%M_x_inv_omega%*%t(G2_1)

a_1<-as.numeric(Q[1])

G4<-kronecker(diag(r_d$W_r),diag(sigma_ej), FUN = "*")

MCPE_BLUP<-G1+G2
MCPE_EBLUP<-MCPE_BLUP

  MCPE_EBLUP<-MCPE_BLUP+2*G3_fin


g1<-diag(G1)
g1<-matrix(g1,nrow=nrow(Z_piu),ncol=n_y,byrow=TRUE)

g2<-diag(G2)
g2<-matrix(g2,nrow=nrow(Z_piu),ncol=n_y,byrow=TRUE)

g3<-0
  g3<-diag(G3_fin)
  g3<-matrix(g3,nrow=nrow(Z_piu),ncol=n_y,byrow=TRUE)

g4<-diag(G4)
g4<-matrix(g4,nrow=nrow(Z_piu),ncol=n_y,byrow=TRUE)

mse_BLUP1<-diag(MCPE_BLUP)
mse_BLUP11<-matrix(mse_BLUP1,nrow=nrow(Z_piu),ncol=n_y,byrow=TRUE)
mse_BLUP[[ba]]<-matrix(mse_BLUP1,nrow=nrow(Z_piu),ncol=n_y,byrow=TRUE)

mse_EBLUP<-diag(MCPE_EBLUP)
mse_EBLUP<-matrix(mse_EBLUP,nrow=nrow(Z_piu),ncol=n_y,byrow=TRUE)

Z_piuu[[ba]]<-Z_piu
cv_BLUP[[ba]]<-100*mse_BLUP11^0.5/stima_omega_XZ_eblup
cv_EBLUP<-100*mse_EBLUP^0.5/stima_omega_XZ_eblup

cv_g1<-100*g1^0.5/stima_omega_XZ_eblup
cv_g2<-100*g2^0.5/stima_omega_XZ_eblup
cv_g3<-100*g3^0.5/stima_omega_XZ_eblup
cv_g4<-100*g4^0.5/stima_omega_XZ_eblup

n_d[[ba]]<-aggregate(as.formula(paste(weights,"~",dom,sep="")),data,sum)
colnames(n_d[[ba]])<-c("dom","nd")

appo<-list()
for (i in 1:length(y_y))
{
    formula_new<-update(formula,as.formula(paste(y_y[i],"~",
                                                 paste("factor(",x_x[x_i!=1],")",collapse="+")
                                                 ,sep="")))

  appo[[i]]<-lm(formula_new,data=data)
}
appo<-lapply(appo,modelPerformance)
names(appo)<-y_y
mod_perf[[ba]]<-appo

rm(list=setdiff(ls(), c("stima_omega_XZ_eblup_dom","stima_omega_XZ_proj_dom","stima_omega_X_proj_dom",
                        "data","data_all","macro","universe","univ_all","ba","bba","broadarea","dom","formula","id","weights",
                        "max_iter","max_diff","phi_u0","rho_u0","REML","domains2","myfun",
                        "mse_BLUP","Z_piuu","cv_BLUP","u_omegaa","n_z","z_i","n_y","n_d","y_y","z_z","x_x","n_x","x_i","beta_omega_broad",
                        "r_effect","sigma_e","sigma_u_fin","ICC","z_i_list","mod_lmer","mod_perf","int1")))
  }

stima_omega_XZ_eblup<-do.call(rbind,stima_omega_XZ_eblup_dom)
colnames(stima_omega_XZ_eblup)[2:(1+n_y)]<-y_y
stima_omega_XZ_proj<-do.call(rbind,stima_omega_XZ_proj_dom)
colnames(stima_omega_XZ_proj)[2:(1+n_y)]<-y_y
stima_omega_X_proj<-do.call(rbind,stima_omega_X_proj_dom)
colnames(stima_omega_X_proj)[2:(1+n_y)]<-y_y
mse_BLUP<-do.call(rbind,mse_BLUP)
mse_BLUP<-cbind(dom=stima_omega_XZ_eblup[,1],mse_BLUP)
colnames(mse_BLUP)[2:(1+n_y)]<-y_y

Z_piuua1<-aggregate(as.formula(paste(weights,"~",dom,sep="")),universe,sum)
cv_BLUP<-do.call(rbind,cv_BLUP)
myfun2<-function(x,y){matrix(x,nrow=length(x)/y,ncol=y,byrow=T)}
u_omegaa<-lapply(u_omegaa,myfun2,y=n_y)

cv_BLUP<-cbind(dom=stima_omega_XZ_eblup[,1],cv_BLUP)
colnames(cv_BLUP)[2:ncol(cv_BLUP)]<-paste("CV_",y_y,sep="")
Nd<-cbind(dom=Z_piuua1[,1],Nd=Z_piuua1[,2])
nd<-do.call(rbind,n_d)

for (ba in 1:length(unique(bba)))
{
appo<-list()
appo1<-z_i_list[[ba]]
z_ii<-c(1,appo1)
for (i in 1:n_z){
  appo[[i]]<-cbind(sort(unique(univ_all[univ_all[,broadarea]==bba[ba],names(z_i)[i]])),u_omegaa[[ba]][sum(z_ii[1:i]):sum(z_ii[2:(i+1)]),])
  colnames(appo[[i]]) <- c(names(z_i)[i],y_y)
}
r_effect[[ba]]<-appo
}

if(any(sort(x_i)==1)){nomi_beta<-sort(x_i)
     nomi_beta[2:length(nomi_beta)]<-nomi_beta[2:length(nomi_beta)]-1
}else{
  nomi_beta<-sort(x_i)
  nomi_beta[2:length(nomi_beta)]<-nomi_beta[2:length(nomi_beta)]-1
  }
beta_omega_broad<-do.call(cbind,beta_omega_broad)
beta_omega_broad<-data.frame(cbind(y=rep(y_y,time=length(names(fixef(mod_lmer[[1]])))),

                                   fixed=paste(rep(rep(names(nomi_beta),times=nomi_beta),each=n_y),
                                               rep(c(seq(nomi_beta[1]),unlist(apply(data.frame(nomi_beta[2:length(nomi_beta)]),1,seq))+1),each=n_y),sep=""),

                                   beta=beta_omega_broad))
colnames(beta_omega_broad)[3:ncol(beta_omega_broad)]<-paste("beta_ba_",1:ba,sep="")

names(mod_perf)<-bba

sigma_e<-do.call(rbind,sigma_e)
colnames(sigma_e)<-paste("sigma_e_",y_y)
sigma_e[,broadarea]<-1
sigma_e[,broadarea]<-unique(macro[,broadarea])

sigma_u_fin<-as.data.frame(do.call(rbind,sigma_u_fin))
colnames(sigma_u_fin)<-paste("sigma_u_",y_y)
sigma_u_fin[,broadarea]<-unique(macro[,broadarea])

ICC<-as.data.frame(do.call(rbind,ICC))
colnames(ICC)<-paste("ICC_",y_y)
ICC[,broadarea]<-unique(macro[,broadarea])

out<-list(EBLUP=stima_omega_XZ_eblup,PROJ=stima_omega_XZ_proj,SYNTH=stima_omega_X_proj,
          mse_EBLUP=as.data.frame(mse_BLUP),cv_EBLUP=as.data.frame(cv_BLUP),Nd=as.data.frame(Nd),nd=nd,r_effect=r_effect,
          beta=beta_omega_broad,mod_performance=mod_perf,sigma_e=sigma_e,sigma_u=sigma_u_fin,ICC=ICC)

return(out)

}
