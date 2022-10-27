benchSAE<-function(estim_sae,benchmark_area,area,name_dom,estimator,Nest)
{
  a1<-which(colnames(benchmark_area)==area)
  a2<-which(colnames(estim_sae)==area)
  
  sae_bench_area<-list()
  
  for(j in estimator)
  {  
    estim_sae_bench<-NULL
    for (i in unique(benchmark_area[,a1])){
      
      estim<-estim_sae[estim_sae[,a2]==i,j]
      ratio_Nd<-estim_sae[estim_sae[,a2]==i,Nest]/sum(estim_sae[estim_sae[,a2]==i,Nest])
      phi<- ratio_Nd/estim
      
      sae_bench <- estim + (1/(sum(ratio_Nd^2/phi))) * (benchmark_area[benchmark_area[,a1]==i,2] - sum(ratio_Nd *estim)) * (ratio_Nd/phi)
      
      if (j==estimator[1])
      {
        sae_bench<-data.frame(domain=estim_sae[estim_sae[,a2]==i,name_dom],  x=sae_bench)
        estim_sae_bench<-rbind(estim_sae_bench,sae_bench)
      }else{
        estim_sae_bench<-c(estim_sae_bench,sae_bench)
      }
      
    }
    sae_bench_area[[j]]<-estim_sae_bench
  }
  sae_bench_area<-do.call(cbind,sae_bench_area)
  colnames(sae_bench_area)<-c(name_dom,paste(estimator,"_bench","_",area,sep=""))
  
  return(sae_bench_area[order(sae_bench_area$d),])
}