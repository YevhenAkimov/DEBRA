
#' @importFrom methods setClass setGeneric setMethod setRefClass
NULL
#' @import locfit
NULL





usethis::use_package("locfit",type="Depends")
usethis::use_package("stats",type="Imports")
usethis::use_package("SummarizedExperiment",type="Imports")
usethis::use_package("DESeq2",type="Imports")
usethis::use_package("DESeq",type="Imports")
usethis::use_package("drc",type="Imports")
usethis::use_package("reshape2",type="Imports")
usethis::use_package("dgof",type="Imports")
usethis::use_package("apeglm",type="Imports")

#' @export
setClass(
     # Set the name for the class
     "DEBRA_class",
     
     # Define the slots
     slots = c(
       counts = "data.frame",
       control_names   = "character",
       condition_names   = "character",
       alpha="numeric",
       method="character",
       trended="logical",
       shrunkLFC="logical",
       modified="logical",
       overlap_data="data.frame",
       sigmoid_fit="list",
       KS_stat="data.frame",
       results="data.frame",
       results_filtered="data.frame",
       filtering_fit="list",
       filter_quantile="numeric"
       
     ),
     
     # Set the default values for the slots. (optional)
     
     prototype=list(
        alpha=-Inf,
        method="DESeq",
        trended=T,
        shrunkLFC=F,
        modified=T,
        filter_quantile=0

     ),
     
     # Make a function that can test to see if the data is consistent.
     # This is not called if you have an initialize function defined!
     validity=function(object)
     {
       
       if (sum(!( object@control_names %in% colnames(object@counts)) ) != 0 ) {
         return("control_names must be present in colnames(counts)")
       }
       
       if (sum(!( object@condition_names %in% colnames(object@counts)) ) != 0 ) {
         return("condition_names must be present in colnames(counts)")
       }
       
       if ( sum( !apply( object@counts[,c(object@control_names,object@condition_names)], 2, is.numeric) ) != 0 ) {
         return("counts must be numeric")
       }
       
      
       
       return(TRUE)
       
     }
     
     
   )
   





setGeneric("testDRB", valueClass = "DEBRA_class", function(theObject) {
  standardGeneric("testDRB")
})

setGeneric("independentFilteringDRB", valueClass = "DEBRA_class", function(theObject) {
  standardGeneric("independentFilteringDRB")
})



setGeneric("KS_plot", valueClass = "NULL", function(theObject) {
  standardGeneric("KS_plot")
})

#' @export
setMethod(f="KS_plot",signature="DEBRA_class", definition=function(theObject) {
  
  if (length(theObject@KS_stat)==0) {stop("run estimateAlpha first")}

  plot(theObject@KS_stat$mean,theObject@KS_stat$KS.p_theor,log="x",col="#343434",cex=0.6,lwd=0,pch=19,ylim=c(  min( c( theObject@KS_stat$KS.p_theor, theObject@KS_stat$KS.p ) ) ,
                                                                                               max( c( theObject@KS_stat$KS.p_theor, theObject@KS_stat$KS.p ))),
       xlab="Mean", ylab="KS statistics")
  points(theObject@KS_stat$mean,theObject@KS_stat$KS.p,col="#EF4136",pch=21,cex=0.6,lwd=1)
  legend("bottomright", legend=c("Theoretical", "Empirical"),
         col=c("#343434", "#EF4136"), cex=0.8, inset=.02,pch=c(19,21)) 
  
  invisible(NULL)  

})




setGeneric("overlapFit", valueClass = "NULL", function(theObject) {
  standardGeneric("overlapFit")
})

### overlapFit ####

#' @export
setMethod(f="overlapFit",signature="DEBRA_class", definition=function(theObject) {
  
  if (length(theObject@overlap_data)==0) {stop("run estimateAlpha first")}
  
  plot(theObject@overlap_data$mean,
       theObject@overlap_data$overlap,
       log="x",col="#343434",cex=0.6,lwd=0,pch=19,
       xlab="Mean", ylab="KS statistics")
  
   lines(theObject@sigmoid_fit$mean_fit,
         theObject@sigmoid_fit$overlap_fit,
          col="#EF4136",cex=1.6,lwd=3)

  
  invisible(NULL)  
  
})




setGeneric("filterPlot", valueClass = "NULL", function(theObject) {
  standardGeneric("filterPlot")
})



#' @export
setMethod(f="filterPlot",signature="DEBRA_class", definition=function(theObject) {

  plot(theObject@filtering_fit[[1]]$data$x,theObject@filtering_fit[[1]]$data$y,cex=0.6,pch=19,xlab="filter quantile", ylab="number of rejections")
  points(theObject@filtering_fit[[1]]$x,theObject@filtering_fit[[1]]$y,cex=0.6,pch=19,col="#EF4136")
  abline(v = theObject@filter_quantile,lwd=2)

  invisible(NULL)
})





setGeneric("resultsDRB", valueClass = "data.frame", function(theObject,filtered) {
  standardGeneric("resultsDRB")
})

#' @export
setMethod(f="resultsDRB",signature="DEBRA_class", definition=function(theObject,filtered=T) {
  
  if (filtered) {return(theObject@results_filtered)}
  else          {return(theObject@results)}
  
})





setGeneric("MAplot", valueClass = "NULL", function(theObject,FDRthr=0.25) {
  standardGeneric("MAplot")
})

#' @export
setMethod(f="MAplot",signature="DEBRA_class", definition=function(theObject,FDRthr=0.25) {
  

  if (toupper(theObject@method)=="DESEQ2(WALD)") {
    logFC="log2FoldChange"
    FDR="padj"
    pval="pvalue"
  }
  if (toupper(theObject@method)=="DESEQ2(LRT)") {
    logFC="log2FoldChange"
    FDR="padj"
    pval="pvalue"
  }
  
  if (toupper(theObject@method)=="DESEQ") {
    logFC="log2FoldChange"
    FDR="padj"
    pval="pval"
  }
  
  
  res_df_clean=theObject@results_filtered[theObject@results_filtered[,logFC]!=0 | is.finite(theObject@results_filtered[,logFC]),]
  res_df_clean=res_df_clean[!res_df_clean$baseMean==0,]
  res_df_clean=res_df_clean[is.finite(res_df_clean$baseMean),]
  
  res_df_clean$log2_mean=log2(res_df_clean$baseMean)
  
  plot(res_df_clean$log2_mean,res_df_clean[,logFC],col="#545454",cex=0.6,pch=19,
       xlab="log2(mean)", ylab="logFC")
  points(res_df_clean$log2_mean[res_df_clean[,FDR]<FDRthr],
         res_df_clean[res_df_clean[,FDR]<FDRthr,logFC],
         col="#EF4136",cex=0.4,pch=19)
  abline(h=0,col="#212121",lwd=2)
  invisible(NULL)  
  
  
})







  
setGeneric("estimateAlpha", valueClass = "DEBRA_class", function(theObject,times=2000) {
  standardGeneric("estimateAlpha")
})

#' @export
setMethod(f="estimateAlpha",signature="DEBRA_class", definition=function(theObject,times=2000) {
  
  validObject(theObject)
  
  if (theObject@alpha!= -Inf) { print("found alpha value, replacing it") }
  
  alpha=estimateAlphaFromDF(
                    counts=theObject@counts[,c(
                                               theObject@control_names,
                                               theObject@condition_names
                                               )],
                    condition_names=theObject@condition_names,
                    sampling_size=30,
                    max_window=200,
                    min_count=0,
                    min_alpha=theObject@alpha,
                    times=times,
                    cv=2,
                    n_obs=30# overlap test sample size
  )

  
  theObject@alpha<-alpha$alpha
  theObject@KS_stat<-alpha$KS_stat
  theObject@overlap_data<-alpha$overlap_data
  theObject@sigmoid_fit<-alpha$sigmoid_fit
  
  return(theObject)
      } )



#' @export
estimateAlphaFromDF=function( counts,
  condition_names, # 
  sampling_size=30,
  max_window=200,
  min_count=0,
  min_alpha=10,
  times=3000,
  cv=2,
  n_obs=30#  test sample size for overlap
) {
  
  ### Checking if mean col is present!!!!!!
  
  mean_colname=paste0(sample(letters,20,replace=T),collapse="")
  while ((mean_colname %in% colnames(counts))==T) { mean_colname=paste0(sample(letters,20,replace=T),collapse="") }


  counts=median_norm(counts)
  counts = round( data.matrix( counts[,condition_names] ) )
  counts = cbind( counts , mean=rowMeans(counts) )
  colnames(counts)[length(colnames(counts))] = mean_colname
 
  counts = counts[ order(counts[,mean_colname]), ]
  counts = counts[rowMeans(counts)>0,]
  
  sample_scale=  as.integer (1  :   (nrow(counts)-sampling_size-1) )
  
  counts=counts[sample(1:nrow(counts)),]
  counts=counts[order(counts[,mean_colname]),]
  
  res=data.frame()
  
  for (i in 1:times) {

    index=sample(sample_scale,1)
    
    if (diff(range(  counts[index:(index+sampling_size),mean_colname]  ))>max_window) {next}
    
    c_fit_data_all= c (counts [ index  :  (index+sampling_size)  , condition_names ]) 
    
    KS.p_theor=c()
    KS.p=c()
    
    for (j in 1:2) {

      size=length(c_fit_data_all)/cv
      c_fit_data=sample( c_fit_data_all,  size)
      c_mean=mean(c_fit_data)
      
      variance=pmax(var(c(c_fit_data)),c_mean*(1+1e-08))
      disp=pmax((variance - c_mean )/c_mean^2,1e-08)
      
      
      
      KS.p[j]= mean(replicate(1,dgof::ks.test( 
        c_fit_data,
        rnbinom(size,mu=c_mean,size=1/disp),exact=F)$statistic
       ))
      
      
      KS.p_theor[j]= mean(replicate(1,dgof::ks.test(
        rnbinom(size,mu=c_mean,size=1/disp),
        rnbinom(size,mu=c_mean,size=1/disp),exact=F)$statistic
       ))
      
    }
    
    res=rbind.data.frame(res,data.frame(mean=c_mean,KS.p=mean(KS.p),KS.p_theor=mean(KS.p_theor)))
    
    
  }
  #
  
  
  print("KS tests complete")
  
  ##### Produce a data.matrix with each row having boundaries of indexes for the sampling range 
  # All possible indexes are probed but only those fitting following criteria are preserved:
  # 1: minimum nuber of observations is equal to n_obs
  # 2: values of the data span the range < W but more than min_window
  # If mean_to_range_ratio is TRUE then W specifies a proportional coefficient.
  # In this case, window is calculated by dividing mean to W

  res1=res[order(res$mean),]
  
  logi_range=logical_sample(res1,sampling_col="mean",
                            n_obs=n_obs,min_range=0,mean_to_range_ratio=T,W=0.5)
  
  res1=cbind(res1,logi_range)
  
  
  cols1="KS.p"
  cols2="KS.p_theor"
  overlap=c()
  mean=c()
  n=1
  for (i in 1:nrow(res1)) {
    
    if (res1[i,"logi_range"]) {
      for (j in 1:5) {
        
        data1= sample(c( res1[i:(i+n_obs),  cols1  ]),n_obs/cv)
        data2= sample(c( res1[i:(i+n_obs),  cols2  ]),n_obs/cv)
        
        par_1=shape_scale_estim_gamma( data1 )
        par_2=shape_scale_estim_gamma( data2 )
        
        lims=c( min( c(data1,data2),na.rm = T )/10 ,  max( c(data1,data2),na.rm = T )*10 )
        
        int=try(integrate(outcome,low=lims[1],upp=lims[2],subdivisions=1000,par_1[1],par_1[2],par_2[1],par_2[2])$value,silent = T) 
        
        overlap[n]=ifelse( sum(is.na(c(par_1,par_2)))>0 | 
                             (class(int)=="try-error"),
                           NA,
                           int  )
        
        mean[n]=ifelse( sum(is.na(c(par_1,par_2)))>0 | 
                          (class(int)=="try-error"),
                        NA,
                        mean(c(c( res1[i:(i+n_obs),  "mean"  ])
                        
                        ),
                        na.rm=T)
        )
        
        n=n+1
      }
    } else {
      overlap[i]=NA
      mean[i]=NA
      n=n+1
    }
  }
  
  print("Overlap complete")
  
  overlap=(overlap[mean>min_count])
  mean=(mean[mean>min_count])
  overlap_data=data.frame(mean=mean,overlap=overlap)

  fm <- try(drc::drm ( overlap ~ (mean) ,  fct = drc::LL.4()), silent = T)
  if (class(fm)=="try-error") { print("Failed fitting LL.4() sigmoid; trying G.3()  ")
    fm <- try(drc::drm ( overlap ~ (mean) ,  fct = drc::G.3()), silent = T)
     if (class(fm)=="try-error") { 
       print("Failed to fit mean-overlap with sigmoid; run overlapFit to assess the NB fit quality and manually set alpha threshold. Using min_alpha value as a final alpha estimate. ") 
       return(list(alpha=min_alpha,overlap_data=overlap_data,sigmoid_fit=list(),KS_stat=res))
       }
  }
  
  max_x=max(mean,na.rm=T)
  min_x=min(mean,na.rm=T)
  

  
  x1=seq(min_x,max_x,length.out = 1000)
  y1=fm$curve[[1]](x1)
  
  sigmoid_fit=data.frame(mean_fit=x1,overlap_fit=y1)
  
  max_over=fm$curve[[1]](max_x)
  min_over=fm$curve[[1]](min_x)
  
  x=seq(min_x,max_x,length.out = 1000)
  fit_direction=sign(fm$curve[[1]](max_x)-fm$curve[[1]](min_x))
  

  if (fit_direction==0) {
    
    if (min_over>0.25) 
          {print("Close to negative binomial, using alpha = 0")
           return(list(alpha=0,overlap_data=overlap_data,sigmoid_fit=sigmoid_fit,KS_stat=res))} 
    else 
          {print("Failed to estimate alpha. Using min_alpha as a final alpha estimate")
           return(list(alpha=min_alpha,overlap_data=overlap_data,sigmoid_fit=sigmoid_fit,KS_stat=res))}
   }
  
  
  
  if ( max_over<0.25) {print ("Warning: Data poorly fits negative binomial, using min_alpha as a final alpha estimate")
    
    return(list(alpha=min_alpha,overlap_data=overlap_data,sigmoid_fit=sigmoid_fit,KS_stat=res))
    
    }
  
  if (min_over>0.25) {
    print("Close to negative binomial, using alpha = 0")
    return(list(alpha=0,overlap_data=overlap_data,sigmoid_fit=sigmoid_fit,KS_stat=res))
    }
 
  
    alpha = approx(x=fm$curve[[1]](x),y=x,xout=max_over*0.8)$y
  
  return( list(alpha=alpha,overlap_data=overlap_data,sigmoid_fit=sigmoid_fit,KS_stat=res))
  
  
  
}







logical_sample=function(data,sampling_col,n_obs=20,min_range=2,mean_to_range_ratio=T,W=0.1) #mean_to_range_ratio : max range value is linked to mean via ratio W=mean*W
{ # min_range minimum value of max_range
  if (is.unsorted(data[,sampling_col])) { warning("Should sampling column be sorted?") }
  
  scale=  as.integer (1  :   (nrow(data)-n_obs-1) )
  
  if (mean_to_range_ratio) 
  { 
    max_ranges=sapply(scale, function(x) 
      mean( data[ x : (x+n_obs) ,sampling_col] )*W )
    
    
  }
  else 
  {max_ranges=W}
  
  
  ranges=c( sapply(scale, function(i) 
    ifelse( (diff(range(  data[ i : (i+n_obs) ,sampling_col]  ))<= pmax(max_ranges[i],min_range) ),
            T,
            F )) 
    , rep(F,nrow(data)-rev(scale)[1])) 
  
  return(ranges)
}  




shape_scale_estim_gamma=function(c_fit_data) {
  scale=var(c_fit_data)/mean(c_fit_data)
  shape=mean(c_fit_data)/scale
  return(c(scale=scale,shape=shape))
}

median_norm=function(data){
  
  size.fac<-DESeq2::estimateSizeFactorsForMatrix(data)
  data<-t(t(data)/size.fac)
  
}






## independentFilteringDRB ####
#' @export
setMethod(f="independentFilteringDRB",signature="DEBRA_class",definition=function(theObject) {
  
  validObject(theObject)
  
  if (length(theObject@results)==0) { stop("no results found, run testDRB first") }
  if (theObject@alpha== -Inf) { print("alpha value not found, using 0") 
    theObject@alpha=0}
  
  res=pCountFilter(theObject@results,method=theObject@method, val.thr=theObject@alpha)
  theObject@results_filtered=res$result
  theObject@filtering_fit=list(res[["fit"]])
  theObject@filter_quantile=res[["filter_quantile"]]
  
  return(theObject)
} )








## testDRB ####



#' @export
setMethod(f="testDRB",signature="DEBRA_class", definition=function(theObject) {
  
  validObject(theObject)
  
  if (theObject@alpha== -Inf) { print(" no alpha value found, using alpha=0 ")
    theObject@alpha=0
  }
  
  
  if ( theObject@shrunkLFC==T & theObject@method=="DESeq") { warning("method 'DESeq' does not support shrunkLFC=T") }
  
  if (theObject@method=="DESeq2(LRT)") {
  res=deseq2_LRT(
                 theObject@counts,
                 theObject@control_names,
                 theObject@condition_names,
                 theObject@modified,
                 theObject@alpha,
                 theObject@trended,
                 theObject@shrunkLFC 
                )
  }

  if (theObject@method=="DESeq2(Wald)") {
    res=DESeq2_Wald(
      theObject@counts,
      theObject@control_names,
      theObject@condition_names,
      theObject@modified,
      theObject@alpha,
      theObject@trended,
      theObject@shrunkLFC 
    )
  }
  
  if (theObject@method=="DESeq") {
    res=DESeq_1(
      theObject@counts,
      theObject@control_names,
      theObject@condition_names,
      theObject@modified,
      theObject@alpha,
      theObject@trended
    )
  }
  
  theObject@results=res
  return(theObject)
} )







deseq2_LRT=function(counts,c_cols,disp_cols,modified,alpha,fitted_only,shrunkLFC) {

  groups=c(rep("A",length(c_cols)),rep("B",length(disp_cols)))
  groups=factor(groups)
  coldata=data.frame(condition=groups)
  design=model.matrix( ~ groups)
  
  dds=DESeq2::DESeqDataSetFromMatrix(countData = counts[,c(c_cols,disp_cols)],colData = coldata,design = ~ condition)
  dds=DESeq2::estimateSizeFactors(dds)
  dds=DESeq2::estimateDispersions(dds,fitType="local")
  
  
  if (modified==F) {
    dds <- DESeq2::nbinomLRT(dds, reduced=~1)
    
  }
  
  if (modified==T) {
    disp_est_dds=DESeq2::DESeqDataSetFromMatrix(countData = counts[,disp_cols],colData=data.frame(condition=c(disp_cols)),design = ~1)
    disp_est_dds@colData@listData[["sizeFactor"]]=dds@colData@listData[["sizeFactor"]][dds@colData@listData[["condition"]]=="B"]
    disp_est_dds=DESeq2::estimateDispersions(disp_est_dds,fitType = "local" )
    
    pseudo_baseMeans=(rowMeans(t(t(counts[,c(c_cols,disp_cols)])/DESeq2::sizeFactors(dds))))
    
    if (fitted_only==T) {
      
      cust_fit=disp_fit_deseq2(counts,disp_cols,sf=disp_est_dds@colData@listData[["sizeFactor"]])
      DESeq2::dispersions(disp_est_dds)=cust_fit$fit(pseudo_baseMeans+1)
      

    } 
    
    DESeq2::dispersions(disp_est_dds)[is.na(DESeq2::dispersions(disp_est_dds))]=max(DESeq2::dispersions(disp_est_dds),na.rm=T)
    DESeq2::dispersions(disp_est_dds)[pseudo_baseMeans<alpha]=max(DESeq2::dispersions(disp_est_dds),na.rm=T)
    
    
    DESeq2::dispersions(dds)=DESeq2::dispersions(disp_est_dds)
    
    dds <- DESeq2::nbinomLRT(dds, reduced=~1)
    
  }
  
  res=DESeq2::results(dds, independentFiltering = F)
  
  if (shrunkLFC) {res <- DESeq2::lfcShrink(dds, coef=2, res=res, type = "apeglm")}
  
  res_df=as.data.frame(res@listData)
  rownames(res_df)=res@rownames
  
  return(res=res_df)
  
  
}




DESeq_1=function(counts,c_cols,disp_cols,modified,alpha,fitted_only) {
  
  
  groups=c(rep("A",length(c_cols)),rep("B",length(disp_cols)))
  
  
  
  #### calc cust disps
  if (modified) {
    
    cds = DESeq::newCountDataSet( counts[,c(c_cols,disp_cols)], condition=groups)
    cds = DESeq::estimateSizeFactors( cds )
    cds = DESeq::estimateDispersions( cds ,fitType="local")
    
    
    sfB=cds@phenoData@data[["sizeFactor"]][cds@phenoData@data[["condition"]]=="B"]
    
    pseudo_baseMeans= base_meanV2(counts,c(c_cols,disp_cols),sf=cds@phenoData@data[["sizeFactor"]])
    
    
    if (fitted_only==T) {
      cust_fit=disp_fit_deseq2(counts,disp_cols,sf=sfB)
      
      cds@featureData@data[["disp_pooled"]]=cust_fit$fit(pseudo_baseMeans+1)
    } 
    
    else{
      
      cust_fit=deseq_disp_shrunk(counts,disp_cols,sf=sfB)
      disp=cust_fit[["disp"]]
      disp[is.na(disp)]=max(disp,na.rm =T)
      
      cds@featureData@data[["disp_pooled"]]<-disp
      
      
    }
    
    cds@featureData@data[["disp_pooled"]][pseudo_baseMeans<alpha]=max(cds@featureData@data[["disp_pooled"]],na.rm=T)
    cds@featureData@data[["disp_pooled"]][is.na(cds@featureData@data[["disp_pooled"]])]=max(cds@featureData@data[["disp_pooled"]],na.rm=T)
    
    
    res1=DESeq::nbinomTest( cds, "A", "B" )
    res1$foldChange=(res1$baseMeanB+0.15)/(res1$baseMeanA+0.15)
    res1$log2FoldChange=log2(res1$foldChange)
    rownames(res1)=res1$id
    res1=res1[,-which(colnames(res1)=="id")]
    
    
  }
  else {
    cds = DESeq::newCountDataSet( counts[,c(c_cols,disp_cols)], condition=groups)
    cds = DESeq::estimateSizeFactors( cds )
    cds = DESeq::estimateDispersions( cds ,fitType="local",method="per-condition")
    
    
    res1=DESeq::nbinomTest( cds, "A", "B" )
    res1$foldChange=(res1$baseMeanB+0.15)/(res1$baseMeanA+0.15)
    res1$log2FoldChange=log2(res1$foldChange)
    rownames(res1)=res1$id
    res1=res1[,-which(colnames(res1)=="id")]
  }
  
  return(res1)
  
  
}





DESeq2_Wald=function(counts,c_cols,disp_cols,modified,alpha,fitted_only,shrunkLFC)  {
  
  groups=c(rep("A",length(c_cols)),rep("B",length(disp_cols)))
  groups=factor(groups)
  coldata=data.frame(condition=groups)
  design=model.matrix( ~ groups)
  
  
  object=DESeq2::DESeqDataSetFromMatrix(countData = counts[,c(c_cols,disp_cols)],colData=coldata,design = design)
  object=DESeq2::estimateSizeFactors(object)
  object=DESeq2::estimateDispersions(object,fitType = "local")
  
  
  if (modified) {
    
    
    disp_est_dds=DESeq2::DESeqDataSetFromMatrix(countData = counts[,disp_cols],colData=data.frame(condition=c(disp_cols)),design = ~1)
    disp_est_dds@colData@listData[["sizeFactor"]]=object@colData@listData[["sizeFactor"]][object@colData@listData[["condition"]]=="B"]
    
    disp_est_dds=DESeq2::estimateDispersions(disp_est_dds,fitType = "local" )
    
    pseudo_baseMeans=rowMeans(t(t(counts[,c(c_cols,disp_cols)])/DESeq2::sizeFactors(object)))
    
    if (fitted_only==T) {
      cust_fit=disp_fit_deseq2(counts,disp_cols,sf=disp_est_dds@colData@listData[["sizeFactor"]])
      DESeq2::dispersions(disp_est_dds)=cust_fit$fit(pseudo_baseMeans+1)
    } 
    
    DESeq2::dispersions(disp_est_dds)[is.na(DESeq2::dispersions(disp_est_dds))]=max(DESeq2::dispersions(disp_est_dds),na.rm=T)
    DESeq2::dispersions(disp_est_dds)[pseudo_baseMeans<alpha]=max(DESeq2::dispersions(disp_est_dds),na.rm=T)
    DESeq2::dispersions(object)=DESeq2::dispersions(disp_est_dds)
    
    
  }
  
  
  object <- DESeq2::nbinomWaldTest(object,modelMatrix = design)
  res=DESeq2::results(object, independentFiltering = F)
  
  if (shrunkLFC) {res <- DESeq2::lfcShrink(object, coef=2, res=res, type = "apeglm")}
  
  res_df=as.data.frame(res@listData)
  rownames(res_df)=res@rownames
  return(res=res_df)
} 



disp_fit_deseq2=function(c_counts,disp_cols=colnames(c_counts),sf=DESeq2::estimateSizeFactorsForMatrix(c_counts)) {
  object=DESeq2::DESeqDataSetFromMatrix( round(c_counts[,disp_cols]), colData =data.frame(condition=c(disp_cols)),design=~1 )
  object@colData@listData[["sizeFactor"]] = sf
  object = DESeq2::estimateDispersions( object ,fitType="local")
  fit=object@dispersionFunction
  
 # dispData=data.frame(means = mcols(object)$baseMean, trended = mcols(object)$dispFit , shrunken=NA  )
  
  return(list(fit=fit))#,dispData=dispData
}



deseq_disp_shrunk=function(c_counts,disp_cols,sf) {
  
  object=DESeq2::DESeqDataSetFromMatrix(countData = c_counts[,c(disp_cols)],colData=data.frame(condition=c(disp_cols)),design = ~1)
  object=DESeq2::estimateSizeFactors(object)
  object@colData@listData[["sizeFactor"]] = sf
  object=DESeq2::estimateDispersions(object,fitType = "local")
  
 # dispData=data.frame(means = mcols(object)$baseMean, trended = SummarizedExperiment::mcols(object)$dispFit , shrunken=DESeq2::dispersions(object)  )
  

  
  return(list(disp=DESeq2::dispersions(object) )) #,dispData=dispData
}


outcome=function(x,scale1,shape1,scale2,shape2){
  first=dgamma(x, shape=shape1 ,  scale = scale1, log = FALSE)
  second=dgamma(x, shape=shape2 ,  scale = scale2, log = FALSE)
  return(first*(first<second)+second*(first>=second))
}


base_meanV2= function(c_counts,disp_cols=colnames(c_counts),sf){
  stopifnot(length(sf)==length(disp_cols))
  x=rowMeans(t(t(c_counts[,disp_cols])/sf))
  return(x)
} 






pCountFilter=function(res,method,val.thr=0,FDR.thr=0.2) {
  
 
  if (toupper(method)=="DESEQ2(WALD)") {
    logFC="log2FoldChange"
    FDR="padj"
    pval="pvalue"
  }
  if (toupper(method)=="DESEQ2(LRT)") {
    logFC="log2FoldChange"
    FDR="padj"
    pval="pvalue"
  }
  
  if (toupper(method)=="DESEQ") {
    logFC="log2FoldChange"
    FDR="padj"
    pval="pval"
  }
  

  theta=seq(0,1,0.01)
  
  temp_pval=res[,pval]
  temp_pval[res[,"baseMean"]<val.thr]=NA
  
  z=genefilter::filtered_p(res[,"baseMean"],temp_pval,theta=theta,method="BH")
  z1=z<FDR.thr
  x=colSums(z1,na.rm=T)
  
  df=data.frame(theta=theta,x=x)
  

  sm=smooth.spline(x=df$theta,y=df$x,df=10)
  

  quantile=df$theta[ which.max(sm$y) ]
  
  
   
  res[,FDR]=z[,which.max(sm$y)]
  
  return(list(  result=res, fit=sm , filter_quantile=quantile ))
}



#' DEBRA - DESeq-based Barcode Representation Analysis
#' 
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' cat_function()
DEBRA=function( counts ,
                control_names  ,
                condition_names ,
                alpha=-Inf,
                method="DESeq",
                trended=T,
                shrunkLFC=F,
                modified=T){
  
  drb = DEBRADataSet(counts=as.data.frame(counts),
                     control_names=control_names,
                     condition_names=condition_names,
                     alpha=alpha,
                     method=method,
                     trended=trended,
                     shrunkLFC=shrunkLFC,
                     modified=modified)
  
  
  if (drb@alpha==-Inf) {
    drb=estimateAlpha(drb)
  }
  
  
  drb=testDRB(drb)
  drb=independentFilteringDRB(drb)
  
  return(drb)
  
}





#' @export
DEBRADataSet=function( counts ,
                       control_names  ,
                       condition_names ,
                       alpha=-Inf,
                       method="DESeq",
                       trended=T,
                       shrunkLFC=F,
                       modified=T){
  
  method = match.arg(method, choices = c("DESeq","DESeq2(Wald)","DESeq2(LRT)"))
  
  drb=new("DEBRA_class",counts=as.data.frame(counts), control_names=control_names,
          condition_names=condition_names, alpha=alpha,method=method,trended=trended,shrunkLFC=shrunkLFC,modified=modified)
  
  
  return(drb)
  
}
  
  

  