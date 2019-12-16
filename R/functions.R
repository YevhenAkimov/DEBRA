#' @importFrom methods setClass setGeneric setMethod setRefClass
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

  "DEBRADataSet",

  slots = c(
    counts = "data.frame",
    control_names   = "character",
    condition_names   = "character",
    beta="numeric",
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
    filter_quantile="numeric",
    filter_FDR= "numeric",
    default_beta = "numeric"),

  prototype=list(
    beta=-Inf,
    method="DESeq",
    trended=T,
    shrunkLFC=F,
    modified=T,
    filter_quantile=0,
    filter_FDR=0.2,
    default_beta=0),

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


#' Statistical significane test for differential barcode representation using parameters provided in DEBRADataSet
#'
#' @param drb DEBRADataSet
#' @examples drb<-testDRB(drb)
#' @export
setGeneric("testDRB", valueClass = "DEBRADataSet", function(drb) {
  standardGeneric("testDRB")
})

#' @export
setMethod(f="testDRB",signature="DEBRADataSet", definition=function(drb) {

  validObject(drb)

  if (drb@beta== -Inf) { print(" no beta value found, using beta = 0 ")
    drb@beta=0
  }

  if ( drb@shrunkLFC==T & drb@method=="DESeq") { warning("method 'DESeq' does not support shrunkLFC=T") }

  if (toupper(drb@method)=="DESEQ2(LRT)") {
    res=deseq2_LRT(
      drb@counts,
      drb@control_names,
      drb@condition_names,
      drb@modified,
      drb@beta,
      drb@trended,
      drb@shrunkLFC
    )
  }

  if (toupper(drb@method)=="DESEQ2(WALD)") {
    res=DESeq2_Wald(
      drb@counts,
      drb@control_names,
      drb@condition_names,
      drb@modified,
      drb@beta,
      drb@trended,
      drb@shrunkLFC
    )
  }

  if (toupper(drb@method)=="DESEQ") {
    res=DESeq_1(
      drb@counts,
      drb@control_names,
      drb@condition_names,
      drb@modified,
      drb@beta,
      drb@trended
    )
  }

  drb@results=res
  return(drb)
} )


#' Plot theoretical and empirical Kolmogorov-Smirnov (KS) statistic values
#'
#' @description  Plot local estimates of theoretical and empirical Kolmogorov-Smirnov test statistic values used for beta threshold estimation
#' @param drb DEBRADataSet object
#' @examples KS_plot(drb)
#' @export
setGeneric("KS_plot",  function(drb) {
  standardGeneric("KS_plot")
})


#' @export
setMethod(f="KS_plot",signature="DEBRADataSet", definition=function(drb) {

  if (length(drb@KS_stat)==0) {stop("run estimateBeta first")}

  plot(drb@KS_stat$mean,drb@KS_stat$KS.p_theor,log="x",col="#343434",cex=0.6,lwd=0,pch=19,ylim=c(  min( c( drb@KS_stat$KS.p_theor, drb@KS_stat$KS.p ) ) ,
                                                                                                   max( c( drb@KS_stat$KS.p_theor, drb@KS_stat$KS.p ))),
       xlab="Mean", ylab="KS statistics")

  points(drb@KS_stat$mean,drb@KS_stat$KS.p,col="#EF4136",pch=21,cex=0.6,lwd=1)
  legend("bottomright", legend=c("Theoretical", "Empirical"),
         col=c("#343434", "#EF4136"), cex=0.8, inset=.02,pch=c(19,21))



})


#' Plot mean count versus overlap values between theoretical and empirical Kolmogorov-Smirnov (KS) test statistics and the sigmoid fit
#'
#' @description This function plots the local overlap between Gamma-modelled empirical and theoretical Kolmogorov-Smirnov statistic values fitted with four parameters sigmoid function.
#' @param drb DEBRADataSet object
#' @examples overlapFit(drb)
#' @export
setGeneric("overlapFit", function(drb) {
  standardGeneric("overlapFit")
})


#' @export
setMethod(f="overlapFit",signature="DEBRADataSet", definition=function(drb) {

  if (length(drb@overlap_data)==0) {stop("run estimateBeta first")}

  plot(drb@overlap_data$mean,
       drb@overlap_data$overlap,
       log="x",col="#343434",cex=0.6,lwd=0,pch=19,
       xlab="Mean", ylab="KS statistics")

  lines(drb@sigmoid_fit$mean_fit,
        drb@sigmoid_fit$overlap_fit,
        col="#EF4136",cex=1.6,lwd=3)

})


#' Visulaize independet filtering results
#'
#' @description Visualize optimization of the number of rejections over the quantiles of the mean counts.
#' @param drb DEBRADataSet object
#' @examples filterPlot(drb)
#' @export
setGeneric("filterPlot", function(drb) {
  standardGeneric("filterPlot")
})


#' @export
setMethod(f="filterPlot",signature="DEBRADataSet", definition=function(drb) {

  plot(drb@filtering_fit[[1]]$data$x,drb@filtering_fit[[1]]$data$y,cex=0.6,pch=19,xlab="filter quantile", ylab="number of rejections")
  points(drb@filtering_fit[[1]]$x,drb@filtering_fit[[1]]$y,cex=0.6,pch=19,col="#EF4136")
  abline(v = drb@filter_quantile,lwd=2)

})


#' Extracts a result table from DEBRADataSet object
#'
#' @description resultsDRB extracts a DEBRA analysis results from a DEBRADataSet object.
#' @param drb DEBRADataSet object
#' @param filtered logical specifying if filtered results should be extracted; default is TRUE
#' @examples resultsDRB(drb, filtered = T)
#' @export
setGeneric("resultsDRB", valueClass = "data.frame", function(drb,filtered=T,sort=T) {
  standardGeneric("resultsDRB")
})

#' @export
setMethod(f="resultsDRB",signature="DEBRADataSet", definition=function(drb,filtered=T,sort=T) {

  if (filtered) {
                   if (sort) {
                    res=drb@results_filtered[order(drb@results_filtered$padj),]
                   }

                   return(res)
                 }

        else   {
                  if (sort) {
                    res=drb@results[order(drb@results$padj),]
                  }

                  return(res)
                }

})




#' MA-plot from count means and log2 fold changes
#' @description A function that creates MA-plot from DEBRA analysis results
#' @param drb DEBRADataSet object
#' @param FDR numeric specifying FDR level. Barcodes above it will be colored in red. Default is 0.25
#' @param filtered logical specifying if filtered results will be used to generate MA-plot
#' @examples MAplot(drb, FDR = 0.25)
#' @export
setGeneric("MAplot", function(drb,FDR=0.25,filtered=T) {
  standardGeneric("MAplot")
})

#' @export
setMethod(f="MAplot",signature="DEBRADataSet", definition=function(drb,FDR=0.25,filtered=T) {

  # set variables specifiyng colnames according to the used method

  if (toupper(drb@method)=="DESEQ2(WALD)") {
    logFC="log2FoldChange"
    FDR_name="padj"
    pval="pvalue"
  }
  if (toupper(drb@method)=="DESEQ2(LRT)") {
    logFC="log2FoldChange"
    FDR_name="padj"
    pval="pvalue"
  }
  if (toupper(drb@method)=="DESEQ") {
    logFC="log2FoldChange"
    FDR_name="padj"
    pval="pval"
  }

  # Select data and remove barcodes with logFC = 0 and non-finite LogFC
  if (filtered) {
    res_df_clean=drb@results_filtered[ drb@results_filtered[,logFC]!=0 | is.finite(drb@results_filtered[,logFC]), ]
  } else {
    res_df_clean=drb@results[ drb@results[,logFC]!=0 | is.finite(drb@results[,logFC]), ]
  }

  # remove barcodes with zero counts and non-finite counts
  res_df_clean=res_df_clean[!res_df_clean$baseMean==0,]
  res_df_clean=res_df_clean[is.finite(res_df_clean$baseMean),]

  # convert counts to log2
  log2_mean=log2(res_df_clean$baseMean)

  plot(log2_mean,res_df_clean[,logFC],col="#545454",cex=0.6,pch=19,
       xlab="log2(mean)", ylab="logFC",main=as.character(FDR))

  points(log2_mean[res_df_clean[,FDR_name]<FDR],
         res_df_clean[res_df_clean[,FDR_name]<FDR,logFC],
         col="#EF4136",cex=0.4,pch=19)

  abline(h=0,col="#212121",lwd=2)

})

#' Independent filtering
#' @description Independent filtering step for DEBRA analysis results
#' @param drb DEBRADataSet object
#' @param filter_FDR user-specified FDR threshold for calculating number of null hypothesis rejections
#' @param beta beta threshold provides a lower read count value for an independent filtering step; if not provided, beta value is extracted from DEBRADataSet
#' @examples drb<-independentFilteringDRB(drb)
#' @export
setGeneric("independentFilteringDRB", valueClass = "DEBRADataSet", function(drb,filter_FDR=NULL,beta = NULL) {
  standardGeneric("independentFilteringDRB")
})


#' @export
setMethod(f="independentFilteringDRB",signature="DEBRADataSet",definition=function(drb,filter_FDR=NULL,beta = NULL) {

  validObject(drb)

  if (length(drb@results)==0) { stop("no results found, run testDRB first") }
  if (drb@beta== -Inf) { print("beta value not found, using 0")
    drb@beta=0}

  res=pCountFilter(drb@results,method=drb@method, val.thr=drb@beta,FDR.thr=filter_FDR)

  drb@results_filtered=res$result
  drb@filtering_fit=list(res[["fit"]])
  drb@filter_quantile=res[["filter_quantile"]]

  return(drb)
} )


#' Estimates the beta threshold for DEBRADataSet
#'
#' @param drb DEBRADataSet object
#' @param s an integer specifying number of replicates for local Kolmogorov-Smirnov statistic estimation. Default is 1500
#' @param obs_window an integer specifying number of observations sampled from ordered mean counts values for theoretical and empirical Kolmogorov-Smirnov statistics estimation. Mean counts are calculated from condition(test) samples
#' @param max_window a numeric specifying maximum size of the sampling window in read counts max(sampled_counts)-min(sampled_counts) < max_window
#' @param KS_test_window an integer specifying number of observations sampled from theoretical and empirical Kolmogorov-Smirnov statistics
#' @param default_beta a numeric specifying beta value used if the estimation has failed
#' @examples drb<-estimateBeta(drb)
#' @export
setGeneric("estimateBeta", valueClass = "DEBRADataSet", function( drb, s = 1500, obs_window = 60, max_window = 150,
                                                                  KS_test_window = 40, default_beta=10 ) {
  standardGeneric("estimateBeta")
})



setMethod(f="estimateBeta",signature="DEBRADataSet", definition=function( drb, s = 1500, obs_window = 60, max_window = 150,
                                                                          KS_test_window = 40, default_beta=10 ) {

  validObject(drb)

  if (drb@beta != -Inf) { print("found beta value, replacing it") }
  if (is.numeric(drb@default_beta)) { drb@default_beta=default_beta }

  beta=estimateBetaFromDF(
    counts=drb@counts[,c(
      drb@condition_names
    )],
    sampling_size=obs_window,
    max_window=max_window,
    default_beta=drb@default_beta,
    times=s,
    n_obs=KS_test_window # overlap test sample size
  )


  drb@beta<-beta$beta
  drb@KS_stat<-beta$KS_stat
  drb@overlap_data<-beta$overlap_data
  drb@sigmoid_fit<-beta$sigmoid_fit

  return(drb)
} )



estimateBetaFromDF=function( counts,

                             sampling_size=60,
                             max_window=150,
                             min_count=1,
                             default_beta=default_beta,
                             times=1500,
                             n_obs=50#  test sample size for overlap
) {


  counts= round( data.matrix(median_norm(counts))) # median normalization

  # remove zero read counts rows
  counts = counts[rowMeans(counts)>0,]



  # generate row indexes for sampling
  sample_scale =  as.integer (1  :   (nrow(counts)-sampling_size-1) )

  counts=counts[sample(1:nrow(counts)),]


  # sort for means
  counts = counts[ order(rowMeans(counts)), ]
  mean=rowMeans(counts)

  res=data.frame()

  for (i in 1:times) { # start sampling

    index=sample(sample_scale,1) # get random index

    # check if the window is less than max_window
    if (diff (range(  mean[ index:(index+sampling_size)]  )) >max_window) {next}

    #extract read counts for testing
    c_fit_data_all= c (counts [ index  :  (index+sampling_size),  ])

    # create empty vectors
    KS.p_theor=c()
    KS.p=c()



    # sampling inside the window
    size=length(c_fit_data_all)/2
    c_fit_data=sample( c_fit_data_all, size)

    # finding mean and dispersion
    c_mean=mean(c_fit_data)
    variance=pmax(var(c(c_fit_data)),c_mean*(1+1e-08))
    disp=pmax((variance - c_mean )/c_mean^2,1e-08)


    # estimate KS statistics
    KS.p= mean(replicate(10,dgof::ks.test(
      c_fit_data,
      rnbinom(size,mu=c_mean,size=1/disp),exact=F)$statistic
    ))


    KS.p_theor=  mean(replicate(10,dgof::ks.test(
      rnbinom(size,mu=c_mean,size=1/disp),
      rnbinom(size,mu=c_mean,size=1/disp),exact=F)$statistic
    ))


    # combine results
    res=rbind.data.frame(res,data.frame(mean=c_mean,KS.p=KS.p,KS.p_theor=KS.p_theor))


  }

  print("KS tests completed",quote = FALSE)

  res1=res[order(res$mean),] # ordering results

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

        data1= sample(c( res1[i:(i+n_obs),  cols1  ]),n_obs/2)
        data2= sample(c( res1[i:(i+n_obs),  cols2  ]),n_obs/2)

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

  print("Overlap estimation completed")

  overlap=(overlap[mean>min_count])
  mean=(mean[mean>min_count])
  overlap_data=data.frame(mean=mean,overlap=overlap)

  fm <- try(drc::drm ( overlap ~ (mean) ,  fct = drc::LL.4()), silent = T)
  if (class(fm)=="try-error") { print("Failed fitting with LL.4() sigmoid; trying G.3() (see drc::drm for details)")
    fm <- try(drc::drm ( overlap ~ (mean) ,  fct = drc::G.3()), silent = T)
    if (class(fm)=="try-error") {
      print("Failed to fit mean-KS overlap with sigmoid; run KS_plot(drb) to assess the NB fit quality and manually set beta threshold. Setting beta to default_beta value. ")
      return(list(beta=default_beta,overlap_data=overlap_data,sigmoid_fit=list(),KS_stat=res))
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
    {print("Close to negative binomial, using beta = 0")
      return(list(beta=0,overlap_data=overlap_data,sigmoid_fit=sigmoid_fit,KS_stat=res))}
    else
    {print("Failed to estimate beta. Using default_beta as a final beta estimate")
      return(list(beta=default_beta,overlap_data=overlap_data,sigmoid_fit=sigmoid_fit,KS_stat=res))}
  }



  if ( max_over<0.25) {print ("Warning: Data poorly fits negative binomial, using default_beta as a final beta estimate")

    return(list(beta=default_beta,overlap_data=overlap_data,sigmoid_fit=sigmoid_fit,KS_stat=res))

  }

  if (min_over>0.25) {
    print("Close to negative binomial, using beta = 0")
    return(list(beta=0,overlap_data=overlap_data,sigmoid_fit=sigmoid_fit,KS_stat=res))
  }


  beta = approx(x=fm$curve[[1]](x),y=x,xout=max_over*0.8)$y

  return( list(beta=beta,overlap_data=overlap_data,sigmoid_fit=sigmoid_fit,KS_stat=res))



}







logical_sample=function(data,sampling_col,n_obs=20,min_range=2,mean_to_range_ratio=T,W=0.1) #mean_to_range_ratio : max range value is linked to mean via ratio W=mean*W
{
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



deseq2_LRT=function(counts,c_cols,disp_cols,modified,beta,fitted_only,shrunkLFC) {

  groups=c(rep("A",length(c_cols)),rep("B",length(disp_cols)))
  groups=factor(groups)
  coldata=data.frame(condition=groups)
  design=model.matrix( ~ groups)

  dds=DESeq2::DESeqDataSetFromMatrix(countData = counts[,c(c_cols,disp_cols)],colData = coldata,design = ~ condition)
  dds=DESeq2::estimateSizeFactors(dds)
  dds=DESeq2::estimateDispersions(dds,fitType="local", quiet=T)


  if (modified==F) {
    dds <- DESeq2::nbinomLRT(dds, reduced=~1)

  }

  if (modified==T) {
    disp_est_dds=DESeq2::DESeqDataSetFromMatrix(countData = counts[,disp_cols],colData=data.frame(condition=c(disp_cols)),design = ~1)
    disp_est_dds@colData@listData[["sizeFactor"]]=dds@colData@listData[["sizeFactor"]][dds@colData@listData[["condition"]]=="B"]
    disp_est_dds=DESeq2::estimateDispersions(disp_est_dds,fitType = "local" , quiet=T)

    pseudo_baseMeans=rowMeans(t(t(counts[,c(c_cols,disp_cols)])/DESeq2::sizeFactors(dds)))

    if (fitted_only==T) {

      cust_fit=disp_fit_deseq2(counts,disp_cols,sf=disp_est_dds@colData@listData[["sizeFactor"]])
      DESeq2::dispersions(disp_est_dds)=cust_fit$fit(pseudo_baseMeans+1)


    }

    DESeq2::dispersions(disp_est_dds)[is.na(DESeq2::dispersions(disp_est_dds))]=max(DESeq2::dispersions(disp_est_dds),na.rm=T)
    DESeq2::dispersions(disp_est_dds)[pseudo_baseMeans<beta]=max(DESeq2::dispersions(disp_est_dds),na.rm=T)


    DESeq2::dispersions(dds)=DESeq2::dispersions(disp_est_dds)

    dds <- DESeq2::nbinomLRT(dds, reduced=~1)

  }

  if(modified==F) {
    res=DESeq2::results(dds)
  } else {
    res=DESeq2::results(dds, independentFiltering = F)
  }
  

  if (shrunkLFC) {res <- DESeq2::lfcShrink(dds, coef=2, res=res, type = "apeglm")}

  res_df=as.data.frame(res@listData)
  rownames(res_df)=res@rownames

  return(res=res_df)


}




DESeq_1=function(counts,c_cols,disp_cols,modified,beta,fitted_only) {


  groups=c(rep("A",length(c_cols)),rep("B",length(disp_cols)))


  if (modified) {

    cds = DESeq::newCountDataSet( counts[,c(c_cols,disp_cols)], condition=groups)
    cds = DESeq::estimateSizeFactors( cds )
    cds = DESeq::estimateDispersions( cds ,fitType="local")

    # sfB=cds@phenoData@data[["sizeFactor"]][cds@phenoData@data[["condition"]]=="B"]
    # sfA=cds@phenoData@data[["sizeFactor"]][cds@phenoData@data[["condition"]]=="A"]
    #
    # pseudo_baseMeans= base_meanV2(counts,c(c_cols), sf=sfA ) #,disp_cols

    sfB=cds@phenoData@data[["sizeFactor"]][cds@phenoData@data[["condition"]]=="B"]
    pseudo_baseMeans= base_meanV2(counts,c(c_cols,disp_cols), sf= cds@phenoData@data[["sizeFactor"]] )


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

    cds@featureData@data[["disp_pooled"]][pseudo_baseMeans<beta]=max(cds@featureData@data[["disp_pooled"]],na.rm=T)
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


DESeq2_Wald=function(counts,c_cols,disp_cols,modified,beta,fitted_only,shrunkLFC)  {

  groups=c(rep("A",length(c_cols)),rep("B",length(disp_cols)))
  groups=factor(groups)
  coldata=data.frame(condition=groups)
  design=model.matrix( ~ groups)


  object=DESeq2::DESeqDataSetFromMatrix(countData = counts[,c(c_cols,disp_cols)],colData=coldata,design = design)
  object=DESeq2::estimateSizeFactors(object)
  object=DESeq2::estimateDispersions(object,fitType = "local" , quiet=T)


  if (modified) {


    disp_est_dds=DESeq2::DESeqDataSetFromMatrix(countData = counts[,disp_cols],colData=data.frame(condition=c(disp_cols)),design = ~1)
    disp_est_dds@colData@listData[["sizeFactor"]]=object@colData@listData[["sizeFactor"]][object@colData@listData[["condition"]]=="B"]

    disp_est_dds=DESeq2::estimateDispersions(disp_est_dds,fitType = "local" , quiet=T)

    pseudo_baseMeans=rowMeans(t(t(counts[,c(c_cols,disp_cols)])/DESeq2::sizeFactors(object)))

    if (fitted_only==T) {
      cust_fit=disp_fit_deseq2(counts,disp_cols,sf=disp_est_dds@colData@listData[["sizeFactor"]])
      DESeq2::dispersions(disp_est_dds)=cust_fit$fit(pseudo_baseMeans+1)
    }

    DESeq2::dispersions(disp_est_dds)[is.na(DESeq2::dispersions(disp_est_dds))]=max(DESeq2::dispersions(disp_est_dds),na.rm=T)
    DESeq2::dispersions(disp_est_dds)[pseudo_baseMeans<beta]=max(DESeq2::dispersions(disp_est_dds),na.rm=T)
    DESeq2::dispersions(object)=DESeq2::dispersions(disp_est_dds)


  }


  object <- DESeq2::nbinomWaldTest(object,modelMatrix = design)
  
  if(modified==F) {
    res=DESeq2::results(object)
  } else {
    res=DESeq2::results(object, independentFiltering = F)
  }
  

  if (shrunkLFC) {res <- DESeq2::lfcShrink(object, coef=2, res=res, type = "apeglm")}

  res_df=as.data.frame(res@listData)
  rownames(res_df)=res@rownames
  return(res=res_df)
}



disp_fit_deseq2=function(c_counts,disp_cols=colnames(c_counts),sf=DESeq2::estimateSizeFactorsForMatrix(c_counts)) {
  object=DESeq2::DESeqDataSetFromMatrix( round(c_counts[,disp_cols]), colData =data.frame(condition=c(disp_cols)),design=~1 )
  object@colData@listData[["sizeFactor"]] = sf
  object = DESeq2::estimateDispersions( object ,fitType="local", quiet=T)
  fit=object@dispersionFunction

  # dispData=data.frame(means = mcols(object)$baseMean, trended = mcols(object)$dispFit , shrunken=NA  )

  return(list(fit=fit))#,dispData=dispData
}



deseq_disp_shrunk=function(c_counts,disp_cols,sf) {

  object=DESeq2::DESeqDataSetFromMatrix(countData = c_counts[,c(disp_cols)],colData=data.frame(condition=c(disp_cols)),design = ~1)
  object=DESeq2::estimateSizeFactors(object)
  object@colData@listData[["sizeFactor"]] = sf
  object=DESeq2::estimateDispersions(object,fitType = "local", quiet=T)


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


  theta=seq(0,1,0.002)

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
#' @param counts a data frame of non-negative read counts with columns of samples and rownames of barcode IDs; samples not included into analysis are allowed
#' @param control_names a character vector specifying the control samples (colnames of the counts data frame)
#' @param condition_names a character vector specifying the condition samples (colnames of the counts data frame)
#' @param beta a numeric specifying beta value providing a lower read count threshold value for an independent filtering step; if beta = -Inf (default), the beta will be estimated from the read counts of condition samples
#' @param method a character specifying the method used for inferring differentially represented barcodes
#' @param trended a logical specifying if the trended dispersion estimates should be used; if trended=FALSE, the shrunken dispersion estimates (as estimated by DESeq2) are used
#' @param default_beta a numeric specifying the beta value used if the beta estimation is failed
#' @param shrinkLFC a logical specifying if the logFC values should be shrunken using "apeglm" shrinkage estimator
#' @param modified a logical, if modified = F then the non-modified version of the correspondig method (DESeq, DESeq2(Wald) or DESeq2(LRT)) will be run; note independent filtering using beta threshold can still be applied
#' @export
DEBRA=function( counts ,
                control_names  ,
                condition_names ,
                beta=-Inf,
                method="DESeq",
                trended=T,
                filter_FDR=0.2,
                default_beta = 0,
                shrunkLFC=F,
                modified=T){

  drb = DEBRADataSet(counts=as.data.frame(counts),
                     control_names=control_names,
                     condition_names=condition_names,
                     beta=beta,
                     method=method,
                     trended=trended,
                     shrunkLFC=shrunkLFC,
                     modified=modified)


  if (drb@beta==-Inf & drb@modified==T) {
    drb=estimateBeta(drb)
  }

  if (is.numeric(filter_FDR)) {drb@filter_FDR=filter_FDR}

  drb=testDRB(drb)
  
  if(drb@modified==T) {
    drb=independentFilteringDRB(drb,filter_FDR=drb@filter_FDR)
  } else {
    drb@results_filtered=drb@results
  }
  
    

  return(drb)

}




#' Construct a DEBRADataSet
#'
#' @param counts a data frame of non-negative read counts with columns of samples and rownames of barcode IDs; samples not included into analysis are allowed
#' @param control_names a character vector specifying the control samples (colnames of the counts data frame)
#' @param condition_names a character vector specifying the condition samples (colnames of the counts data frame)
#' @param beta a numeric specifying beta value providing a lower read count threshold value for an independent filtering step; if beta = -Inf (default), the beta will be estimated from the read counts of condition samples
#' @param method a character specifying the method used for inferring differentially represented barcodes
#' @param trended a logical specifying if the trended dispersion estimates should be used; if trended=FALSE, the shrunken dispersion estimates (as estimated by DESeq2) are used
#' @param default_beta a numeric specifying the beta value used if the beta estimation is failed
#' @param shrinkLFC a logical specifying if the logFC values should be shrunken using "apeglm" shrinkage estimator
#' @param modified a logical, if modified = F then the non-modified version of the correspondig method (DESeq, DESeq2(Wald) or DESeq2(LRT)) will be run; note independent filtering using beta threshold can still be applied
#' @export
DEBRADataSet = function( counts ,
                         control_names  ,
                         condition_names ,
                         method=c("DESeq","DESeq2(Wald)","DESeq2(LRT)"),
                         beta=-Inf,
                         trended=T,
                         shrunkLFC=F,
                         modified=T,
                         default_beta=10) {

  method = match.arg(method, choices = c("DESeq","DESeq2(Wald)","DESeq2(LRT)"))

  drb=new("DEBRADataSet",counts=as.data.frame(counts), control_names=control_names,
          condition_names=condition_names, beta=beta,method=method,trended=trended,
          shrunkLFC=shrunkLFC,modified=modified,default_beta=default_beta)

  return(drb)

}

#' DATA
#'
#' @docType data
#' @usage data(bar_counts)
"bar_counts"
