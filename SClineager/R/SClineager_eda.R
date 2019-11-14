##------- Source from Fun_SClineager_v2.0.R: do not edit by hand


#' @name explore_sclineager
#' @title Exploratory analysis by visualization
#' @description Explore data before applying SClineager.
#' @param folders A path where each subfolder has the raw output from one \link{read_sclineager} call (for cells of one sample), which contains the ``cleaned.RData'' file.
#' @param coverage_cutoff Coverage below this cutoff is considered not observable.
#' @param file A file name (extension pdf) at which outputs (figures) are saved. For example, ``figure_eda.pdf''.
#' @return This returns nothing, but figures in \code{file}.
#' @importFrom vioplot vioplot
explore_sclineager<-function(folders,coverage_cutoff,file)
{
  categories=c("chrM","frameshift substitution","nonframeshift substitution",
               "nonsynonymous SNV","synonymous SNV","splicing","UTR5","UTR3",
               "stopgain","stoploss","other","unknown")
  vaf_all=list() 
  count_all=c()
  sd_all=list()
  missing_all=c()
  coverage_all=c()
  mutations_all=c()
  type_all=c()
  
  for (folder in folders)
  {
    # load data
    load(paste(folder,"/cleaned.RData",sep=""))
    runinfo=results$runinfo
    mutations_mat=results$mutations_mat
    annotation=results$annotation
    coverage_mat=results$coverage_mat
    annotation=results$annotation
    
    # vaf, sd, and count of mutations
    vaf_all[[folder]]=aggregate(apply(mutations_mat,1,function(x) {x=x[!is.na(x)];mean(x)}),
                                by=list(factor(annotation$ExonicFunc.refGene)),mean)
    sd_all[[folder]]=aggregate(apply(mutations_mat,1,function(x) {x=x[!is.na(x)];sd(x)}),
                               by=list(factor(annotation$ExonicFunc.refGene)),mean)
    count_all=rbind(count_all,cbind(annotation,id=folder))
    
    # proportion of missing data points
    a=is.na(coverage_mat[grepl("chrM",rownames(coverage_mat)),])
    b=is.na(coverage_mat[!grepl("chrM",rownames(coverage_mat)),])
    missing_all=rbind(missing_all,c(sum(a)/length(a),sum(b)/length(b)))
    
    # allelic bias
    keep=apply(mutations_mat,1,function(x) max(x,na.rm = T)-min(x,na.rm = T)>0.3)
    tmp=t(t(coverage_mat)/apply(coverage_mat,2,function(x) median(x,na.rm=T)))
    coverage_all=c(coverage_all,as.vector(tmp[keep,]))
    mutations_all=c(mutations_all,as.vector(mutations_mat[keep,]))
    type_all=c(type_all,rep(annotation$ExonicFunc.refGene[keep],dim(mutations_mat)[2]))
  }
  
  # plot
  pdf(file,height=9,width=5)
  
  # vaf, sd, and count of mutations
  par(mfrow=c(3,1),mar=c(8,8,8,8))
  boxplot(as.vector(unlist(sapply(vaf_all,function(x) x[,2])))~
            factor(unlist(sapply(vaf_all,function(x) x[,1])),levels=categories),ylab="Average of VAF",
          xaxt="n",col="cyan",pch=19,ylim=c(0,1))
  text(1:length(categories),par("usr")[3]-0.1,srt=45,adj=1,xpd=TRUE,labels=categories,cex=1)
  
  boxplot(as.vector(unlist(sapply(sd_all,function(x) x[,2])))~
            factor(unlist(sapply(sd_all,function(x) x[,1])),levels=categories),ylab="Variation of VAF",
          xaxt="n",col="red",pch=19)
  text(1:length(categories),par("usr")[3]-0.01,srt=45,adj=1,xpd=TRUE,labels=categories,cex=1)
  
  tmp=as.matrix(table(count_all$id,count_all$ExonicFunc.refGene))
  boxplot(as.vector(tmp)~factor(rep(colnames(tmp),each=dim(tmp)[1]),levels=categories),xaxt="n",
          ylab="Count",col="orange",pch=19)
  text(1:length(categories),par("usr")[3]-2,srt=45,adj=1,xpd=TRUE,labels=categories,cex=1)
  
  # proportion of missing data points
  par(mfrow=c(3,2),mar=c(4,4,4,4))
  plot(density(missing_all[,1]),col="red",xlim=c(0,0.5),lwd=2,xlab="Missing%",ylab="Density")
  lines(density(missing_all[,2]),col="blue",lwd=2)
  
  # allelic bias
  par(mfrow=c(4,1),mar=c(2,2,2,2))
  keep=type_all %in% c("chrM","nonsynonymous SNV","synonymous SNV","UTR5","UTR3","other")
  tmp=data.frame(x=coverage_all[keep],y=abs(0.5-mutations_all)[keep],type=type_all[keep],
                 stringsAsFactors = F)
  tmp=tmp[(!is.na(tmp$x)),]
  tmp$type[tmp$type!="chrM"]="non-chrM"
  
  for (group in c("chrM","non-chrM"))
  {
    tmp1=tmp[tmp$type %in% group,]
    vioplot(tmp1$y[tmp1$x>quantile(tmp1$x,0.9)],
            tmp1$y[tmp1$x<=quantile(tmp1$x,0.9) & tmp1$x>quantile(tmp1$x,0.8)],
            tmp1$y[tmp1$x<=quantile(tmp1$x,0.8) & tmp1$x>quantile(tmp1$x,0.6)],
            tmp1$y[tmp1$x<=quantile(tmp1$x,0.6) & tmp1$x>quantile(tmp1$x,0.3)],
            tmp1$y[tmp1$x<=quantile(tmp1$x,0.3)],names=c(">90%","80%-90%","60%-80%","30%-60%","<30%"))
    title(main=group)
  }
  
  dev.off()
  
  invisible()
}


