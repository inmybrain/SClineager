##------- Source from Fun_SClineager_v2.0.R: do not edit by hand


#' @name read_sclineager
#' @title Read datsets for SClineager
#' @description 
#' @param runinfo A data frame with ``Path'' and ``Cell''. ``Path'' is the path to the raw mutation output folder for each cell. ``Cell'' is the name of each single cell.
#' @param coverage_cutoff Coverage below this cutoff is considered not observable.
#' @param coverage_percentage A numeric value in [0,1]. If a mutation has more than \code{coverage_percentage} of cells whose coverages are below \code{coverage_cutoff}, then it will not considered in modeling.
#' @param cell_percentage A numeric value in [0,1]. If a cell has more than \code{cell_percentage} of mutations not observable, then it will not considered in modeling.
#' @param folder A folder at which outputs are saved.
#' @param artefact_percentage A numeric value in [0,1]. Only keep mutations that appeared in more than \code{artefact_percentage} of cells.
#' @param species A character value that indicates germline_mutations. Default is ``hg38'' and other possibilities are ``hg19'' and ``mm10''.
#' @param keep_other A logical value whether to keep intronic and intergenic mutations. Default is \code{FALSE}.
#' @return ``summary.pdf'' and ``results.RData'' are saved under \code{folder}. The RData file contains a list with components (also returned):
#' \item{mutations_mat}{The variant allele frequency (VAF) matrix, whose rows are mutations and columns are samples.}
#' \item{coverage_mat}{The total sequencing coverage matrix, whose rows are mutations and columns are samples.}
#' \item{runinfo}{Same as \code{runinfo} in the input, but trimmed to contain only cells that appeared in \code{mutations_mat}, in the same order.}
#' \item{annotation}{Annotation information of the variants in \code{mutations_mat}.}
#' @examples
#' # See the github page https://github.com/inmybrain/SClineager.
#' @importFrom gplots heatmap.2
#' @export
read_sclineager <-
  function(runinfo,
           coverage_cutoff,
           coverage_percentage,
           cell_percentage,
           folder,
           artefact_percentage,
           species = "hg38",
           keep_other = FALSE){
    # read mutations
  mutations=c()
  coverages=list()
  
  for (i in 1:dim(runinfo)[1])
  {
    tryCatch({
      print(runinfo$Path[i])
      if (!file.exists(runinfo$Path[i])) {next}
      
      file=paste(runinfo$Path[i],"/coverage.txt",sep="")
      if ((!file.exists(file))) {next}
      tmp=read.table(file,header=F,stringsAsFactors = F,sep="\t",comment="")[,c(1,2,4)]
      tmp=tmp[tmp[,3]>coverage_cutoff,]
      coverages[[runinfo$Cell[i]]]=tmp
      
      file=paste(runinfo$Path[i],"/germline_mutations_",species,".txt",sep="")
      if ((!file.exists(file))) {next}
      mutation=read.table(file,header=T,stringsAsFactors = F,sep="\t")
      mutation$id=runinfo$Cell[i]
      mutation$mutation=paste(mutation$Chr,mutation$Start,mutation$Ref,mutation$Alt)
      
      #file=paste(runinfo$Path[i],"/intermediate/variants.vcf",sep="")
      #if (file.exists(paste(file,".gz",sep=""))) {system(paste("gzip -f -d ",file,".gz",sep=""))}
      #qual=read.table(file,stringsAsFactors = F,sep="\t")
      #rownames(qual)=paste(qual[,1],qual[,2],qual[,4],qual[,5])
      #qual=qual[mutation$mutation,6]
      mutations=rbind(mutations,mutation)
    },error=function(e) {})
  }
  
  runinfo=runinfo[runinfo$Cell %in% unique(mutations$id),]
  coverages=coverages[runinfo$Cell]
  
  # convert to matrix format for mutations
  mutations$vaf=mutations$Tumor_alt/(mutations$Tumor_alt+mutations$Tumor_ref)
  mutations_mat=matrix(0,ncol=length(unique(mutations$id)),nrow=length(unique(mutations$mutation)))
  colnames(mutations_mat)=runinfo$Cell
  rownames(mutations_mat)=unique(mutations$mutation)
  
  for (id in unique(mutations$id))
  {
    a=mutations[mutations$id==id,"mutation"]
    b=mutations[mutations$id==id,"vaf"]
    mutations_mat[a,id]=b 
  }
  
  # annotation
  other=c("SIFT_pred","Polyphen2_HVAR_pred","cosmic70",
          "esp6500siv2_all","ExAC_ALL","X1000g2015aug_all")
  other=other[other %in% colnames(mutations)]
  annotation=mutations[,c("Func.refGene","mutation","ExonicFunc.refGene","Gene.refGene",other)]
  annotation=annotation[!duplicated(annotation),]
  rownames(annotation)=annotation$mutation
  keep1=annotation[rownames(mutations_mat),"Func.refGene"] %in% 
    c("exonic","UTR3","UTR5","exonic;splicing","splicing","UTR5;UTR3")
  keep2=grepl("chrM",rownames(mutations_mat))
  if (!keep_other) {mutations_mat=mutations_mat[keep1 | keep2,,drop=F]}
  print(dim(mutations_mat))
  
  # add coverage info
  coverage_mat=mutations_mat
  coverage_mat[]=NA
  chr_pos=sapply(strsplit(rownames(mutations_mat)," "),function(x) paste(x[1],x[2]))
  
  for (i in 1:dim(mutations_mat)[2])
  {
    cat(paste("Processing coverage data for",i,"th cell\n"))
    coverage=coverages[[colnames(mutations_mat)[i]]]
    coverage_mat[,i]=coverage$V4[match(chr_pos,paste(coverage$V1,coverage$V2))]
  }
  
  mutations_mat[is.na(coverage_mat)]=NA
  
  # remove mutation artefact
  cutoff=max(artefact_percentage,1/dim(mutations_mat)[2])
  mutations_mat=mutations_mat[!apply(mutations_mat,1,function(x) all(is.na(x))),]
  mutations_mat=mutations_mat[apply(mutations_mat,1,function(x) mean(x>0.1,na.rm=T)>cutoff),,drop=F]
  print(dim(mutations_mat))
  
  # remove mutations that are not estimatable in too many cells
  keep=apply(is.na(mutations_mat),1,mean)<coverage_percentage
  mutations_mat=mutations_mat[keep,,drop=F]
  print(table(keep))
  mutations_mat=mutations_mat[apply(mutations_mat,1,function(x) sd(x,na.rm=T))>0,,drop=F]
  mutations_mat=mutations_mat[apply(mutations_mat,1,function(x) sum(x>0,na.rm=T))>1,,drop=F]
  mutations_mat=mutations_mat[apply(mutations_mat,1,function(x) sum(x<1,na.rm=T))>1,,drop=F]
  annotation=annotation[rownames(mutations_mat),]
  coverage_mat=coverage_mat[rownames(mutations_mat),]
  
  # cell percentage
  keep=apply(is.na(mutations_mat),2,mean)<cell_percentage
  mutations_mat=mutations_mat[,keep]
  coverage_mat=coverage_mat[,keep]
  print(table(keep))
  
  # adjust annotation
  annotation$ExonicFunc.refGene[annotation$Func.refGene=="UTR5"]="UTR5"
  annotation$ExonicFunc.refGene[annotation$Func.refGene=="UTR3"]="UTR3"
  annotation$ExonicFunc.refGene[grepl("splicing",annotation$Func.refGene)]="splicing"
  annotation$ExonicFunc.refGene[grepl("chrM",rownames(mutations_mat))]="chrM"
  annotation$ExonicFunc.refGene[annotation$ExonicFunc.refGene=="."]="other"
  
  # wrap up
  if (!file.exists(folder)) {dir.create(folder)}
  results = list(
    mutations_mat = mutations_mat,
    runinfo = runinfo,
    annotation = annotation,
    coverage_mat = coverage_mat
  )
  save(results,file=paste(folder,"/cleaned.RData",sep=""))
  
  pdf(paste(folder,"/summary.pdf",sep=""))
  heatmap.2(mutations_mat,na.color="grey",trace="none",Rowv=F,Colv=F,dendrogram="none")
  par(mar=c(5,15,5,5))
  plot(apply(mutations_mat,1,function(x) {x=x[!is.na(x)];mean(x)})~
         factor(annotation$ExonicFunc.refGene),horizontal=T,las=1,pch=19,col="green")
  dev.off()
  
  return(results)
}

