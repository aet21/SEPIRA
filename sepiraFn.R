### sepiraFn.R
### Description: This file contains two R functions which together implement the SEPIRA (Systens EPigenomics Inference of Regulatory Activity) algorithm.
### Author: Andrew E Teschendorff (a.teschendorff@ucl.ac.uk)
### Date: 21 July 2017
### This software is released under a GPL2 licence

### required libraries
library(corpcor);
library(parallel);
library(limma);

### Auxiliary functions

LimmaFn <- function(pheno.v,data.m){

### construct model matrix
sampletype.f <- as.factor(pheno.v);
design.sample <- model.matrix(~0 + sampletype.f);
colnames(design.sample) <- levels(sampletype.f);
sampletypes.v <- levels(sampletype.f);

### do linear model fit
lmf.o <- lmFit(data.m,design.sample);

### construct contrast matrix
ntypes <- length(levels(sampletype.f));
ncomp <- 0.5*ntypes*(ntypes-1);
cont.m <- matrix(0,nrow=ncol(design.sample),ncol=ncomp);
tmp.v <- vector();
c <- 1;
for(i1 in 1:(ntypes-1)){
 for(i2 in (i1+1):ntypes){
   cont.m[i1,c] <- -1;
   cont.m[i2,c] <- 1;
   tmp.v[c] <- paste(sampletypes.v[i2],"--",sampletypes.v[i1],sep="");
   c <- c+1;
 }
}
rownames(cont.m) <- sampletypes.v; # sampletype.v determined separately
colnames(cont.m) <- tmp.v;

### do linear model to contrasts
lmf2.o <- contrasts.fit(lmf.o,cont.m);

### empirical Bayesian estimation of differentially expressed genes (DEGs)
bay.o <- eBayes(lmf2.o);

### build ranked list of DEGs for each comparison
top.lm <- list();
for(c in 1:ncol(cont.m)){
top.lm[[c]] <- topTable(bay.o,coef=c,adjust="fdr",number=nrow(data.m));
}

return(list(top=top.lm,cont=cont.m));

} ## end of limma function

ComputePCOR <- function(idx,mapTG.idx,mapTF.idx,selbinNET.m,exp.m){
    g <- idx;
    reg.idx <- which(selbinNET.m[g,]==1);
    if(length(reg.idx)>=2){
        tmp.idx <- c(mapTG.idx[g],mapTF.idx[reg.idx]);
        cor.m <- cor(t(exp.m[tmp.idx,]));
        pcor.m <- cor2pcor(cor.m);
    }
    else{
        pcor.m <- NULL;
    }
    return(pcor.m);
} ## end of ComputePCOR function

InferTFact <- function(exp.v,regnet.m){
  act.v <- apply(regnet.m,2,function(tmp.v){lm.o <- lm(exp.v ~ tmp.v); act <- summary(lm.o)$coeff[2,3];return(act);})
  return(act.v);
} ### end of InferTFact function

InferTFactPRL <- function(idx,tmp.m,regnet.m){
  exp.v <- tmp.m[,idx];
  act.v <- apply(regnet.m,2,function(tmp.v){lm.o <- lm(exp.v ~ tmp.v); act <- summary(lm.o)$coeff[2,3];return(act);})
  return(act.v);
} ### end of InferTFactPRL function

### MAIN SEPIRA FUNCTIONS
### sepiraInfNet main function
### DESCRIPTION: this function takes as input a large normalized gene expression data matrix (data.m) with rownames annotated to unique Entrez Gene IDs, encompassing multiple tissue types specified by a vector (tissue.v) of the same length as the number of columns in data.m. The user needs to specify the tissue of interest (toi) for which the regulatory network is desired. toi must be one of the tissues listed in tissue.v.
### Further INPUT parameters include:
### cft: a character or character vector containing the tissues which are deemed to be confounding. For instance, blood can be seen as confounder, since immune-cells infiltrate epithelial tissues. If specified, must be one of the tissues in tissue.v
### regEID.v: a user-specified vector of Entrez gene IDs specifying the regulator genes e.g. a vector of human transcription factors. However, the regulators need not be transcription factors.
### sdth: a threshold on the standard deviation of gene expression in data.m, to select genes with standard deviation larger than sdth.
### sigth: the significance P-value threshold for correlations. If not specified it uses the Bonferroni threshold
### pcorth: a significance threshold on the partial correlation values. By default this is set to 0.2. User needs to select or estimate this for each dataset and choose an appropriate threshold design to find the appropriate number of target genes per regulator
### degth.v: P-value thresholds for calling significant differential expression between the tissue of interest and all other tissues (1st element in deth.v) and optionally confounding tissues as specified in cft argument. There should be at least as many elements in degth.v as 1 + number of confounding tissues.
### lfcth.v: As degth.v, but specifying thresholds on log fold-change. This assumes data.m is in a log2-basis.
### minNtgts: Desired minimum number of target genes per regulator/TF in regulatory network
### ncores: number of cores to use with parallel package.

### OUTPUT:
### netTOI: the regulatory network for the tissue of interest. This is a high-confidence set of regulators highly expressed in the tissue of interest plus predicted downstream targets. Entries are 1 for activating interactions, 0 no interactions, and -1 for inhibiting interactions.
### sumnet: summary statistics for the regulators, including number of targets, plus number of positive and negative interactions
### topTOI.lm: list of top-ranked differentially expressed regulators for each comparison of interest. Entry-1 is for toi against all other tissues, further entries for comparisons of toi to confounding tissues as specified in cft

sepiraInfNet <- function(data.m,tissue.v,toi,cft=NULL,regEID.v,sdth=0.25,sigth=NULL,pcorth=0.2,degth.v=rep(0.05,3),lfcth.v=c(1,log2(1.5),log2(1.5)),minNtgts=10,ncores=4){

   tt.v <- levels(factor(tissue.v));
   if( length(intersect(toi,tt.v))==0){
       print("Your tissue of interest is not in tissue.v, or you have mispelled toi");
       stop;
   }
   
   ### remove genes with no or little variance
   sd.v <- apply(data.m,1,sd);
   selG.idx <- which(sd.v>sdth);
   exp.m <- data.m[selG.idx,];
   ### find representation of regulators in data, and define regulatees/targets
   tfEID.v <- regEID.v;
   repTF.v <- intersect(tfEID.v,rownames(exp.m));
   tgtsEID.v <- setdiff(rownames(exp.m),tfEID.v);
   match(repTF.v,rownames(exp.m)) -> mapTF.idx;
   match(tgtsEID.v,rownames(exp.m)) -> mapTGTS.idx;
   ### compute correlations and estimate P-values
   corNET.m <- cor(t(exp.m[mapTGTS.idx,]),t(exp.m[mapTF.idx,]));
   zNET.m <- 0.5*log( (1+corNET.m)/(1-corNET.m) );
   stdev <- 1/sqrt(length(tt.v)-3); ### this is not the number of independent samples but the number of independent tissues. the latter is used because of the strong dependence between samples from the same tissue, but also because it leads to a more stringent significance threshold
   pvNET.m <- 2*pnorm(abs(zNET.m),0,stdev,lower.tail=FALSE)
   ### for each gene, now identify the TFs which are correlated univariately- for these then run multivariate regression
   if(is.null(sigth)){
    sigth <- 0.05/prod(dim(pvNET.m));
   }
   binNET.m <- pvNET.m;
   binNET.m[pvNET.m < sigth] <- 1;
   binNET.m[pvNET.m >= sigth] <- 0;
   ### number of targets per tf
   ntgTF.v <- apply(binNET.m,2,sum)
   ### number of regulators per gene
   nregG.v <- apply(binNET.m,1,sum)

   ### select TFs with at least minNtgts targets
   selTF.idx <- which(ntgTF.v>=minNtgts);
   selbinNET.m <- binNET.m[,selTF.idx];
   selpvNET.m <- pvNET.m[,selTF.idx];
   selzNET.m <- zNET.m[,selTF.idx];
   selcorNET.m <- corNET.m[,selTF.idx];

   mapTG.idx <- match(rownames(selbinNET.m),rownames(exp.m));
   mapTF.idx <- match(colnames(selbinNET.m),rownames(exp.m));

   idx.l <- as.list(1:nrow(selbinNET.m));
   print("Computing Partial Correlations");
   pcor.l <- mclapply(idx.l,ComputePCOR,mapTG.idx,mapTF.idx,selbinNET.m,exp.m,mc.cores=ncores);

   pcorNET.m <- matrix(0,nrow=nrow(selbinNET.m),ncol=ncol(selbinNET.m));
   rownames(pcorNET.m) <- rownames(selbinNET.m);
   colnames(pcorNET.m) <- colnames(selbinNET.m);
   for(g in 1:length(pcor.l)){
    reg.idx <- which(selbinNET.m[g,]==1);
    if(length(reg.idx)>=2){
        pcorNET.m[g,reg.idx] <- pcor.l[[g]][1,-1];
    }
    else if (length(reg.idx)==1){
        pcorNET.m[g,reg.idx] <- selcorNET.m[g,reg.idx];
    }
    print(g);
   }

   gtexNET.m <- sign(pcorNET.m);
   gtexNET.m[abs(pcorNET.m) < pcorth] <- 0;

   sumnetTF.m <- matrix(nrow=4,ncol=ncol(gtexNET.m));
   rownames(sumnetTF.m) <- c("nTG","nUP","nDN","P");
   colnames(sumnetTF.m) <- colnames(gtexNET.m);
   sumnetTF.m[1,] <- apply(abs(gtexNET.m),2,sum);

   for(tf in 1:ncol(gtexNET.m)){
    sumnetTF.m[2,tf] <- length(which(gtexNET.m[,tf]==1));
    sumnetTF.m[3,tf] <- length(which(gtexNET.m[,tf]==-1));    
   }
   pv.v <- pbinom(apply(sumnetTF.m[2:3,],2,max),size=sumnetTF.m[1,],prob=0.5,lower.tail=FALSE);
   sumnetTF.m[4,] <- pv.v;

   gtexNETf.m <- gtexNET.m[,which(sumnetTF.m[1,]>= minNtgts)];

   ### now which TFs are overexpressed in toi?
   topTOI.lm <- list();
   ### compare toi to all other tissues
   pheno.v <- rep(0,ncol(exp.m));
   pheno.v[which(tissue.v==toi)] <- 1;
   lim.o <- LimmaFn(pheno.v,exp.m);
   topTOI.lm[[1]] <- lim.o$top[[1]];

   ### now compare toi to other tissues (in order to avoid confounding by immune or stromal cell infiltrates)
   if(!is.null(cft)){
    ti <- 2;
    for(t in cft){
     sel.idx <- which(tissue.v %in% c(toi,t));
     tmp.v <- tissue.v[sel.idx];
     tmpPH.v <- rep(0,length(sel.idx));
     tmpPH.v[which(tmp.v==toi)] <- 1;
     lim.o <- LimmaFn(tmpPH.v,exp.m[,sel.idx]);
     topTOI.lm[[ti]] <- lim.o$top[[1]];
     ti <- ti+1;
    }
   }

   ### now find tissue-specific TFs
   toiTF.lv <- list();
   for(i in 1:length(topTOI.lm)){
    statTF.m <- topTOI.lm[[i]][match(colnames(gtexNETf.m),rownames(topTOI.lm[[i]])),c(1,3,4,5)];
    toiTF.idx <- intersect(which(statTF.m[,4] < degth.v[i]),which(statTF.m[,1]>lfcth.v[i]));
    toiTF.lv[[i]] <- rownames(statTF.m[toiTF.idx,]);
   }

   toiTF.v <- toiTF.lv[[1]];
   if(length(toiTF.lv)>1){
     for(i in 2:length(toiTF.lv)){
       toiTF.v <- intersect(toiTF.v,toiTF.lv[[i]]);
     }
   }
   map.idx <- match(toiTF.v,colnames(gtexNETf.m));
   #### tissue-specific regulatory network is:
   netTOI.m <- gtexNETf.m[,map.idx];

   distNet.m <- matrix(nrow=ncol(netTOI.m),ncol=3);
   rownames(distNet.m) <- colnames(netTOI.m);
   colnames(distNet.m) <- c("nTGTS","Act","Rep");
   distNet.m[,1] <- apply(abs(netTOI.m),2,sum);
   for(c in 1:ncol(netTOI.m)){
    act.idx <- which(netTOI.m[,c]==1)
    inact.idx <- which(netTOI.m[,c]==-1)
    distNet.m[c,2:3] <- c(length(act.idx),length(inact.idx));
   }

   return(list(netTOI=netTOI.m,sumnet=distNet.m,top=topTOI.lm));
}  ### end of function sepiraInfNet

### MAIN SEPIRA function: sepiraRegAct
### DESCRIPTION: estimates regulatory activity in independent data
### INPUT ARGUMENTS:
### data: either a vector or a data matrix with entries/rows labeling genes and columns labeling samples. This data is assumed normalized.
### type: type of data, expression or DNAm
### regnet.m: the regulatory network to be used, i.e. the netTOI output argument of sepiraInfNet
### norm: further normalization of the data (if a matrix). Option "c" will center each gene to mean zero. Option "z" will z-score normalise each gene. Beware of genes with zero variance in data which should be removed if present.

### OUTPUT:
### actTF: a vector or matrix of estimated TF-activity levels. If a matrix, rows label the regulators/TFs, columns the samples.

sepiraRegAct <- function(data,type=c("mRNA","DNAm"),regnet.m,norm=c("c","z"),ncores=4){

 if(type=="DNAm"){
    regnet.m <- -regnet.m;
 }
 if(is.vector(data)){
   actTF <- InferTFact(data,regnet.m);
   names(actTF) <- colnames(regnet.m);
 }
 else if (is.matrix(data)){
  ndata <- data - rowMeans(data);
  if(norm=="z"){
   ndata <- (data - rowMeans(data))/apply(data,1,sd);
  }
  idx.l <- as.list(1:ncol(data));
  prl.o <- mclapply(idx.l,InferTFactPRL,ndata,regnet.m,mc.cores=ncores);
  actTF <- matrix(unlist(prl.o),nrow=ncol(regnet.m),ncol=length(prl.o),byrow=FALSE)
  rownames(actTF) <- colnames(regnet.m);
  colnames(actTF) <- colnames(data);
 }

 return(actTF);
}
