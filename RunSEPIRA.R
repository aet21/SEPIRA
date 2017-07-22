### RunSEPIRA.R

### This is a template R-script for running SEPIRA

source("sepiraFn.R");

### prepare your data matrix as data.m, with rows labeling genes (Entrez gene IDs) and columns labeling samples (tissues)
### prepare your tissue.v vector, specifying the tissue type of each sample in data.m
### decide on a tissue of interest (e.g lung)
### does your tissue of interest have a lot of immune cell infiltration? (you can check this with the ESTIMATE algorithm. If yes, set cft=c("blood").
### specify a vector of human transcription factors annotated to Entrez gene IDs: regEID.v

### now run SEPIRA to infer regulatory network
net.o <- sepiraInfNet(data.m,tissue.v,toi="lung",cft=c("blood"),regEID.v,sdth=0.25,sigth=NULL,pcorth=0.2,degth.v=rep(0.05,3),lfcth.v=c(1,log2(1.5),log2(1.5)),minNtgts=10,ncores=4);
### now estimate TF-activity in the same samples (as a sanity check)
actTF.m <- sepiraRegAct(data.m,type="mRNA",regnet.m=net.o$netTOI);
### check that mean TF-activity is highest in your tissue of interest
boxplot(colMeans(actTF.m) ~ tissue.v);
