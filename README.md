# SEPIRA
Systems EPigenomics Inference of Regulatory Activity
SEPIRA is a novel systems-epigenomics algorithm which estimates activity of a regulator (e.g. a transcription factor) in a sample for which
there is a transcriptome or DNA methylome profile available. Briefly, it takes as input a large gene expression dataset encompassing many
different tissue types (user must provide this) and for a tissue of interest, it construct a tissue-specific regulatory network. From this
network it then estimates regulatory activity for a sample with an input expression or DNA methylation profile.

SEPIRA will be developed into an R and BioC package. For the time being, we provide a file containing all the necessary R-scripts (sepiraFn.R) and an R-script that can be used as a template for running SEPIRA (RunSEPIRA.R). At this stage, the user needs to provide all the input data, including the list of desired regulators (transcription factors).
