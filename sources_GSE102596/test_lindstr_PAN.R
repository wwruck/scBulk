################################################################################################
# R code for analysis of sc-RNAseq kidney dataset GSE102596
# and deconvolution of one kidney organoid microarray datasets bundled with this software
# R version 4.3.2
#
# input:
# sc-data already preprocessed (counts reduced to markers from NCBI GEO GSE202109): lindstr_counts_CmarkersU.txt
# (unzip lindstr_counts_CmarkersU.zip)
# k-means clusters and associated Lindstroem clusters: clusters_kmeans.txt
# bulk kidney organoid microarray dataset: PAN_kidneyorgs_sym.txt
# celltype cluster lookup table: clusterdesc_lindstr_jasn.txt
#
# author: Wasco Wruck, 2024
# institute: Institute for Stem Cell Research and Regenerative Medicine
# University Clinic Duesseldorf, Germany
################################################################################################

library(preprocessCore)
library(sva)
library(gplots)
library(InstaPrism)

# unzip file lindstr_counts_CmarkersU.zip
# read lindstroem count table reduced to unique cluster markers directly 
t2 <- read.table("lindstr_counts_CmarkersU.txt", header=TRUE, sep="\t",row.names=1,as.is=TRUE)
tkm <- read.table("clusters_kmeans.txt", header=TRUE, sep="\t",row.names=1,as.is=TRUE) # read clusters from file
ac=tkm$lindstr_JASN_cl
ct=read.table("clusterdesc_lindstr_jasn.txt", header=T) # read cluster descriptions
cs=ct[ac,"desc"]

# read bulk microarray data of kidney organoids treated with PAN to be tested
tdata <- read.table("PAN_kidneyorgs_sym.txt", header=TRUE, sep="\t",row.names=1,as.is=TRUE)
m=intersect(rownames(t2),rownames(tdata))
tdata_red=tdata[m,]
t2_red=t2[m,]

mat=as.matrix(cbind(tdata_red,t2_red))
mat2=ComBat(dat=mat,batch=c(rep(1,dim(tdata_red)[2]),rep(2,dim(t2_red)[2])),prior.plots =TRUE)
matn=normalize.quantiles(mat2)

tdata_redn=matn[,1:dim(tdata_red)[2]]
t2_redn=matn[,(dim(tdata_red)[2]+1):dim(mat)[2]]
colnames(tdata_redn)=colnames(tdata_red)
rownames(tdata_redn)=rownames(tdata_red)
colnames(t2_redn)=colnames(t2_red)
rownames(t2_redn)=rownames(t2_red)

# InstaPrism deconvolution
InstaPrism.res = InstaPrism(input_type = 'raw',sc_Expr = t2_redn,bulk_Expr = tdata_redn,
                    cell.type.labels = cs,cell.state.labels = cs)
prediction.matrix=InstaPrism.res@Post.ini.cs@theta[,1:4]
write.table(prediction.matrix,file="prediction.matrix_GSE102596_PAN.txt",col.names=T,row.names=T,sep="\t",quote=F)
pdf("predictionmapInstaPrism_cs_GSE102596_PAN.pdf")
rgb.palette <- colorRampPalette(c("#000000", "white", "green"), space = "rgb")
par(cex.main=1)
heatmap.2(prediction.matrix,  margins=c(13, 13),
          main="sc-Type - InstaPrism Prediction",
          Rowv=F, Colv=F, symm=FALSE, dendrogram="none",
          cexRow=0.7, cexCol=0.6, sepwidth=c(0.0001,0.0001),
          colsep=1:ncol(prediction.matrix), rowsep=1:nrow(prediction.matrix),
          breaks=seq(0,1,1/100),
          keysize=1,
          col=rgb.palette(100),
          trace="none", density.info="none")
dev.off()


