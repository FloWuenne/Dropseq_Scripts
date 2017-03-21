## Dropseq functions for data analysis and plotting
## Functions used from other sources are generally marked as such

#### Beautiful dot plot with legends function and S4 class definition taken from 
# S4 class file for Shekhar et al., 
# "Comprehensive classification of retinal bipolar cells using single-cell transcriptomics", Cell, 2016

# Required packages
require(Matrix)
require(igraph)
require(gmodels)
require(ggplot2)
require(sva)
require(RANN)
require(reshape)

# Define slots in S4 object
scDrop <- setClass("scDrop", slots = 
                     c(count.data="data.frame", data="data.frame",batch.data="data.frame", scale.data="matrix", 
                       group="vector", pca.load="data.frame",pca.scores="data.frame",
                       meta="data.frame",tsne.y="data.frame", cell.names="vector"))

## Function to call dot.plot
setGeneric("dot.plot", function(object,features.use=NULL, group.use=NULL, group.names=NULL, thresh.use=0,do.transpose=FALSE,max.val.perc=NULL, max.val.exp=NULL,max.size=10,min.perc=0,...) standardGeneric("dot.plot"))
setMethod("dot.plot", "scDrop", 
          function(object,features.use=NULL,group.use=NULL, group.names = NULL,thresh.use=0,do.transpose=FALSE,max.val.perc=NULL, max.val.exp=NULL,max.size=10,min.perc=0,...) {
            
            
            features.use=features.use[features.use %in% rownames(object@data)]
            if (is.null(group.use)) group.use = levels(object@group)
            if (is.null(group.names)) group.names = group.use
            
            if (length(group.names) != length(group.use)){
              print("Error : group.names must be of the same length as the groups.use/ number of clusters. Using cluster numbers as labels ...")
              group.names = group.use
            }
            
            #Initialize matrix of percent expressing cells
            PercMat = matrix(0, nrow=length(features.use), ncol = 0)
            rownames(PercMat) = features.use; 
            
            #Initialize matrix of average transcript levels
            ExpMat = PercMat;
            
            #Count mat
            Count.mat = object@count.data[features.use, colnames(object@data)]
            
            
            for (i in group.use){
              cells.in.cluster = names(object@group)[which(object@group== i)]
              print(cells.in.cluster)
              vec.exp = apply(object@data[features.use, cells.in.cluster], 1, function(x) sum(x>thresh.use)/length(x)) 
              PercMat = cbind(PercMat,vec.exp)
              
              
              vec.exp = apply(Count.mat[features.use, cells.in.cluster], 1, function(x) if (sum(x>0) > 1){ mean(x[x>0]) } else {sum(x)})
              
              ExpMat = cbind(ExpMat, vec.exp)
              
            }
            colnames(ExpMat) = group.names
            colnames(PercMat) = group.names
            
            
            
            rows.use = rownames(PercMat)[apply(PercMat, 1, function(x) max(x) >= min.perc)];
            
            PercMat = PercMat[rows.use,]
            ExpMat = ExpMat[rows.use,]
            features.use = rows.use
            if (!is.null(max.val.perc)) PercMat[PercMat > max.val.perc] = max.val.perc
            if (!is.null(max.val.exp)) ExpMat[ExpMat > max.val.exp] = max.val.exp
            
            
            ExpVal = melt(ExpMat)
            PercVal = melt(PercMat)
            colnames(ExpVal) = c("gene","cluster","nTrans")
            ExpVal$percExp = PercVal$value*100
            
            if (!do.transpose){
              ExpVal$gene = factor(ExpVal$gene, levels=features.use)
              ExpVal$cluster = factor(ExpVal$cluster, levels= rev(group.names))
              p=ggplot(ExpVal, aes(y = factor(cluster),  x = factor(gene))) + geom_point(aes(colour = nTrans,  size =percExp)) + 
                scale_color_gradient(low ="blue",   high = "red", limits=c( 1, max(ExpVal$nTrans) ))+scale_size(range = c(0, max.size))+   theme_light() +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
              p = p + ylab("Cluster") + xlab("Gene") + theme(axis.text.x=element_text(size=12, face="italic", angle=45, hjust=1)) + 
                theme(axis.text.y=element_text(size=12, face="italic"))
              print(p)
            } else {
              ExpVal$gene = factor(ExpVal$gene, levels=rev(features.use))
              ExpVal$cluster = factor(ExpVal$cluster, levels= group.names)
              p=ggplot(ExpVal, aes(y = factor(gene),  x = factor(cluster))) + geom_point(aes(colour = nTrans,  size =percExp)) + 
                scale_color_gradient(low ="blue",   high = "red", limits=c(1, max(ExpVal$nTrans) ))+scale_size(range = c(0, max.size))+   theme_bw() +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
              p = p + xlab("Cluster") + ylab("Gene") + theme(axis.text.x=element_text(size=12, face="italic", hjust=1)) + 
                theme(axis.text.y=element_text(size=12, face="italic")) 
              
              print(p)
              
              
              
            }
            
          }
)
