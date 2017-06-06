rm(list=ls(all=TRUE))
getwd()
#library(ggplot2) # for ggplot
library(data.table) # for fread function
#library(reshape2) # for cast
library(RColorBrewer) # for colour palette
library(Rtsne)
#library('e1071')

setwd("/Users/Zireael/Desktop/hackathon2") # replace with your working directory

expr<-fread("data/matched_normal_samples.log_transformed.csv")
meta<-fread("data/metadata.txt", header = F)

# drop samples without pair
id<-gsub( "-.*", "", expr$Sample)
expr<-expr[duplicated(id) | duplicated(id, fromLast=TRUE),]

# subset random N samples
#id<-as.numeric(sample(seq(1,nrow(expr),by=2), size = 300))
id<-as.numeric(sample(seq(1,nrow(expr),by=2), size = 702))
subs<-expr[sort(c(id, id+1)),]

# add 1st column with sample tissue types
id<-gsub( "-.*", "", subs$Sample)
subs<-cbind(meta$V2[as.numeric(id)+1], subs)

# compute PCA
pca<-prcomp(subs[,-c(1,2)])
#str(expr)
summary(pca)

# plot PCA
tmp<-as.data.frame(pca$x[,1:2])             # get first two comp
tmp$group<-gsub( ".*-", "", subs$Sample)    # group by normal/not
tmp$type<-gsub( ".*-", "", subs$V1)         # group by tissue type
tmp$names<-gsub( ".*\\.", "", rownames(subs))

# colours by tissue type
my_palette <- colorRampPalette(rev(brewer.pal(11,"Spectral")))(length(unique(meta$V2)))
tmp$color<-""
k<-0
for(i in sort(unique(meta$V2))){
  k<-k+1
  tmp$color[which(tmp$type==i)]<-my_palette[k]
}

cairo_pdf('graphs/set_all.pdf', width = 11, height = 10) # save plot

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE) # change graph params for legend 

plot(x=tmp$PC2, y=tmp$PC1, col=tmp$color, pch=16, cex=1)#+  # check - 1,2 PC -15% variance , label=tmp$names
arrows(x0=tmp$PC2[which(tmp$group=="0")], y0=tmp$PC1[which(tmp$group=="0")], 
       x1=tmp$PC2[which(tmp$group=="1")], y1=tmp$PC1[which(tmp$group=="1")], 
       lty =1, lwd =0.5, length=0.1)  
legend("topright", inset=c(-0.18,0), bty="n", legend = sort(unique(tmp$type)), 
       col = my_palette[which(my_palette%in%sort(unique(tmp$color)))], lty= 0, pch = 16)

dev.off()


# tmp$color<-"red" # colour by two types: tumor/normal
# tmp$color[which(tmp$group=="1")]<-"royal blue" 

library(ecodist)
#tt<-bcdist(subs[,-c(1,2)])
rowSums(subs[,-c(1,2)])

### try tSNE
tsne <- Rtsne(pca$x, dims = 2, 
              perplexity=30, verbose=TRUE, max_iter = 500) # change perplexity?
#tsne <- Rtsne(tt, dims = 2, 
#              perplexity=30, verbose=TRUE, max_iter = 500)


# plot tsne
tmp<-as.data.frame(tsne$Y)             # get first two comp
tmp$group<-gsub( ".*-", "", subs$Sample)    # group by normal/not
tmp$type<-gsub( ".*-", "", subs$V1)         # group by tissue type
#tmp$names<-gsub( ".*\\.", "", rownames(subs))
# colours by tissue type
my_palette <- colorRampPalette(rev(brewer.pal(11,"Spectral")))(length(unique(meta$V2)))
tmp$color<-""
k<-0
for(i in sort(unique(meta$V2))){
  k<-k+1
  tmp$color[which(tmp$type==i)]<-my_palette[k]
}

#cairo_pdf('graphs/set_all_tsne_1.pdf', width = 11, height = 10)

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(x=tmp$V2, y=tmp$V1, col=tmp$color, pch=16, cex=1)#+  # check - 1,2 PC -15% variance , label=tmp$names
arrows(x0=tmp$V2[which(tmp$group=="0")], y0=tmp$V1[which(tmp$group=="0")],
       x1=tmp$V2[which(tmp$group=="1")], y1=tmp$V1[which(tmp$group=="1")],
       lty =1, lwd =0.3, length=0.1)
legend("topright", inset=c(-0.18,0), bty="n", legend = sort(unique(tmp$type)), 
       col = my_palette, lty= 0, pch = 16)

#dev.off()

### plot new samples on this graph -> from PCA comp and tsne comp
news<-fread("data/combined.c10.n10.csv")
# do the same PCA transformation
pcs<-as.matrix(news[,-c(1,2)])%*%as.matrix(pca$rotation)
pcs<-cbind(news[,c(1,2)], pcs)

x<-as.numeric(pcs[1,-c(1,2)])
y<-as.numeric(pcs[2,-c(1,2)])

head(pca$rotation[,1:10])

sim<-function(x,y){
  1/(1+as.numeric(x-y)%*%as.numeric(x-y))
}

part<-10

dim(pcs)
dim(pca$x)

rowSums(pcs)
rowSums(pca$x)

# why there is such difference?
rowSums(news[,-c(1,2)])
rowSums(subs[,-c(1,2)])

rowSds

all<-rbind(pcs[1:part,], pca$x)
dt<-matrix(0, nrow(subs)+part, nrow(subs)+part)
for (i in 1:part) {
  sim(all[i, -c(1,2)])
}


















# which(subs$V1=="Adrenal Gland")
# which(tmp$type=="Adrenal Gland")
# 
# 
# which(subs$V1=="Bile Duct")
# 
# 
# which(subs$V1=="Thyroid")
# which(tmp$type=="Adrenal Gland")



### by types
# subset random 702 samples
id<-as.numeric(sample(seq(1,nrow(expr),by=2), size = 702))
subs<-expr[sort(c(id, id+1)),]

# add 1st column with sample tissue types
id<-gsub( "-.*", "", subs$Sample)
subs<-cbind(meta$V2[as.numeric(id)+1], subs)

onet<-subs[which(subs$V1=="Kidney")]

# compute PCA
pca<-prcomp(onet[,-c(1,2)])
str(expr)
View(summary(pca))

# plot PCA
tmp<-as.data.frame(pca$x[,1:2])             # get first two comp
tmp$group<-gsub( ".*-", "", onet$Sample)    # group by normal/not
tmp$type<-gsub( ".*-", "", onet$V1)         # group by tissue type
tmp$names<-gsub( ".*\\.", "", rownames(onet))
# colours by tissue type
my_palette <- colorRampPalette(rev(brewer.pal(11,"Spectral")))(length(unique(meta$V2)))
tmp$color<-""
k<-0
for(i in unique(meta$V2)){
  k<-k+1
  tmp$color[which(tmp$type==i)]<-my_palette[k]
}
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)

plot(x=tmp$PC2, y=tmp$PC1, col=tmp$color, pch=16, cex=1)#+  # check - 1,2 PC -15% variance , label=tmp$names
arrows(x0=tmp$PC2[which(tmp$group=="0")], y0=tmp$PC1[which(tmp$group=="0")], 
       x1=tmp$PC2[which(tmp$group=="1")], y1=tmp$PC1[which(tmp$group=="1")], 
       lty =1, lwd =0.5)  
legend("topright", inset=c(-0.18,0), bty="n", legend = sort(unique(tmp$type)), 
       col = my_palette[which(my_palette%in%sort(unique(tmp$color)))], lty= 0, pch = 16)



########### MDS is bad
library(ecodist)
ttt<-subs[,-c(1,2)]/rowSums(subs[,-c(1,2)])
tt<-bcdist(ttt)

library(MASS)
mds<-isoMDS(tt)
# plot MDS
tmp<-as.data.frame(mds$points)   # get first two comp
tmp$group<-gsub( ".*-", "", subs$Sample)
tmp$type<-gsub( ".*-", "", subs$V1)         # group by tissue type
my_palette <- colorRampPalette(rev(brewer.pal(11,"Spectral")))(length(unique(meta$V2)))
tmp$color<-""
k<-0
for(i in sort(unique(meta$V2))){
  k<-k+1
  tmp$color[which(tmp$type==i)]<-my_palette[k]
}

cairo_pdf('graphs/set_all_mds.pdf', width = 11, height = 10)

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(x=tmp$V2, y=tmp$V1, col=tmp$color, pch=16, cex=1)#+  # check - 1,2 PC -15% variance , label=tmp$names
arrows(x0=tmp$V2[which(tmp$group=="0")], y0=tmp$V1[which(tmp$group=="0")],
       x1=tmp$V2[which(tmp$group=="1")], y1=tmp$V1[which(tmp$group=="1")],
       lty =1, lwd =0.3, length=0.1)
legend("topright", inset=c(-0.18,0), bty="n", legend = sort(unique(tmp$type)), 
       col = my_palette, lty= 0, pch = 16)
dev.off()


library(vegan)
mds<-metaMDS(tt)
# plot MDS
tmp<-as.data.frame(mds$points)   # get first two comp
tmp$group<-gsub( ".*-", "", subs$Sample)    # group by normal/not
tmp$type<-gsub( ".*-", "", subs$V1)         # group by tissue type
#tmp$names<-gsub( ".*\\.", "", rownames(subs))
# colours by tissue type


length(unique(gsub(".*-","",expr$Sample)))
xx<-gsub(".*-","",expr$Sample)
length(which(xx=="0"))
length(which(xx=="1"))
