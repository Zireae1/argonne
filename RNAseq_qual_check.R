rm(list=ls(all=TRUE))
getwd()

library(data.table) # for fread function
library(RColorBrewer) # for colour palette
library(scales) # for points transparency on plots
library(Rtsne)  # for tSNE

setwd("/Users/Zireael/Desktop/hackathon2") # replace with your working directory

# load data
expr<-fread("data/matched_normal_samples.log_transformed.csv", stringsAsFactors = F)
#expr<-fread("http://bioseed.mcs.anl.gov/~fangfang/matched_normal/matched_normal_samples.log_transformed.csv", stringsAsFactors = F)
meta<-fread("data/metadata.txt", header = F)
#meta<-fread("http://bioseed.mcs.anl.gov/~fangfang/matched_normal/metadata", header = F)

# drop samples without pair
id<-gsub( "-.*", "", expr$Sample)
expr<-expr[duplicated(id) | duplicated(id, fromLast=TRUE),]

# subset random N samples
# id<-as.numeric(sample(seq(1,nrow(expr),by=2), size = 300))
# subs<-expr[sort(c(id, id+1)),]
subs<-expr# all paired samples

# add tissue types as 1st column
id<-gsub( "-.*", "", subs$Sample)
subs<-cbind(meta$V2[as.numeric(id)+1], subs)

# load samples from Rick 
simulated<-fread("data/decoded_out.csv", stringsAsFactors = F)

# assuming order of genes is the same as the other GDC:
labs<-as.vector(fread("data/feature.order.txt", header = F, stringsAsFactors = F))
#labs<-as.vector(fread("http://bioseed.mcs.anl.gov/~fangfang/tmp/feature.order", header = F, stringsAsFactors = F))
colnames(simulated)<-as.matrix(labs)[1:ncol(simulated)]

# reorder genes alphabetically as in GDC
simulated<-setcolorder(simulated, colnames(simulated)[order(colnames(simulated))])

# renormalize
simulated_resc<-t(apply(simulated, 1, function(x) (x-min(x))/(max(x)-min(x))))

# for GDC data table: subset genes which are present in simulated data 
subs<-as.data.frame(subs)
subs_lg<-subs[,which(colnames(subs)%in%colnames(simulated_resc))]

# renormalize
#subs_lg_resc<-subs_lg/rowSums(subs_lg)
subs_lg_resc<-t(apply(subs_lg, 1, function(x) (x-min(x))/(max(x)-min(x))))
subs_lg_resc<-cbind(subs[,1:2], subs_lg_resc) # first two columns - type and sample id

# compute PCA
pca<-prcomp(subs_lg_resc[,-c(1,2)])
summary(pca)$importance[,1:10] # check - 1,2 PC 10-15% of variance 

# calc tSNE for rescaled data
tsne <- Rtsne(pca$x, dims = 2, initial_dims = ncol(pca$x),
              perplexity=30, verbose=TRUE, max_iter = 500, pca=F) # change perplexity?

# plot tsne for rescaled data -> looks highly similar to the previous one
tmp<-as.data.frame(tsne$Y)                  # get first two comp
tmp$group<-gsub( ".*-", "", subs$Sample)    # add info to group by normal/not
tmp$type<-gsub( ".*-", "", subs$V1)         # add info to group by tissue type

# colours by tissue type
my_palette <- colorRampPalette(rev(brewer.pal(11,"Spectral")))(length(unique(meta$V2)))
tmp$color<-"" # create column for custom colors
k<-0
for(i in sort(unique(meta$V2))){
  k<-k+1
  tmp$color[which(tmp$type==i)]<-my_palette[k]
}

cairo_pdf('graphs/tsne_rescaled.pdf', width = 11, height = 10) # save plot for GDC
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(x=tmp$V2, y=tmp$V1, col=tmp$color, pch=16, cex=1)#+  # check - 1,2 PC -15% variance , label=tmp$names
arrows(x0=tmp$V2[which(tmp$group=="0")], y0=tmp$V1[which(tmp$group=="0")],
       x1=tmp$V2[which(tmp$group=="1")], y1=tmp$V1[which(tmp$group=="1")],
       lty =1, lwd =0.3, length=0.1)
legend("topright", inset=c(-0.18,0), bty="n", legend = sort(unique(tmp$type)), 
       col = my_palette, lty= 0, pch = 16)
dev.off()

# coord of GDC in same format:
pctr<-cbind(subs_lg_resc[,c(1,2)], pca_resc$x)

# calculate coord of simulated samples in PC space
simulated_resc<-as.data.frame(simulated_resc)
pcsim<-scale(simulated_resc, pca_resc$center, pca_resc$scale) %*% pca_resc$rotation 

#  create the same structure as for GDC data with type "sim"
pcsim<-cbind(as.data.frame(rep("sim", nrow(pcsim)), stringsAsFactors = F), 
             as.data.frame(rownames(simulated_resc), stringsAsFactors = F), pcsim)

# train and predict x coord
colnames(pctr)[1]<-"type"
data<-cbind(tmp$V1,pctr[,-c(1,2)])
colnames(data)[1]<-"V1"
fit<-lm(V1~., data=data)
predx<-predict(fit, as.data.frame(pcsim[,-c(1,2)]))

# train and predict y coord
data<-cbind(tmp$V2,pctr[,-c(1,2)])
colnames(data)[1]<-"V2"
fit<-lm(V2~., data=data)
predy<-predict(fit, pcsim[,-c(1,2)])

new<-cbind(predx, predy, rep(2, length(predx)), pcsim$`rep("sim", nrow(pcsim))`, rep("black", length(predx)))
new<-as.data.frame(new, stringsAsFactors =F)
colnames(new)<-colnames(tmp)[1:5]

# plot simulated samples on tSNE graph
all_plot<-rbind(tmp, new)

cairo_pdf('graphs/tsne_resc_with_pred_sim.pdf', width = 11, height = 10) # save plot
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(x=all_plot$V2, y=all_plot$V1, col=alpha(all_plot$color, 0.0), pch=all_plot$symb, cex=1)# label=tmp$names
#arrows(x0=all_plot$V2[which(all_plot$group=="0")], y0=all_plot$V1[which(all_plot$group=="0")],
#       x1=all_plot$V2[which(all_plot$group=="1")], y1=all_plot$V1[which(all_plot$group=="1")],
#       lty =1, lwd =0.3, length=0.1)
points(x=all_plot$V2[which(all_plot$group==2)], y=all_plot$V1[which(all_plot$group==2)], 
       col=alpha(all_plot$color[which(all_plot$group==2)],1), pch=17, cex=0.7)
points(x=all_plot$V2[which(all_plot$group!=2)], y=all_plot$V1[which(all_plot$group!=2)], 
       col=alpha(all_plot$color[which(all_plot$group!=2)],0.8), pch=16, cex=0.7)

legend("topright", inset=c(-0.18,0), bty="n", legend = sort(unique(tmp$type)), 
       col = my_palette, lty= 0, pch = 16)

dev.off()

# recalc tsne for combined datamatrix
tsne <- Rtsne(rbind(pctr[,-c(1,2)],pcsim[,-c(1,2)]), dims = 2, initial_dims = ncol(pca_resc$x),
              perplexity=30, verbose=TRUE, max_iter = 500, pca=F) # change perplexity?

# plot tsne
tmp<-as.data.frame(tsne$Y)   # get first two comp
tmp$group<-all_plot$group    # add info to group by normal/not
tmp$type<-all_plot$type      # add info to group by tissue type

# colours by tissue type
my_palette <- colorRampPalette(rev(brewer.pal(11,"Spectral")))(length(unique(meta$V2)))
tmp$color<-"black"
k<-0
for(i in sort(unique(meta$V2))){
  k<-k+1
  tmp$color[which(tmp$type==i)]<-my_palette[k]
}

cairo_pdf('graphs/tsne_with_new_recalc.pdf', width = 11, height = 10) # save plot
par(mar=c(5.1, 4.1, 4.1, 12.1), xpd=TRUE)
plot(x=tmp$V2, y=tmp$V1, col=alpha(tmp$color, 0.0), pch=tmp$symb, cex=1)
points(x=tmp$V2[which(tmp$group==2)], y=tmp$V1[which(tmp$group==2)], 
       col=alpha(tmp$color[which(tmp$group==2)],1), pch=17, cex=0.7)
points(x=tmp$V2[which(tmp$group!=2)], y=tmp$V1[which(tmp$group!=2)], 
       col=alpha(tmp$color[which(tmp$group!=2)],0.8), pch=16, cex=0.7)
legend("topright", inset=c(-0.31,0), bty="n", legend = sort(unique(tmp$type)), 
       col = my_palette, lty= 0, pch = 16)
dev.off()




