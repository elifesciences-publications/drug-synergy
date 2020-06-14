#april 23 2019
#combine MCF7 combos, self-synergy, and lncap combos and remove batch effects for pca

library(limma)
setwd('~/Google_Drive/ibm_synergy/gene_expression/PCA')

cpm = read.csv('mcf7_log_cpm.csv',row.names=1)
batch = substr(colnames(cpm),1,4)
cpm_new = removeBatchEffect(cpm,batch)
write.csv(cpm_new,'mcf7_log_cpm_removedbatch.csv')

rmeans = read.csv('mcf7_log_cpm_means.csv',row.names=1)
batch = substr(colnames(rmeans),1,4)
means_new = removeBatchEffect(rmeans,batch)
write.csv(means_new,'mcf7_log_cpm_means_removedbatch.csv')

means = read.csv('mcf7_batchremoved_means.csv',row.names=1)

fittt <- prcomp(t(means),
                center = TRUE,
                scale. = TRUE)

xxx=summary(fittt)
aaa=xxx$importance
y=aaa[3,]
x=1:length(y)

plot(x,y,main='Cumulative Contribution to Variance',
     xlab = 'Principal components',ylab='Percentage Contribution to Variance',
     type='l')

x<- fittt$x[,1]
y<- fittt$x[,2]
plot(x,-y)
text(x,-y,labels=row.names(fittt[['x']]),cex=0.5)
