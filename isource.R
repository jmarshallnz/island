sc <- c("CHICKEN","CATTLE","SHEEP", "WATER/ENVIRONMENT")
ng = length(sc)

#########################
### READ IN THE FILES ###
#########################
### SET THE DIRECTORY
mcmc_dir = "" #~/Documents/C++/Campy/source/Distribute/XP/"
### LIST THE FILENAME(S)
fnames = c("out1")
### READ IN THE FILES (MAKE TAKE A WHILE)
mcmc = NULL; fmcmc = NULL;
for(i in 1:length(fnames)) {
  fname = fnames[i]
  mcmc_file = paste(mcmc_dir,fname,".txt",sep="")
  mcmc = rbind(mcmc,read.table(mcmc_file,header=T,comment.char=""))
  fmcmc_file = paste(mcmc_dir,"f_",fname,".txt",sep="")
  fmcmc = rbind(fmcmc,read.table(fmcmc_file,header=T,comment.char=""))
  g_file = paste(mcmc_dir,"g_",fname,".txt",sep="")
  if(i==1) g = matrix(scan(g_file,what=double(0),sep="\t"),nrow=ng)
  else g = g + matrix(scan(g_file,what=double(0),sep="\t"),nrow=ng)
}
g = t(g)/length(fnames)
### SET THE BURN-IN
gd = mcmc$iter>=5000
fd = fmcmc$iter>=1000

################################################################
### PLOT 1                                                   ###
### VISUALISE DIRECTLY THE MCMC OUTPUT FOR PARAMETER F       ###
### (THE PROPORTION OF ISOLATES ATTRIBUTABLE TO EACH SOURCE) ###
################################################################
pdf(file=paste(fname,".pdf",sep=""), width=11, height=8)

plot(fmcmc$f0[fd],type="l",ylim=c(0,1),col=2,ylab="Proportion")
for(i in 2:ng) lines(fmcmc[fd,(1+i)],col=rainbow(ng)[i])

#############################################################
### PLOT 2                                                ###
### HISTOGRAMS OF THE MARGINAL DISTRIBUTIONS OF F[i]      ###
### (THE PROPORTION OF ISOLATES ATTRIBUTABLE TO SOURCE i) ###
#############################################################
par(mfrow=c(2,2))
### PLOT THE HISTOGRAMS
for(i in 1:ng) hist(fmcmc[fd,(1+i)],30,col=rainbow(ng)[i],main=sc[i],prob=T,xlim=c(0,1),xlab="Proportion")

########################################################
### TABLE 1                                          ###
### SUMMARIES OF THE POSTERIOR DISTRIBUTIONS OF F[i] ###
########################################################
df = fmcmc[fd,2:(ng+1)]; names(df) <- sc;
pe = apply(df,2,function(x)c("mean"=mean(x),"median"=median(x),"sd"=sd(x),quantile(x,c(.025,.975))))
print(pe)

#################################################################################
### PLOT 3                                                                    ###
### BARCHART OF THE ESTIMATED PROPORTION OF CASES ATTRIBUTABLE TO EACH SOURCE ###
#################################################################################
par(mfrow=c(1,1))
mp = barplot(pe[1,],col=rainbow(ng),ylim=c(0,1),ylab="Proportion of human cases")
segments(mp,pe[4,],mp,pe[5,],lwd=2)

#################################################
### PLOT 4                                    ###
### MCMC TRACE OF THE EVOLUTIONARY PARAMETERS ###
#################################################
### 3 BY 3 PANE
par(mfrow=c(2,2))
COL = c(rainbow(ng),"black");
for(i in 0:(ng-1)) {
  plot(mcmc$iter[gd],mcmc[[paste("A",i,0,"",sep=".")]][gd],type="l",col=COL[1],ylim=c(0,1),xlab="iter",ylab="M,R",main=sc[i+1])
  for(j in 2:(ng+1)) lines(mcmc$iter[gd],mcmc[[paste("A",i,j-1,"",sep=".")]][gd],col=COL[j])
  lines(mcmc$iter[gd],mcmc[[paste("r",i,sep="")]][gd],col="grey")
}

#################################################################
### TABLE 2                                                   ###
### PROBABILITY THAT EACH SOURCE IS RESPONSIBLE FOR THE MOST, ###
### SECOND MOST, THIRD MOST, ETC, NUMBER OF CASES             ###
#################################################################
od = t(apply(df,1,order,decreasing=T))
### THE TABLE
print(apply(od,2,function(x)table(factor(x,levels=1:ng,labels=sc)))/nrow(df))
enc = apply(od,1,function(x) sum(x*((0:(ng-1))^ng)))
encmode = as.numeric(levels(factor(enc)))[which.max(table(enc))]
### MOST LIKELY ORDER
print(sc[od[which(enc==encmode)[1],]])
### POSTERIOR PROBABILITY OF THIS, THE MOST LIKELY ORDER
print(max(table(enc))/nrow(df))

#################################################################
### TABLE 3                                                   ###
### MOST PROBABLE ORDERINGS AND THEIR POSTERIOR PROBABILITIES ###
#################################################################
enc_od = as.numeric(levels(factor(enc)))[order(table(enc),decreasing=TRUE)]
enc_pr = table(enc)[order(table(enc),decreasing=TRUE)]/nrow(df)
unenc_od = t(sapply(enc_od,function(i) sc[od[which(enc==i)[1],]]))
print(cbind(unenc_od,"Posterior probability"=enc_pr)[enc_pr>.05,])

### Ordering of top 3 no. cases: full order
good = apply(od[,1:3],1,function(x)all(sort(x)==c(1,2,5)))
odg = od[good,1:3]
encg = apply(odg,1,function(x) sum(x*((0:2)^3)))
encg_od = as.numeric(levels(factor(encg)))[order(table(encg),decreasing=TRUE)]
encg_pr = table(encg)[order(table(encg),decreasing=TRUE)]/nrow(df)
unencg_od = t(sapply(encg_od,function(i) sc[odg[which(encg==i)[1],]]))
print(cbind(unencg_od,encg_pr)[1:min(4,nrow(unencg_od)),])

#################################################
### PLOT 5                                    ###
### PIE CHARTS OF THE EVOLUTIONARY PARAMETERS ###
#################################################
par(mfrow=c(2,2))
COL = c(rainbow(ng),"black");
for(i in 0:(ng-1)) {
  wh0 = which(names(mcmc)==paste("A",i,0,"",sep="."))
  whng = which(names(mcmc)==paste("A",i,ng,"",sep="."))
  pie(apply(mcmc[gd,wh0:whng],2,mean),col=COL,labels="",main=sc[i+1],col.main=COL[i+1],radius=1.,border="white")
}

########################################################
### PLOT 6                                           ###
### POSTERIOR PROBABILITY OF SOURCE FOR EACH ISOLATE ###
########################################################
cod = c(1,2,3,4)
G = g[,cod]
od = order(G[,1],G[,2],G[,3],G[,4])
res = 1000
tp = apply(G,1,function(x)sort(sample(1:ncol(G),res,replace=TRUE,prob=x)))
stretch = function(x,res) {
  y = floor(x*res)
  while(sum(y)<res) {
    i = which.max( abs(y/sum(y)-x) )
    y[i] = y[i]+1
  }
  rep(1:length(x),times=y)
}
tp2 = apply(G,1,stretch,res)

### DO THE PLOT
par(mfrow=c(1,2))
image(1:nrow(G),seq(0,1,len=res),t(tp2[,od]),col=COL[cod],ylab="Source probability",xlab="Human cases",bty="n")

### DO THE PLOT, BUT RE-ORDERED
wm = apply(G,1,which.max)
od = rev(order(wm!=4,wm!=5,G[,1],G[,2],G[,3],G[,4]))
image(1:nrow(G),seq(0,1,len=res),t(tp2[,od]),col=COL[cod],ylab="Source probability",xlab="Human cases",bty="n")

##########################################################
### TABLE 4                                            ###
### PROBABILITIES OF SOURCE ATTRIBUTION FOR EACH HUMAN ###
### (ORDERED AS IN THE ORIGINAL FILE)                  ###
### NB: HUMANS WITH THE SAME ST WILL HAVE THE SAME     ###
### ATTRIBUTION PROBABILITIES, SO CROSS-REFERRING THIS ###
### LIST BACK TO THE ORIGINAL LIST OF STs WOULD ALLOW  ###
### YOU TO CONSTRUCT A CONDENSED LIST OF ST-SPECIFIC   ###
### SOURCE ATTRIBUTION PROBABILITIES                   ###
##########################################################
colnames(g) <- sc
print(g)

dev.off()
