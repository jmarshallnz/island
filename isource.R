sc <- c("Poultry","Cattle","Sheep", "Water-Environment", "Ruminants")
o <- c(1,2,3,4)
col = c("#FF7F00","#CF0000","#004FCF", "#009F9F","#8F006F","#9F5F3F","#FFAFAF")
#col = brewer.pal(8, "Set1")[c(5,1,2,3,4,7,8)]

ng = length(o)  # TODO: Automate this
nt = 32 # TODO: Automate this

#########################
### READ IN THE FILES ###
#########################
### SET THE DIRECTORY
mcmc_dir = "" #~/Documents/C++/Campy/source/Distribute/XP/"
### LIST THE FILENAME(S)
fnames = c("out2")
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
fd = fmcmc$iter>=1000 #10000

# read in our humans and count times...
humans <- read.table("humans.txt", header=T)
num_times <- length(unique(humans[,ncol(humans)]));
human_counts <- rep(0,num_times)
for (i in 1:num_times)
	human_counts[i] <- sum(humans[,ncol(humans)]==i)

# read in our design matrix for names
design <- read.table("design.txt", header=T)

################################################################
### PLOT 1                                                   ###
### VISUALISE DIRECTLY THE MCMC OUTPUT FOR PARAMETER F       ###
### (THE PROPORTION OF ISOLATES ATTRIBUTABLE TO EACH SOURCE) ###
################################################################
pdf(file=paste(fname,".pdf",sep=""), width=11, height=8)

par(mfrow=c(3,2))
for (t in 1:nt) {
plot(fmcmc[fd,((t-1)*ng+2)],type="l",ylim=c(0,1),col=col[o[1]],ylab="Proportion",main=paste("F, t=",t,sep=""))
for(i in 2:ng) lines(fmcmc[fd,((t-1)*ng+i+1)],col=col[o[i]])
}
for (t in 0:nt) {
plot(fmcmc[fd,nt*ng+t*(ng-1)+2],type="l",ylim=c(-3,3),col=col[o[1]],ylab="Proportion",main=paste("f, t=",t,sep=""))
for(i in 2:ng-1) lines(fmcmc[fd,nt*ng+t*(ng-1)+i+1],col=col[o[i]])
}

off <- nt*ng + (nt+1)*(ng-1)
plot(fmcmc[fd,off+2],type="l",ylim=range(fmcmc[fd,off+3*1:(ng-1)-1]),col=2,ylab="Value",main="mu")
for(i in 2:(ng-1)) lines(fmcmc[fd,off+3*i-1],col=i+1)
plot(fmcmc[fd,off+3],type="l",ylim=range(fmcmc[fd,off+3*1:(ng-1)+0]),col=2,ylab="Value",main="rho")
for(i in 2:(ng-1)) lines(fmcmc[fd,off+3*i+0],col=i+1)
plot(fmcmc[fd,off+4],type="l",ylim=range(fmcmc[fd,off+3*1:(ng-1)+1]),col=2,ylab="Value",main="tau")
for(i in 2:(ng-1)) lines(fmcmc[fd,off+3*i+1],col=i+1)

if (1) {
	pdf("mean_correlation.pdf", width=11, height=8)
	par(mfrow=c(2,2))
	l <- list()
	xlim <- NULL
	ylim <- NULL
	n_alpha <- ncol(design)+2
	for (i in 1:(ng-1)) {
		l[[i]] <- density(fmcmc[fd,off+n_alpha*(i-1)+2])
		xlim <- range(xlim, l[[i]]$x);
		ylim <- range(ylim, l[[i]]$y);
	}
	plot(l[[1]], xlim=xlim, ylim=ylim, col=col[o[1]],lwd=2,main="rho")
	for(i in 2:(ng-1)) lines(l[[i]],lwd=2,col=col[o[i]])

	xlim <- NULL
	ylim <- NULL
	for (i in 1:(ng-1)) {
		l[[i]] <- density(fmcmc[fd,off+n_alpha*(i-1)+3])
		xlim <- range(xlim, l[[i]]$x);
		ylim <- range(ylim, l[[i]]$y);
	}
	plot(l[[1]], xlim=xlim, ylim=ylim, col=col[o[1]],lwd=2,main="tau")
	for(i in 2:(ng-1)) lines(l[[i]],col=col[o[i]],lwd=2)

	for (k in 1:ncol(design)) {
		xlim <- NULL
		ylim <- NULL
		for (i in 1:(ng-1)) {
			l[[i]] <- density(fmcmc[fd,off+n_alpha*(i-1)+k+3])
			xlim <- range(xlim, l[[i]]$x);
			ylim <- range(ylim, l[[i]]$y);
		}
		plot(l[[1]], xlim=xlim, ylim=ylim, col=col[o[1]],lwd=2,main=names(design)[k])
		for(i in 2:(ng-1)) lines(l[[i]],col=col[o[i]],lwd=2)
	}
	dev.off()
}

plot(density(fmcmc[fd,off+2]), xlim=range(fmcmc[fd,off+3*1:(ng-1)-1]), col=2,main="mu")
for(i in 2:(ng-1)) lines(density(fmcmc[fd,off+3*i-1]),col=i+1)
plot(density(fmcmc[fd,off+3]), xlim=range(fmcmc[fd,off+3*1:(ng-1)+0]), col=2,main="rho")
for(i in 2:(ng-1)) lines(density(fmcmc[fd,off+3*i+0]),col=i+1)
plot(density(fmcmc[fd,off+4]), xlim=range(fmcmc[fd,off+3*1:(ng-1)+1]), col=2,main="tau")
for(i in 2:(ng-1)) lines(density(fmcmc[fd,off+3*i+1]),col=i+1)


#############################################################
### PLOT 2                                                ###
### HISTOGRAMS OF THE MARGINAL DISTRIBUTIONS OF F[i]      ###
### (THE PROPORTION OF ISOLATES ATTRIBUTABLE TO SOURCE i) ###
#############################################################
par(mfrow=c(2,2))
### PLOT THE HISTOGRAMS
#for (t in 1:nt) {
#	for(i in 1:ng)
#		hist(fmcmc[fd,(t-1)*ng+(1+i)],30,col=rainbow(ng)[i],main=paste(sc[i]," t=",t,sep=""),prob=T,xlim=c(0,1),xlab="Proportion")
#}

########################################################
### TABLE 1                                          ###
### SUMMARIES OF THE POSTERIOR DISTRIBUTIONS OF F[i] ###
########################################################
for (t in 1:nt) {
	df = fmcmc[fd,((t-1)*ng+2):(t*ng+1)]; names(df) <- sc[o];
	pe = apply(df,2,function(x)c("mean"=mean(x),"median"=median(x),"sd"=sd(x),quantile(x,c(.025,.975))))
	print(pe)
}
#################################################################
### TABLE 2                                                   ###
### PROBABILITY THAT EACH SOURCE IS RESPONSIBLE FOR THE MOST, ###
### SECOND MOST, THIRD MOST, ETC, NUMBER OF CASES             ###
#################################################################
if (0) {
	od = t(apply(df,1,order,decreasing=T))
### THE TABLE
	print(apply(od,2,function(x)table(factor(x,levels=1:ng,labels=sc)))/nrow(df))
	enc = apply(od,1,function(x) sum(x*((0:(ng-1))^ng)))
	encmode = as.numeric(levels(factor(enc)))[which.max(table(enc))]
### MOST LIKELY ORDER
	print(sc[od[which(enc==encmode)[1],]])
### POSTERIOR PROBABILITY OF THIS, THE MOST LIKELY ORDER
	print(max(table(enc))/nrow(df))
}


#################################################################################
### PLOT 3                                                                    ###
### BARCHART OF THE ESTIMATED PROPORTION OF CASES ATTRIBUTABLE TO EACH SOURCE ###
#################################################################################

# proportions
if (1) {
	pdf("proportions.pdf", width=8, height=5)
	x = c(0,seq(1.5,nt-1.5),nt,nt,seq(nt-1.5,1.5),0)
	g <- matrix(0,nt,ng)
	for (t in 1:nt) {
		df = fmcmc[fd,((t-1)*ng+2):(t*ng+1)]; names(df) <- sc[o];
		pe = apply(df,2,function(x)c("mean"=mean(x),"median"=median(x),"sd"=sd(x),quantile(x,c(.025,.975))))
		g[t,] = pe[1,];
	}

	prop = matrix(0, nrow = nt, ncol = ng+1)
	for(j in 1:ng)
		prop[,j+1] = prop[,j] + g[,j]
	print(prop)
	ymax <- 1

	alpha <- "7F"
	plot(NULL, xlim=c(0,nt), ylim=c(0,ymax),ylab="Proportion of human cases",xlab="", col.axis="transparent", xaxp=c(0,32,8), xaxs="i", yaxs="i")
	for(j in 1:ng)
	{  
	  y = c(prop[,j],rev(prop[,j+1]))
	  polygon(x,y,col=paste(col[o[j]],alpha,sep=""))
	}

	for (i in 1:8)
	   	mtext(2004+i, 1, line=0.5, at=4*i-2)
	for(i in 0:5)
		mtext(i/5,2,line=1,at=i/5)

	legend("bottomleft", legend=sc[o], fill=col[o])
	dev.off()
}

# totals
if (1) {
	pdf("totals.pdf", width=8, height=5)
	x = c(0,seq(1.5,nt-1.5),nt,nt,seq(nt-1.5,1.5),0)
	g <- matrix(0,nt,ng)
	for (t in 1:nt) {
		df = fmcmc[fd,((t-1)*ng+2):(t*ng+1)]*human_counts[t]; names(df) <- sc[o];
		pe = apply(df,2,function(x)c("mean"=mean(x),"median"=median(x),"sd"=sd(x),quantile(x,c(.025,.975))))
		g[t,] = pe[1,];
	}

	prop = matrix(0, nrow = nt, ncol = ng+1)
	for(i in 1:ng)
		prop[,i+1] = prop[,i] + g[,i]
	ymax <- 80

	alpha <- "7F"
	plot(NULL, xlim=c(0,nt), ylim=c(0,ymax),ylab="Cases",xlab="", col.axis="transparent", xaxp=c(0,32,8), xaxs="i", yaxs="i")
	for(j in 1:ng)
	{  
	  y = c(prop[,j],rev(prop[,j+1]))
	  polygon(x,y,col=paste(col[o[j]],alpha,sep=""))
	}

	for (i in 1:8)
	   	mtext(2004+i, 1, line=0.5, at=4*i-2)
	for(i in 0:(ymax/20))
		mtext(i*20,2,line=1,at=i*20)

	legend("topright", legend=sc[o], fill=col[o])
	dev.off()
}

# totals per source
for (j in 1:ng)
{
	pdf(paste("totals_",sc[o[j]],".pdf",sep=""), width=8, height=5)
	x = c(0,seq(1.5,nt-1.5),nt)

	pee <- matrix(0,nt,3)
	for (t in 1:nt) {
		df = fmcmc[fd,((t-1)*ng+2):(t*ng+1)]*human_counts[t]; names(df) <- sc[o];
		pe = apply(df,2,function(x)c("mean"=mean(x),quantile(x,c(.025,.975))))
		pee[t,] <- pe[,j]
  	}

	ymax <- 80

	plot(NULL, xlim=c(0,nt), ylim=c(0,ymax),ylab="Cases",xlab="", col.axis="transparent", xaxp=c(0,32,8), xaxs="i", yaxs="i", main="")
	y = c(pee[,2],rev(pee[,3]))
	polygon(c(x,rev(x)),y,col="grey80", border="grey80")
	lines(x, pee[,1], col="black", lwd="2")

	for (i in 1:8)
		mtext(2004+i, 1, line=0.5, at=4*i-2)
	for(i in 0:(ymax/20))
		mtext(i*20,2,line=1,at=i*20)
	text(2,ymax-2,sc[o[j]], adj=c(0,0))
	box()
	dev.off()
}
# ruminants...
if (0) {

	x = c(0,seq(1.5,nt-1.5),nt)

	for (j in 2:3) {

	pdf(paste("ruminant_cis_",sc[o[j]],".pdf",sep=""), width=8, height=3)
	ymax <- 30

	plot(NULL, xlim=c(0,nt), ylim=c(0,ymax),ylab="Cases",xlab="", col.axis="transparent", xaxp=c(0,32,8), xaxs="i", yaxs="i", main="")

	pee <- matrix(0,nt,3)
	for (t in 1:nt) {
		df = fmcmc[fd,((t-1)*ng+2):(t*ng+1)]*human_counts[t]; names(df) <- sc[o];
		pe = apply(df,2,function(x)c("mean"=mean(x),quantile(x,c(.025,.975))))
		pee[t,] <- pe[,j]
  	}

	y = c(pee[,2],rev(pee[,3]))
	polygon(c(x,rev(x)),y,col="grey80", border="grey80")
	y <- pee[,1]
	lines(x, y, col="black", lwd="2")
	for (i in 1:8)
		mtext(2004+i, 1, line=0.5, at=4*i-2)
	for(i in 0:(ymax/10))
		mtext(i*10,2,line=1,at=i*10)
	text(2,ymax-3,sc[o[j]])
	box()
	dev.off()
	}
}


#################################################################
### TABLE 3                                                   ###
### MOST PROBABLE ORDERINGS AND THEIR POSTERIOR PROBABILITIES ###
#################################################################
if (0) {
	enc_od = as.numeric(levels(factor(enc)))[order(table(enc),decreasing=TRUE)]
	enc_pr = table(enc)[order(table(enc),decreasing=TRUE)]/nrow(df)
	unenc_od = t(sapply(enc_od,function(i) sc[od[which(enc==i)[1],]]))
	print(cbind(unenc_od,"Posterior probability"=enc_pr)[enc_pr>.05,])

## Ordering of top 3 no. cases: full order
	good = apply(od[,1:3],1,function(x)all(sort(x)==c(1,2,5)))
	odg = od[good,1:3]
	encg = apply(odg,1,function(x) sum(x*((0:2)^3)))
	encg_od = as.numeric(levels(factor(encg)))[order(table(encg),decreasing=TRUE)]
	encg_pr = table(encg)[order(table(encg),decreasing=TRUE)]/nrow(df)
	unencg_od = t(sapply(encg_od,function(i) sc[odg[which(encg==i)[1],]]))
	print(cbind(unencg_od,encg_pr)[1:min(4,nrow(unencg_od)),])
}

#################################################
### PLOT 4                                    ###
### MCMC TRACE OF THE EVOLUTIONARY PARAMETERS ###
#################################################
### 3 BY 3 PANE
par(mfrow=c(2,2))
COL = col[c(o,6)];
for(i in 0:(ng-1)) {
  plot(mcmc$iter[gd],mcmc[[paste("A",i,0,"",sep=".")]][gd],type="l",col=COL[1],ylim=c(0,1),xlab="iter",ylab="M,R",main=sc[o[i+1]])
  for(j in 2:(ng+1)) lines(mcmc$iter[gd],mcmc[[paste("A",i,j-1,"",sep=".")]][gd],col=COL[j])
  lines(mcmc$iter[gd],mcmc[[paste("r",i,sep="")]][gd],col=col[7])
}

#################################################
### PLOT 5                                    ###
### PIE CHARTS OF THE EVOLUTIONARY PARAMETERS ###
#################################################
if (1) {
	COL = paste(col[c(o,6)],"7F",sep="");
	for(i in 0:(ng-1)) {
	pdf(paste("migration_mutation",i+1,".pdf",sep=""), width=3, height=3)
	par(omi=rep(0,4))
	  wh0 = which(names(mcmc)==paste("A",i,0,"",sep="."))
	  whng = which(names(mcmc)==paste("A",i,ng,"",sep="."))
	  pie(apply(mcmc[gd,wh0:whng],2,mean),col=COL,labels="",radius=1.05)
	dev.off()
	}
}

dev.off()

########################################################
### PLOT 6                                           ###
### POSTERIOR PROBABILITY OF SOURCE FOR EACH ISOLATE ###
########################################################
if (0) {
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
	od = rev(order(wm!=4,wm!=5,G[,1],G[,2],G[,3]))
	image(1:nrow(G),seq(0,1,len=res),t(tp2[,od]),col=COL[cod],ylab="Source probability",xlab="Human cases",bty="n")
}

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
colnames(g) <- sc[o]
print(g)

dev.off()
