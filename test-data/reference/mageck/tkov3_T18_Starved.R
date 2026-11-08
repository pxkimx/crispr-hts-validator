pdf(file='tkov3_T18_Starved.pdf',width=4.5,height=4.5);
gstable=read.table('tkov3_T18_Starved.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("POLR3H","RPE","POLR2L","RPL14","ELP5","COPA","UBE2S","RPP21","PRPF38A","RPL12")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='T18_A_Starved,T18_B_Starved,T18_C_Starved_vs_T0_- neg.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(282.73132040687995,0.0,0.0,0.0),c(250.06515191365477,1.2940582783324055,3.351029143882061,0.0),c(235.98490687347152,1.2940582783324055,1.1170097146273537,6.3157757266731585),c(189.2384933400631,0.0,5.585048573136769,0.0))
targetgene="POLR3H"
collabel=c("T0_-","T18_A_Starved","T18_B_Starved","T18_C_Starved")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(177.41108750630914,12.940582783324055,1.1170097146273537,0.0),c(290.05304782777523,7.7643496699944325,0.0,13.894706598680948),c(217.96219322203694,0.0,2.2340194292547073,0.0),c(306.9493418759952,16.822757618321273,10.053087431646183,6.3157757266731585))
targetgene="RPE"
collabel=c("T0_-","T18_A_Starved","T18_B_Starved","T18_C_Starved")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(321.59279671778575,0.0,0.0,15.15786174401558),c(406.07426695888535,6.470291391662028,31.2762720095659,3.789465436003895),c(157.69874445005257,0.0,0.0,6.3157757266731585),c(252.88120092169143,3.8821748349972163,7.819068002391475,0.0))
targetgene="POLR2L"
collabel=c("T0_-","T18_A_Starved","T18_B_Starved","T18_C_Starved")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(193.1809619513144,15.528699339988865,2.2340194292547073,0.0),c(482.6707999774823,68.58508875161749,14.521126290155598,13.894706598680948),c(232.60564806382754,2.588116556664811,16.755145719410304,13.894706598680948),c(242.7434244927595,6.470291391662028,10.053087431646183,50.52620581338527))
targetgene="RPL14"
collabel=c("T0_-","T18_A_Starved","T18_B_Starved","T18_C_Starved")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(337.362671162791,15.528699339988865,113.93499089199007,2.5263102906692634),c(284.98415961330926,0.0,0.0,1.2631551453346317),c(209.51404619792697,3.8821748349972163,68.13759259226858,11.368396308011686),c(146.9977582195133,0.0,7.819068002391475,0.0))
targetgene="ELP5"
collabel=c("T0_-","T18_A_Starved","T18_B_Starved","T18_C_Starved")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(323.84563592421506,49.17421457663141,3.351029143882061,15.15786174401558),c(363.83353183833555,58.232622524958245,22.340194292547075,40.420964650708214),c(296.81156544706323,45.29203974163419,39.095340011957376,10.105241162677054),c(256.82366953294274,2.588116556664811,13.404116575528244,29.052568342696528))
targetgene="COPA"
collabel=c("T0_-","T18_A_Starved","T18_B_Starved","T18_C_Starved")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(247.81231270722546,0.0,6.702058287764122,1.2631551453346317),c(215.7093540156076,6.470291391662028,20.106174863292367,10.105241162677054),c(483.23400977908966,1.2940582783324055,17.87215543403766,1.2631551453346317),c(183.60639532398977,3.8821748349972163,8.93607771701883,101.05241162677054))
targetgene="UBE2S"
collabel=c("T0_-","T18_A_Starved","T18_B_Starved","T18_C_Starved")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(273.7199635811627,266.57600533647553,30.15926229493855,34.105188924035055),c(126.15899556004206,0.0,0.0,23.999947761358),c(207.82441679310497,40.11580662830457,10.053087431646183,0.0),c(152.6298562355866,0.0,0.0,0.0))
targetgene="RPP21"
collabel=c("T0_-","T18_A_Starved","T18_B_Starved","T18_C_Starved")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(138.54961119540332,0.0,3.351029143882061,2.5263102906692634),c(189.2384933400631,19.410874174986084,1.1170097146273537,2.5263102906692634),c(318.7767477097491,2.588116556664811,35.74431086807532,27.789413197361895),c(192.05454234809972,5.176233113329622,11.170097146273537,137.68391084147484))
targetgene="PRPF38A"
collabel=c("T0_-","T18_A_Starved","T18_B_Starved","T18_C_Starved")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(117.14763873432476,2.588116556664811,5.585048573136769,0.0),c(122.77973675039807,1.2940582783324055,0.0,0.0),c(185.29602472881177,6.470291391662028,0.0,10.105241162677054),c(83.35505063788493,0.0,2.2340194292547073,0.0))
targetgene="RPL12"
collabel=c("T0_-","T18_A_Starved","T18_B_Starved","T18_C_Starved")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=9
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("DOT1L","LacZ","KCTD5","ARID4A","INPPL1","MIB1","PMAIP1","ITM2A","MYH9","KMT2E")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='T18_A_Starved,T18_B_Starved,T18_C_Starved_vs_T0_- pos.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(41.114315517335136,247.16513116148946,156.3813600478295,354.9465958390315),c(254.57083032651343,1023.6000981609327,1053.3401608935944,1024.4188228663863),c(137.986401393796,262.6938305014783,576.3770127477145,635.3670381033197),c(197.68664036417303,700.0855285778314,776.3217516660108,804.6298275781604))
targetgene="DOT1L"
collabel=c("T0_-","T18_A_Starved","T18_B_Starved","T18_C_Starved")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(41.114315517335136,129.40582783324055,48.03141772897621,66.94722270273547),c(254.57083032651343,376.57095899473,371.96423497090876,160.42070345749823),c(56.884189962340386,178.58004240987196,159.73238919171158,200.84166810820645),c(96.87208587646086,284.6928212331292,201.06174863292367,44.21043008671211),c(114.33158972628812,219.98990731650895,166.4344474794757,64.42091241206622),c(176.2846679030945,201.87309141985526,121.75405889438154,260.2099599389341),c(106.44665250378549,370.10066760306796,195.47670005978688,229.89423645090298),c(524.3483252964247,557.7391179612667,537.2816727357571,361.26237156570465),c(97.43529567806819,64.70291391662028,151.9133211893201,16.421016889350213),c(347.50044759172295,208.34338281151727,399.8894778365926,151.5786174401558),c(157.69874445005257,155.28699339988867,258.0292440789187,266.5257356656073),c(115.45800932950277,186.3443920798664,272.5503703690743,627.788107231312),c(269.77749496991134,335.161094088093,198.82772920366895,217.26268499755665),c(238.23774607990083,328.690802696431,443.4528567070594,443.3674560124557),c(72.65406440734564,85.40784636993877,84.89273831167888,61.89460212139695),c(34.91900769965449,15.528699339988865,27.925242865683842,21.47363747068874),c(171.2157796886285,164.3454013482155,306.0606618078949,318.31509662432717),c(112.64196032146612,177.28598413153955,235.6890497863716,106.10503220810907),c(219.0886128252516,241.98889804815983,195.47670005978688,118.73658366145538),c(171.2157796886285,309.2799285214449,584.1960807501059,597.4723837432808),c(212.33009520596363,232.93049009983298,148.56229204543803,183.1574960735216),c(157.13553464844523,194.10874174986083,206.64679720606043,75.7893087200779),c(206.69799718989032,234.2245483781654,226.75297206935278,227.3679261602337),c(318.2135379081418,520.211427889627,349.62404067836167,325.89402749633496),c(172.9054090934505,156.58105167822106,269.19934122519226,319.5782517696618),c(21.965182262685893,67.29103047328509,77.07367030928741,13.894706598680948),c(361.58069263190623,395.9818331697161,621.0574013328086,466.10424862847907),c(435.9243864440739,652.2053722795324,336.21992410283343,530.5251610405453),c(237.67453627829352,341.63138547975507,284.8374772299752,416.8411979604285),c(129.53825436968603,223.87208215150616,444.56986642168675,306.9467003163155),c(215.7093540156076,432.21546496302346,350.74105039298905,447.1569214484596),c(59.70023897037704,209.6374410898497,251.32718579115456,95.999791045432),c(149.2505974259426,279.51658811979956,184.30660291351336,106.10503220810907),c(73.7804840105603,169.52163446154512,350.74105039298905,190.7364269455294),c(113.76837992468077,219.98990731650895,510.4734395847006,209.68375412554886),c(326.09847513064443,194.10874174986083,148.56229204543803,402.9464913617475),c(152.6298562355866,165.6394596265479,332.8688949589514,248.84156363092245),c(225.28392064293223,267.8700636148079,500.4203521530544,132.63129026013632),c(117.14763873432476,151.40481856489146,96.06283545795242,202.10482325354107),c(171.2157796886285,288.57499606812644,224.5189526400981,385.26231932706264),c(50.68888214465975,133.28800266823777,99.41386460183448,95.999791045432),c(59.70023897037704,76.34943842161192,102.76489374571653,160.42070345749823),c(92.36640746360221,172.10975101820992,175.37052519649453,97.26294619076664),c(295.1219360422412,278.2225298414672,312.762720095659,452.2095420297982),c(68.148385994487,104.81872054492484,100.53087431646183,319.5782517696618),c(198.24985016578037,137.17017750323498,156.3813600478295,125.05235938812854),c(43.930364525371786,91.8781377616008,177.60454462574924,165.47332403883675),c(177.97429730791646,139.75829405989978,180.9555737696313,296.8414591536384),c(108.13628190860747,59.52668080329065,116.16901032124478,197.05220267220255),c(105.32023290057082,238.1067232131626,150.79631147469274,281.68359740962285),c(56.884189962340386,59.52668080329065,77.07367030928741,49.26305066805064),c(26.470860675544536,64.70291391662028,50.26543715823092,0.0),c(56.32098016073306,1.2940582783324055,1.1170097146273537,18.947327180019474),c(235.98490687347152,586.2084000845797,524.9945658748562,471.1568692098176),c(153.19306603719392,245.87107288315704,72.60563145077799,125.05235938812854),c(390.86760231548743,245.87107288315704,428.9317304169038,113.68396308011685),c(3.9424686112513143,77.64349669994434,42.44636915583944,75.7893087200779),c(162.2044228629112,138.4642357815674,256.91223436429135,106.10503220810907),c(181.35355611756046,390.8056000563865,340.68796296134286,198.31535781753718),c(58.57381936716238,72.4672635866147,65.90357316301386,16.421016889350213),c(61.953078176806365,100.93654570992763,25.691223436429134,53.05251610405453),c(56.884189962340386,113.87712849325169,145.21126290155598,324.63087235100033),c(142.49207980665463,191.52062519319603,120.6370491797542,381.4728538910588),c(208.95083639631966,163.0513430698831,136.27518518453715,251.3678739215917),c(121.6533171471834,315.75021991310695,406.5915361243567,185.68380636419087),c(69.27480559770166,121.64147816324612,144.09425318692863,185.68380636419087),c(142.49207980665463,157.87510995655347,338.45394353208815,166.73647918417137),c(283.2945302084873,288.57499606812644,351.8580601076164,285.47306284562677),c(382.9826650929848,517.6233113329622,555.1538281697948,553.2619536565687),c(153.75627583880126,372.6887841597328,201.06174863292367,373.893923019051),c(101.94097409092683,148.81670200822663,73.72264116540534,109.89449764411296),c(97.99850547967552,159.1691682348859,170.9024863379851,122.52604909745928),c(136.85998179058134,117.7593033282489,68.13759259226858,93.47348075476275),c(119.96368774236142,135.87611922490257,136.27518518453715,113.68396308011685),c(85.04468004270691,51.76233113329622,59.20151487524974,32.842033778700426),c(101.94097409092683,181.16815896653677,262.4972829374281,263.999425374938),c(156.00911504523057,385.6293669430568,168.6684669087304,458.5253177564713),c(174.03182869666514,146.22858545156183,147.44528233081067,108.63134249877832),c(41.677525318942465,256.2235391098163,146.32827261618334,65.68406755740085),c(526.0379547012468,491.7421457663141,711.5351882176243,671.9985373180241),c(28.16049008036653,251.04730599648667,73.72264116540534,34.105188924035055),c(210.64046580114163,952.4268928526504,414.4106041267482,357.4729061297008),c(208.38762659471232,401.1580662830457,264.7313023666828,358.7360612750354),c(274.8463831843773,151.40481856489146,158.61537947708422,287.999373136296),c(215.1461442140003,166.93351790488032,390.95340011957376,367.5781472923778),c(431.98191783282255,601.7370994245686,373.08124468553615,199.5785129628718),c(129.53825436968603,260.1057139448135,177.60454462574924,183.1574960735216),c(360.45427302869155,236.8126649348302,330.63487552969667,341.05188924035053),c(47.8728331366231,80.23161325660914,97.17984517257977,157.89439316682896),c(40.5511057157278,124.22959471991092,135.1581754699098,194.52589238153328),c(72.09085460573831,318.3383364697718,120.6370491797542,269.0520459562766),c(223.59429123811023,557.7391179612667,368.6132058270267,434.5253699951133),c(248.3755225088328,772.552792164446,804.2469945316947,1067.3660978077637),c(273.15675377955534,367.5125510464032,195.47670005978688,327.15718264166964))
targetgene="LacZ"
collabel=c("T0_-","T18_A_Starved","T18_B_Starved","T18_C_Starved")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(211.7668854043563,500.80055371464096,890.2567425580008,708.6300365327284),c(250.62836171526212,542.2104186212779,464.67604128497914,1393.2601253040987),c(166.14689147416252,657.381605392862,573.0259836038324,346.1045098216891),c(135.73356218736666,465.86098019966596,527.2285853041109,635.3670381033197))
targetgene="KCTD5"
collabel=c("T0_-","T18_A_Starved","T18_B_Starved","T18_C_Starved")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(317.0871183049271,1004.1892239859467,567.4409350306956,1525.891415564235),c(105.88344270217814,328.690802696431,423.34668184376704,574.7355911272574),c(155.44590524362323,669.0281298978537,623.2914207620634,275.3678216829497),c(474.22265295337235,900.6645617193542,693.6630327835866,495.1568169711756))
targetgene="ARID4A"
collabel=c("T0_-","T18_A_Starved","T18_B_Starved","T18_C_Starved")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(156.57232484683792,900.6645617193542,988.553597445208,762.9457077821176),c(184.1696051255971,392.09965833471887,434.51677899004056,630.3144175219812),c(140.23924060022532,653.4994305578648,562.9728961721862,976.4189273436702),c(336.7994613611837,749.2597431544627,936.0541408577224,1504.4177780935463))
targetgene="INPPL1"
collabel=c("T0_-","T18_A_Starved","T18_B_Starved","T18_C_Starved")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(100.25134468610484,217.40179075984412,370.84722525628143,287.999373136296),c(73.7804840105603,153.99293512155626,250.2101760765272,303.1572348803116),c(221.34145203168092,735.0251020928064,598.7172070402615,464.84109348314445),c(131.22788377450803,703.9677034128285,602.0682361841436,319.5782517696618))
targetgene="MIB1"
collabel=c("T0_-","T18_A_Starved","T18_B_Starved","T18_C_Starved")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(266.39823616026734,1017.1298067692708,1054.4571706082218,501.47259269784877),c(252.88120092169143,1024.894156439265,1714.6099119529879,813.4719135955028),c(122.77973675039807,446.4501060246799,440.1018275631773,363.78868185637396),c(137.986401393796,204.46120797652006,220.05091378158866,188.21011665486012))
targetgene="PMAIP1"
collabel=c("T0_-","T18_A_Starved","T18_B_Starved","T18_C_Starved")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(142.49207980665463,320.92645302643655,359.67712811000786,385.26231932706264),c(253.44441072329877,762.2003259377868,644.5146053399831,486.3147309538332),c(117.14763873432476,569.3856424662584,390.95340011957376,271.5783562469458),c(134.607142584152,389.51154177805404,249.09316636189988,485.05157580849857))
targetgene="ITM2A"
collabel=c("T0_-","T18_A_Starved","T18_B_Starved","T18_C_Starved")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(77.72295262181161,751.8478597111276,416.6446235560029,197.05220267220255),c(165.5836816725552,410.21647423137256,671.3228384910395,712.4195019687322),c(216.8357736188223,631.5004398262139,629.9934790498274,629.0512623766466),c(198.24985016578037,87.99596292660357,269.19934122519226,135.1576005508056))
targetgene="MYH9"
collabel=c("T0_-","T18_A_Starved","T18_B_Starved","T18_C_Starved")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(139.11282099701066,465.86098019966596,158.61537947708422,335.99926865901205),c(391.4308121170948,535.7401272296158,1005.3087431646182,759.1562423461137),c(86.73430944752891,380.4531338297272,245.7421372180178,274.10466653761506),c(176.8478777047018,767.3765590511165,443.4528567070594,309.47301060698476))
targetgene="KMT2E"
collabel=c("T0_-","T18_A_Starved","T18_B_Starved","T18_C_Starved")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



dev.off()
Sweave("tkov3_T18_Starved_summary.Rnw");
library(tools);

texi2dvi("tkov3_T18_Starved_summary.tex",pdf=TRUE);

