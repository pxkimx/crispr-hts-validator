pdf(file='demo3.pdf',width=4.5,height=4.5);
gstable=read.table('demo3.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("C1orf109","ANAPC4","C7orf26","C3orf17","C14orf80","C16orf80","ANAPC2","C9orf78","C9orf41","BUB3")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='HL60.final,KBM7.final_vs_HL60.initial,KBM7.initial neg.'


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
targetmat=list(c(138.45596176113995,80.98170955333637,21.389263381426552,0.0),c(1042.7464620135852,1037.2550883214574,173.70674503703987,90.87354691914342),c(3386.4020647412144,3607.9936128656673,886.0340315881847,387.3776198796819),c(905.7327498541238,606.501314101583,95.27944597180918,0.0),c(1711.9502771924283,1025.1939826433008,405.7478447506976,387.3776198796819),c(711.0290536275207,654.745736814209,211.94815532504492,178.2519574183198),c(584.1110886798092,563.4259366795956,31.111655827529532,0.0),c(1110.5321932924767,1538.652481513391,154.91011964124078,34.368841463009375),c(237.9711842769593,399.73950247604336,73.89018259038264,190.4849348882045),c(1256.1994030620092,1611.01911558233,75.83466107960324,390.2902335629878))
targetgene="C1orf109"
collabel=c("HL60.initial","KBM7.initial","HL60.final","KBM7.final")

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
targetmat=list(c(866.7920106088033,1399.088258666152,493.24937676562445,413.5911430294348),c(181.72344981149618,175.74753988170872,284.5420189226138,174.75682099835274),c(340.3709059961357,439.3688497042718,40.18588877722564,32.62127325302585),c(764.3922888896268,1144.0820243279861,147.78036518076527,48.93190987953877),c(503.34511098581083,356.66412505405594,1060.3889361216316,369.90193777984666),c(475.9423685539186,380.7863364103689,34.352453309563856,0.0),c(517.7676070025963,716.7742803018708,217.7815907927067,29.708659569719966),c(1051.3999596236565,982.1186052213135,249.54140611664312,27.96109135973644),c(771.6035368980195,1087.2225261309627,662.4190053278162,91.4560696558046),c(1.442249601678541,218.82291730369616,41.48220777003937,28.543614096397615))
targetgene="ANAPC4"
collabel=c("HL60.initial","KBM7.initial","HL60.final","KBM7.final")

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
targetmat=list(c(1468.2100945087548,1910.8237424393624,367.5064344626926,170.09663910506333),c(278.35417312395845,752.9575973363403,176.29938302266734,266.7954133908185),c(426.9058820968482,1061.3772996777702,262.5045960447804,81.55318313256461),c(90.86172490574809,468.66010635122325,57.686195180211,32.62127325302585),c(262.48942750549446,606.501314101583,312.412877268109,5.242704629950582),c(711.0290536275207,603.055283907824,232.04109971365776,70.48525113600228),c(223.54868826017386,890.7988050867001,240.467173166947,0.0),c(137.0137121594614,268.7903551132016,82.96441554007875,15.145591153190571),c(545.1703494344886,263.6213098225631,152.31748165561334,40.19406882962113),c(760.0655400845911,704.7131746237144,284.5420189226138,930.8713331845589))
targetgene="C7orf26"
collabel=c("HL60.initial","KBM7.initial","HL60.final","KBM7.final")

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
targetmat=list(c(1012.4592203783359,1088.9455412278423,281.3012214405795,75.72795576595286),c(451.4241253253834,904.582925861736,610.566245615267,545.2412815148606),c(359.12015081795676,248.1141739506476,117.3168688496426,230.09648098116446),c(787.4682825164834,597.8862386171855,280.00490244776574,53.59209177282818),c(1129.2814381142978,1542.09851170715,653.3447723781202,797.4736264891498),c(315.8526627676005,341.1569891821405,123.79846381371125,18.640727573157626),c(807.659776939983,1166.4812205874196,126.39110179933871,244.07702666103268),c(15.864745618463953,335.98794389150197,9.074232949696112,1.1650454733223516),c(786.0260329148049,759.8496577238583,555.4726884206835,321.55255063696904),c(253.83592989542325,279.12844569447856,39.53772928081878,57.08722819279523))
targetgene="C3orf17"
collabel=c("HL60.initial","KBM7.initial","HL60.final","KBM7.final")

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
targetmat=list(c(1025.4394667934428,1023.4709675464214,88.79785100774053,76.31047850261403),c(246.62468188703053,136.1181926534803,0.0,0.0),c(2783.5417312395844,2908.4494835325913,694.8269801481595,139.8054567986822),c(641.8010727469508,1109.6217223903961,473.8045918734185,222.52368540456916),c(1909.5384726223883,1862.5793197267365,1130.390161733573,550.4839861448112),c(126.91796494771162,460.04503086682575,246.30060863460878,387.3776198796819),c(338.92865639445716,206.76181162553968,60.27883316583846,22.13586399312468),c(199.03044503163866,196.4237210442627,58.98251417302473,72.2328193459858),c(1185.5291725797608,863.2305635366281,358.43220151299647,361.164096729929),c(1039.861962810228,1007.9638316745059,797.2361805804442,941.9392651811213))
targetgene="C14orf80"
collabel=c("HL60.initial","KBM7.initial","HL60.final","KBM7.final")

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
targetmat=list(c(1828.7724949283902,1493.854088994524,1099.9266654024502,402.5232110328725),c(722.5670504409491,1207.8335829125276,777.1432361918314,506.2122581585618),c(360.56240041963525,124.0570869753238,51.20460021614235,215.53341256463506),c(582.6688390781306,921.813076830531,472.50827288060475,477.66864406216416),c(1759.5445140478203,994.1797108994699,696.1232991409732,191.64998036152684),c(535.0746022227388,692.6520689455579,478.98986784467337,163.10636626512922),c(1651.3757939219295,1304.3224283377795,1221.132491230534,2020.1888507409578),c(437.00162930859796,790.8639294676892,125.09478280652499,4.077659156628231),c(861.023012202089,1209.556598009407,221.02238827474105,81.55318313256461),c(654.7813191620577,530.6886498388851,81.66809654726502,0.0))
targetgene="C16orf80"
collabel=c("HL60.initial","KBM7.initial","HL60.final","KBM7.final")

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
targetmat=list(c(480.2691173589542,472.10613654498223,297.5052088507511,717.0854888299075),c(872.5610090155174,415.2466383479588,1018.2585688551852,551.6490316181334),c(1419.1736080516844,1581.7278589353784,1803.179719003899,1287.3752480211986),c(1639.8377971085013,2708.5797322945696,1140.112554179676,344.27093736675494),c(790.3527817198406,1280.2002169814664,776.4950766954245,895.3374462482273),c(51.920985660427476,80.98170955333637,81.01993705085815,0.0),c(842.273767380268,787.4178992739303,166.57699057656436,136.8928431153763),c(1797.0430036914622,1249.1859452376355,672.7895572703261,218.44602624794092),c(572.5730918663808,441.0918648011513,86.85337251851993,0.0),c(406.7143876733486,608.2243291984626,53.797238201769815,0.0))
targetgene="ANAPC2"
collabel=c("HL60.initial","KBM7.initial","HL60.final","KBM7.final")

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
targetmat=list(c(1140.819434927726,1076.8844355496858,459.54508295246745,410.678529346129),c(956.2114859128727,1125.1288582623117,887.3303505809985,838.250218055432),c(1070.1492044454776,1347.397805759767,701.956734608635,552.2315543547946),c(1882.1357301904961,1402.5342888599107,1083.7226779922787,401.3581655595501),c(1195.6249197915106,1433.5485606037416,499.73097172969307,226.6013445611974),c(1396.0976144248277,1207.8335829125276,734.3647094289782,376.3096878831196),c(497.5761125790967,484.1672422231387,553.5282099314629,298.251641170522),c(128.36021454939015,118.8880416846853,77.13098007241696,69.9027283993411),c(641.8010727469508,199.86975123802168,186.02177546877033,681.5516018935757),c(837.9470185752324,963.165439155639,119.2613473388632,11.06793199656234))
targetgene="C9orf78"
collabel=c("HL60.initial","KBM7.initial","HL60.final","KBM7.final")

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
targetmat=list(c(344.6976548011713,198.14673614114218,17.50030640298536,53.59209177282818),c(483.15361656231124,508.2894535794517,568.4358783488208,286.6011864372985),c(868.2342602104817,1779.8745950765208,792.6990641055961,782.3280353359592),c(829.2935209651612,616.83940468286,263.80091503759417,231.84404919114797),c(219.22193945513825,329.095883503984,23.33374187064715,92.62111512912695),c(363.44689962299236,223.99196259433464,182.780977986736,460.77548469899006),c(1015.343719581693,1393.9192133755132,389.543857340526,233.00909466447033),c(1460.998846500362,1356.0128812441644,603.4364911547915,75.72795576595286),c(203.3571938366743,806.3710653396047,193.15152992924584,131.65013848542574),c(390.8496420548846,260.1752796288041,171.11410705141242,148.54329784859982))
targetgene="C9orf41"
collabel=c("HL60.initial","KBM7.initial","HL60.final","KBM7.final")

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
targetmat=list(c(413.9256356817413,742.6195067550633,814.0883274870228,601.7459869709946),c(493.24936377406107,451.42995538242826,555.4726884206835,925.0461058179473),c(876.887757820553,628.9005103610165,676.0303547523604,109.51427449230106),c(1547.5338226010747,897.6908654742181,200.9294438861282,0.0),c(562.477344654631,585.8251329390291,29.167177338308935,0.0),c(865.3497610071247,721.9433255925094,251.4858846058637,187.57232120489863),c(1296.5823919090085,932.151167411808,563.8987618739727,676.8914200002863),c(1080.2449516572274,697.8211142361964,849.7370997894003,1366.5983402071186),c(591.3223366882019,823.6012163083997,289.72729489386876,76.89300123927521),c(744.2007944661273,687.4830236549194,89.4460105041474,77.47552397593638))
targetgene="BUB3"
collabel=c("HL60.initial","KBM7.initial","HL60.final","KBM7.final")

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
targetgenelist=c("BIN1","C7orf34","CAMK1G","CASP8","ANKK1","C9orf69","CAMKV","BRSK2","BRE","BNC2")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='HL60.final,KBM7.final_vs_HL60.initial,KBM7.initial pos.'


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
targetmat=list(c(315.8526627676005,199.86975123802168,427.785267628531,358.8340057832843),c(552.3815974428812,180.91658517234723,624.825754536218,301.7467775904891),c(480.2691173589542,218.82291730369616,405.7478447506976,461.35800743565125),c(569.6885926630238,584.1021178421496,1064.2778931000728,674.5613290536417),c(1768.1980116578914,1483.5159984132472,2085.129099940885,2054.5576922039672),c(666.319315975486,304.97367214767104,595.0104177015022,894.1724007749049),c(100.95747211749787,184.3626153661062,338.9874166207905,341.94084642011023),c(356.23565161459965,261.89829472568357,276.76410496573146,523.687940258397),c(473.0578693505615,170.57849459107024,352.5987660453347,498.6394625819665),c(415.3678852834198,501.3973931919337,727.2349549685027,710.0952159899733))
targetgene="BIN1"
collabel=c("HL60.initial","KBM7.initial","HL60.final","KBM7.final")

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
targetmat=list(c(256.7204290987803,484.1672422231387,633.2518279895073,496.3093716353218),c(585.5533382814876,1056.2082543871318,1101.871143891671,1142.909609329227),c(1283.6021454939016,1338.7827302753694,1681.3257336794084,1359.6080673671843),c(693.7220584073783,885.6297597960616,1199.7432278491076,963.4926064375848),c(641.8010727469508,161.96341910667275,692.882501658939,644.2701467472605),c(1364.3681231879,1616.1881608729684,2162.2600800133023,2010.2859642177177),c(514.8831077992392,516.9045290638492,663.0671648242231,666.4060107403851),c(4.326748805035623,34.460301937589946,64.81594964068653,230.09648098116446),c(987.9409771498007,802.9250351458458,763.5318867672872,1432.4234094498313),c(449.98187572370483,361.83317034469445,679.9193117308016,299.4166866438444))
targetgene="C7orf34"
collabel=c("HL60.initial","KBM7.initial","HL60.final","KBM7.final")

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
targetmat=list(c(432.67488050356235,323.9268382133455,409.6368017291388,505.04721268523946),c(856.6962633970534,587.5481480359085,1264.5591774897941,779.4154216526532),c(1367.2526223912569,1464.5628323475728,1794.1054860542029,1699.801345577311),c(888.4257546339813,718.4972953987503,1754.5677567733842,901.7451963515002),c(586.9955878831662,410.07759305732037,784.9211501487138,630.2896010673923),c(364.8891492246709,535.8576951295237,769.365322234949,419.9988931327078),c(1044.1887116152639,580.6560876483906,1105.1119413737051,957.6673790709731),c(556.7083462479169,616.83940468286,552.2318909386491,1224.4627924617917),c(1374.4638703996497,1611.01911558233,1709.1965920249036,1896.6940305687885),c(1273.5063982821518,1051.0392090964933,1693.6407641111389,1716.694504940485))
targetgene="CAMK1G"
collabel=c("HL60.initial","KBM7.initial","HL60.final","KBM7.final")

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
targetmat=list(c(902.8482506507668,477.27518183562074,591.7696202194679,1702.713959260617),c(826.4090217618041,882.1837296023026,1363.7275804400444,2459.4109941834845),c(584.1110886798092,520.3505592576082,754.4576538175911,1187.1813373154764),c(561.0350950529524,1361.1819265348029,797.2361805804442,1360.7731128405067),c(516.3253574009177,410.07759305732037,222.96686676396163,998.4439706372553),c(542.2858502311315,508.2894535794517,469.9156348949773,808.541558485712),c(1383.1173680097208,933.8741825086876,820.5699224510913,2387.7606975741596),c(755.7387912795556,1011.4098618682649,578.1582707949237,1843.6844615326215),c(522.0943558076319,351.49507976341744,869.8300441780132,939.6091742344767),c(312.9681635642434,180.91658517234723,285.1901784190207,318.057414217002))
targetgene="CASP8"
collabel=c("HL60.initial","KBM7.initial","HL60.final","KBM7.final")

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
targetmat=list(c(105.2842209225335,406.6315628635614,486.1196223051489,355.92139209997845),c(1058.6112076320492,789.1409143708097,1078.5374020210238,2276.498854871875),c(894.1947530406954,963.165439155639,596.9548961907228,682.1341246302369),c(37.49848964364207,31.01427174383095,508.1570451829823,774.7552397593638),c(1022.5549675900857,892.5218201835796,1739.6600883560263,1325.8217486408362),c(1191.298170986475,1952.1761047644704,2503.84013461972,2313.197787281529),c(1032.6507148018354,534.1346800326442,707.14201057989,650.6778968505334),c(892.7525034390169,970.057499543157,1655.3993538231337,1939.2181903450544),c(187.49244821821034,544.4727706139212,650.7521343924926,371.6495059898302),c(563.9195942563096,372.1712609259714,462.78588043450173,324.46516432027494))
targetgene="ANKK1"
collabel=c("HL60.initial","KBM7.initial","HL60.final","KBM7.final")

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
targetmat=list(c(1042.7464620135852,785.6948841770508,1168.631572021578,1036.3079485202318),c(2294.619116270559,1380.1350926004773,1826.5134608745461,2665.6240429615405),c(242.2979330819949,408.35457796044085,599.5475341763503,566.2121000346629),c(331.71740838606445,480.7212120293797,1070.7594880641414,294.1739820138938),c(627.3785767301654,549.6418159045596,431.67422460697225,501.5520762652724),c(409.59888687670565,696.0980991393169,514.638640147051,877.861764148392),c(60.574483270498725,172.30150968794973,233.98557820287834,262.1352314975291),c(930.250993082659,580.6560876483906,1261.9665395041666,1050.2884942001),c(5.768998406714164,1.7230150968794973,603.4364911547915,387.3776198796819),c(2487.8805628954833,2126.2006295492997,2864.2168146219374,2246.207672565494))
targetgene="C9orf69"
collabel=c("HL60.initial","KBM7.initial","HL60.final","KBM7.final")

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
targetmat=list(c(138.45596176113995,168.85547949419075,256.6711605771186,572.0373274012746),c(1572.0520658296098,1144.0820243279861,1564.6570243261726,1798.2476880730499),c(1534.5535761859678,1652.3714779074378,1829.1060988601737,3086.2054588309097),c(1120.6279405042264,990.733680705711,869.8300441780132,1149.8998821691612),c(556.7083462479169,597.8862386171855,648.1594964068652,639.6099648539711),c(413.9256356817413,184.3626153661062,305.9312823040404,272.0381180207691),c(245.18243228535198,525.5196045482467,976.7763610851458,1165.0454733223517),c(372.10039723306363,199.86975123802168,1217.2435342520928,548.1538951981664),c(895.637002642374,639.2386009422935,1157.6128605826614,756.1145121862062),c(631.705325535201,1109.6217223903961,1339.7456790729905,848.7356273153332))
targetgene="CAMKV"
collabel=c("HL60.initial","KBM7.initial","HL60.final","KBM7.final")

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
targetmat=list(c(256.7204290987803,503.1204082888132,672.1413977739193,640.1924875906323),c(1384.5596176113995,828.7702615990382,1858.2732761984826,1511.06397889909),c(1215.8164142150101,761.5726728207378,1479.0999708004665,1251.8413610848668),c(445.6551269186692,454.8759855761873,434.9150220890066,376.89221061978077),c(1106.2054444874411,1847.072183854821,2117.5370747612287,1457.471887126262),c(927.366493879302,873.5686541179051,970.9429256174841,1693.3935954740382),c(907.1749994558023,685.7600085580399,1399.376352742422,1186.016291842154),c(533.6323526210602,623.731465070378,886.0340315881847,353.59130115333375),c(1018.22821878505,1648.925447713679,1508.2671481387754,1707.9566638905676),c(1026.8817163951212,1171.6502658780582,1418.172978138221,1384.6565450436149))
targetgene="BRSK2"
collabel=c("HL60.initial","KBM7.initial","HL60.final","KBM7.final")

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
targetmat=list(c(561.0350950529524,1080.3304657434448,1199.0950683527005,1589.1220256116876),c(230.75993626856658,513.4584988700902,353.2469255417415,611.6488734942346),c(532.1901030193817,727.1123708831478,963.1650116606017,1053.201107883406),c(385.08064364817045,270.51337021008106,412.87759921117316,186.40727573157625),c(478.82686775727564,318.757792922707,499.73097172969307,621.5517600174746),c(553.8238470445598,777.0798086926533,968.3502876318566,823.1046269022414),c(83.65047689735539,141.28723794411877,109.53895489276022,279.0283908607032),c(1045.6309612169423,572.0410121639931,933.9978343222928,1106.793199656234),c(1008.1324715733002,1331.8906698878513,1298.911630799358,1689.31593631741),c(1100.436446080727,585.8251329390291,1029.277280294102,929.7062877112367))
targetgene="BRE"
collabel=c("HL60.initial","KBM7.initial","HL60.final","KBM7.final")

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
targetmat=list(c(1058.6112076320492,1287.0922773689845,1667.0662247584573,1364.8507719971349),c(1319.6583855358651,833.9393068896767,1097.334027416823,1165.0454733223517),c(343.25540519949277,630.623525457896,1284.003962382,860.9686047852178),c(383.63839404649195,501.3973931919337,809.5512110121747,669.3186244236911),c(549.4970982395241,682.313978364281,551.5837314422423,520.19280383843),c(279.796422725637,330.81889860086346,288.43097590105504,233.00909466447033),c(620.1673287217727,930.4281523149285,1070.1113285677345,2024.8490326342471),c(575.4575910697379,584.1021178421496,1053.259181661156,956.5023335976507),c(449.98187572370483,1014.8558920620239,923.627282379783,1247.7637019282386),c(184.60794901485326,399.73950247604336,615.7515215865219,581.9402139245146))
targetgene="BNC2"
collabel=c("HL60.initial","KBM7.initial","HL60.final","KBM7.final")

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
Sweave("demo3_summary.Rnw");
library(tools);

texi2dvi("demo3_summary.tex",pdf=TRUE);

