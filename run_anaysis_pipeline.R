
repoDir = "C:/Users/user/pCloud local/boulot/Github/TrophicNetEncoder/"
#repoDir = "/home/christophe/pCloud local/boulot/Github/EcoGraph Encoder/Simulate networks/trophicSBM9/"
source(paste(repoDir,'functions_to_source.R',sep=""))

masterDir="C:/Users/user/pCloud local/boulot/data/Simu_networks/trophicSBM10/"
#masterDir='/home/christophe/pCloud local/boulot/data/Simu_Networks/trophicSBM9/'
graphsDir = paste(masterDir,'Graphs/',sep="")
ebdDir =paste(masterDir,'Embeddings/',sep="")

# Load property table and embeddings
setwd(ebdDir)
load(file = 'ebds_and_Graphproperties')

ebds = ebds[c('Groups2Vec','Metrics','Motifs2Vec','Graph2Vec_dp2','Graph2Vec_lab_dp2','ShortPaths2Vec','ShortPaths2Vec_lab')]
#ebds = ebds[c('Graph2Vec_dp2','Graph2Vec_lab_dp2','Graph2Vec_dp1','Graph2Vec_lab_dp1','Graph2Vec_dp4','Graph2Vec_lab_dp4')]
# Append Graph2Vec recto versions
#initSize = length(ebds)
#for(i in 1:initSize){ebds[[i+initSize]] = ebds[[i]][,1:30]}
#names(ebds)[(initSize+1):(2*initSize)] = paste(names(ebds)[1:initSize],'_recto')

# Order of rows in embeddings equal to properties table
orderedNetNames = rownames(ebds[[1]])

# Old names -> new names 
if(T){
  methNames = data.frame(
    meth=c('Groups2Vec',
           'Metrics',
           'Motifs2Vec',
           'Graph2Vec_dp2',
           'Graph2Vec_lab_dp2',
           'ShortPaths2Vec',
           'ShortPaths2Vec_lab'),
    print=c('Groups2Vec',
            'Metrics2Vec',
            'Motifs2Vec',
            'Graph2Vec',
            'Graph2Vec_lab',
            'ShortPaths2Vec',
            'ShortPaths2Vec_lab'),stringsAsFactors = F)
  propNames = data.frame(
    prop=c('maxTropLen',
           'TropLens',
           'nchains',
           'omni',
           'generalism',
           'loop'),
    print=c('maxTrophLen',
            'TrophLens',
            'nModules',
            'omni',
            'generalism',
            'loop'),stringsAsFactors = F)
}

#####
# Plot example networks
#####

require(ggnet)
selec = paste("site_",c(1696,284,1956,1194),".graphml",sep="")

ToSample = exp[exp$fileName%in%selec,]

palette = c( colorRampPalette(c("palegreen", "green4"))(5),
             colorRampPalette(c("gold", "darkorange3"))(5))
tmp = expand.grid(1:5,c('c1_','c2_'))
names(palette) = paste(tmp[,2],tmp[,1],sep="")

plot.list.nets = function(filesDf){
  return(lapply(filesDf$fileName,function(file){
  setwd(graphsDir)
  g = read.graph(file,format="graphml")
  gnet = intergraph::asNetwork(g)
  
  # Define node position based on SBM group 
  SBMgroups = network::get.vertex.attribute(gnet,attrname="feature") 
  gnet %v% "x" = as.numeric(sapply(SBMgroups,function(el) strsplit( strsplit(el,'_')[[1]][1]  ,'c')[[1]][2] ))
  gnet %v% "y" = as.numeric(sapply(SBMgroups,function(el) as.numeric(strsplit(el,'_')[[1]][2])))
  nodeX = network::get.vertex.attribute(gnet,attrname="x")
  nodeY = network::get.vertex.attribute(gnet,attrname="y")
  
  newNodeX = nodeX
  newNodeY = nodeY
  for(posX in unique(nodeX)){
    for(posY in unique(nodeY)){
      cd = nodeX==posX & nodeY==posY
      if(sum(cd)>0){
        newNodeX[cd] = posX + (1:sum(cd))*.07
        newNodeY[cd] = posY + (round(1:sum(cd)/2)==1:sum(cd)/2)*.2
      }
    }
  }
  gnet %v% "x" = newNodeX
  gnet %v% "y" = newNodeY
  
  cd = ToSample$fileName==file
  
  trophLens = eval(parse(text=paste('c(',ToSample$TropLens[cd],')')))
  trophLens = paste(trophLens - 1,collapse=",")
  
  p = ggnet2(gnet, 
             mode=c('x','y'),
             palette=palette,
             node.color="feature",
             arrow.gap = .05,
             arrow.size = 7,
             edge.size = 1,
             node.size = 10)
  #p = p + geom_text(data=textToAdd,aes(x=textToAdd$x,y=textToAdd$y,label=as.character(textToAdd$lab)),size=40)
  p = p + ggtitle(paste('nModules=',ToSample$nchains[cd],'; loop=',ToSample$loop[cd],
                        '\n trophLens=',trophLens,
                        '; omni=',ToSample$omni[cd],
                        '; generalism=',ToSample$generalism[cd],sep=""))
  p = p +theme(title = element_text(size=35),legend.position = "None")
  print(p)
  return(p)}))
}

plotList = plot.list.nets(filesDf=ToSample)

setwd(ebdDir)
png('net_examples.png',height=2400,width=2200)
multiplot(plots=plotList,cols = 2)
dev.off()

#####
# Compute predictive accuracy
# AND variables importance
#####

caracToPredict = c("maxTropLen","TropLens","nchains","omni","generalism","loop","n")

results = expand.grid(method=names(ebds),
                      carac=caracToPredict,oob_accuracy=NA)
results$method = as.character(results$method)
results$carac = as.character(results$carac)

for(i in 1:length(ebds)){
  method = names(ebds)[i]
  print(i)
  print(method)
  mat = ebds[method][[1]]
  
  if(method%in%c('ShortPaths2Vec_lab')){
    ebdSvd = svd(mat)$u[,1:64]
    Xpre = as.matrix(ebdSvd)
  }else{
    Xpre = mat
  }
  
  X= Xpre
  methodName = as.character(method) 
  
  for(col in caracToPredict){
    print(caracToPredict[caracToPredict==col])
    Y = factor(exp[,col])
    
    # Select mtry
    nbvars <- seq(1, min(ncol(X) - 1, 13) ,  2) 
    oobsMtry <- sapply(nbvars, function(nbv) {
      RF <- randomForest(Y ~ ., data = X, ntree = 100, mtry = nbv)
      return(RF$err.rate[RF$ntree, "OOB"])})
    
    # RF with optimal mtry
    forest = randomForest(Y ~ ., 
                          data = X,
                          ntree=300,
                          mtry=nbvars[which.min(oobsMtry)],
                          importance=T)
    
    # Store RF OOB accuracy
    results$oob_accuracy[results$method==methodName & results$carac==col] = 1- as.numeric(forest$err.rate[forest$ntree, "OOB"])
    
  }
  setwd(ebdDir)
  write.table(results,"Segregation_And_VarImp_EbdsAndRandom.csv",sep=";",row.names=F,col.names=T)
}

#####
# Table of predictive accuracy (Table 3)
#####

setwd(ebdDir)
results = read.csv("Segregation_And_VarImp_EbdsAndRandom.csv",sep=";",header=T)

#results_tmp = results[regexpr('RF_rdVars',results$method)<0,]
#tmp = data.frame(method=unique(results_tmp$method))

tmp = data.frame(method=unique(results$method))
tmp$method = as.character(tmp$method)
for(car in unique(results$carac)){
  eval(parse(text=paste('tmp$',car,'=results$oob_accuracy[results$cara=="',car,'"]',sep="")))
}
setwd(ebdDir)
write.table(tmp,'segregationTable.csv',sep=";",row.names=F,col.names=T)

rownames(tmp)=tmp$method
tmp = tmp[,colnames(tmp)!="method"]
print(100*tmp,digits=2)

#####
# Table of pseudo R2 on full embeddings (R2-ebd)
#####

ldists <- lapply(ebds, dist)
props=c("maxTropLen","TropLens","nchains","omni","generalism","loop","n")
res_R2 <- matrix (0, nrow = length(ebds), ncol = length(props))
rownames(res_R2) <- names(ebds)
colnames(res_R2) <- colnames(exp[,props])
for(i in 1:nrow(res_R2)){
  for(j in 1:ncol(res_R2)){
    ado = vegan::adonis(ldists[[i]] ~ as.factor(exp[,props[j]]), permutations = 0)
    res_R2[i, j] <- ado$aov.tab$R2[1]
  }
}

print(res_R2[,colnames(res_R2)!='name'],digits=2)

res_R2= as.data.frame(res_R2)
res_R2$name = rownames(res_R2)

setwd(ebdDir)
write.table(res_R2,'R2anderson_Ebds.csv',sep=";",row.names=F,col.names=T)

#####
# Compute UMAP plans 
#####

umaps = list()
for(i in 1:length(ebds)){
  method = names(ebds)[i]
  print(i)
  print(method)
  mat = ebds[method][[1]]
  
  # reduce dimension to 2 with Umap using default settings
  ebdUmap = as.data.frame(umap(mat,n_neighbors=min(150,dim(mat)[1]))$layout)
  ebdUmap = cbind(ebdUmap,rownames(ebds[[i]]))
  colnames(ebdUmap) = c('umap1','umap2','name')
  umaps[[i]] = ebdUmap
  setwd(ebdDir)
  write.table(ebdUmap,paste("umap_",method,".csv",sep=""),sep=";",row.names=F,col.names=T)
  gc(reset=T)
}
names(umaps) = names(ebds)

umaps = lapply(umaps,function(el)el[,colnames(el)!="name"])

setwd(ebdDir)
save(umaps,file="umaps")

#####
# Table pseudo R2 on UMAP plan : R2-umap
#####

setwd(ebdDir)
load(file="umaps")

ldists <- lapply(umaps, dist)
props=c("maxTropLen","TropLens","nchains","omni","generalism","loop","n")
res_R2 <- matrix (0, nrow = length(umaps), ncol = length(props))
rownames(res_R2) <- names(umaps)
colnames(res_R2) <- props
for(i in 1:nrow(res_R2)){
  for(j in 1:ncol(res_R2)){
    ado = vegan::adonis(ldists[[i]] ~ as.factor(exp[,props[j]]), permutations = 0)
    res_R2[i, j] <- ado$aov.tab$R2[1]
  }
}

res_R2= as.data.frame(res_R2)
res_R2$name = rownames(res_R2)

setwd(ebdDir)
write.table(res_R2,'R2anderson_Umaps.csv',sep=";",row.names=F,col.names=T)

#####
# Latex table R2-ebd + R2-umap (Table 4)
#####

setwd(ebdDir)
r2_umap = read.csv('R2anderson_Umaps.csv',sep=";",header=T)
r2_ebd = read.csv('R2anderson_Ebds.csv',sep=";",header=T)
rownames(r2_umap) = r2_umap$name
rownames(r2_ebd) = r2_ebd$name

methods = c("Groups2Vec","Metrics","Motifs2Vec","Graph2Vec_dp2","Graph2Vec_lab_dp2",
            "ShortPaths2Vec","ShortPaths2Vec_lab")

tmp = expand.grid(spa=c('embedding','umap'),prop=c('maxTropLen','TropLens','nchains','omni','generalism','loop','n'))
newCols = paste(tmp$spa,'_',tmp$prop,sep="") 

r2_all = matrix(NA,nrow = length(methods),ncol=length(newCols))
for(i in 1:dim(r2_all)[1]){
  for(j in 1:dim(r2_all)[2]){
    if(tmp$spa[j]=="embedding"){
      r2_all[i,j] = r2_ebd[methods[i],tmp$prop[j]]
    }else{
      r2_all[i,j] = r2_umap[methods[i],tmp$prop[j]]
    }
  }
}

r2_all = round(100*r2_all)/100
rownames(r2_all) = methods
colnames(r2_all) = newCols

print(xtable(r2_all))

#####
# Figure UMAPs (Figure 2 and 3) 
#####

setwd(ebdDir)
load(file="umaps")
oob = read.csv('segregationTable.csv',sep=";",header=T,stringsAsFactors = F)
r2 = read.csv('R2anderson_Umaps.csv',sep=";",header=T,stringsAsFactors = F)

exp$maxTropLen = exp$maxTropLen - 1

colfunc =colorRampPalette(c("blue", "red"))

samp = sample(1:dim(ebds[[1]])[1],min(600,dim(ebds[[1]])[1]))

# function to get list of plots
GetPlotList = function(){return(lapply(1:dim(selec)[1],function(i){
  meth = as.character(selec$meth[i])
  meName = as.character(methNames$print[methNames$meth==meth])
  prop = as.character(selec$prop[i])
  prName = as.character(propNames$print[propNames$prop==prop])
  toPlot = umaps[meth][[1]]
  toPlot[,prop] = factor(exp[,prop])
  toPlot$n = factor(exp$n)
  toPlot_tmp=toPlot[samp,]
  
  dotSize = 3
  scoreSize = 7
  textSize = 24
  
  textToAdd= data.frame(lab=c(paste('oob accur.:',round(100*oob[oob$method==meth,prop])),
                              paste('R2 umap:',round(100*r2[r2$name==meth,prop])/100)),
                        x=rep( min(toPlot_tmp$umap1)+0.12*(max(toPlot_tmp$umap1)-min(toPlot_tmp$umap1)) ,2),
                        y=c( max(toPlot_tmp$umap2)-0.00*(max(toPlot_tmp$umap2)-min(toPlot_tmp$umap2)),max(toPlot_tmp$umap2)-0.03*(max(toPlot_tmp$umap2)-min(toPlot_tmp$umap2)))) 
  
  preP = ggplot()+
    geom_point(data=toPlot_tmp,
               aes(x=toPlot_tmp$umap1,
                   y=toPlot_tmp$umap2,
                   color=toPlot_tmp[,prop],
                   shape=toPlot_tmp$n),size=dotSize)+labs(color=prop)+
    geom_text(data=textToAdd,aes(x=textToAdd$x,y=textToAdd$y,label=textToAdd$lab),size=scoreSize)+
    theme_bw()+xlab('umap1')+ylab('umap2')+
    ggtitle(paste(selec$let[i],'. Method:',meName,'; color:',prName))+
    theme(plot.title = element_text(size=textSize),
          axis.title.x=element_text(size=textSize),
          axis.title.y=element_text(size=textSize),
          axis.text.x = element_text(size=textSize),
          axis.text.y = element_text(size=textSize),
          legend.title = element_text(size=textSize),
          legend.text = element_text(size=textSize))+
    scale_shape_manual(values=c(4,16))+labs(color=prName,shape="size")
  if(prop=="maxTropLen"){
    return(preP + scale_color_manual(values=colfunc(4)))
  }else{
    return(preP + scale_color_manual(values=colfunc(2)))
  }
}))}

# Figure_Umaps_1 (Figure 2)
selec = data.frame(meth=c('Motifs2Vec',
                          'Motifs2Vec',
                          'Graph2Vec_dp2',
                          'Graph2Vec_dp2',
                          'Metrics',
                          'Metrics'),
                   prop=c('omni',
                          'nchains',
                          'omni',
                          'loop',
                          'loop',
                          'nchains'))
selec$let = LETTERS[1:dim(selec)[1]]
selec = selec[c(1,3,5,2,4,6),]

plotlist = GetPlotList()
setwd(ebdDir)
png('Figure_Umaps_1.png',height=2000,width=1385)
multiplot(plots=plotlist,cols=2)
dev.off()

# Figure_Umaps_2
selec = data.frame(meth=c('Graph2Vec_lab_dp2',
                          'Graph2Vec_lab_dp2',
                          'ShortPaths2Vec',
                          'ShortPaths2Vec',
                          'ShortPaths2Vec_lab',
                          'ShortPaths2Vec_lab'),
                   prop=c('omni',
                          'loop',
                          'loop',
                          'maxTropLen',
                          'nchains',
                          'generalism'))
selec$let = LETTERS[1:dim(selec)[1]]
selec = selec[c(1,3,5,2,4,6),]

plotlist = GetPlotList()
setwd(ebdDir)
png('Figure_Umaps_2.png',height=2000,width=1385)
multiplot(plots=plotlist,cols=2)
dev.off()

#####
# Plot UMAP multiplots for Appendix
#####

setwd(ebdDir)
load(file="umaps")
oob = read.csv('segregationTable.csv',sep=";",header=T,stringsAsFactors = F)
r2 = read.csv('R2anderson_Umaps.csv',sep=";",header=T,stringsAsFactors = F)

colfunc =colorRampPalette(c("blue", "red"))

samp = sample(1:dim(ebds[[1]])[1],min(600,dim(ebds[[1]])[1]))

props =  c('maxTropLen','nchains','omni','generalism','loop')
methodSelec = c("Groups2Vec","Metrics",           
                "Motifs2Vec","Graph2Vec_dp2",     
                "Graph2Vec_lab_dp2","ShortPaths2Vec",
                "ShortPaths2Vec_lab")

for(prop in props){
  print(prop)
  plotlist = lapply(1:length(methodSelec),function(i){
    meth = methodSelec[i]
    toPlot = umaps[meth][[1]]
    toPlot[,prop] = factor(exp[,prop])
    toPlot$n = factor(exp$n)
    toPlot_tmp=toPlot[samp,]
    
    textToAdd= data.frame(lab=c(paste('oob accur.:',round(100*oob[oob$method==meth,prop])),
                                paste('R2 umap:',round(100*r2[r2$name==meth,prop])/100)),
                          x=rep( min(toPlot_tmp$umap1)+0.12*(max(toPlot_tmp$umap1)-min(toPlot_tmp$umap1)) ,2),
                          y=c( max(toPlot_tmp$umap2)-0.00*(max(toPlot_tmp$umap2)-min(toPlot_tmp$umap2)),max(toPlot_tmp$umap2)-0.03*(max(toPlot_tmp$umap2)-min(toPlot_tmp$umap2)))) 
    if(prop=="maxTropLen"){
      return( 
        ggplot()+
          geom_point(data=toPlot_tmp,
                     aes(x=toPlot_tmp$umap1,
                         y=toPlot_tmp$umap2,
                         color=toPlot_tmp[,prop],
                         shape=toPlot_tmp$n),size=2)+labs(color=prop)+
          geom_text(data=textToAdd,aes(x=textToAdd$x,y=textToAdd$y,label=textToAdd$lab))+
          theme_bw()+xlab('umap1')+ylab('umap2')+
          ggtitle(paste('Method:',methNames$print[methNames$meth==meth],'; color:',propNames$print[propNames$prop==prop]))+
          scale_color_manual(values=colfunc(4))+scale_shape_manual(values=c(4,16))
      )
    }else{
      return(ggplot()+
               geom_point(data=toPlot_tmp,
                          aes(x=toPlot_tmp$umap1,
                              y=toPlot_tmp$umap2,
                              color=toPlot_tmp[,prop],
                              shape=toPlot_tmp$n),size=2)+labs(color=prop)+
               geom_text(data=textToAdd,aes(x=textToAdd$x,y=textToAdd$y,label=textToAdd$lab))+
               theme_bw()+xlab('umap1')+ylab('umap2')+
               ggtitle(paste('Method:',methNames$print[methNames$meth==meth],'; color:',propNames$print[propNames$prop==prop]))+
               scale_color_manual(values=colfunc(2))+scale_shape_manual(values=c(4,16)))
      
    }
  })
  
  setwd(ebdDir)
  png(paste(which(props==prop),'_Figure_Umaps_',prop,'.png',sep=""),height=1500,width=900)
  multiplot(plots=plotlist,cols=2)
  dev.off()
}

#####
# Figure Appendice: Density vs n_b, generalism, nModules
#####

gList = load.graphs.list(graphsDir)

setwd(ebdDir)
load("ebds_and_Graphproperties")
splitted = strsplit(exp$fileName,'.g')
orderedNetNames = sapply( 1:dim(exp)[1] ,function(el) splitted[[el]][1] ) 

coco = exp
coco$density = sapply( 1:dim(coco)[1] , function(i){ length(E(gList[orderedNetNames[i]][[1]]))/length(V(gList[orderedNetNames[i]][[1]]))^2 } )
coco$nb = NA
coco$nb[coco$nchains==1] = as.numeric( coco$TropLens[coco$nchains==1])
coco$nb[coco$nchains==2] = as.numeric(substr(coco$TropLens[coco$nchains==2],1,1))+as.numeric(substr(coco$TropLens[coco$nchains==2],3,3))
coco$generalism = factor(coco$generalism)
coco$nchains = factor(coco$nchains)
coco$n = factor(coco$n)
coco$n_b = factor(coco$nb) 

# With size
#coco$combi = factor(paste(coco$generalism,coco$nchains,coco$n))
#p = ggplot()+geom_boxplot(data=coco,aes(x=nb,y=density,fill=combi))
#p = p +scale_y_continuous(limits=c(0,.3))+theme_bw()
#print(p)  

# With loop
coco$combi3 = factor(paste(coco$generalism,coco$loop))
p = ggplot()+geom_boxplot(data=coco,
  aes(x=n_b,y=density,fill=combi3))+labs(fill="generalism - loop")
p = p +scale_y_continuous(limits=c(0,.3))+theme_bw()+xlab('number of trophic groups')
print(p)  

# FIGURE
coco$combi2 = factor(paste(coco$generalism,coco$nchains))
p = ggplot()+geom_boxplot(data=coco,
        aes(x=n_b,y=density,fill=combi2))+xlab('Number of trophic groups')+ylab('directed density')
p = p +scale_y_continuous(limits=c(0,.3))+theme_bw()+labs(fill="generalism - nModules")
p = p + theme(axis.text.x = element_text(size=25),
              axis.text.y = element_text(size=25),
              axis.title.x = element_text(size=35),
              axis.title.y = element_text(size=35),
              legend.title = element_text(size=35),
              legend.text = element_text(size=25))
setwd(ebdDir)
png('Connectance.png',height=1000,width=1500)
print(p)
dev.off()

#####
# Table of R2-loss in embedding
#####

Ns = c(60,120)
resList = list()
for(N in Ns){
  SAMPLE = exp$n==N
  ldists <- lapply(1:length(ebds),function(i) dist(ebds[[i]][SAMPLE,])  )
  props=c("maxTropLen","TropLens","nchains","omni","generalism","loop")
  res_R2 <- matrix (0, nrow = length(ebds), ncol = length(props))
  rownames(res_R2) <- names(ebds)
  colnames(res_R2) <- colnames(exp[,props])
  for(i in 1:nrow(res_R2)){
    for(j in 1:ncol(res_R2)){
      ado = vegan::adonis(ldists[[i]] ~ as.factor(exp[SAMPLE,props[j]]), permutations = 0)
      res_R2[i, j] <- ado$aov.tab$R2[1]
    }
  }
  resList[[which(Ns==N)]] = res_R2
}

results = (resList[[1]] + resList[[2]])/2
results= as.data.frame(results)
results$name = rownames(results)

setwd(ebdDir)
write.table(results,'MeanR2anderson_Ebds.csv',sep=";",row.names=F,col.names=T)


# Print R2-loss table
setwd(ebdDir)
R2mean = read.csv('MeanR2anderson_Ebds.csv',sep=";",header=T)
R2mixed = read.csv('R2anderson_Ebds.csv',sep=";",header=T) # Output from FullAnalysisPipeline.R

meths = as.character(R2mixed$name)

R2mean= R2mean[,colnames(R2mean)!="name"]
R2mixed = R2mixed[,!colnames(R2mixed)%in%c('name','n')]

print('R2_loss full embedding')
res = (R2mean - R2mixed)/R2mean
rownames(res)=meths
print( round(100*res),digits=2)

#####
# Table of R2-loss on Umap
#####

Ns = c(60,120)
resList = list()
for(N in Ns){
  SAMPLE = exp$n==N
  ldists <- lapply(1:length(umaps),function(i) dist(umaps[[i]][SAMPLE,])  )
  props=c("maxTropLen","TropLens","nchains","omni","generalism","loop")
  res_R2 <- matrix (0, nrow = length(umaps), ncol = length(props))
  rownames(res_R2) <- names(umaps)
  colnames(res_R2) <- colnames(exp[,props])
  for(i in 1:nrow(res_R2)){
    for(j in 1:ncol(res_R2)){
      ado = vegan::adonis(ldists[[i]] ~ as.factor(exp[SAMPLE,props[j]]), permutations = 0)
      res_R2[i, j] <- ado$aov.tab$R2[1]
    }
  }
  resList[[which(Ns==N)]] = res_R2
}

results = (resList[[1]] + resList[[2]])/2
results= as.data.frame(results)
results$name = rownames(results)

setwd(ebdDir)
write.table(results,'MeanR2anderson_Umaps.csv',sep=";",row.names=F,col.names=T)

# Print R2-loss in Umap plan
setwd(ebdDir)
R2mean = read.csv('MeanR2anderson_slicedEbds.csv',sep=";",header=T)
R2mixed = read.csv('R2anderson_slicedEbds.csv',sep=";",header=T)
meths = as.character(R2mixed$name)
R2mean= R2mean[,colnames(R2mean)!="name"]
R2mixed = R2mixed[,!colnames(R2mixed)%in%c('name','n')]
res = (R2mean - R2mixed)/R2mean
rownames(res)=meths
print( round(100*res),digits=2)


#####
# Euclidian distance Scatter plot inter-methods (investigate complementarity)
#####

nSamp = 1000
nSampledDists = 10000

samp = sample(1:dim(ebds[[1]]),nSamp)
NetNames = orderedNetNames[samp]
distances = expand.grid(j=as.integer(1:nSamp),i=as.integer(1:nSamp))
distances = distances[distances$j<distances$i,]
for(i in 1:length(ebds)){
  method = names(ebds)[i]
  print(i)
  print(method)
  mat = ebds[method][[1]][samp,]
  
  # For embeddings, Distance = Euclidean distance between points in the embedding
  Ds = dist(mat,method = "euclidian",diag=F,upper = F)  
  vec = as.numeric(Ds)
  vec = vec/max(vec)
  print('dists computed')
  
  eval(parse(text= paste('distances$',method,'=vec',sep="" ) ))  
}

# We add a row full of 0 to insure that point (0,0) is in all plots windows
toAdd = distances[1,,drop=F]
toAdd[1,] = 0
distances = rbind(toAdd,distances)

sampledDists = c(1,sample(1:dim(distances)[1],nSampledDists))

# Order of methods to plot
refNames = c('Groups2Vec','ShortPaths2Vec_lab','Motifs2Vec','Graph2Vec_lab_dp2','Metrics','ShortPaths2Vec','Graph2Vec_dp2')
toPlotNames = c('Groups2Vec','ShortPaths2Vec_lab','Motifs2Vec','Graph2Vec_lab','Metrics2Vec','ShortPaths2Vec','Graph2Vec')

formu = as.formula( paste('~ ', paste(refNames,collapse="+"))) 
print(formu)

# Plot and save multi scatter plot
setwd(ebdDir)
png('methods_distances_pairs_final.png',height=1700,width=1700)
pairs( formu ,data=distances[sampledDists,] , 
       diag.panel = panel.hist, 
       upper.panel = NULL , 
       labels= toPlotNames,
       pch=3 ,
       col=  alpha(rep('black',dim(distances)[1]),0.05) , 
       cex.labels = 3.5 , 
       cex.axis = 3)
dev.off()

#####
# Distance correlations across embeddings
#####
require(corrplot)

mcor = matrix(NA,length(refNames),length(refNames))
for(i in 1:length(refNames)){
  for(j in 1:length(refNames)){
    mcor[i,j] = energy::dcor( ebds[refNames[i]][[1]] , ebds[refNames[j]][[1]])
  }
}
rownames(mcor) = toPlotNames
colnames(mcor) = toPlotNames

setwd(ebdDir)
saveRDS(mcor,'mat_dCor')

#print(mcor,digits=2)
setwd(ebdDir)
mcor = readRDS('mat_dCor')

setwd(ebdDir)
png('methods_distances_correlations_final.png',height=2000,width=2000)
print(corrplot(mcor, type="upper", tl.col="white", tl.srt=45, diag=F,cl.cex=3,cl.lim=c(0.25,1),addCoef.col = "grey80",number.cex=4))
dev.off()

# Old agreement metric
if(F){
  library(corrplot)
  
  mcor = matrix(NA,length(toPlot),length(toPlot))
  for(i in 1:length(toPlot)){
    for(j in 1:length(toPlot)){
      mcor[i,j] = cor(distances[,as.character(names(ebds)[toPlot[i]])],distances[,as.character(names(ebds)[toPlot[j]])])
    }
  }
  rownames(mcor) = names(ebds)[toPlot]
  colnames(mcor) = names(ebds)[toPlot]
  
  print(mcor,digits=2)
  
  setwd(ebdDir)
  png('methods_distances_pairs_correlations_OLD.png',height=2000,width=2000)
  print(corrplot(mcor[reord,reord], type="upper", tl.col="black", tl.srt=45, diag=F,cl.cex=3))
  dev.off()
  
}


#####
# correlations across metrics
#####

met = ebds['Metrics'][[1]]
str(met)

ordered = c("density","dirCon",
            "predPerPrey","vulnerabilitySD","generalitySD",
            "modul",
            "omniLev","omniProp","mean_tropLevel","tropLen",
            "diameter",
            "interProp","transitivity",
            "mean_distance" ,
            "basalProp","topProp","assortativity")

formu = as.formula(paste("~",paste(ordered,collapse="+")))

samp = sample(1:dim(met)[1],800)

setwd(ebdDir)
png("metrics_pairs.png",height=2000,width=2000)
pairs( formu ,data=met[samp,] , 
       diag.panel = panel.hist, 
       upper.panel = NULL , 
       labels= ordered,
       pch=3 ,
       col=  alpha(rep('black',length(samp)),0.05) , 
       cex.labels = 1.6 , 
       cex.axis = 2)
dev.off()

# Correlations across metrics
mcor = matrix(NA,length(ordered),length(ordered))
for(i in 1:length(ordered)){
  for(j in 1:length(ordered)){
    mcor[i,j] = cor(met[,ordered[i]] , met[,ordered[j]] )
  }
}
rownames(mcor) = ordered
colnames(mcor) = ordered

require(corrplot)

setwd(ebdDir)
png('metrics_correlations.png',height=1000,width=1200)
print(corrplot(mcor, type="upper", tl.col="black", tl.srt=45, diag=F,cl.cex=1.5,cl.lim=c(-1,1),addCoef.col = "grey80",number.cex=1))
dev.off()

