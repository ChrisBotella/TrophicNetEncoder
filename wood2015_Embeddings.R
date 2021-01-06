repoDir = getwd() #OR REPLACE BY DIRECTORY CONTAINING SCRIPTS
source(paste(repoDir,'functions_to_source.R',sep=""))

masterDir = getwd() #OR REPLACE BY DIRECTORY CONTAINING SCRIPTS

graphsDir = paste(masterDir,"Graphs/",sep="")
ebdDir = paste(masterDir,"Embeddings/",sep="")

# Web scales
setwd(ebdDir)
load("woodData")
orderedNetNames = as.character(webScales$name)

# Assumes .graphml files are already created in .../Graphs/
sched = data.frame(ebd_method=c('Groups2Vec',
                                'Metrics',
                                "Motifs2Vec",
                                "Graph2Vec_dp2",
                                "Graph2Vec_lab_dp2",
                                "Graph2Vec",
                                "Graph2Vec_lab",
                                "ShortPaths2Vec",
                                "ShortPaths2Vec_lab"),
                   ebdType=c('Groups2Vec',
                             'Metrics',
                             "Motifs2Vec",
                             "Graph2Vec",
                             "Graph2Vec",
                             "Graph2Vec",
                             "Graph2Vec",
                             "ShortPaths2Vec",
                             "ShortPaths2Vec"),
                   params=c('list(groupsAttributeName="SBMgroupId",groupsIds=1:6)',
                            "list(metrics=c('density','dirCon','modul',
                                        'omniLev','omniProp','predPerPrey',
                                        'mean_tropLevel','basalProp','interProp','topProp',
                                        'vulnerabilitySD','generalitySD','transitivity','diameter',
                                        'mean_distance','assortativity','tropLen'))",
                            "list()",
                            "list(graphsDir=graphsDir,depth=2,dim=30,lab=F)",
                            "list(graphsDir=graphsDir,depth=2,dim=30,lab=T)",
                            "list(graphsDir=graphsDir,depth=6,dim=30,lab=F)",
                            "list(graphsDir=graphsDir,depth=6,dim=30,lab=T)",
                            "list(nodeLabelAttribName=NULL)",
                            "list(nodeLabelAttribName='SBMgroupId')"))

toRun = 1:dim(sched)[1]
for(i in toRun){
  print(paste('operation',i))
  method = as.character(sched$ebd_method[i])
  ebdType = as.character(sched$ebdType[i])
  print(as.character(method))
  # Create save directory if doesnt exist yet
  setwd(masterDir)
  if(!'Embeddings'%in%list.files()){dir.create('Embeddings')}
  setwd(paste(masterDir,'Embeddings/',sep=""))
  if(!method%in%list.files()){dir.create(as.character(method))}
  
  # Extract params
  eval(parse(text=paste('params= ',sched$params[i]))) 
  # Compute embedding
  ebd= get.embedding(gList,ebdType = ebdType,params = params)
  # add networks names 
  if(ebdType!="Graph2Vec"){ebd$name = names(gList)}
  # Order networks consistently for all embeddings
  rownames(ebd) = ebd$name
  ebd = ebd[orderedNetNames,,drop=F]
  
  ### save embedding
  if(!is.null(ebd)){
    outputName = as.character(method)
    write.table(ebd,paste(masterDir,'Embeddings/',method,'/',outputName,".csv",sep=""),sep=";",row.names=F,col.names=T)
  }
}

# Make all embeddings along with properties table in R format
ebds = list()
for(i in 1:length(sched$ebd_method)){
  file = paste(masterDir,'Embeddings/',sched$ebd_method[i],'/',sched$ebd_method[i],".csv",sep="")
  print(as.character(sched$ebd_method[i]))
  ebd = read.csv(file,sep=";",header=T,stringsAsFactors = F)
  rownames(ebd)= ebd$name
  ebds[[i]] = as.matrix(ebd[orderedNetNames,colnames(ebd)!="name"])
}
names(ebds) = as.character(sched$ebd_method)

setwd(paste(masterDir,'Embeddings/',sep=""))
save(ebds,file="WoodEbds")

#####
# Umap plans
#####

toRun = 1:dim(methods)[1]
NetIds = sapply( strsplit(orderedNetNames,'te_'), function(el) el[[2]])  

toPlot = list()
for(i in 1:length(toRun)){
  print(names(ebds)[toRun[i]])
  
  mat = as.matrix(ebds[[toRun[i]]])
  
  # reduce dimension to 2 with Umap using default settings
  ebdUmap = as.data.frame(umap(mat,n_neighbors=350)$layout)
  
  tmp = cbind(ebdUmap,rownames(ebd))
  colnames(tmp)=c('umap1','umap2','name')
  toPlot[[i]] = tmp
  
  setwd(ebdDir)
  write.table(toPlot[[i]],paste('toPlot_umap_',names(ebds)[toRun[i]],'.csv',sep=""),sep=";",row.names=F,col.names=T)
}
names(toPlot) = names(ebds)[toRun]
setwd(ebdDir)
save(toPlot,file = 'umaps')

# MultiPlot
selec = names(toPlot)
plotlist = list()
for(i in 1:length(selec)){
  method = as.character(selec[i])
  tmp = toPlot[selec[i]][[1]]
  tmp$webScale = webScales$WebScale
  
  p = ggplot(tmp,aes(x=umap1,y=umap2,colour=webScale))+
    geom_point(size=1.5)+
    scale_colour_manual(values=rainbow(length(levels(webScales$WebScale))))+
    theme_bw()+
    theme(legend.position = "none")+
    ggtitle(as.character(method))
  
  plotlist[[i]] = p
  #png(paste(ebdDir,method,'_Umap_alongWith_webScale.png',sep=""),height=400,width=500)
  #print(p)
  #dev.off()
}

multiplot <- function(plots=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  numPlots = length(plots)
  print(numPlots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

setwd(ebdDir)
png('Figure_wood_raw.png',height=700,width=900)
multiplot(plots=plotlist,cols=3)
dev.off()
