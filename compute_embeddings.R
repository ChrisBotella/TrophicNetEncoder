
repoDir = "C:/Users/user/pCloud local/boulot/Github/EcoGraph Encoder/Simulate networks/trophicSBM9/"
#repoDir = "/home/christophe/pCloud local/boulot/Github/EcoGraph Encoder/Simulate networks/trophicSBM9/"
source(paste(repoDir,'functions_to_source.R',sep=""))

masterDir="C:/Users/user/pCloud local/boulot/data/Simu_networks/trophicSBM9/"
#masterDir='/home/christophe/pCloud local/boulot/data/Simu_Networks/trophicSBM9/'

graphsDir = paste(masterDir,'Graphs/',sep="")

# Order of rows in embeddings equal to properties table
exp = read.csv(paste(masterDir,'graphsParameters.csv',sep=""),sep=";",header=T,stringsAsFactors = F)
splitted = strsplit(exp$fileName,'.g')
orderedNetNames = sapply( 1:dim(exp)[1] ,function(el) splitted[[el]][1] ) 

# Assumes .graphml files are already created in .../Graphs/
sched = data.frame(ebd_method=c('Groups2Vec',
                                 'Metrics',
                                 "Motifs2Vec",
                                "Graph2Vec_dp1",
                                "Graph2Vec_lab_dp1",
                                "Graph2Vec_dp2",
                                "Graph2Vec_lab_dp2",
                                "Graph2Vec_dp4",
                                "Graph2Vec_lab_dp4",
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
                              "Graph2Vec",
                              "Graph2Vec",
                                 "Graph2Vec",
                                 "Graph2Vec",
                                 "ShortPaths2Vec",
                                 "ShortPaths2Vec"),
                    params=c('list(groupsAttributeName="SBMgroupId",groupsIds=1:10)',
                                "list(metrics=c('density','dirCon','modul',
                                        'omniLev','omniProp','predPerPrey',
                                        'mean_tropLevel','basalProp','interProp','topProp',
                                        'vulnerabilitySD','generalitySD','transitivity','diameter',
                                        'mean_distance','assortativity','tropLen'))",
                                "list()",
                             "list(graphsDir=graphsDir,depth=1,dim=30,lab=F)",
                             "list(graphsDir=graphsDir,depth=1,dim=30,lab=T)",
                             "list(graphsDir=graphsDir,depth=2,dim=30,lab=F)",
                             "list(graphsDir=graphsDir,depth=2,dim=30,lab=T)",
                             "list(graphsDir=graphsDir,depth=4,dim=30,lab=F)",
                             "list(graphsDir=graphsDir,depth=4,dim=30,lab=T)",
                                "list(graphsDir=graphsDir,depth=6,dim=30,lab=F)",
                                "list(graphsDir=graphsDir,depth=6,dim=30,lab=T)",
                                "list(nodeLabelAttribName=NULL)",
                                "list(nodeLabelAttribName='SBMgroupId')"))
print(as.character(sched$ebd_method))
gList = load.graphs.list(graphsDir)

toRun = c(4:9,12,13)
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
save(ebds,exp,file="ebds_and_Graphproperties")
