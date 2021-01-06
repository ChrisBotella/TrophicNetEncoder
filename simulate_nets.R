repoDir = getwd() #OR REPLACE BY DIRECTORY CONTAINING SCRIPTS
source(paste(repoDir,'functions_to_source.R',sep=""))
masterDir = getwd() #OR REPLACE BY DIRECTORY CONTAINING SCRIPTS
if(masterDir!=getwd()){setwd(masterDir)}
create.dir('Graphs')
graphsDir = paste(masterDir,'Graphs/',sep="")

set.seed(32)

nPerSize = 2500

exp = data.frame(maxTropLen=NA,
                 TropLens = NA,
                 nchains = NA,
                 omni = NA ,
                 loop = NA ,
                 generalism = NA,
                 nbasals=NA,
                 ntop=NA,
                 dirCon=NA,
                 density=NA,
                 n=c(rep(60,nPerSize),rep(120,nPerSize)))

maxChains = 2
minChainLength = 2
maxChainLength = 5
nGraphs = dim(exp)[1]
pBase = .7

# Associate SBMgroup label with integer id
chains = paste('c',1:maxChains,sep="")
chainsLength = rep(maxChainLength,maxChains)
grNames = NULL
for(ch in 1:length(chains)){grNames = c(grNames, paste(chains[ch],'_',1:chainsLength[ch],sep=""))}
groupsIDs = data.frame(SBMgroup=grNames,SBMgroupId = 1:length(grNames))

for(j in 1:nGraphs){
  nNodes = exp$n[j]
  
  # Firstly we build the SBM groups
  nc = round(runif(1,.5001,maxChains+.4999)) # Draw number of chains
  chains = paste('c',1:nc,sep="")
  Ls = round(runif(nc,.5001+minChainLength-1,maxChainLength+.4999))
  grNames = unlist(lapply(1:length(chains),function(ch) paste(chains[ch],'_',1:Ls[ch],sep="")))
  
  # Draw nodes membership uniformly across groups
  SBMgrMembership = c(grNames , grNames[round(runif(nNodes-length(grNames),.5001,length(grNames)+.4999))] )
  
  # Define the architectural properties of the network
  exp$maxTropLen[j] = max(Ls)
  exp$TropLens[j] =paste(Ls,collapse = ",")
  exp$nchains[j]=nc
  exp$omni[j]=as.logical(round(runif(1,0,1)))
  exp$generalism[j] = as.logical(round(runif(1,0,1)))
  exp$loop[j]=as.logical(round(runif(1,0,1)))

  # set the block connections probabilities
  nb = sum(Ls)
  fLoop = nNodes*(nNodes-nb)*nb^2/(2*nNodes^2*nb*(nb-nc))
  fOmni = sum(sapply(1:nc,function(j)  choose(Ls[j]-1, 2) ))/(nb-nc)
  pOmni = exp$omni[j] * .2
  pLoop = exp$loop[j] * .15
  pGeneralism = exp$generalism[j]* .2
  pPred =  pBase + pGeneralism - pOmni*fOmni  - pLoop*fLoop
  pInter = .1
  
  # Build the SBM template
  adj = matrix(0,length(grNames),length(grNames))
  rownames(adj) = grNames
  colnames(adj) = grNames
  grs = strsplit(grNames,'_')
  chainMemb = sapply(1:length(grs),function(i)grs[[i]][1])
  tlMemb = as.numeric(sapply(1:length(grs),function(i)grs[[i]][2]))
  for(gr1 in 1:length(grNames)){
    for(gr2 in 1:length(grNames)){
      if(gr1==gr2){
        adj[gr1,gr2] = pLoop
      }else if(tlMemb[gr1]+1==tlMemb[gr2]){
        if(chainMemb[gr1]==chainMemb[gr2]){
          adj[gr1,gr2] = pPred
        }else{
          adj[gr1,gr2] = pInter
        }
      }else if(tlMemb[gr1]+1<tlMemb[gr2] & chainMemb[gr1]==chainMemb[gr2]){
        adj[gr1,gr2] = pOmni
      }
    }
  }
  
  g = simu.sbm(SBMgrMembership,adj=adj,selfLink=F)
  # SBMgroup as character attribute 
  g = set_vertex_attr(g, 'feature', index = V(g), value = vertex_attr(g, "SBMgroup") )
  # SBMgroup attribute as integer
  toJoin = data.frame(SBMgroup=vertex_attr(g, "SBMgroup"),row.id=1:length(vertex_attr(g, "SBMgroup")))
  tmp = merge(toJoin,groupsIDs,by="SBMgroup",sort=F) 
  g = set_vertex_attr(g, 'SBMgroupId', index = V(g), value =  tmp[order(tmp$row.id),'SBMgroupId'])
  # Change order of attributes
  tmp = vertex_attr(g)[[1]]
  vertex_attr(g)[[1]] = vertex_attr(g)[[3]]
  vertex_attr(g)[[3]] = tmp
  names(vertex_attr(g)) = c('SBMgroupId','feature','SBMgroup')
  
  # Get basic metrics
  Adj = as.matrix(as_adjacency_matrix(g))
  exp$nbasals[j] = sum( rowSums(as.matrix(Adj) ) == 0 )
  exp$ntop[j] = sum( colSums(as.matrix(Adj))==0)
  exp$dirCon[j] = sum(as.vector(as.matrix(Adj)))/nNodes^2
  exp$density[j] = sum( (as.vector(Adj[lower.tri(Adj)]) + as.vector(t(Adj)[lower.tri(t(Adj))])) >0 ) / (nNodes*(nNodes-1)/2)
  
  if(j/100==round(j/100)){
    flush.console()
    cat('    \r     Process...',100*j/nGraphs,'%       \r    ')}
  
  setwd(graphsDir)
  write_graph(g,file=paste('site_',j,".graphml",sep=""),format="graphml")

}

exp$fileName = paste('site_',1:nGraphs,'.graphml',sep="")
setwd(masterDir)
write.table(exp,'graphsParameters.csv',sep=";",row.names=T,col.names=T)

