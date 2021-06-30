# This script contains the functions to source 
# to run all the embedding, reduction and analysis pipeline
require(igraph)
require(network)
require(intergraph)

require(ggplot2)
require(sna)
require(GGally)

require(reticulate)
require(umap)
require(moments)
require(randomForest)

library(energy)
require(xtable)

pcName =as.character(Sys.info()["nodename"]) 
if(pcName=="DESKTOP-RUARS8N"){
  use_python("C:/Users/user/miniconda3/python.exe",required = T)
  source_python('C:/Users/user/pCloud local/boulot/Github/EcoGraph Encoder/embed_reduce_analyse_wrappers/py_get_ebd.py')
}

#####
# Plot functions
#####


panel.hist <- function(x, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE,breaks="fd")
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, ...)
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


#####
# Simulate graphs
#####

# Simulate a graph according to a SBM
simu.sbm = function(nodesGroups,adj=NULL,selfLink=F){
  if(is.null(rownames(adj)) | is.null(colnames(adj))){
    print('Error: rownames and colnames of adj must be the SBM groups names / identifiers')
    return(NULL)
  }else{
    rownames(adj) = unique(nodesGroups)
    colnames(adj) = unique(nodesGroups)
    
    ADJ = matrix(0,length(nodesGroups),length(nodesGroups))
    
    draw = runif(length(nodesGroups)^2,0,1)
    
    k = 1
    for(i in 1:length(nodesGroups)){
      for(j in 1:length(nodesGroups)){
        if(draw[k]<adj[nodesGroups[i],nodesGroups[j]]){
          ADJ[i,j] = 1
        }
        k=k+1
        if(k/1000==round(k/1000)){
          flush.console()
          cat('     \r       Process...',100*k/length(nodesGroups)^2,'%         \r  ')
        }
      }
    }
    if(!selfLink){
      diag(ADJ)=rep(0,length(nodesGroups))
    }
    G = graph_from_adjacency_matrix(ADJ,mode="directed")
    G = set_vertex_attr(G, "SBMgroup", index = V(G), nodesGroups)
    return(G)
  }
}

# Simulate a trophic network from a trophic SBM template
Graph.from.TrophicSBM.template = function( N,nTLs = 3, GrAdj = NULL,SpPerGroup = NULL,canib=F){
  if(is.null(GrAdj)){
    GrAdj = matrix(0,nTLs,nTLs)
    for(i in 1:nTLs){GrAdj[i,i:nTLs] = 1}
  }
  
  if(is.null(SpPerGroup)){
    nodesGroups = round(runif(N,.5,nTLs+.5))
  }else{
    nodesGroups = NULL
    for(tl in 1:nTLs){
      nodesGroups = c(nodesGroups,rep(tl,SpPerGroup[tl]) )
    }
  }
  G = simu.sbm(nodesGroups,adj=GrAdj,selfLink= canib)
  
  return(G)
}




#####
# graphs list loader
#####

load.graphs.list = function(graphsDir){
  setwd(graphsDir)
  files = list.files()
  gList = list()
  for(i in 1:length(files)){
    gList[[i]] = read.graph(files[i],format="graphml")
  }
  split = strsplit(files,'.graphm')
  names = sapply(1:length(split),function(i) split[[i]][1])
  names(gList)=names
  return(gList)
}

#####
# Foodweb metrics functions
#####

clusteringCoef = function(adjmat){
  tmp = data.frame(n=rep(0,2*dim(adjmat)[1]),ones=0)
  for(k in 1:dim(adjmat)[1]){
    ids = which(adjmat[k,]>0)
    if(length(ids)>1){
      pairs = expand.grid(p1=ids,p2=ids)
      pairs = pairs[pairs$p1<pairs$p2,,drop=F]
      vTmp = sapply(1:dim(pairs)[1],function(j){adjmat[pairs$p1[j],pairs$p2[j]]})
      tmp$ones[k] = sum(vTmp)
      tmp$n[k] = length(vTmp)
    }
  }
  for(l in 1:dim(adjmat)[1]){
    ids = which(adjmat[,l]>0)
    if(length(ids)>1){
      pairs = expand.grid(p1=ids,p2=ids)
      pairs = pairs[pairs$p1<pairs$p2,,drop=F]
      vTmp = sapply(1:dim(pairs)[1],function(j)adjmat[pairs$p1[j],pairs$p2[j]])
      tmp$ones[k+l] = sum(vTmp)
      tmp$n[k+l] = length(vTmp)
    }
  }
  return(sum(tmp$ones)/sum(tmp$n))
}

short.weighted.trophic.levels = function(g){
  # The short weighted trophic level (SWTL) 
  # is the average (directed) distance to the basal species (who have no prey)
  
  # /!\ By convention 
  # i predates j iif as_adjacency_matrix(g)[i,j]=1
  
  # NB: The SWTL of a species is computed only over 
  # basal species that can be reached from it 
  # following edges directions. 
  # It can be thus computed on non-connected
  # graphs
  
  SWtrophicLevels = NULL
  
  Adj = as.matrix(as_adj(g))
  # Remove Cannibalism to integrate basal species that are cannibals
  diag(Adj)=0
  
  # ids of the basal species (who have no prey)
  basals = which(rowSums(Adj)==0)
  # Compute shortest distances between species 
  D = igraph::distances(g,mode='out')
  
  if(length(basals)>0){
    for(i in 1:dim(Adj)[1]){
      # Short Weighted Trophic Level of species i
      # toDisp=which(is.finite(D[i,]))
      #D[i,toDisp]
      #Ord2 = which(is.finite(D[toDisp[toDisp!=i],]))
      DirectedDistToBasals = D[i,basals]
      DirectedDistToBasals = DirectedDistToBasals[is.finite(DirectedDistToBasals)]
      SWtrophicLevels[i] = mean(DirectedDistToBasals)
    }
    return(SWtrophicLevels)
  }else{
    print('Short Weighted Trophic Level can not be computed, there is no basal species')
    return(NULL)
  }
}

omnivory.levels = function(g){
  # Standard Deviation of the Trophic Level (Kefi) of the preys
  
  # By convention it is 0 for basal species
  # and species with only one prey.
  
  trophicLevels = trophiclevel(g)
  Adj =  as.matrix(as_adj(g))
  omnivoryLevels = NULL
  for(i in 1:dim(Adj)[1]){
    # Omnivory level of each species
    if(sum(Adj[i,]>0)>1){
      omnivoryLevels[i] = sd(trophicLevels[Adj[i,]>0])
    }else{
      omnivoryLevels[i] = 0 # For basals and those with only one prey
    }
  }
  return(omnivoryLevels)
}

mean.shortest.path.length = function(g){
  # Mean shortest directed path length between all oriented pairs of nodes 
  # discarding non-connected oriented pairs and self distance
  
  D = igraph::distances(g,mode="out")
  diag(D)=Inf
  tmp = as.vector(D)
  return(mean(tmp[is.finite(tmp)]))
}

basal.species.proportion=function(g){
  Adj = as.matrix(as_adj(g))
  basals = which(rowSums(Adj)==0)
  return(length(basals)/dim(Adj)[1])
}

intermediary.species.proportion=function(g){
  Adj = as.matrix(as_adj(g))
  inter = which(rowSums(Adj)>0 & colSums(Adj)>0)
  return(length(inter)/dim(Adj)[1])
}

top.species.proportion=function(g){
  Adj = as.matrix(as_adj(g))
  top = which(colSums(Adj)==0)
  return(length(top)/dim(Adj)[1])
}


## From Sonia Kefi's code
## C++ to R translation
## Warning: works only when basal species exist (i.e. some nodes with out-degree=0)
### NOUVEAU CODE TL corrige par Christophe
# Ce code met TL = 1 aux noeuds qui n'ont pas de proies, y compris les isol?s
trophiclevel <- function(g){
  adjmat <- as.matrix(as_adj(g))
  N <- nrow(adjmat)
  TLvec <- rep(1,N)
  for(rep in 1:10){
    TLtemp <- rep(0,N)
    for(i in 1:N){
      temp1 <- sum(adjmat[i,])  ## temp1 contains the number of prey of species i
      if(temp1>0){ ## If species i has at least one prey
        ## calculate the mean TL of the preys of species i
        TLtemp[i] = sum( adjmat[i,] * TLvec ) / temp1
      }
    }
    TLvec <- 1 + TLtemp ## TL = 1 + av TL of the prey
  }
  TLvec = TLvec - min(TLvec) +1
  return(TLvec)
}

### NOUVEAU CODE SPARSE MacKay (Marc)
trophicLevelMackay <- function(G) {
  A = as.matrix(get.adjacency(G))
  names_loc = rownames(A)
  u  = igraph::degree(G)
  v =  igraph::degree(G,mode ='in') -  igraph::degree(G,mode = 'out')  
  A = A[-1,-1]
  u = u[-1]
  v = v[-1]
  L = diag(u) - A - t(A)  
  L_sparse = as(L, "sparseMatrix")
  v_sparse = as(v, "sparseMatrix")  
  TL_vec = Matrix::solve(L,v)
  TL_vec = c(0,TL_vec)
  TL_vec = TL_vec - min(TL_vec)
  names(TL_vec) = names_loc
  return(TL_vec)
}

robustness = function(adj,n=100,extinctionOrder=c('random','rarest'),SpPresences=NULL){
  adj = as.matrix(adj)
  nSp = dim(adj)[1]
  
  # Get basal species 
  basals =  rownames(adj)[rowSums(adj)==0]
  
  Fragilite = matrix(0,nSp,n)
  Fragilite[1,] = 1
  for(j in 1:n){
    # we initialize the bag of species with all species of the network
    bagOfSpecies = rownames(adj)
    k=1
    extinction=NULL
    while(length(bagOfSpecies)>0 & k<nSp){
      # We draw the species for primary extinction
      if(extinctionOrder=="random"){
        extinction[k] = sample(bagOfSpecies,1)
      }else if(extinctionOrder=="rarest"){
        remainingPresences = SpPresences[bagOfSpecies]
        rarest = names(remainingPresences)[which(remainingPresences==min(remainingPresences))]
        extinction[k] = sample(  rarest , 1 )
      }
      # we remove it from the bag of species
      bagOfSpecies = setdiff(bagOfSpecies,extinction[k])
      if(length(bagOfSpecies)>0){
        
        adjTmp = adj[setdiff(bagOfSpecies,basals),bagOfSpecies,drop=F]
        # secondary extinctions: species that have no preys left except basals
        secondary =  rownames(adjTmp)[rowSums(adjTmp)==0]
        
        bagOfSpecies = setdiff(bagOfSpecies,secondary)
        Fragilite[k+1,j] = length(secondary)/(dim(adj)[1]-length(basals))
      }
      k=k+1
    }
  }
  robustness = sapply(1:n,function(j) 1 - sum(Fragilite[,j]) )
  
  return(list(meanRobustness=mean(robustness),robustness=robustness,extinctionChain=extinction))
}


# Robustness over a region of networks

robustnessInRegion = function(adjs,n=100,extinction=c('random','rarest')){
  adjs = lapply(1:length(adjs),function(j) as.matrix(adjs[[j]]))
  speciesPerSite = lapply(1:length(adjs),function(j) rownames(adjs[[j]]))
  SpPresences = rep(0,length(unique(unlist(speciesPerSite))))
  names(SpPresences)= unique(unlist(speciesPerSite))
  # sum presences of each species across all sites
  for(sp in names(SpPresences)){
    SpPresences[sp] = sum(sapply(speciesPerSite,function(localPool) sp%in%localPool))
  }
  
  # Get basal species in each site 
  basals = list()
  for(i in 1:length(adjs)){
    basals[[i]] = rownames(adjs[[i]])[rowSums(adjs[[i]])==0] 
  }
  
  nSp = length(SpPresences)
  Fragilite = array(0,dim = c(nSp,n,length(adjs)))
  localExtinction = 0. * Fragilite
  extinctionChains = list()
  for(j in 1:n){
    print(paste('path ',j))
    bagOfSpecies = names(SpPresences)
    k=1
    extinctionChains[[j]] = rep(NA,nSp)
    
    while(length(bagOfSpecies)>0 & k<nSp){
      if(extinction=="random"){
        extinctionChains[[j]][k] = sample(bagOfSpecies,1)
      }else if(extinction=="rarest"){
        remainingPresences = SpPresences[bagOfSpecies]
        rarest = names(remainingPresences)[which(remainingPresences==min(remainingPresences))]
        extinctionChains[[j]][k] = sample(  rarest , 1 )
      }
      print(paste('Primary extinction n?',k))
      print(SpPresences[extinctionChains[[j]][k]])
      
      bagOfSpecies = setdiff(bagOfSpecies,extinctionChains[[j]][k])
      if(length(bagOfSpecies)>0){
        Left = NULL
        for(i in 1:length(adjs)){
          adj = adjs[[i]]
          # Report local extinction event
          if(extinctionChains[[j]][k]%in%rownames(adj)){localExtinction[k,j,i] =1}
          localLeft = intersect(bagOfSpecies,rownames(adj))
          # consumers adjacency (without the original basal species)
          adjTmp = adj[setdiff(localLeft,basals[[i]]),localLeft,drop=F]
          #print(paste('Left after extinction:',length(localLeft)))
          # consumers who doesn't have any prey left + non-extinct basals
          localLeftAfterSecondaryExt = setdiff( c( rownames(adjTmp)[rowSums(adjTmp)>0] , basals[[i]] ) , extinctionChains[[j]] )
          #print(paste('Left after 2ndary extinction:', length(localLeftAfterSecondaryExt)))
          Left = union(Left,localLeftAfterSecondaryExt)
          if(i/50 == round(i/50)){
            flush.console()
            cat('   \r    Processed ',100*i/length(adjs),'%...    \r   ')
          }
        }
        secondary = setdiff(bagOfSpecies,Left)
        bagOfSpecies = Left
        print(paste(length(bagOfSpecies),' left species globally after secondary extinctions'))
        Fragilite[k+1,j,] = sapply(1:length(adjs),function(i){
          adj = adjs[[i]]
          length(intersect(secondary,rownames(adj))) / (dim(adj)[1]-length(basals[[i]])) })
      }
      k=k+1
    }
  }
  
  robustness = matrix(NA,n,length(adjs)) 
  for(j in 1:n){
    for(i in 1:length(adjs)){
      robustness[j,i] = 1 - sum( localExtinction[,j,i] * Fragilite[,j,i] )    
    }
  }
  
  #plot(1:sum(localExtinction[,j,i]),ratioSpLeft[localExtinction[,j,i]==1,j,i],type="l",col="red")
  #lines(1:sum(localExtinction[,j,i]),1-(1:sum(localExtinction[,j,i]))/sum(localExtinction[,j,i]))
      
  meanRobust = apply(robustness,2,mean)
  return(list(meanAUCs=meanRobust,AUCs=robustness,fragilite=Fragilite,extinctionChains=extinctionChains))
}

#####
# get_foodweb_metrics
#####

get_foodweb_metrics = function(gList,
                               tmpSavePath=NULL,
                               metrics=c('nV','nE','density','dirCon','modul','clust',
                                         'omniLev','omniProp','caniProp','predPerPrey','preyPerPred',
                                         'meanSWtrophLevel','meanShortestPathLength','basalProp','interProp','topProp',
                                         'vulnerabilitySD','generalitySD','sd_predPerPrey','skew_predPerPrey',
                                         'sd_preyPerPred','skew_preyPerPred','transitivity','diameter',
                                         'mean_distance','assortativity','tropLen',
                                         'mean_tropLevel','medi_tropLevel')){
  
  
  df=data.frame(nV = rep(NA,length(gList)),
                nE = NA,
                density = NA,
                dirCon = NA,
                modul = NA,
                clust = NA,
                omniLev = NA,
                omniProp = NA,
                caniProp = NA,
                predPerPrey = NA,
                sd_predPerPrey = NA,
                skew_predPerPrey = NA,
                preyPerPred = NA,
                sd_preyPerPred = NA,
                skew_preyPerPred = NA,
                meanSWtrophLevel = NA,
                meanShortestPathLength = NA,
                basalProp = NA,
                interProp = NA,
                topProp = NA,
                vulnerabilitySD = NA,
                generalitySD = NA,
                transitivity = NA,
                diameter = NA,
                mean_distance = NA,
                assortativity = NA,
                tropLen = NA,
                mean_tropLevel = NA,
                medi_tropLevel = NA,
                tropLen_MK = NA,
                mean_tropLevel_MK = NA,
                medi_tropLevel_MK = NA)
  df = df[,metrics,drop=F]
  df$name = names(gList)
  
  for(i in 1:dim(df)[1]){
    print(i)
    g = gList[[i]]
    
    n = length(V(g))
    if("nV"%in%metrics){df$nV[i] = n}
    
    if("nE"%in%metrics){df$nE[i] = length(E(g))}
    
    adjmat = t(as.matrix(as_adj(g)));rownames(adjmat)=NULL;colnames(adjmat)=NULL # Fixed transpose Adj
    
    if('density'%in%metrics){
      df$density[i] = sum( (as.vector(adjmat[lower.tri(adjmat)]) + as.vector(t(adjmat)[lower.tri(t(adjmat))])) >0 ) / (n*(n-1)/2)
    }
    
    if('dirCon'%in%metrics){df$dirCon[i] = sum(as.vector(adjmat))/(n^2)}
    
    if('modul'%in%metrics){
      adj = as.matrix(as_adjacency_matrix(g))
      adj_undir = 1. * matrix(as.vector(adj + t(adj))>0,dim(adj)[1],dim(adj)[1])
      g_undir = graph_from_adjacency_matrix(adj_undir, mode = "undirected")
      degs = igraph::degree(g_undir)
      membership = igraph::components(g_undir,mode = "weak")$membership
      globMembership = as.character(membership)
      for(m in unique(membership)){
        g_sub = igraph::induced.subgraph(g_undir,vids=which(membership==m))
        spinGlassClust = cluster_spinglass(g_sub)
        globMembership[which(membership==m)] = paste(m,spinGlassClust$membership)
      }
      nLinks = sum(adj_undir[lower.tri(adj_undir)])
      comMatrix = 0. * adj_undir
      for(m in unique(globMembership)){
        comMatrix[which(globMembership==m),which(globMembership==m)] = 1
      }
      expectedRandomLinkProbaMatrix = matrix(degs,length(degs),1) %*% matrix(degs,1,length(degs)) / (2*nLinks-1)
      termsMat = comMatrix * (adj_undir - expectedRandomLinkProbaMatrix)
      df$modul[i] = sum(termsMat[lower.tri(termsMat)])/(2*nLinks)
    }
    
    if('clust'%in%metrics){
      df$clust[i] = clusteringCoef(adjmat) # 0.988 correlation with Kortsch version
    }
    if(sum(c('omniLev','omniProp','caniProp')%in%metrics)>0){
      omniLevels = omnivory.levels(g)
      if('omniLev'%in%metrics){df$omniLev[i] = mean(omniLevels)}
      if('omniProp'%in%metrics){df$omniProp[i] = sum(omniLevels>0)/n}
      if('caniProp'%in%metrics){df$caniProp[i] = sum(diag(adjmat))/n}
    }
    
    # In Degrees distribution moments
    if(sum(regexpr('predPerPrey',metrics)>0)>0){
      npreds = rowSums(adjmat)
      # average count of predators
      if('predPerPrey'%in%metrics){df$predPerPrey[i] = mean(npreds)}
      if('sd_predPerPrey'%in%metrics){df$sd_predPerPrey[i] = sd(npreds)}
      if('skew_predPerPrey'%in%metrics){df$skew_predPerPrey[i] = skewness(npreds)}
    }
    
    # Out Degrees distribution moments
    if(sum(regexpr('preyPerPred',metrics)>0)>0){
      npreys = colSums(adjmat)
      # average count of preys
      if('preyPerPred'%in%metrics){df$preyPerPred[i] = mean(npreys)}
      if('sd_preyPerPred'%in%metrics){df$sd_preyPerPred[i] = sd(npreys)}
      if('skew_preyPerPred'%in%metrics){df$skew_preyPerPred[i] = skewness(npreys)}
    }
    
    # mean short weighted trophic level
    if('meanSWtrophLevel'%in%metrics){df$meanSWtrophLevel[i] = mean(short.weighted.trophic.levels(g))}
    
    # 
    if('meanShortestPathLength'%in%metrics){df$meanShortestPathLength[i] = mean.shortest.path.length(g)}
    
    
    # proportion of basal species
    if('basalProp'%in%metrics){df$basalProp[i] = basal.species.proportion(g)}
    # proportion of top species
    if('topProp'%in%metrics){df$topProp[i] = top.species.proportion(g)}
    # proportion of intermediary species (non basal and non top)
    if('interProp'%in%metrics){df$interProp[i] = intermediary.species.proportion(g)}
    
    if('vulnerabilitySD'%in%metrics){df$vulnerabilitySD[i] = sd(igraph::degree(g, mode = c("in")))  }# nodes in-degrees
    if('generalitySD'%in%metrics){df$generalitySD[i] = sd(igraph::degree(g, mode = c("out")))}# nodes out-degrees
    
    ## transitivity or clustering coefficient
    if('transitivity'%in%metrics){df$transitivity[i] = transitivity(g)}
    ## diameter or longuest shortest path
    if('diameter'%in%metrics){df$diameter[i]= diameter(g, directed = TRUE)}
    ## average shortest paths
    if('mean_distance'%in%metrics){df$mean_distance[i] = mean_distance(g, directed=TRUE)}
    ## assortativity or degree correlation
    if('assortativity'%in%metrics){df$assortativity[i] = assortativity.degree(g, directed=TRUE)}
    
    # Trophic levels metrics
    if(sum(c('tropLen','mean_tropLevel','medi_tropLevel')%in%metrics)>0){
      trophicLevels = trophiclevel(g)
      if('tropLen'%in%metrics){df$tropLen[i] = max(trophicLevels) - min(trophicLevels)}
      if('mean_tropLevel'%in%metrics){df$mean_tropLevel[i] = mean(trophicLevels)}
      if('medi_tropLevel'%in%metrics){df$medi_tropLevel[i] = median(trophicLevels)}
    }
    
    
    if(sum(c('tropLen_MK','mean_tropLevel_MK','medi_tropLevel_MK')%in%metrics)>0){
      trophicLevels = tryCatch(trophicLevelMackay(g),
               error=function(cond) {
                 message(cond)
                 return(rep(0, length(V(g)) )  )
               })
      if('tropLen_MK'%in%metrics){df$tropLen_MK[i] = max(trophicLevels) - min(trophicLevels)}
      if('mean_tropLevel_MK'%in%metrics){df$mean_tropLevel_MK[i] = mean(trophicLevels)}
      if('medi_tropLevel_MK'%in%metrics){df$medi_tropLevel_MK[i] = median(trophicLevels)}
    }
    
    
    if(!is.null(tmpSavePath)){
      write.table(df,tmpSavePath,sep=";",row.names=F,col.names=T)
    }
  }
  return(df)
}

#####
# Shortest-Paths embedding
#####

# Embedding equivalent au SPKernel:
CalculateShortestPathEmbedding = function(G,nodeLabelAttribName=NULL,directedPath=T){
  G.floyd <- as.list(rep(NA, length(G)))
  pathsList = as.list(rep(NA, length(G)))
  for (i in 1:length(G)){
    
    # On calcule la transformee de Floyd du graphe (graphe des plus courts chemins) 
    if(directedPath){
      D = igraph::distances(G[[i]],mode="out")
    }else{
      D = igraph::distances(G[[i]],mode="all")
    }
    # On marque les chemins statiques 
    diag(D) = NA
    
    df = expand.grid(from=1:length(V(G[[i]])),to=1:length(V(G[[i]])),weight=NA)
    df$weight = as.vector(D)
    # On retire les chemins statiques
    df = df[!is.na(df$weight),]
    
    if(!is.null(nodeLabelAttribName)){
      # On associe chaque chemin a un type 
      # au format: LabelDepart-Longueur-LabelArrivee
      nodeLabel = igraph::get.vertex.attribute(G[[i]],nodeLabelAttribName)
      tmp=data.frame(from=1:length(nodeLabel),fromLabel= nodeLabel)   
      df = merge(df,tmp,by="from")
      colnames(tmp) = c('to','toLabel')
      df = merge(df,tmp,by="to")
      paths = paste(df$fromLabel,'-',df$weight,'-',df$toLabel,sep="")
    }else{
      paths = df$weight
    }
    
    # On compte chaque type
    pathsList[[i]] = table(paths)
    
    if(i/10==round(i/10)){
      flush.console()
      cat('    \r        Process...',100*i/length(G),'%           \r ')
    }
  }
  
  # On fait l'union des types de chemins (futures colonnes de l embedding)
  pathsVocab = unique(names(unlist(pathsList)))
  
  # On cree l embedding
  ebd = matrix(0,length(pathsList),length(pathsVocab))
  colnames(ebd) = pathsVocab
  for(i in 1:length(G)){
    n = length(V(G[[i]]))
    # Size-Normalize each counts vector by the number of distinct nodes pairs in its graph
    ebd[i,names(pathsList[[i]])] = pathsList[[i]]/(n*(n-1))
  }
  return(ebd)
}

#####
# Triangular motifs count
#####

Count.Triangular.Directed.Motifs = function(Adj){
  counts = motifs(graph_from_adjacency_matrix(Adj),size=3)
  counts = counts[!is.na(counts)]
  counts = counts[c(9,5,10,6,2,1,3,11,12,4,7,8,13)]
  return(counts)
}

#####
# Generic function to get an embedding
#####

get.embedding = function(gList,ebdType=NULL,params=NULL,scaleCols=F){
  # gList: list of igraph objects, each element is a directed unweighted graph, possibly having nodes attributes
  # ebdType: character, available types are "Metrics" , "Motifs2Vec" , "Groups2Vec" , "ShortPaths2Vec , "Graph2Vec"
  # params: named list, see each method possible parameters
  if(ebdType=="Metrics"){
    # params
    # metrics : character vector, the list of metrics to compute, see get_foodweb_metrics for all available metrics.
    if(!is.null(params) && "metrics"%in%names(params)){
      ebd = get_foodweb_metrics(gList=gList,metrics=params$metrics)
    }else{
      ebd = get_foodweb_metrics(gList=gList)
      print("Info: no parameters found for applying the method")
    }
    if(scaleCols){
      ebd[,!colnames(ebd)=="name"] = scale(ebd[,!colnames(ebd)=="name"],center=T,scale=T)
    }
    
  }else if(ebdType=="Motifs2Vec"){
    # No parameters passed to this method
    ebd = as.data.frame(matrix(NA,length(gList),13))
    for(j in 1:length(gList)){
      Adj=  as.matrix(as_adjacency_matrix(gList[[j]]))
      n = dim(Adj)[1]
      ebd[j,]=Count.Triangular.Directed.Motifs(Adj)/(n*(n-1)*(n-2)/6)
      if(j/100==round(j/100)){
        flush.console()
        cat('   \r      Process...',100*j/length(gList),'%        \r ')}
    }
    
  }else if(ebdType=="Groups2Vec"){
    # params
    # groupsAttributeName: character
    # groupsIds: 
    ebd = as.data.frame(matrix(NA,length(gList),length(params$groupsIds)))
    colnames(ebd) = params$groupsIds
    for(j in 1:length(gList)){
      ebd[j,] = table( c(params$groupsIds , igraph::get.vertex.attribute(gList[[j]],params$groupsAttributeName)) ) - rep(1,length(params$groupsIds))  
      # To proportions
      ebd[j,] = ebd[j,] / sum(ebd[j,])
    }
    
  }else if(ebdType=="ShortPaths2Vec"){
    # params
    # nodeLabelAttribName: character
    if(!is.null(params$nodeLabelAttribName)){
      ebd = as.data.frame(CalculateShortestPathEmbedding(gList,nodeLabelAttribName = params$nodeLabelAttribName,directedPath = T))
    }else{
      ebd = as.data.frame(CalculateShortestPathEmbedding(gList,directedPath = T))
    }
  }else if(ebdType=="Graph2Vec"){
    # params
    # graphsDir: character
    # depth: integer
    # dim: integer
    # lab: logical
    ebd = py_get_graph2Vec_rectoverso_embedding(graphsPath=params$graphsDir,
                                     wl_it=as.integer(params$depth),
                                     ebd_dim=as.integer(params$dim),
                                     attrib=params$lab,
                                     down_samp=0.0001,
                                     n_epochs=as.integer(100),
                                     lr=0.015)  
    
  }
  return(ebd)
}

#####
# get_RFseparationPower
#####

get_RFseparationPower = function(X,Y){
  # Check NAs
  NAs = sapply(1:dim(X)[1],function(i)sum(is.na(X[i,])))
  nLinesWithNA = sum(NAs>0)

  if(nLinesWithNA>0){
    print(paste('Error: X has ',nLinesWithNA,'row containing with NA(s)'))
    return(NULL)
  }else if(is.numeric(Y)){
    # REGRESSION
    
    # Select mtry
    nbvars <- 1:(ncol(X) - 1)
    oobsMtry <- sapply(nbvars, function(nbv) {
      RF <- randomForest(Y ~ ., data = X, ntree = 100, mtry = nbv)
      return(RF$rsq[length(RF$rsq)])})
    
    # RF with optimal mtry
    RFDef = randomForest(Y ~ ., 
                         data = X,
                         ntree=300,
                         mtry=nbvars[which.max(oobsMtry)])
    
    # Compute OOB pseudo R2 
    pseudoR2 = as.numeric(RFDef$rsq[RFDef$ntree])
    
    sepPowers = list(pseudoR2,RFDef$predicted)
    names(sepPowers) = c('OobPseudoR2','OobPrediction')
    return(sepPowers)
    
  }else{
    # CLASSIFICATION
    # Select mtry
    nbvars <- 1:(ncol(X) - 1)
    oobsMtry <- sapply(nbvars, function(nbv) {
      RF <- randomForest(Y ~ ., data = X, ntree = 50, mtry = nbv)
      return(as.numeric(RF$err.rate[RF$ntree, "OOB"]))})
    
    # RF with optimal mtry
    RFDef = randomForest(Y ~ ., 
                         data = X,
                         ntree=300,
                         mtry=nbvars[which.min(oobsMtry)])
    
    # Compute uniform OOB classification error rate
    errRate = as.numeric(RFDef$err.rate[RFDef$ntree, "OOB"])
    
    # Compute Proba-weighted OOB classification error rate
    OOB_votes = RFDef$votes
    w = apply(OOB_votes,1,max)
    OobPred = sapply(1:dim(OOB_votes)[1],function(j) colnames(OOB_votes)[which.max(OOB_votes[j,])] )
    errors = OobPred!=as.character(Y)
    wErrRate = sum(w * errors) /sum(w)
    
    # Confusions des classes en OOB
    confusionProbaMatrix = RFDef$confusion[,colnames(RFDef$confusion)!="class.error"]
    for(i in 1:dim(confusionProbaMatrix)[1]){
      confusionProbaMatrix[i,] = confusionProbaMatrix[i,] / sum(confusionProbaMatrix[i,])
    }
    
    sepPowers = list(errRate,wErrRate,OobPred,w,OOB_votes,confusionProbaMatrix)
    names(sepPowers) = c('OobErrRate','OobWeightedErrRate','OobPrediction','PredictionProba','OobVotes','confusionProbaMatrix')
    return(sepPowers)
  }
}

#####
# Plot points 2D
#####

writePlot = function(df,outputLocAndName,
                     isLabel=NULL,
                     isColor=NULL,
                     isSize=NULL,alpha=.6,
                     legend=F,
                     xlab="x",ylab="y"){
  
  p = ggplot()
  
  if(!is.null(isColor) & is.null(isSize)){
    p = p + geom_point(aes(x=df[,1],y=df[,2],colour=df[,isColor]),alpha=alpha)
  }else if(is.null(isColor) & !is.null(isSize)){
    if(is.numeric(df[,isSize])){
      p = p + geom_point(aes(x=df[,1],y=df[,2],size=df[,isSize]),alpha=alpha)
    }else{
      print('Error: column indicated by isSize should be numeric')
      return(NULL)
    }
  }else if(!is.null(isColor) & !is.null(isSize)){
    p = p + geom_point(aes(x=df[,1],y=df[,2],colour=df[,isColor],size=df[,isSize]),alpha=alpha)
  }else if(is.null(isColor) & is.null(isSize)){
    p = p + geom_point(aes(x=df[,1],y=df[,2]),alpha=alpha)
  }
  
  if(!is.null(isLabel)){
    if(is.character(df[,isLabel])){
      p = p+ geom_text(data=df,aes(x=df[,1],y=df[,2],label=label,fontface=2),size=7)
    }else{
      print('Error: column indicated by isLabel should be a character vector')
      return(NULL)
    }
  }
  
  # Define points colors mode
  if(!is.null(isColor)){
    if(is.numeric(df[,isColor])){
      p = p + scale_colour_gradient(low = "darkorchid4",high="goldenrod")
    }
  }
  
  if(!legend){p=p+theme(legend.position = "None")}
  
  bandX = max(df[,1])-min(df[,1])
  bandY = max(df[,2])-min(df[,2])
  p = p+scale_x_continuous(limits=c(min(df[,1])-bandX/8,max(df[,1])+bandX/8))
  p = p+scale_y_continuous(limits=c(min(df[,2])-bandY/8,max(df[,2])+bandY/8))
  p = p + xlab(xlab) + ylab(ylab) +theme_bw()
  
  png(outputLocAndName,height=600,width=600)
  print(p)
  dev.off()
} 

#####
# Plot .graphml
#####

# Function to create a single graph plot with ggnet
graphml.to.single.plot.ggnet = function(fileLoc,plotName){
  G = read.graph(fileLoc,format="graphml")
  Gnet = intergraph::asNetwork(G)
  
  if("weight"%in%list.edge.attributes(Gnet)){
    bidirectional = T
    trophicIntake = get.edge.attribute(Gnet,attrname = "weight")==1
    edgesColors = ifelse(trophicIntake,"chartreuse4",'red2')  
    set.edge.attribute(Gnet, attrname="color", value= edgesColors, e=1:length(trophicIntake))
    set.edge.attribute(Gnet, attrname="lty", value= rep("dotted",length(trophicIntake)), e=1:length(trophicIntake))
  }else{
    bidirectional = F
    set.edge.attribute(Gnet, attrname="color", value= "red2")
  }
  
  
  if(is.directed(Gnet)){
    if(bidirectional){
      p = ggnet2(Gnet,label="vertex.names",arrow.size = 8, arrow.gap = 0.025,edge.color = "color",edge.lty="lty",edge.size=1)
    }
    p = ggnet2(Gnet,label="vertex.names",arrow.size = 8, arrow.gap = 0.025,edge.color = "color")
  }else{
    p = ggnet2(Gnet,label="vertex.names")
  }
  
  p = p + ggtitle(plotName)+ theme(plot.title = element_text(size=22),panel.border = element_rect(colour = "black", fill=NA, size=3))
  return(p)
}


# Function to create a single graph plot with ggnet
plot.graph.with.ggnet = function(G,coos=NULL){
  Gnet = intergraph::asNetwork(G)
  
  if("weight"%in%network::list.edge.attributes(Gnet)){
    bidirectional = T
    trophicIntake = network::get.edge.attribute(Gnet,attrname = "weight")==1
    edgesColors = ifelse(trophicIntake,"chartreuse4",'red2')  
    network::set.edge.attribute(Gnet, attrname="color", value= edgesColors, e=1:length(trophicIntake))
    network::set.edge.attribute(Gnet, attrname="lty", value= rep("dotted",length(trophicIntake)), e=1:length(trophicIntake))
  }else{
    bidirectional = F
    network::set.edge.attribute(Gnet, attrname="color", value= "red2")
  }
  
  if(network::is.directed(Gnet)){
    if(!is.null(coos)){
      Gnet = network::set.vertex.attribute(Gnet,'x',coos[,1])
      Gnet = network::set.vertex.attribute(Gnet,'y',coos[,2])
      p = GGally::ggnet2(Gnet,label="vertex.names",arrow.size = 8, arrow.gap = 0.025,edge.color = "color",mode=c('x','y'))
    }else{
      if(bidirectional){
        p = GGally::ggnet2(Gnet,label="vertex.names",arrow.size = 8, arrow.gap = 0.025,edge.color = "color",edge.lty="lty",edge.size=1)
      }
      p = GGally::ggnet2(Gnet,label="vertex.names",arrow.size = 8, arrow.gap = 0.025,edge.color = "color")
    }
  }else{
    p = GGally::ggnet2(Gnet,label="vertex.names")
  }
  
  print(p)
  return(NULL)
}

