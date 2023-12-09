rm(list=ls())
cat("\014")  #is the code to send ctrl+L to the console and therefore will clear the screen
graphics.off()
#dev.off()
plot.new()

# The purpose of this script is the inference of a candidate signaling network
 # via a "time-slice approach" from a underlying protein-protein interaction 
 # network and a phosphoproteomics time course, and subsequent extraction of 
 # candidate signaling paths.
# Steps:
 # Prepare the underlying network and phosphoproteomics data (A - D)
 # Map phosphoproteomics data onto the PPI network (E)
 # Generate Steiner Trees on slices of two adjacent timepoints (F)
 # Combine time slices and divide into subgraphs (G)
 # Extract paths and plot choronopaths (H, I)
 # Analyse and evaluate the derived network and paths (J, K)

library(ggpubr)
library(igraph)
library(stringr)
library(tidyverse)
library(graphlayouts) #(dplyr, ggplot2, tidyr, tibble)  https://datacarpentry.org/R-ecology-lesson/03-dplyr.html#:~:text=The%20tidyverse%20package%20is%20an,as%20ggplot2%20%2C%20tibble%20%2C%20etc.
library(ggraph)
library(reshape2)
library(animation)
library(patchwork)
library(gridGraphics)

source(paste('./20221024_PCST_ILP','20230327_PCST_ILP_DIR_FUN.R',sep='/'))
source(paste('./20221024_PCST_ILP/20221209_testing/testcases.R',sep=''))
source(paste('.','20230128_SteinerNet_mod.R',sep='/'))

loc= './'

testmodeAns= readline("Should I run in testmode (y/n)?:")

if(testmodeAns == 'n'){
  simulAns = readline('Use simulated data? (y/n):')
  gSol_exist= readline('In case of PCST: Do Gurobi solutions already exist (y/n)?:')
  if(simulAns == 'n'){
    # Phosphoproteomics data
    # Kanshin, 2015 (5s resolution) data:
    TCfile = './TabS1_mmc2_mod.csv'  #include location
  } else if (simulAns == 'y'){
    TCfile = readline("Give location/name of simulated phosphoproteomics dataset:\n")
  } else {
    stop('Incorrect simulation answer!')
  }
  PC_ans_vec= c('n','y')
  TCfile_split = strsplit(TCfile,'\\/|\\\\')[[1]]
  TCfile_loc = paste(
    paste(TCfile_split[1:length(TCfile_split)-1], collapse = '/'),
    '/', sep = '')
  TCfile_base = str_replace((tail(str_split(TCfile,'/')[[1]], n = 1)), '.csv', '')
  data_orig= read.csv(TCfile, stringsAsFactors=FALSE)
  #data_orig contains column "Gene", "ORF" and "pSite" with gene name; columns with "log2FC_..." (16-28), 
  #and with "log10int_..." (55-67)
  data_orig$Gene= stringr::str_to_title(data_orig$Gene)
} else if (testmodeAns== 'y'){
  PC_ans_vec= 'y'
  testcase_fun= testcase19_dir  ########### modify
  testString= '19_dir'        ########### modify
  gSol_exist= 'y'
} else {
  stop("Incorrect testmode answer!")
}

if (gSol_exist != 'y' & gSol_exist != 'n') {
  stop("Incorrect Gurobi-answer!")
}

options(warnPartialMatchArgs=FALSE) #prevent 'partial argument match' warnings in
#igraph plotting
options(warnPartialMatchDollar = TRUE) #https://stackoverflow.com/questions/43876289/selecting-non-existent-columns-in-data-table

TP_interval= 5
TPs= seq(5,60,TP_interval)
TPslices= c()
for (i in 1:(length(TPs)-1)){
  TPslices[i]= paste(TPs[i],TPs[i+1],sep='_')
}
FC_cutoff= 3   #avoid that first-shell expansion is greater than what STRING can do (OPTION 2: 'kansh')
expEvid_cut= 0.7
prizeToCost= 1 
max_pathLenSt= 10  ########### max. considered depth in extracting paths
######## max. Steiner nodes between terminals in extracting paths
#(2 matches the 1-node expansion from signif. prot. in STRING)
max_St_betwSliceTerm_class= 2 #for classical Steiner tree
max_St_betwSliceTerm_PCST= 2 #for PCST; #5 is too slow for real data
  # in testmode, max_St_betwSliceTerm is set ~line 2017

today= Sys.Date()
today=gsub("-","",today)
plotDir= paste(loc,'plots/',today,'/',sep='')
if (! dir.exists(plotDir)){
  dir.create(plotDir)
}


################ findSubgraphs ###################
findSubgraphs= function(G) {   #G can be matrix or iGraph graph
  
  #findSubgraphs: find connected subgraph in graph: start at one node,
  #go to next node you can reach; continue for all nodes until all connected nodes are visited
  #then do the same from first unconnected node
  
  #allAllocNodes= c() #all nodes that have been allocated to subnetworks
  G_orig= G #original G, before removing subgraphs
  G_undirMtx= as_adjacency_matrix(as.undirected(G)) #undirected to find both up-
  #and downstream neighbors
  subgraphs= list()
  
  last_node=c()
  while (dim(G_undirMtx)[1]>0) {  #while there are still unallocated nodes in network
    #'as.matrix' in case graph is reduced to one element
    startIdx= 1
    #print(paste('\nstart node:',startIdx,V(G_orig)$name[startIdx]))
    inSet_loc_idxes= c(startIdx)  #local indices: reflecting changing size of G
    inSet_names= c(rownames(G_undirMtx)[startIdx]) #add idx of starting node to set
    toAdd= which(as.matrix(G_undirMtx)[startIdx,]>=1) #connected nodes
    #print(paste('connected:',paste(toAdd,collapse=','),
    #         paste(V(G_orig)$name[toAdd],collapse=',')))
    while(length(toAdd)>0) {  #while there are nodes to be added
      toAddCand= which(G_undirMtx[toAdd[1],]>=1)
      toAdd= c(toAdd,toAddCand[is.na(match(toAddCand,inSet_loc_idxes))]) #candidate nodes that are not already in set
      inSet_loc_idxes= c(inSet_loc_idxes,toAdd[1])  #add evaluated candidate from set...
      toAdd= toAdd[-1] #and remove it from toAdd
      toAdd= unique(toAdd)
      #print(paste('to add:',paste(toAdd,collapse=', ')))
    }
    inSet_names= rownames(G_undirMtx)[inSet_loc_idxes]
    #print(paste(c('in set indices:',inSet_loc_idxes),sep=' '))
    #print(paste(c('in set names:',inSet_names),sep=' '))
    inSet_glob_idxes= match(inSet_names,rownames(as.matrix(G_orig))) #indices in original G
    subgraphs[[length(subgraphs)+1]]= inSet_glob_idxes #node indices
    
    #allAllocNodes= c(allAllocNodes,inSet_names)
    if(dim(G_undirMtx)[1]-length(inSet_loc_idxes)==1){ #last node reached
      last_node= row.names(G_undirMtx)[-inSet_loc_idxes]
      #print('last node in set index: 1')
      #print(paste(c('last node in set name:',last_node),sep=' '))
      inSet_glob_idxes= match(last_node,rownames(as.matrix(G_orig))) #indices in original G
      subgraphs[[length(subgraphs)+1]]= inSet_glob_idxes #node indices
      break()
    }
    G_undirMtx= G_undirMtx[-inSet_loc_idxes,-inSet_loc_idxes] #remove allocated nodes
    #print('remaining graph:')
    #print(G_undirMtx)
  }
  return(subgraphs)
}
#####################################################

################ get Steiner trees ####################
getSteiner= function(underlNetw,subgr,t1,t2,subgr_idx,design,
                     FC_thres,expEvid_cut,prizeToCost,max_St_betwTerm,gSol_exist,testmodeAns) {
  #find Steiner trees for given nodes in subgraph
  #use directed edges (both directions for undirected edge)
  
  #design choices:
  #fullCon: fully connected: connect all nodes via minimum-spanning tree
  #visitOne: one terminal per timepoint is visited (ensured by making flow go
  #out of terminals at first and into terminals at second TP)
  #visitAll: all terminals of given TP can be visited by making edges at each
  #TP undirected => allows signaling between signif. proteins at same TP; this is more consistent with the way PCST is implemented
  
  #input: for PCST, undelNetwork needs to have 'prizes' and 'costs'
  
  if(design== 'fullCon_visitOne' | design=='fullCon_visitAll'){
    algo= "KB"   #Kruskal-Based Heuristic (min. span. tree)
  } else if (design== '2node_shortestP') {
    algo= "SP"
  }
  
  #all nodes not at t1 or t2 are potential Steiner nodes and are assigned intermediate time
  V(underlNetw)$t[is.na(V(underlNetw)$t)]= t1+(t2-t1)/2
  V(underlNetw)$t[V(underlNetw)$t!= t1 & V(underlNetw)$t!= t2]= t1+(t2-t1)/2
  V(underlNetw)$type= rep('Steiner',length(V(underlNetw)))
  V(underlNetw)$type[V(underlNetw)$t==t1]= 'out-term'
  V(underlNetw)$type[V(underlNetw)$t==t2]= 'in-term'
  
  subgr_nodes= V(underlNetw)$name[subgr]
  curr_subgr= subgraph(underlNetw,subgr)
  #plot(curr_subgr,vertex.size=8,vertex.label=NA)
  out_term= V(curr_subgr)$name[V(curr_subgr)$t== t1]
  in_term= V(curr_subgr)$name[V(curr_subgr)$t== t2]
  curr_terminals= c(out_term, in_term)
  curr_terminals_t= c(V(curr_subgr)$t[V(curr_subgr)$t== t1],V(curr_subgr)$t[V(curr_subgr)$t== t2])
  noTermIn= length(curr_terminals) #number of terminals in input
  
  if (length(curr_terminals) <= 1) {  # >1 protein required for tree
    print(paste("No/Single node in subgraph",subgr_idx,':',curr_terminals))
    curr_tree0= make_empty_graph(n=0, directed=TRUE)
    unallocNodes= cbind.data.frame(curr_terminals,curr_terminals_t)
    return(list(curr_tree0, unallocNodes, noTermIn, 0))
  } else {
    if(design != 'prize_collecting') {
      print(paste('Calculating Classical Steiner Tree for ',t1,' - ',t2,' (',length(curr_terminals),' terminals) ...',sep=''))
      curr_tree0= steinertree_mod(type= algo, terminals= curr_terminals, 
                                  graph= curr_subgr, color=FALSE,optimize=FALSE,merge=FALSE)[[1]]
      print('Done!')
      V(curr_tree0)$color= NA
      V(curr_tree0)$color[V(curr_tree0)$name %in% curr_terminals &
                            V(curr_tree0)$t ==t1  ]= 'yellow'
      V(curr_tree0)$color[V(curr_tree0)$name %in% curr_terminals &
                            V(curr_tree0)$t ==t2  ]= 'orange'
      
      #Trim to remove paths with number of Steiner nodes above limit
      print('  Start trimming.. (class. Steiner)')
      #find all nodes on directed paths between out- and in-termini
      #"simple path": path on which each vertex is visited only once
      V_onPaths_idxes= c()  #will include terminal nodes
      for(o in out_term) {
        #print('Trimming to paths between terminals: Please wait!')
        V_onPaths_idxes= c(V_onPaths_idxes, unlist(
          all_simple_paths(graph=curr_tree0,from=o,to=in_term,mode="out",cutoff=max_St_betwTerm+1) ))  
        #cutoff necessary for speed
        print(V_onPaths_idxes)
      }
      V_onPaths_idxes= unique(V_onPaths_idxes)
      #remove nodes that are not on these paths
      #this also removes termini that are not on path from in-terminus to out-terminus
      curr_tree= subgraph(graph= curr_tree0,vids= V_onPaths_idxes) 
      print('  Trimming complete (class. Steiner)!')
    } else if(design== 'prize_collecting') {
      ### PRIZE COLLECTING
      #make graph directed
      curr_subgr_dir= as.directed(graph=curr_subgr,mode="mutual")
      #remove outgoing edges from in-terminals and vice versa
      length(E(curr_subgr_dir))
      for(i in out_term){
        curr_subgr_dir= delete.edges(curr_subgr_dir,
                                     incident(graph=curr_subgr_dir,v=i,mode='in'))
      }
      for(i in in_term){
        curr_subgr_dir= delete.edges(curr_subgr_dir,
                                     incident(graph=curr_subgr_dir,v=i,mode='out'))
      }
      length(E(curr_subgr_dir))
      
      print('Calculating tree---')
      PCST_OUT= PCST_DIR(dir_graph= curr_subgr_dir,out_term,in_term,t1,t2,FC_thres=FC_cutoff,expEvid_cut=expEvid_cut, 
                         prizeToCost,max_St_betwTerm=max_St_betwTerm,gSol_exist,
                         testmodeAns,testcase_fun=testcase_fun,testString=testString)
      curr_tree= PCST_OUT[[1]]
      print(paste('Runtime',PCST_OUT[[2]],'min',sep=' '))
      if(length(V(curr_tree))==0) {
        unallocNodes= cbind.data.frame(curr_terminals,curr_terminals_t)
        return(list(curr_tree, unallocNodes, noTermIn, 0))
      }
    } #end of design== 'prize-collecting'
    
    # distances(graph=curr_tree,v=V(curr_tree)[V(curr_tree)$type=='out-term'],
    #          to=V(curr_tree)[V(curr_tree)$type=='in-term'],mode='out')
    
    unallocNodes= V(curr_tree)$name[! V(curr_tree)$name %in% curr_terminals]
  } #end of if-else(number of terminals)
  
  
  if(length(E(curr_tree))> 0) {
    #if using classical Steiner tree (not prize collecting):
    #as edges are undirected, create reverse edge for each edge, except:
    #-they are from diff. TPs, in this case earlier TP is 'from'
    
    if(design != 'prize_collecting'){
      #terminal edges from the same TP are treated according to design choice:
      #if 'visitAll': they are connected via undirected edges
      #if 'visitOne': terminal is 'from' if it is at t1 and 'to' if it is at t2
      
      print('Generate edges in both directions between Steiner nodes...')
      
      #a, all Steiner nodes with bi-directional edges in between
      sS= as.directed(      #igraph object
        subgraph(graph=curr_tree,vids=V(curr_tree)[V(curr_tree)$type=='Steiner']),
        mode='mutual')##
      #b, add terminal nodes w/o edges
      V_toAdd= V(curr_tree)[V(curr_tree)$type=='in-term' | V(curr_tree)$type=='out-term']  #igraph vertices
      g_toAdd= as.directed(subgraph(curr_tree,
                                    v=V(curr_tree)[V(curr_tree)$type=='in-term' | V(curr_tree)$type=='out-term'])) #igraph object
      g_toAdd= delete_edges(g_toAdd,edges=E(g_toAdd)) #delete all edges
      
      stS1= sS + g_toAdd  #merge
      stS1= clean_union_attributes(stS1)
      # par(mar=c(0,0,3,0))
      # plot(stS1,vertex.label.cex=0.5,arrow.width=0.1,layout=layout_with_fr(stS1),
      #      main='Graph with edge deletions') 
      #only edges between pairs of Steiner nodes exist so far
      #c, add edges to/from terminals in correct orientation
      stS2= stS1
      if(design== 'fullCon_visitAll'){
        for(v_idx in seq_along(V(curr_tree)$name[V(curr_tree)$type=='out-term' | V(curr_tree)$type=='in-term'])) {
          i= V(curr_tree)$name[V(curr_tree)$type=='out-term' | V(curr_tree)$type=='in-term'][v_idx]
          aV= neighbors(graph=curr_tree,v=i,mode='all')$name  #'aV' = adjacent vertices
          i_time= V(curr_tree)$t[V(curr_tree)$type=='out-term' | V(curr_tree)$type=='in-term'][v_idx]
          stS2= add.edges(graph= stS2, edges=as.character(   #both directions
            matrix(c(rep(i,length(aV)),aV,aV,rep(i,length(aV))),nrow=2,byrow = TRUE)),
            attr=list(costs=rep(NA,2*length(aV)), t= rep(i_time,2*length(aV)),type=rep('terminal',2*length(aV)) ) )
        }
      } else if(design== 'fullCon_visitOne') {
        #out-edges from out-term, in-edges to in-term
        for(v_idx in seq_along(V(curr_tree)$name[V(curr_tree)$type=='out-term'])){
          i= V(curr_tree)$name[V(curr_tree)$type=='out-term'][v_idx]
          aV= neighbors(graph=curr_tree,v=i,mode='all')$name
          i_time= V(curr_tree)$t[V(curr_tree)$type=='out-term'][v_idx]
          stS2= add.edges(graph= stS2, 
                          edges=as.character(matrix(c(rep(i,length(aV)),aV),nrow=2,byrow = TRUE)),
                          attr=list(costs=rep(NA,length(aV)), t= rep(i_time,length(aV)),type=rep('terminal',length(aV))) )
        }
        for(v_idx in seq_along(V(curr_tree)$name[V(curr_tree)$type=='in-term'])){
          i= V(curr_tree)$name[V(curr_tree)$type=='in-term'][v_idx]
          aV= neighbors(graph=curr_tree,v=i,mode='all')$name
          i_time= V(curr_tree)$t[V(curr_tree)$type=='in-term'][v_idx]
          stS2= add.edges(graph= stS2, 
                          edges=as.character(matrix(c(aV,rep(i,length(aV))),nrow=2,byrow = TRUE)),
                          attr=list(costs=rep(NA,length(aV)), t= rep(i_time,length(aV)),type=rep('terminal',length(aV)) ) )
        }
        print('Done!')
      }
    } else { #PCST
      stS2= curr_tree    
      plot(stS2,vertex.label.cex=0.5,arrow.width=0.1)
    }
  } else { #no edges
    stS2= curr_tree ##
  } 
  
  noTermOut= sum(V(stS2)$type=='in-term')+sum(V(stS2)$type=='out-term')
  return(list(stS2,unallocNodes,noTermIn,noTermOut))  #new tree, unalloc nodes, #of terminals in input, #in output
}
########################################

depthFirstSearch
# Aim: find all paths over consecutive TPs, without skipping TP
# Store three lists: 1, stack, 2, recorded paths 3, visited nodes
# Start from a node significant at the first TP; add it to stack
# Move on to adjacent node: any adjacent node is allowed, but if it is not of
 # the subsequent TP, count it as Steiner node
 #record first node as visited; add new node on top of stack
# Keep moving and recording visited notdes until no further consecutive node or 
 #max. Steiner limit is exceeded;
 # then record the path (immediately trim Steiner nodes from end)
 # remove the last node from stack and try new adjacent node to node that is now
 # on top of the stack if adjacent node has not yet been visited
 # record all paths in this way; if node connects to already visited node, search
  # all recorded paths for this node; consider their downstream protions, but
  # re-score all nodes to the subsequent TP as Steiner nodes

 #As you allow other TPs than the first as starting TP, move on to next TP
  # repeat the algorithm keeping the recorded paths and visited nodes you already
  # have
  # before starting, eliminate all start nodes that have already be recorded on
   # paths with their TP of significance


################################
################################
extract_paths= function(currEdges_exp,max_St_betwSliceTerm,max_pathLenSt,end_prot){
  #function to extract paths from tree using given source terminals
  #input: currEdges_exp: DF of DIRECTED graph with columns 'node_from/to','TP_node_from/to', 
  #'type_from/to','idx_from/to';
  #end_prot (optional) is a single node where paths should end
  #core algorithm: add nodes to path as long as possible (if multiple proteins are
  #to be added to node n, fork into multiple paths at n+1); once no new proteins to be
  #added, or max. depth is reached, or previous protein in path is visited, store the 
  #path as complete
  #this function does NOT ensure temporal order of nodes
  
  end_prot[is.na(end_prot)]= ''  #use empty value, rather than NA for 'no end_prot'
  
  idxes_from= which(currEdges_exp$type_from=='out-term')
  ##idxes_from= 3 ####### for testing!
  termNames_from= currEdges_exp$node_from[idxes_from]
  TP_from= currEdges_exp$TP_from[idxes_from]
  types_from= currEdges_exp$type_from[idxes_from]
  uniqIdx= which(!duplicated(termNames_from))
  sourceIdxes= idxes_from[uniqIdx]
  sourceNodes= termNames_from[uniqIdx]
  sourceTPs= TP_from[uniqIdx]
  sourceTypes= types_from[uniqIdx]
  cat('\n\nSource Nodes:')
  print(sourceNodes)
  
  #"paths" is the list at start of going deeper into paths, while paths2 is modified during stepping in;
  #at the end of the round, "paths" is created from updated "paths2"
  #as lists
  pathsN= as.list(sourceNodes)  #nodes on path
  pathsTPs= as.list(sourceTPs) #TPs on path
  pathsType= as.list(sourceTypes)
  pathsIdxes= as.list(sourceIdxes)
  depth= 1
  
  compl_pathsN= list()
  compl_pathsTPs= list()
  compl_pathsType= list()
  compl_pathsIdxes= list()
  
  while(length(pathsN)>0) {
    #print(paste('Depth:',depth))
    cat('\nDepth:',depth)
    cat('number pathsN: ',length(pathsN))
    #at each steps the paths lists are created new
    pathsN2= list()
    pathsTPs2= list()
    pathsType2= list()
    pathsIdxes2= list()
    for(pathIdx in seq_along(pathsN)) {
      upIdx= which(currEdges_exp$node_from== pathsN[[pathIdx]][depth] &  #indicates node that has edge from current node
                     currEdges_exp$TP_from== pathsTPs[[pathIdx]][depth])  #timepoint from which edge originates is current TP
      currEdges_exp$node_from[upIdx]
      
      print(paste('path (before):', paste(pathsN[[pathIdx]],collapse=' - ')))
      cat('Candidates:',currEdges_exp$node_to[upIdx],'(t: ',currEdges_exp$TP_to[upIdx],')\n')
      if (length(upIdx)>0  & depth< max_pathLenSt) { #node(s) to be added to path
        candCount=0
        for(nN in upIdx) {  #nN= 'new node'
          nextNode= currEdges_exp$node_to[nN]
          candCount= candCount+1
          print(paste(' Candidate:', candCount,':',nextNode))
          
          #do not allow addition of nodes that are already in path (removed 20230321b)
          #(loops not allowed)
          # if(nextNode %in% pathsN[[pathIdx]]){
          #   if(end_prot=='' & any(pathsType[[pathIdx]]=='in-term')) { #only if no end_prot is required; last protein 
          #     #is not end_prot, otherwise would have terminated last round
          #     # & only save path if it contains in-term
          #     compl_pathsN[[length(compl_pathsN)+1]]= pathsN[[pathIdx]] #path is complete
          #     compl_pathsTPs[[length(compl_pathsTPs)+1]]= pathsTPs[[pathIdx]]
          #     compl_pathsType[[length(compl_pathsType)+1]]= pathsType[[pathIdx]]
          #     compl_pathsIdxes[[length(compl_pathsIdxes)+1]]= pathsIdxes[[pathIdx]]
          #     print(paste(' Completed (loop):',
          #                 paste(compl_pathsN[[length(compl_pathsN)]],collapse=' - ')))
          #     print(paste(compl_pathsType[[length(compl_pathsN)]],collapse=' - '))
          #   }
          #   next
          # }
          
          #if added protein is end-protein, add path WITH new protein to complete paths
          if(nextNode== end_prot){
            compl_pathsN[[length(compl_pathsN)+1]]= c(pathsN[[pathIdx]],currEdges_exp$node_to[nN])
            compl_pathsTPs[[length(compl_pathsTPs)+1]]= c(pathsTPs[[pathIdx]],currEdges_exp$TP_to[nN])
            compl_pathsType[[length(compl_pathsType)+1]]= c(pathsType[[pathIdx]],currEdges_exp$type_to[nN])
            compl_pathsIdxes[[length(compl_pathsIdxes)+1]]= c(pathsIdxes[[pathIdx]],currEdges_exp$idx_to[nN])
            print(paste(' Completed (end prot):',
                        paste(compl_pathsN[[length(compl_pathsN)]],collapse=' - ')))
            print(paste(compl_pathsType[[length(compl_pathsN)]],collapse=' - '))
            next
          }
          
          #if Steiner count is exceeded, save the path as it is (Steiner nodes later need to be trimmed from end)
          if(length(pathsType[[pathIdx]]) >= max_St_betwSliceTerm){
            finalIndixes= (length(pathsType[[pathIdx]])-max_St_betwSliceTerm+1):length(pathsType[[pathIdx]]) #indices to check for Steiner nodes; covering max. number
            if( all( (pathsType[[pathIdx]][finalIndixes])=='Steiner') &   #if they are all Steiner nodes and..
                currEdges_exp$type_to[nN]=='Steiner' ) { #...the new node is also a Steiner node -> save and exit
              compl_pathsN[[length(compl_pathsN)+1]]= pathsN[[pathIdx]] #path is complete
              compl_pathsTPs[[length(compl_pathsTPs)+1]]= pathsTPs[[pathIdx]]
              compl_pathsType[[length(compl_pathsType)+1]]= pathsType[[pathIdx]]
              compl_pathsIdxes[[length(compl_pathsIdxes)+1]]= pathsIdxes[[pathIdx]]
              print(paste(' Too many Steiner nodes:',
                          paste(compl_pathsN[[length(compl_pathsN)]],collapse=' - ')))
              next
            }
          }
          
          #otherwise just add the protein
          pathsN2[[length(pathsN2)+1]]= 
            c(pathsN[[pathIdx]],currEdges_exp$node_to[nN])
          pathsTPs2[[length(pathsTPs2)+1]]= 
            c(pathsTPs[[pathIdx]],currEdges_exp$TP_to[nN])
          pathsType2[[length(pathsType2)+1]]= 
            c(pathsType[[pathIdx]],currEdges_exp$type_to[nN])
          pathsIdxes2[[length(pathsIdxes2)+1]]= 
            c(pathsIdxes[[pathIdx]],currEdges_exp$idx_to[nN])
        }
      } else if (length(upIdx)==0 | depth== max_pathLenSt){ #no new node found or max path-length reached
        if(end_prot=='' & any(pathsType[[pathIdx]]=='in-term')) { #only if no end_prot is required; last protein 
          #is not end_prot, otherwise would have terminated last round
          # &there is in-term in path
          compl_pathsN[[length(compl_pathsN)+1]]= pathsN[[pathIdx]] #path is complete
          compl_pathsTPs[[length(compl_pathsTPs)+1]]= pathsTPs[[pathIdx]]
          compl_pathsType[[length(compl_pathsType)+1]]= pathsType[[pathIdx]]
          compl_pathsIdxes[[length(compl_pathsIdxes)+1]]= pathsIdxes[[pathIdx]]
          print(paste(' Completed (no more edges or max. depth):',
                      paste(compl_pathsN[[length(compl_pathsN)]],collapse=' - ')))
          print(paste(compl_pathsType[[length(compl_pathsN)]],collapse=' - '))
        }
      } else {
        stop ('ERROR: Unexpected case!')
      }
    }  
    pathsN= pathsN2  #update paths with new node
    pathsTPs= pathsTPs2
    pathsType= pathsType2
    pathsIdxes= pathsIdxes2
    depth= depth+1
  }
  return(list(compl_pathsN,compl_pathsTPs,compl_pathsType,compl_pathsIdxes))
}
##############################

##############################
chronoPaths= function(nodePaths_list, TPsPaths_list, typesPaths_list, sigDF, maxLen, plotmode,
                      plotSuffix,altColor_l= NULL) {
  #function to plot paths over time
  #input: nodePaths_list/TPsPaths_list/typesPaths_list: list with each entry corresponding to vector 
  #of nodes/TPs/types on path; : sigDF: DF that lists significant proteins and TP of 
  #significance; maxLen= maximum length in extracting paths; plotmode: indicates if all paths into one plot
  #(combPlot) or individual plots (indivPlots)
  #altColor_l is an optional list (matching the nodePaths_list) of colors to be applied to nodes
  
  if(plotmode!= 'indivPlots' & plotmode!= 'combPlot'){
    stop('Plotmode not recognized!')
  } else if (plotmode== 'combPlot') {
    plotSuffix= paste(plotSuffix,'_comb',sep='')
  } else if (plotmode== 'indivPlots') {
    plotSuffix= paste(plotSuffix,'_indiv',sep='')
  }
  
  TPs= sort(unique(unlist(TPsPaths_list)))
  
  #generate indices of trimmed paths
  uniq_pathNodes= unique(unlist(nodePaths_list))
  
  #manual selection:
  ########## regular #########
  node_idx_trim= cbind.data.frame(node=uniq_pathNodes, idx=seq_along(uniq_pathNodes),
                                  sigfTP= sigDF$firstSigTP[match(uniq_pathNodes,sigDF$Gene)])
  
  ## !!!!! you need to write node_idx_trim of CST and PCST to files and combine them manually !!!###
  #write.csv(node_idx_trim,file=paste(loc,today,'_nodeIdx_',plotSuffix,'.csv',sep=''),row.names=FALSE)
  #pdf(file=paste(plotDir,today,'_chronopaths_',plotSuffix,'.pdf',sep=''),width=7,height=14)  #inches
  ##### open file
  ########## if you instead want order of combination of CST and PCST
  node_idx_trim= read.csv(paste(loc,'20230328_CSTcombPCSTnodeIdx_FC3_maxTree2_maxPath10.csv',sep=''),
                          header=TRUE,stringsAsFactors = FALSE)
  if(plotmode=='indivPlots') {
    pdf(file=paste(plotDir,today,'_CSTcombPCSTchronopaths_',plotSuffix,'.pdf',sep=''),width=7,height=14)  #inches
  }
  ##########################
  node_idx_trim[1:5,]
  dim(node_idx_trim)
  
  if (plotmode== 'combPlot') {  #all paths into one plot
    dev.new(width=3.5,height=10.5,units='in',noRStudioGD = TRUE)
    par(mar= c(2.5,3,1,1))  #default mai: c(0.8, 1.0, 0.2, 0.2)
    print(
      plot(type='n',x=0,y=0,xlim= c(0,60),ylim=c(1,max(c(1,node_idx_trim$idx))),
           xlab= 't/s',ylab='',yaxt='n') )
    axis(side=2, at= node_idx_trim$idx[order(node_idx_trim$idx)], 
         labels = node_idx_trim$node[order(node_idx_trim$idx)],las=2,
         cex.axis=0.4,tck=1,col='grey',lty=2)  #cex.axis=1 #not working: col.ticks=axCol, col.labels= axCol
    axis(side=1,at=seq(0,60,5),tck=1,col='grey',lty=2) #tck=1 adds grid
  }

  #lineColors= palette(distinctColorPalette(length(nodePaths_list)))
  for (i in seq_along(nodePaths_list)){
    if(plotmode== 'indivPlots') {
      sum(duplicated(node_idx_trim$node))  #should be 0; each only at time of significance, or NA
      startNodeOnPath= nodePaths_list[[i]][1]
      endNodeOnPath= nodePaths_list[[i]][length(nodePaths_list[[i]])]
      plot(type='n',x=0,y=0,xlim= c(TPs[1],TPs[length(TPs)]),ylim=c(1,max(c(1,node_idx_trim$idx))),
           xlab= 't/s',ylab='',yaxt='n', main= paste('Path ',i,' (',startNodeOnPath,' - ',
                                                     endNodeOnPath,')',sep=''))
      axis(2, at= node_idx_trim$idx[order(node_idx_trim$idx)], 
           labels = node_idx_trim$node[order(node_idx_trim$idx)],las=2,
           cex.axis=0.4,tck=1,col='grey',lty=2)  #not working: col.ticks=axCol, col.labels= axCol
      axis(1,tck=1,col='grey',lty=2)  #tck=1 adds grid
    }
    
    if(!is.null(altColor_l)) {  #alternative color is specified
      colr= altColor_l[[i]]
    } else {     #alternative color is not specified
      # color by terminals
      colr= rep('black',length(nodePaths_list[[i]]))
      colr[(typesPaths_list[[i]]!='in-term' & typesPaths_list[[i]]!='out-term') &
             nodePaths_list[[i]] %in% sigDF$Gene]='green' #significant but not terminal (i.e. sigf. at diff. TP)
      colr[which(typesPaths_list[[i]]=='in-term' | typesPaths_list[[i]]=='out-term')]=
        'red' #terminal at it`s sigf. TP
    }
    
    TPs_onPath= TPsPaths_list[[i]]
    nodeIDs_onPath= node_idx_trim$idx[match(nodePaths_list[[i]],node_idx_trim$node)]
    points(TPs_onPath,nodeIDs_onPath,pch=21,bg=colr,cex=2)   #pch=16,col=color
    lines(TPs_onPath,nodeIDs_onPath,lwd=2)
  }
  
  if(plotmode== 'combPlot'){
    recPlot= recordPlot()
    dev.off()
    pdf(file=paste(plotDir,today,'_CSTcombPCSTchronopaths_',plotSuffix,'.pdf',sep=''),width=7,height=14)  #inches
    #png(file=paste(plotDir,today,'_CSTcombPCSTchronopaths_',plotSuffix,'.png',sep=''),width=7,height=14,
    #   unit='in',res=500)  #inches
    print(recPlot)
  } else {
    recPlot= plot.new()
  }
  dev.off()  #close both combPlot and indiv plots
  
  plot(x=0,y=0,type='n',axes=FALSE,bty='n',xlab='',ylab='')
  legend(-0.5,0,c('terminal','sigf., but not term.'),pch=21,pt.bg=c('red','green'),
         pt.cex=2, cex=.8, bty="n", ncol=1)
  return(list(node_idx_trim,recPlot))
}
#################################

clean_union_attributes= function(tree){
  #function to restore vertex attributes of tree to 'sigf','prizes','t','type','name' and 'time_slice'
  #and edge attributes to "type", "t", "costs",
  #after creating a graph.union, i.e. attributes will be 'sigf_1', 'sigf_2', ...
  #for each value of a attribute, the function will choose the first entry which is not NA
  
  tree_out= tree
  
  vertex.attributes(tree_out)= list(name= vertex.attributes(tree)$name) #keep 'name' as only attribute
  if(length(V(tree))>0) {
    #keywords= c('sigf','prizes','t','type')
    
    tree_attr= names(vertex.attributes(tree)) #e.g. "sigf_1" "sigf_2" "t_1""t_2" "color_1" "color_2", etc 
    attr_bases= unique(sapply(tree_attr,function(x) {str_split(x,'_')[[1]][1]}))  #e.g. sigf, t, color
    keywords= attr_bases[attr_bases!= 'name']
    
    out_list= list()
    for (k in keywords) {
      idxes= which(str_detect(names(vertex.attributes(tree)),
                              paste('^',k,'_*[0-9]*$',sep='')))
      DF= do.call('rbind',vertex.attributes(tree)[idxes])
      #DF= cbind.data.frame(vertex.attributes(tree)[idxes]) #list to DF
      if(k != 'time.slice') {
        first_Vattr_vec= unlist(lapply(data.frame(DF,stringsAsFactors = FALSE),
                                       FUN=function(x) {
                                         if(all(is.na(x))){
                                           return(NA)
                                         } else {
                                           return(x[min(which(!is.na(x)))])
                                         } 
                                       } ))
      } else {  #for time.slice make a list for each node of the time.slices in which it is found
        allElements= lapply(data.frame(DF,stringsAsFactors = FALSE),
                            function(x){ unlist(sapply(x,function(y){str_split(y,'_')}) )})
        first_Vattr_vec= lapply(allElements,function(x){     #e.g. 5_10 and 15 -> 5_10_15
          if(all(is.na(x))){
            return(NA)
          } else {
            return(x[!is.na(x)])
          } 
        } )
      } 
      names(first_Vattr_vec)= NULL
      if(length(first_Vattr_vec)>0){
        out_list[[length(out_list)+1]]= first_Vattr_vec
        names(out_list)[length(out_list)]= k
      }
    }
    vertex.attributes(tree_out)= c(vertex.attributes(tree_out),out_list)
    
    keywords= c('type','t','costs')
    out_list= list()
    for (k in keywords) {
      idxes= which(str_detect(names(edge.attributes(tree)),
                              paste('^',k,'_*[0-9]*$',sep='')))
      if (length(idxes)>0){
        DF= do.call('rbind',edge.attributes(tree)[idxes])
        firstAttr_vec= unlist(lapply(data.frame(DF,stringsAsFactors = FALSE),
                                     FUN=function(x) {
                                       if(all(is.na(x))) {
                                         return(NA)
                                       } else {
                                         return(x[min(which(!is.na(x)))])
                                       }
                                     }
        ))
        names(firstAttr_vec)= NULL
        out_list[[length(out_list)+1]]= firstAttr_vec
      } else { #keyword not present in input graph
        out_list[[length(out_list)+1]]= rep(NA,length(E(tree)))
      }
      names(out_list)[length(out_list)]= k
    }
    edge.attributes(tree_out)= out_list
  }
  return(tree_out)
}
#################################

plot_nodes_edges_atTP= function(nodeDF, edgeDF){
  #plot number of nodes and edges at each TP / timeslice
  #nodeDF needs column 't' and edgeDF needs columns 'TP_node1' and 'TP_node2'
  #even though columns are called '_from' and '_to', the DF is expected to be from
  #UNDIRECTED graph
  
  #Steiner nodes may be NA or 'Steiner', etc
  if( sum(is.na(nodeDF$t)) + sum(nodeDF$t== 'Steiner', na.rm=TRUE) > 0 ){
    print('Non-numeric nodes in nodeDF: Possibly Steiner nodes?')
    nonSteinerInd= (nodeDF$t != 'Steiner')
    nonSteinerInd[is.na(nonSteinerInd)]= FALSE
    nodeDF= nodeDF[nonSteinerInd,]
  }
  if( sum(is.na(edgeDF$TP_node1)) + sum(is.na(edgeDF$TP_node2)) + 
      sum(edgeDF$TP_node1== 'Steiner', na.rm=TRUE) + sum(edgeDF$TP_node2== 'Steiner', na.rm=TRUE) > 0 ){
    print('Non-numeric nodes in edgeDF: Possibly Steiner nodes?')
    nonSteinerInd= (edgeDF$TP_node1!= 'Steiner' & edgeDF$TP_node2!= 'Steiner') #indicator if node is not Steiner,
    #may contain NA
    nonSteinerInd[is.na(nonSteinerInd)]= FALSE
    edgeDF= edgeDF[nonSteinerInd,]
  }
  
  par(mar=c(3,3,3,0))
  nodeDF %>%
    count(t) %>%
    ggplot(aes(x=t,y=n))+
    geom_bar(stat='identity')+
    ggtitle('Signif. node counts in network')
  
  
  #sort as edges are undirected
  from_to_colNos= c(which(colnames(edgeDF)=='TP_node1'),which(colnames(edgeDF)=='TP_node2'))
  for(i in 1:dim(edgeDF)[1]){
    edgeDF[i,from_to_colNos]= sort(c(edgeDF$TP_node1[i],edgeDF$TP_node2[i]))
    edgeDF$TP_from_to[i]= paste(edgeDF[i,from_to_colNos],collapse='_')
  }
  uniqTP_DF= t(sapply(unique(edgeDF$TP_from_to), function(x){as.numeric(str_split(x,'_')[[1]])}))
  uniqTP_order= order(uniqTP_DF[,1],uniqTP_DF[,2])
  edgeDF$TP_from_to= factor(edgeDF$TP_from_to,levels= 
                              unique(edgeDF$TP_from_to)[uniqTP_order])
  
  edgeDF %>%
    count(TP_from_to) %>%
    ggplot(aes(x=TP_from_to,y=n))+
    geom_bar(stat='identity')+
    theme(axis.text= element_text(angle=90))+
    ggtitle('Signif. edge counts in network')
}
##################################

reduce_paths= function(paths,TPs,Types) {
  #function to remove duplicated and sub-paths
  
  ### First remove duplicated items
  #convert paths and TPs to vector and bind into DF
  paths_vec= unlist(lapply(paths,paste,collapse='_'))
  TPs_vec= unlist(lapply(TPs,paste,collapse='_'))
  DF= cbind.data.frame(paths_vec,TPs_vec)
  
  #find duplicated items in the DF
  dup_idx= which(duplicated(DF))
  
  #and remove from original lists
  if(length(dup_idx)>0){ #if dup_idx is length 0, the "-dup_idx" would create empty vector
    pathsOut= paths[-dup_idx]
    TPsOut= TPs[-dup_idx]
    TypesOut= Types[-dup_idx]
  } else if (length(dup_idx)==0){
    pathsOut= paths
    TPsOut= TPs
    TypesOut= Types
  }
  
  #### Then remove subpaths
  #convert paths and TPs to vector and bind into DF
  pathsOut_vec= unlist(lapply(pathsOut,paste,collapse='_'))
  TPsOut_vec= unlist(lapply(TPsOut,paste,collapse='_'))
  DFOut= cbind.data.frame(pathsOut_vec,TPsOut_vec)
  
  #detect which strings are substrings of other string in list
  substrIdx= which(
    sapply(pathsOut_vec,function(x){sum(str_detect(pathsOut_vec,x))>1}) &  # '>1', because needs to detect one more than itself
    sapply(TPsOut_vec,function(x){sum(str_detect(TPsOut_vec,x))>1})   #strictly you would need to check that TP match is in the same place as name match
  )
  if (length(substrIdx)==0){
    return(list(pathsOut,TPsOut,TypesOut))
  } else {
    pathsOut2= pathsOut[-substrIdx]
    TPsOut2= TPsOut[-substrIdx]
    TypesOut2= TypesOut[-substrIdx]
  }
  
  return(list(pathsOut2,TPsOut2,TypesOut2))
}
########################

########################
calcPathsLen= function(nodeDF, graph, t1, t2) {
  #summarize the number and lengths of paths between termini in the underlying network
  #nodeDF lists the terminals from/to which to calculate path lengths; needs to contain 'node' and 't' columns
  print(' Calculating path lengths ...')
  pathNoDF= data.frame(matrix(NA,ncol=3,nrow=0)) #lists number of simple path BELOW CUTOFF between termini
  pathLengthDist= data.frame(matrix(NA,ncol=3,nrow=0)) #lists all simple path length BELOW CUTOFF between termini
  for (o in nodeDF$nodes[nodeDF$t==t1]) {  #out out
    for (i in nodeDF$nodes[nodeDF$t==t2]) { #in terminals
      asp= all_simple_paths(graph=STRING_net,from=o,
                            to=i,mode="out",cutoff= max_St_betwSliceTerm+1 )
      pathNoDF= rbind.data.frame(pathNoDF,c(o,i,length(asp)))
      pathLengthDist= rbind.data.frame(pathLengthDist,
                                       cbind.data.frame(rep(o,length(asp)),rep(i,length(asp)),unlist(lapply(asp,
                                                                                                            FUN= function(x) {length(x)-1} ))) )  #-1 to go from nodes to edges
    }
  }
  colnames(pathNoDF)= c('out_term','in_term','noPaths')
  colnames(pathLengthDist)=c('out_term','in_term','pathLength')
  ggplot(data=pathLengthDist,aes(x=pathLength))+
    geom_histogram(position='dodge',binwidth=0.5)+  #aes(fill=paste0(out_term,in_term) ),
    #theme(legend.position="none")+
    geom_vline(xintercept= max_St_betwSliceTerm+1.2)+
    xlab('path length')+
    ggtitle(paste(t1,' to ',t2,'min'))
  return(pathLengthDist)
  print('Done!')
}
##########################

#############
plot_graph_overTime= function(all_t_graph_l,plotSuffix=plotSuffix,
                              layout=NULL,vertexColors=NULL, vertexSizes=NULL,
                              vertexLabCol=NULL, vertexLabCex= NULL,
                              vertexLabDist=NULL, vertexLabDegree=NULL,vertexSize_scaleFac=1) { 
  #plot graph with fill indicating time; ONLY FOR DISPLAY; only subgraph 1 plotted!!
  #all_t_graph_l needs to be a list of graphs over time (the first one of which is plotted)
  #optional parameters for layout, vertexColors and vertexSizes, etc. to override default
  
  all_t_graph_l
  if(length(V(all_t_graph_l[[1]]))>0) {
    vertex.attributes(all_t_graph_l[[1]])
    V(all_t_graph_l[[1]])$color= rep('white',length(V(all_t_graph_l[[1]])))
    V(all_t_graph_l[[1]])$frameCol= rep('black',length(V(all_t_graph_l[[1]])))
    V(all_t_graph_l[[1]])$frameCol[V(all_t_graph_l[[1]])$sigf]= 'red'
    greenLevel= V(all_t_graph_l[[1]])$t/60
    V(all_t_graph_l[[1]])$color= rgb(0.3,greenLevel,0.6)
    V(all_t_graph_l[[1]])$labCol= rep('white',length(greenLevel))
    V(all_t_graph_l[[1]])$labCol[greenLevel>0.3]='black' 
    set.seed(1)
    l_orig= layout.fruchterman.reingold(all_t_graph_l[[1]])
    l_nodes= V(all_t_graph_l[[1]])$name
    row.names(l_orig)= l_nodes
    
    ###
    #par(mar=c(0,0,0,0))
    plot_WH= 7  #plot width/height
    #rescale_mtx= matrix( c(rep(1/20,dim(l_orig)[1]),rep(1/20,dim(l_orig)[1])),ncol=2,byrow = FALSE ) #rescale manually to fit plot
    #rescale_mtx= matrix( c(rep(1/6,dim(l_orig)[1]),rep(1/6,dim(l_orig)[1])),ncol=2,byrow = FALSE ) #rescale manually to fit plot
    #rescale_mtx= matrix( c(rep(1/20,dim(l_orig)[1]),rep(1/20,dim(l_orig)[1])),ncol=2,byrow = FALSE ) #rescale manually to fit plot
    maxPlotDim= 3.8*max(apply(l_orig,2,function(x){max(x)-min(x)})) #max. dimension of the plot
    #factor 3.8 is empirical; not clear where it comes from
    rescale_mtx= matrix( c(rep(plot_WH/maxPlotDim,dim(l_orig)[1]),
                           rep(plot_WH/maxPlotDim,dim(l_orig)[1])),ncol=2,byrow = FALSE )
    l_resc= l_orig*rescale_mtx
    
    if (PC_ans== 'y' & testmodeAns=='n'){   # TODO: any shift necessary if simulation
      l= l_resc
      l[,1]= l[,1]-0.1  #shift in x-direction
      l[,2]= l[,2]-0.2  #shift in y-direction
    } else if (PC_ans== 'n' & testmodeAns=='n') {
      l= l_resc
      l[,1]= l[,1]-0.1  #shift in x-direction
    }
    
    #alternative user choices
    if(! is.null(layout)){
      l= layout
    }
    if( is.null(vertexSizes) ){
      vertexSizes= 10
    }
    if (is.null(vertexColors) ){
      vertexColors= V(all_t_graph_l[[1]])$color
    }
    if(is.null(vertexLabCol)){
      vertexLabCol= V(all_t_graph_l[[1]])$labCol
    }
    if(is.null(vertexLabDist)){
      vertexLabDist= 0
    }
    if(is.null(vertexLabDegree)){
      vertexLabDegree= 0
    }
    if( is.null(vertexLabCex) ){
      vertexLabCex= 0.4
    }
    
    range(l[,1])
    range(l[,2])
    
    xyMax= 1
    dev.new(width=4,height=4,unit='in',noRStudioGD = TRUE)
    par(mar=c(0,0.5,0,1))
    plot(all_t_graph_l[[1]],layout=l,vertex.color=vertexColors,vertex.size=vertexSizes,
         edge.arrow.width=.7, edge.arrow.size=.5,
         vertex.label.cex=vertexLabCex,vertex.frame.color=V(all_t_graph_l[[1]])$frameCol,
         vertex.frame.width=2,
         vertex.label.dist=vertexLabDist,vertex.label.degree=vertexLabDegree,vertex.label.color=vertexLabCol,
         rescale=FALSE,xlim=c(-xyMax,xyMax),ylim=c(-xyMax,xyMax))  #, axes=TRUE
    #main= paste('Subgraph','1','over time1'))  #####?,rescale=FALSE
    ## hardcode legend timepoints to avoid problems with missing trees at certain TPs
    leg_TPs_withIntermed= seq(5,60,2.5) #HARDCODED
    fullTP_idx= match(seq(5,60,5),leg_TPs_withIntermed)  #indices of "non-between" times
    leg_TPs= TPs
    
    #leg_TPs_withIntermed= unique(sort(V(all_t_graph_l[[1]])$t))  #includes times between timepoints
    #fullTP_idx= match(seq(5,60,5),leg_TPs_withIntermed)  #indices of "non-between" times
    #leg_TPs= leg_TPs_withIntermed[fullTP_idx]
    leg_colors= rgb(0.3,leg_TPs_withIntermed/60,0.6)[fullTP_idx]

    if (PC_ans== 'n'){
      # legend(0.9,-0.1, leg_TPs,pt.bg=leg_colors,  title='t/s',       #TP legend for PC_ans=='y
      #        pch=21,col="#777777", pt.cex=2, cex=.8, bty="n", ncol=1)
      legend(0.8,1.1,c('not sig.','sig.'),  title='significance',
             pch=1,col=c('black','red'), pt.cex=1, cex=.6, bty="n", ncol=1)
      if(str_detect(plotSuffix,'btwns')) {
        legend(0.8,0.6,c(1,3,4)*vertexSize_scaleFac,  title='betw. cent.',
               pch=1,col=c('black','black','black'), pt.cex=c(1,3,4)*0.375, cex=.8, bty="n", ncol=1)
        #apparently 0.375 is a 'magic multiplier for pch=1 symbols to convert from cex to plot units
        #https://stat.ethz.ch/pipermail/r-help/2003-July/035796.html
      }
    } else {
      legend(0.9,1.1, leg_TPs,pt.bg=leg_colors,  title='t/s',
             pch=21,col="#777777", pt.cex=1, cex=.6, bty="n", ncol=1)
      # legend(0.3,1.1,c('not sig.','sig.'),  title='significance',    #sigf. legend for PC_ans=='n'
      #        pch=1,col=c('black','red'), pt.cex=2, cex=.8, bty="n", ncol=1)
      if(str_detect(plotSuffix,'btwns')) {
        legend(0.5,0.6,c(1,3,4)*vertexSize_scaleFac,  title='betw. cent.',
               pch=1,col=c('black','black','black'), pt.cex=c(1,3,4)*0.375, cex=.8, bty="n", ncol=1)
      }
    }
    all_t_plot= recordPlot()
    dev.off()
    svg(paste(plotDir,today,'_netwOvrTime_',plotSuffix,'.svg',sep=''),height=plot_WH,width=plot_WH,
        pointsize=12)  #width and height in inches; cannot be set
    print(all_t_plot)
    dev.off()
    ###
    
    STRING_net_blank= STRING_net
    V(STRING_net_blank)$color= NA
    if(!is.null(V(STRING_net_blank)$commColWSigf)){ #you may not have communities due to testmode
      idx_in_all_t_graph= match(V(STRING_net_blank)$name,V(all_t_graph_l[[1]])$name )
      idx_in_all_t_graph_noNA= idx_in_all_t_graph[!is.na(idx_in_all_t_graph)]
      V(STRING_net_blank)$color[idx_in_all_t_graph_noNA]= 
        V(all_t_graph_l[[1]])$commColWSigf [idx_in_all_t_graph_noNA]
      svg(paste(plotDir,today,'_netwWithCommOnPaths_',plotSuffix,'.svg',sep=''),width=7,height=7)
      par(mar=c(0,0,0,0))
      plot(STRING_net_blank,layout=lay_STRING,vertex.size=3,arrow.width=0.1,
           vertex.label= NA) #vertex.label.cex=0.5,
      dev.off()
    }
    
    #indicate in graph which nodes are part of time-slices over time
    #the time.slice attribute contains a list for each vertex at which slice number(s) it is found
    ##https://yihui.org/animation/example/ani-record/
    ##https://www.r-bloggers.com/2015/03/animation-in-r/
    par(mar=c(0,0,1,0))
    ani.record(reset = TRUE)  # clear history before recording
    for(sliceNo in 1:(length(TPs)-1)){
      videoColor= rep('white',length(all_t_graph_l[[1]]))
      sliceNo_ind= unlist(lapply( vertex.attributes(all_t_graph_l[[1]])$time.slice, 
                                  function(x){sliceNo %in% x}))  #indicator for each vertex if it is found at given slice
      videoColor[sliceNo_ind]= 'green'
      plot(all_t_graph_l[[1]],layout=l,vertex.size=vertexSizes,arrow.width=0.1,vertex.color=videoColor,
           vertex.frame.color=V(all_t_graph_l[[1]])$frameCol, vertex.frame.width=2,
           vertex.label.dist=vertexLabDist,vertex.label.degree=vertexLabDegree,    #vertex.label.color=vertexLabCol,
           main=paste(TPs[sliceNo],' to ',TPs[sliceNo+1],' s',sep=''))
      ani.record() #record current frame
    }
    oopt= ani.options(interval = 2,loop=1)
    #saveGIF(ani.replay(),movie.name=paste(today,'NetwVideo_',plotSuffix,'.gif',sep='')) # movie.name = "testVideo.gif",outdir=paste(loc,'plots',sep='')
    saveGIF(ani.replay(),movie.name=paste(today,'_netwVideo_',plotSuffix,'.gif',sep=''),outdir=plotDir) #,movie.name = "testVideo.gif",outdir=paste(loc,'plots',sep='')
    #saveSWF(ani.replay(),movie.name=paste(today,'_netwVideo_',plotSuffix,'.swf',sep=''),outdir=plotDir) #Flash animation
    dev.off()
  }
  return(list(all_t_graph_l,l,all_t_plot))
}
#######################

make_underlNet= function(edgeDF, sigfDF, maxAbslogDF) {
  #function to produce undirected network from edge DF, significance and log2FC information
  #adds sigf, TP_node1, TP_node2, prizes
  # ARGS:
  #edgeDF needs columns 'node1' and 'node2'
  #sigfDF needs columns 'Gene' and 'firstSigTP'
  #maxAbslogDF needs columns 'Gene' and 'maxAbslog3FC'
  
  STRING_links= cbind.data.frame(edgeDF,
                                 costs= 1-edgeDF$experimentally_determined_interaction)
  STRING_nodes= tibble( name= unique(c(STRING_links[,1],STRING_links[,2])),
                        sigf= rep(FALSE,length(unique(c(STRING_links[,1],STRING_links[,2])))))
  STRING_nodes$sigf[STRING_nodes$name %in% sigfDF$Gene]= TRUE
  STRING_nodes$t= sigfDF$firstSigTP[ match(STRING_nodes$name, sigfDF$Gene) ]
  STRING_nodes$prizes= maxAbslogDF$maxAbslog3FC[ match(STRING_nodes$name, maxAbslogDF$Gene) ] #not just for significant proteins
  #from 27/3/23 based on log3
  STRING_nodes$prizes[is.na(STRING_nodes$prizes)]= 0 # prize of 0 if no FC (e.g. no phosphosite detected by Kanshin)
  STRING_nodes[1:4,]
  summary(STRING_nodes) #should not have NA in prizes
  
  STRING_links$TP_node1= sigfDF$firstSigTP[match(STRING_links$node1,sigfDF$Gene)]
  STRING_links$TP_node2= sigfDF$firstSigTP[match(STRING_links$node2,sigfDF$Gene)]
  
  STRING_net_allEdges= graph_from_data_frame(d= STRING_links, vertices= STRING_nodes,
                                             directed= FALSE)  #contains sigf....
  #V(STRING_net_allEdges)$t= sigfDF$firstSigTP[match(V(STRING_net_allEdges)$name,sigfDF$Gene)]
  par(mfrow=c(1,1),mar=c(0,0,2,0))
  plot(STRING_net_allEdges,main= 'Underl. expand. netw. with all edges')
  
  plot_nodes_edges_atTP(nodeDF= STRING_nodes, edgeDF= STRING_links)
  
  STRING_net= graph_from_data_frame(d= STRING_links, vertices= STRING_nodes,
                                    directed= FALSE)
  return(list(STRING_net, STRING_nodes,STRING_links))
}
##############################

mark_communities= function(STRING_net, plot_L) {   
  # function to mark communities in STRING_net and plot them, colored according
   # to manually defined GO terms
   # Additionally, analyze within-cluster sum-of-squares and plot
  #ARGV:
   # STRING_net (igraph object)
   # plot_L: list of recorded plots
  # RETURNS:
   # STRING_net: with added community membership
   # plot_L: updated list of recorded plots
  
  #cluster_optimal(STRING_net)  #this runs out of memory
  g_comm= cluster_edge_betweenness(STRING_net)
  
  lay_STRING1= layout_with_fr(STRING_net)
  rescale_mtx= matrix( c(rep(1/20,dim(lay_STRING1)[1]),rep(1/20,dim(lay_STRING1)[1])),ncol=2,byrow = FALSE ) #rescale manually to fit plot
  lay_STRING1= lay_STRING1*rescale_mtx
  plot(g_comm,STRING_net,layout= lay_STRING1, vertex.size=3,vertex.label= NA,
       vertex.color= nodCol,mark.groups = communities(g_comm),rescale=FALSE)
  STRING_net_comm= STRING_net
  comm1= V(STRING_net_comm)$name[membership(g_comm)==1]
  #I didn`t find a good way for GO-analysis in R, so I export and import
  comm_listList= list()
  for (l1 in seq_along(1:length(g_comm))){
    comm_listList[[l1]]= V(STRING_net_comm)$name[membership(g_comm)==l1]
  }
  #lapply(comm_listList, write, paste(loc,today,'_communities.txt',sep=''), append=TRUE, ncolumns=1000) #specify columns, otherwise truncated
  #write.csv(V(STRING_net)$name,paste(loc,today,'_protSTRINGnet.txt',sep='')) #background
  
  #comm1DF= cbind.data.frame(name=V(STRING_net_comm)$name,commInd= rep(0,length(V(STRING_net_comm)$name)))
  #comm1DF$commInd[membership(g_comm)==1]= 1
  #I select 'Function' and 'Process' at https://yeastgenome.org/goTermFinder at default settings and just choose the term that sounds most representative
  manuGOterms= c('actin\ncytoskeleton','ribosome/\ntranslation',NA, 'chromatin\norganization',
                 NA,'signal\ntransduction (I)',NA,'TOR signaling','signal\ntransduction (II)',
                 'transcription',NA,NA,'cell cycle',NA,'RNA splicing',NA,'proteasomal\ncatabolism',
                 NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,'vesicle\ntransport',NA,NA)
  
  #plot individual communities to find out which one on plot corresponds to which
  # STRING_net_comm1= subgraph(STRING_net,vids=V(STRING_net)[V(STRING_net_comm)$name %in% communities(g_comm)[[1]]])
  # lay_STRING_comm1= lay_STRING1[match(V(STRING_net_comm1)$name,V(STRING_net)$name),]
  newColsUniq= c( 'pink','orange','white','yellow',
                  'white','blue','white','violet','lightblue',
                  'darkgreen',rep('white',2),'brown','white','green','white','black',
                  rep('white',10),'grey',rep('white',2)
  )
  newCols= newColsUniq[membership(g_comm)]
  commInfo= cbind.data.frame(names= g_comm$names, commNo= as.numeric(membership(g_comm)),
                             commColor= newColsUniq[as.numeric(membership(g_comm))])  #g_comm as global variable
  
  V(STRING_net)$comm= membership(g_comm)
  V(STRING_net)$commCol= newCols   #color according to communities
  V(STRING_net)$commColWSigf= newCols   #color according to communities with significant nodes in red
  V(STRING_net)$commColWSigf[STRING_nodes$sigf]='red'
  #firstSigDF_prot$comm= V(STRING_net)$comm[match(firstSigDF_prot$Gene,V(STRING_net)$name)] #6/2/23 removed again
  #firstSigDF_prot$commCol= V(STRING_net)$commCol[match(firstSigDF_prot$Gene,V(STRING_net)$name)] #6/2/23 removed again
  
  dev.new(width=6,height=6)
  par(mar=c(0,1,0,0))
  plot(STRING_net, layout= lay_STRING, vertex.size=3,vertex.label= NA,
       vertex.color= V(STRING_net)$commColWSigf)
  text(x=c(-.25, 0.8 , -.77,   -.90,  .23, -.905 ,    .10, -.05,  .8, .71,  .66  ),
       y=c(  .8, 0.07, 0.22,   -.50, -.78, -.05,     .17, -.39, -.3, .44,  .83),
       manuGOterms[c(1,2,4,6,8,9,10,13,15,17,28)],
       col=c('pink','orange','yellow',
             'blue','violet','lightblue',
             'darkgreen','brown','green','black',
             'grey') )
  #col= newColsUniq[!is.na(manuGOterms)] ) #not using this simpler form, because I like the colors in the code for manual adjustment
  
  C1= recordPlot()
  dev.off()
  png(paste(plotDir,today,'_commun_STRING_expEvid',expEvid_cut,'.png',sep=''),width=7,height=7,unit='in',res=300)
  print(C1)
  dev.off()
  plot_L= c(plot_L,'C1'=list(C1))
  
  dev.new(width=3.5,height=3.5)
  par(mar=c(1,1,1,1))
  plot_dendrogram(g_comm,mode="phylo",cex=0.5,use.edge.length=FALSE,
                  palette=newColsUniq)  #xlim
  C2= recordPlot()
  dev.off()
  png(paste(plotDir,today,'_commDendro_STRING_expEvid',expEvid_cut,'.png',sep=''))
  print(C2)
  dev.off()
  plot_L= c(plot_L,'C2'=list(C2))
  
  # Get within-cluster sum of squares as measure for cluster compactness
  communityWSS=
    sapply(unique(g_comm$membership), function(memb) {  #for each community..
      sub_g1= induced.subgraph(STRING_net, V(STRING_net)[g_comm$membership== memb])
      n= vcount(sub_g1)
      distMtx= distances(sub_g1,v=V(sub_g1),to=V(sub_g1)) #matrix of pairwise distances
      totalDist_vec= apply(distMtx,1,sum) #shortest total distance to all other nodes 
      center_idx= which(totalDist_vec== min(totalDist_vec))[1] #center as node with shortest total distance; choose the first if multiple
      sum( distMtx[,center_idx]^2 ) /(n-1) #within-cluster-SSQ: squared distance from center normalized by n-1 (I suppose '-1', because distance from center to itself is 0)
    })
  
  dev.new(width=4,height=4)
  par(mar=c(7,4,1,.6))
  barplot(communityWSS[!is.na(manuGOterms)],names.arg=manuGOterms[!is.na(manuGOterms)],las=2,
          ylab='Within-cluster SSQ',col=newColsUniq[!is.na(manuGOterms)],cex.names = .8,lheight=.8) #
  #title(ylab='Within-cluster SSQ',mgp=c(1,1,7))   #label distance from axis
  C3= recordPlot()
  dev.off()
  png(paste(plotDir,today,'_communDensity_STRING_expEvid',expEvid_cut,'.png',sep=''),
      width=7,height=7,unit='in',res=300)
  print(C3)
  dev.off()
  plot_L= c(plot_L,'C3'=list(C3))
  
  return(list(STRING_net, plot_L))
}
##############################

addToTreePlotList= function(redrawTree, treePlot_list){
  # Function to draw a tree and add it to a list of tree plots
  # ARGV:
   # redrawTree (igraph object): tree to be drawn
   # treePlot_list: list of recorded plots (class: "recordedplot")
  # RETURNS:
   # treePlot_list: updated list of recorded plots
  
  set.seed(10) 
  l=layout_with_fr(redrawTree) 
  dev.new(width=2,height=2)
  par(mar=c(0,0,0,0))
  print( plot(redrawTree,vertex.size=24,arrow.width=0.1,vertex.label.cex=.5) ) #,vertex.label.font=2 (bold),main=paste(t1,'to',t2,sep=' ')
  text(x=-.8,y=1.1,paste(t1,'-',t2,'s',sep=' '),cex=1)
  box(lwd=1,lty=3)
  #vertex.shape='rectangle',
  #vertex.size=(strwidth(V(redrawTree)$name) + strwidth("oo")) * 500,
  #vertex.size2=strheight("I") * 2 * 100,
  #main=paste('Steiner Tree at',t1,'to',t2,'s (min. span. tree)')))
  #main=paste('Steiner Tree at',t1,'to',t2,'s; Subgraph',subgr_idx)) )
  #https://stackoverflow.com/questions/14472079/match-vertex-size-to-label-size-in-igraph
  D_curr= recordPlot() 
  dev.off()
  png(file=paste(plotDir,today,'_Tree_sG',sG_idx,'_',t1,'_',t2,'_',plotSuffix,'.png',sep=''),width=7,height=7,
      units='in',res=300)  ###
  print(D_curr)
  dev.off()   ###
  treePlot_list[[t_idx]]= D_curr
  
  return(treePlot_list)
}
##############################

plot_timeOrderedSlice= function(slice_graph, t1, t2) {
  # plot graph of time slice with t1-nodes (blue) as network on top,
  # t2-nodes (red) as network on bottom and nodes from intermediate
  # time points (green) in between
  # ARGV: 
   #slice_graph: igraph object of slice-graph (two TPs plus 
    #intermediary time) to plot
   #t1/t2: first/second time point
  
  V(sG_ST_sG_filt)$lvl= rep(3,length(V(sG_ST_sG_filt)$name))
  V(sG_ST_sG_filt)$lvl[V(sG_ST_sG_filt)$t==(t1 + t2)/2]= 2
  V(sG_ST_sG_filt)$lvl[V(sG_ST_sG_filt)$t==t2]= 1 #nodes at t2 are level 1
  
  if(length(V(sG_ST_sG_filt))>3) { #no proper plot if less than four nodes
    # generate layout
    ST_sG_filt_xy_orig= tryCatch(    #avoid error from 'layout_as_multilevel' that I don`t understand
      {
        layout_as_multilevel(sG_ST_sG_filt,type = "all", alpha = 25, beta = 45)
      }, error= function(err) {
        print(paste("Caught error:",err,"subgraph number:",subgraphNo))
        ltree= layout_as_tree(sG_ST_sG_filt)
        return(ltree)
      }, warning= function(war) {
        print(paste("Caught warning:",war))
        lmult= layout_as_multilevel(sG_ST_sG_filt,type = "all", alpha = 25, beta = 45)
        return(lmult)
      }
    )
    
    #improve separation between nodes by multiplying the distance of level-2 nodes from median in
    #y-direction 
    ST_sG_filt_xyDF1= data.frame(ST_sG_filt_xy_orig)
    ST_sG_filt_xyDF1$yDistFromMed= rep(0,dim(ST_sG_filt_xyDF1)[1])
    ST_sG_ymed= median(ST_sG_filt_xyDF1[,2])
    ST_sG_filt_xyDF1$yDistFromMed[V(sG_ST_sG_filt)$lvl==2]=
      ST_sG_filt_xyDF1[V(sG_ST_sG_filt)$lvl==2,2]-ST_sG_ymed
    ST_sG_filt_xyDF= ST_sG_filt_xyDF1
    ST_sG_filt_xyDF[V(sG_ST_sG_filt)$lvl==2,2]= rep(ST_sG_ymed,sum(V(sG_ST_sG_filt)$lvl==2))+ 
      3.3*ST_sG_filt_xyDF$yDistFromMed[V(sG_ST_sG_filt)$lvl==2]
    ST_sG_filt_xyDF= ST_sG_filt_xyDF[,-3]
    
    #print(paste(today,'_p2c',prizeToCost,'_',t1,'_to_',t2,'.svg',sep=''))

    print( 
      ggraph(sG_ST_sG_filt, "manual", x = ST_sG_filt_xyDF[, 1], y = ST_sG_filt_xyDF[, 2]) +
        geom_edge_link(aes(filter = (node1.lvl == 1 & node2.lvl == 1)),
                       edge_colour = "red",alpha = 0.5,edge_width = 0.5,
                       arrow= arrow(length=unit(3,'mm')),end_cap = circle(1, 'mm') ) + #do not end in center of node
        geom_edge_link(aes(filter = (node1.lvl != node2.lvl)),
                       alpha = 0.3,edge_width = 0.5,edge_colour = "black",
                       arrow= arrow(length=unit(3,'mm')),end_cap = circle(1, 'mm') ) + #do not end in center of node
        geom_edge_link(aes(filter = (node1.lvl == 2 & node2.lvl == 2)),
                       edge_colour = "green",edge_width = 0.5,alpha = 0.5,
                       arrow= arrow(length=unit(3,'mm')),end_cap = circle(1, 'mm')) +
        geom_edge_link(aes(filter = (node1.lvl == 3 & node2.lvl == 3)),
                       edge_colour = "blue",edge_width = 0.5,alpha = 0.5,
                       arrow= arrow(length=unit(3,'mm')),end_cap = circle(1, 'mm')) +
        geom_node_point(aes(color = as.factor(lvl)), fill = "grey25", size = 3) +
        geom_node_text(aes(label = name,color = as.factor(lvl)), 
                       repel = FALSE,cex=2.5,nudge_y=0.04 ) +   #,nudge_x=0.1
        #scale_shape_manual(values = c(21, 22,23)) +
        coord_cartesian(clip = "off", expand = TRUE) +
        ggtitle( paste(t1,' to ',t2,' min',sep=''))+
        #ggtitle( paste('t1: ',t1,', t2: ',t2,'(subgr.: ',subgraphNo,')'))+
        theme_void()+
        theme(legend.position = "none",plot.title = element_text(hjust = 0.5))   )
  }
}
#####################

trimPaths= function(compl_pathsN, compl_pathsTPs, compl_pathsType,max_St_betwSliceTerm){
  #trim paths to begin at out-term and end on in-term
  compl_pathsN2= list()
  compl_pathsTPs2= list()
  compl_pathsType2= list()
  
  for (i in seq_along(compl_pathsType)) {
    cat(' #: ',i,'\n')
    outTermIdxes= sort(which(compl_pathsType[[i]]=='out-term'))  #in order
    inTermIdxes= sort(which(compl_pathsType[[i]]=='in-term'))    #in order
    
    #only continue if at least one out-term and one in-term in path
    if(length(outTermIdxes)>0 & length(inTermIdxes)>0){
      cat('Path with out- and in-term: ',compl_pathsN[[i]],'\n')
      #get distance between first out-term and last in-term
      longestPathDist= max(inTermIdxes)-min(outTermIdxes)
      cat('i: ',i,'longest: ',longestPathDist,'\n')      ##expect: class Steiner: i=104,105,108
      
      #don`t continue if max(inTermIdxes) < min(outTermIdxes)
      if(max(inTermIdxes) < min(outTermIdxes)) {
        next
      }
      
      #do not remove paths exceeding allowed number of Steiner nodes:
      #Steiner limit is taken care of while extracting paths
      
      #record paths that pass
      compl_pathsN2[[length(compl_pathsN2)+1]]=       compl_pathsN[[i]][min(outTermIdxes) : max(inTermIdxes)]
      compl_pathsTPs2[[length(compl_pathsTPs2)+1]]=   compl_pathsTPs[[i]][min(outTermIdxes) : max(inTermIdxes)]
      compl_pathsType2[[length(compl_pathsType2)+1]]= compl_pathsType[[i]][min(outTermIdxes) : max(inTermIdxes)]
      cat('accepted: ',compl_pathsN2[[length(compl_pathsN2)]],'\n\n')
    }
  }
  return(list(compl_pathsN2,compl_pathsTPs2,compl_pathsType2))
}
##########################

plotIndivTermPaths= function(terminals, nodeType_paths_list,
                             nodeIdxDF, leg_TPs) {
  # plot paths for individually termini
  # mark the selected terminus in orange, terminal nodes in red and nodes
  # that are significant, but not termini in green
  # ARGS:
  # terminals: character vector of terminal nodes
  # nodeType_paths_list: list of paths, where each element is vector of node
  # type ("out-term", "Steiner", "in-term")
  # nodeIdxDF: DF with columns 'node', 'idx', 'sigfTP'
  # leg_TPs: numeric vector of time points to be plotted
  # OUTPUT:
  # plots to standard graphics device
  
  for (tm in terminals) {
    plot(type='n',x=0,y=0,xlim= c(leg_TPs[2],leg_TPs[length(leg_TPs)]),
         ylim=c(0,max(1,node_idx_trim$idx)),
         xlab= 't/s',ylab='',yaxt='n',
         main=paste(tm,': terminal to terminal (max. path length:',max_pathLenSt,')') )
    axis(2, at= node_idx_trim$idx[order(node_idx_trim$idx)], cex.axis=0.4,
         labels = node_idx_trim$node[order(node_idx_trim$idx)],las=2)  #,col.ticks=axCol
    for (i in 1:length(compl_pathsN3)){
      if(any(compl_pathsN3[[i]] == tm)){
        color= rep('black',length(compl_pathsN3[[i]]))
        color[ (compl_pathsType3[[i]]!='out-term' & compl_pathsType3[[i]]!='in-term') &
                 compl_pathsN3[[i]] %in% firstSigDF_prot$Gene]='green'  #not terminal, but significant (at diff. TP)
        color[which(compl_pathsType3[[i]]=='out-term' | compl_pathsType3[[i]]=='in-term')]='red'
        color[which(compl_pathsN3[[i]]==tm)]='orange'   #terminal in focus
        points(compl_pathsTPs3[[i]], node_idx_trim$idx[
          match(compl_pathsN3[[i]],node_idx_trim$node)],pch=16,col=color)
        lines(compl_pathsTPs3[[i]], node_idx_trim$idx[
          match(compl_pathsN3[[i]],node_idx_trim$node)])
      }
    }
    legend('topleft',legend=c('term. in focus','term.','sigf. but not term'),pt.bg=c('orange','red','green'),
           pch=21)
  }
}
###########################
#---------------------------- End of functions ---------------------------------#

plot_LoL= list()  #list with lists for each PC_ans, each containing plots to be arranged
for (PC_ans_idx in seq_along(PC_ans_vec)){
  
  # Continue setup and generation of file names
  PC_ans= PC_ans_vec[PC_ans_idx]
  plot_L= list()

  #node prize and edge cost function is specified elsewhere
  if(PC_ans=='y'){
    max_St_betwSliceTerm= max_St_betwSliceTerm_PCST
  } else if (PC_ans=='n'){
    max_St_betwSliceTerm= max_St_betwSliceTerm_class
  }
  
  #create a string to identify algorithm choice in plots
  if(testmodeAns== 'y'){
    plotSuffix= paste('test',testString,'_',sep='')
  } else if (testmodeAns=='n'){
    plotSuffix= paste('FC',FC_cutoff,'_',sep='')
  } else if (testmodeAns == 'simul'){
    plotSuffix= paste('simul_FC',FC_cutoff,'_',sep='')
  }
  if (PC_ans== 'n'){
    plotSuffix= paste(plotSuffix,'classSteiner_maxTree',max_St_betwSliceTerm_class,sep='')
  } else if (PC_ans== 'y') {
    plotSuffix= paste(plotSuffix,'PCST_p2c',prizeToCost,'_maxTree',max_St_betwSliceTerm_PCST,sep='')
  }
  plotSuffix= paste(plotSuffix,'_maxPath',max_pathLenSt,sep='')
  
  
  if (testmodeAns== 'n') {
    
    # A, Process phosphoproteomics data --------------------------------------------
    colIdxToSelect= c(which(colnames(data_orig) == 'Gene'),
                      which(colnames(data_orig) == 'ORF'),
                      which(colnames(data_orig) == 'pSite'),
                      which(str_detect(colnames(data_orig), '^log2_')),
                      which(str_detect(colnames(data_orig), '^log10int_'))
    )
    if(any(is.na(colIdxToSelect))){
      stop("Error: Expected columns missing in phosphoproteomics data file!")
    }
    data2= cbind.data.frame(data_orig[, colIdxToSelect])
    data2[1:2,]
  
    data2
    str(data2)
    summary(data2)
  
    #some genes may have since been named
    IDmap= read.csv(paste(loc,'20150727_UniProt_yeast_IDmapping.csv',sep=''), stringsAsFactors = FALSE)
    IDmap$Gene.Name= stringr::str_to_title(IDmap$Gene.Name)
    IDmap[1:3,]
    data2$Gene[data2$Gene %in% IDmap$ORF.Name]=
      IDmap$ORF.Name[match(data2$Gene,IDmap$ORF.Name)[data2$Gene %in% IDmap$ORF.Name]]
    
    #I assume data2 lists each pSite only once, even though it is a ppep table; here is the proof
    sum(duplicated(paste(data2$ORF,data2$pSite,sep='_')))   # 0 -> i.e. no duplicates in sites
  
    #You need to calculate max(abs(logFC)) already for full data, because you need it for prizes also for Steiner nodes
    d2_log2colNos_not0= which(str_detect(colnames(data2),'^log2_') &  ##########
                               colnames(data2)!='log2_T00')
    data2$maxAbslog2FC=
      apply(data2[,d2_log2colNos_not0],1,function(x){max(abs(x),na.rm=TRUE)})
    #summarize to proteins so I can later use it to assign node prizes
    data2_prot= data.frame( data2 %>%
           group_by(ORF,Gene) %>%
           summarize(maxAbslog2FC= max(maxAbslog2FC,na.rm=TRUE)) )
    data2_prot$maxAbslog3FC= log( 2^data2_prot$maxAbslog2FC,3)
    dim(data2)
    dim(data2_prot)
    summary(data2_prot)
    
    # B, Extract significant proteins ------------------------------------------
    ################## OPTION 1 ####################
    #### DEFINE SIGNIFICANT SITES YOURSELF #########
    #determine significant values (at least ...x change)
    data_signif= data2[, d2_log2colNos_not0]> log(FC_cutoff,2) | 
                 data2[, d2_log2colNos_not0] < log(1/FC_cutoff,2)   #don`t consider changes at t=0
    colnames(data_signif)= str_replace(colnames(data_signif),"log2_","sigf_")
    sigfIndDF= cbind.data.frame(Gene= data2$Gene,ORF=data2$ORF,pSite= data2$pSite,data_signif)
  
    #find for each site the first timepoint at which it is significant
    firstSigfIdx= apply(as.matrix(sigfIndDF[,4:dim(sigfIndDF)[2]]),1,
                        function(x){min(99,which(x))})
    TPs= as.numeric(str_replace( colnames(sigfIndDF)[4:dim(sigfIndDF)[2]],'sigf_T',''))
    firstSigTP= TPs[firstSigfIdx]
    log2_colNos= which(str_detect(colnames(data2),'log2_T[0-9]+$'))
    sigfIndDF= cbind.data.frame(sigfIndDF,firstSigTP,abs(data2[,log2_colNos]))
    colnames(sigfIndDF)= str_replace(colnames(sigfIndDF),'^log2_','abs_log2_')
    firstSigDF_sites= sigfIndDF[!is.na(sigfIndDF$firstSigTP),]
    firstSigDF_sites= firstSigDF_sites[firstSigDF_sites$firstSigTP!=0,]
    dim(firstSigDF_sites) #number of sigf. sites (before collapsing on proteins)
    firstSigDF_sites[,c("Gene","ORF","firstSigTP")]
  
    abs_log2_colNos= which(str_detect(colnames(firstSigDF_sites),'abs_log2_T[0-9]+$'))
    abs_log2_0subs= firstSigDF_sites[,abs_log2_colNos]
    abs_log2_0subs[is.na(abs_log2_0subs)]= 0
    firstSigDF_sites_0subs= cbind.data.frame( firstSigDF_sites[,1:(min(abs_log2_colNos)-1)],abs_log2_0subs)
  
    firstSigDF_prot0= data.frame( firstSigDF_sites_0subs %>%
      # note that this approach may use log2-ratios from different sites across samples
          group_by(Gene,ORF) %>%
          summarise(
            firstSigTP= min(firstSigTP,na.rm=TRUE),
            across(starts_with("abs_log2_T"), \(x) max(x, na.rm = TRUE)) 
          )
    )
    
    #~~~~ checkpoint ~~~#
    chck = 0
    while(chck != 1){
      if( sum(firstSigDF_prot0$Gene == 'Msn2') == 1){
        chck = 1
      } else {
        chck= readline("Incorrect count of Msn2 among signif. proteins.
                        Press '1' to acknowledge?")
      }
    }
    #~~~~~~~~~~~~~~~~~~~#
    
    firstSigDF_protA= firstSigDF_prot0[!is.na(firstSigDF_prot0$firstSigTP),]
    firstSigDF_protA= firstSigDF_protA[firstSigDF_protA$firstSigTP != 0,] #don`t use sites changing before salt addition
    dim(firstSigDF_protA) #number of significant proteins
    firstSigDF_protA[1:4,]
    dim(firstSigDF_protA)[1]   ## Number of signif. diff. phos. proteins
    if(simulAns == 'y'){
    write.csv(firstSigDF_protA,
              paste(TCfile_loc, TCfile_base,
                    '_allSignif', FC_cutoff, 'Xcutoff.csv', sep=''),
              row.names = FALSE)
    } else if (simulAns == 'n'){
      write.csv(firstSigDF_protA,
      paste(loc,today,'_allSignif',FC_cutoff,'Xcutoff.csv',sep=''),
      row.names = FALSE)
    }
    #######################################
    #######################################
    
    ############ OPTION 2 #################
    #### TAKE THE SIGNIFICANCE AS DETERMINED BY KANSHIN 2015 ######
    #the problem with this approach is that Kanshin, 2015 does not specify TP of significance
    # 
    # #merge Tab1 (all quantifications) and Tab2 based on sites
    # Kanshin_ppeps= read.csv(paste(loc,'20230227_sigfSites_Kanshin_fromS2.csv',sep=''))
    # Kanshin_ppeps$Gene= stringr::str_to_title(Kanshin_ppeps$Gene)
    # dim(Kanshin_ppeps)
    # 
    # #Tab2 contains 'T723/S726', while Tab1 only lists individ sites
    # Kanshin_sites= data.frame()
    # for (i in 1:dim(Kanshin_ppeps)[1]) { 
    #   currSites= str_split(Kanshin_ppeps$pSite[i],'/')[[1]]
    #   currPpep_DF= data.frame(matrix(rep(as.character(Kanshin_ppeps[i,]),times=length(currSites)), 
    #                 nrow=length(currSites),byrow=TRUE))
    #   currPpep_DF[,which(str_detect(colnames(Kanshin_ppeps),'pSite'))]= currSites
    #   Kanshin_sites= rbind.data.frame(Kanshin_sites,currPpep_DF)
    # }
    # colnames(Kanshin_sites)= colnames(Kanshin_ppeps)
    # dim(Kanshin_sites)
    # 
    # dim(data2)
    # data2_KanshSigf= merge.data.frame(x=Kanshin_sites,y=data2,all=FALSE,by.x=c('ORF','pSite'),by.y=c('ORF','pSite'))
    # dim(data2_KanshSigf)  #should have same dim(Kanshin_sites)[1], as exactly one entry should be found in data2
    # head(data2_KanshSigf)
    # 
    # #find maximum absolute FC
    # d2K_log2colNos= which(str_detect(colnames(data2_KanshSigf),'^log2_'))  #includes t=0
    # 
    # #find TPs (among log2 columns, where half maximum is reached)
    # maxAbsMtx= matrix( rep(data2_KanshSigf$maxAbslog2FC,times= length(d2K_log2colNos)),
    #                   ncol=length(d2K_log2colNos),byrow=FALSE)
    # halfMaxInd= abs(data2_KanshSigf[,d2K_log2colNos]) > maxAbsMtx/2  #indicates if value > half max
    # colnames(halfMaxInd)= paste('halfMaxInd',colnames(halfMaxInd),sep='_')
    # data2_KanshSigf2= cbind.data.frame(data2_KanshSigf[,1:max(d2K_log2colNos)],
    #                           maxAbslog2FC= data2_KanshSigf$maxAbslog2FC,halfMaxInd)
    # 
    # #get first sigf TP
    # halfMaxInd_colNos= which(str_detect(colnames(data2_KanshSigf2),'halfMaxInd'))
    # TPs_with0= c(0,TPs)
    # firstHalfMaxInd_idx= apply(data2_KanshSigf2[,halfMaxInd_colNos],1,
    #                       function(x) min(which(x),na.rm=TRUE)) #indices of first TP where half max. is reached
    # data2_KanshSigf2$firstHalfMax= TPs_with0[firstHalfMaxInd_idx]
    # 
    # #summarize to proteins
    # Kanshin_prot1= data.frame( data2_KanshSigf2 %>%
    #   group_by(ORF,Gene.x) %>%
    #   summarize(firstHalfMax= min(firstHalfMax,na.rm=TRUE),maxAbslog2FC= max(maxAbslog2FC,na.rm=TRUE)) )
    #      #it may not be the site that reaches the half max. first, which has the highest max(abs(log2FC))
    # dim(Kanshin_prot1)
    # Kanshin_prot= Kanshin_prot1[Kanshin_prot1$firstHalfMax!=0,]  #remove all PROTEINS with any SITE changing at 0
    # dim(Kanshin_prot)  # 298, compared to 332 reported by Kanshin, 2015
    # 
    # #which proteins significant in Kanshin, 2015 were not found
    # length(unique(Kanshin_sites$ORF))
    # sum(unique(Kanshin_sites$ORF) %in% data2$ORF)  #326 (some sites of Tab S2 are not in Tab S1)
    # unique(Kanshin_sites$ORF)[unique(Kanshin_sites$ORF) %in% data2$ORF == FALSE]  # 6 of the ORFs are not in Tab. S1
    # #proteins for which merge did not work
    # mergeORFs= unique(Kanshin_sites$ORF[
    #   (paste(Kanshin_sites$ORF,Kanshin_sites$pSite,sep='_') %in% paste(data2$ORF,data2$pSite,sep='_'))] )
    # nonMergeORFs= unique(Kanshin_sites$ORF[ Kanshin_sites$ORF %in% mergeORFs == FALSE])
    # cbind.data.frame(data2$ORF[match(nonMergeORFs,data2$ORF)],data2$pSite[match(nonMergeORFs,data2$ORF)])
    # #Example
    # cbind.data.frame(ORF=data2$ORF[data2$ORF=='YOR014W'],pSite=data2$pSite[data2$ORF=='YOR014W'])
    # cbind.data.frame(ORF=Kanshin_sites$ORF[Kanshin_sites$ORF=='YOR014W'],pSite=Kanshin_sites$pSite[Kanshin_sites$ORF=='YOR014W'])
    #   #some sites do not match between tab S1 and S2; e.g. YOR014W is T257 in S2 and S263/264 in S1
    # sum(unique(Kanshin_sites$ORF) %in% data2_KanshSigf$ORF)
    # sum(unique(Kanshin_sites$ORF) %in% data2_KanshSigf2$ORF)
    # sum(unique(Kanshin_sites$ORF) %in% Kanshin_prot1$ORF)  #318
    # sum(mergeORFs %in% data2_KanshSigf$ORF)
    # sum(mergeORFs %in% data2_KanshSigf2$ORF)
    # sum(mergeORFs %in% Kanshin_prot1$ORF)
    # sum(mergeORFs %in% Kanshin_prot$ORF) #further loss because of proteins with sites that are already changed at t=0
    # 
    # #write output
    # #write.csv(Kanshin_prot,paste(loc,'20230228_sigfProt_Kanshin_fromS2.csv',sep=''),row.names = FALSE)
    # dim(Kanshin_prot)  #298
    # 
    # firstSigDF_protA= Kanshin_prot #for consistency with OPTION 1, I call it 'firstSigDF_protA'
    # colnames(firstSigDF_protA)[c(2,3)]= c('Gene','firstSigTP')
    
    ######################################
    ######################################
    
    ############
    #Generate STRING networks w/o expansion and with 1 shell expansion (400 nodes)
     #to annotate proteins from phosphoproteomics correctly
     #use exp.evid., DB, textmining; minimum total score of 'expEvid_cut'
  
    #STRING was searched with ORFs and may have replaced some gene names
    if (simulAns == 'n'){
    STRING_1shell= read.table(paste(loc,'20230313_allSignif_STRING1exp_',FC_cutoff,
                                    'Xcutoff_TexExDB',expEvid_cut,'.tsv',sep=''),sep='\t',
                              header=TRUE, stringsAsFactors = FALSE)
    # STRING_1shell= read.table(paste(loc,'20230228_STRING_Kanshin_fromS2_1shell_texExDB_0.7.tsv',sep=''),sep='\t',
    #                           header=TRUE, stringsAsFactors = FALSE)
    } else if (simulAns == 'y'){
      STRING_1shell= read.table(
        paste(TCfile_loc, TCfile_base, '_STRING1exp_', FC_cutoff,
              'Xcutoff_TexExDB',expEvid_cut,'.tsv',sep=''),
        sep='\t', header=TRUE, stringsAsFactors = FALSE)
    }
       #use the maximally expanded STRING network I am going to use in this analysis
    
    #convert from gene names to protein names
    STRING_1shell$node1= sapply(STRING_1shell$node1,str_to_title)
    STRING_1shell$node2= sapply(STRING_1shell$node2,str_to_title)
    match(STRING_1shell$node1,firstSigDF_protA$Gene)
    STRING_1shell$sigf_node1= STRING_1shell$node1 %in% firstSigDF_protA$Gene
    STRING_1shell$sigf_node2= STRING_1shell$node2 %in% firstSigDF_protA$Gene
    
    head(STRING_1shell)
    nodes_1shell= cbind.data.frame(
      nodes_names= c(STRING_1shell$node1,STRING_1shell$node2),
      nodes_ORFs= c(STRING_1shell$node1_string_id,STRING_1shell$node2_string_id) #e.g. 4932.YCR088W
    )
    nodes_1shell= nodes_1shell[!duplicated(nodes_1shell$nodes_names),]
    nodes_1shell$ORFs= sapply(nodes_1shell$nodes_ORFs,function(x) {str_split(x,'\\.')[[1]][2]})
    head(nodes_1shell)
    firstSigDF_protA$Gene= nodes_1shell$nodes_names[match(firstSigDF_protA$ORF,nodes_1shell$ORFs)]
    dim(firstSigDF_protA)
    
    #not all nodes may have been found by STRING
    firstSigDF_prot= firstSigDF_protA[ !is.na(firstSigDF_protA$Gene), ]
    dim(firstSigDF_prot)  ### Number of proteins mapped to STRING
    
    ###########################################
    ######## ANALYSIS of the unexpanded network
    load_intConfcut= 0.7 #cutoff for file to load
    STRunexpand_allEvid= read.table(
      paste(TCfile_loc, TCfile_base,
            '_STRING1exp_', FC_cutoff,'Xcutoff_TexExDB',expEvid_cut,'.tsv',sep=''),
      sep='\t', header=TRUE, stringsAsFactors = FALSE)
    
    #convert from gene names to protein names
    STRunexpand_allEvid$node1= sapply(STRunexpand_allEvid$node1,str_to_title)
    STRunexpand_allEvid$node2= sapply(STRunexpand_allEvid$node2,str_to_title)
    
    dim(STRunexpand_allEvid)[1] #number of edges
    dim(STRunexpand_allEvid[STRunexpand_allEvid$experimentally_determined_interaction>0.7,])[1]
    dim(STRunexpand_allEvid[STRunexpand_allEvid$experimentally_determined_interaction>0.4,])[1]
    dim(STRunexpand_allEvid[STRunexpand_allEvid$experimentally_determined_interaction>0.15,])[1]
    
    #filter by 'experimentally_determined_interaction' score
    STRunexpand= STRunexpand_allEvid[
      STRunexpand_allEvid$experimentally_determined_interaction>expEvid_cut,]
    dim(STRunexpand)
    
    #add information which proteins are significant  (all should be, since unexpanded)
    STRunexpand$sigf_node1= STRunexpand$node1 %in% firstSigDF_prot$Gene
    STRunexpand$sigf_node2= STRunexpand$node2 %in% firstSigDF_prot$Gene
    head(STRunexpand)
    
    # convert to graph (mainly for compatibility with analysis of expanded network)
    unexNet_out = make_underlNet(edgeDF = STRunexpand, sigfDF = firstSigDF_prot,
                                maxAbslogDF = data2_prot)
    STRunexpand_net= unexNet_out[[1]]
    STRINGunexpand_nodes= unexNet_out[[2]]
    STRunexpand_new= unexNet_out[[3]]
    vertex.attributes(STRunexpand_net)
    
    plot_nodes_edges_atTP(nodeDF= STRINGunexpand_nodes, edgeDF= STRunexpand_new)
    
    #plot
    g_unexpand= graph_from_data_frame(d=STRunexpand_new[,1:2],directed = FALSE,
                                      vertices= unique(c(STRunexpand_new[,1],STRunexpand_new[,2])))
    V(g_unexpand)$t= firstSigDF_prot$firstSigTP[match(V(g_unexpand)$name,firstSigDF_prot$Gene)]
    V(g_unexpand)$label= paste(V(g_unexpand)$name,V(g_unexpand)$t,sep='_')
    length(V(g_unexpand))  #nodes
  
    l= layout_with_fr(g_unexpand)
    dev.new(width=3.5,height=3.5,units='in',noRStudioGD = TRUE)
    par(mar=c(1,1,1,1),mfrow=c(1,1))
    plot(g_unexpand,vertex.label=V(g_unexpand)$label,vertex.size=0.5,edge.size=5)
        # main='Underl. netw. with all edges')  #layout= l*0.2,,rescale=FALSE
    A1= recordPlot()
    dev.off()
    png(paste(plotDir,today,'_unexpand_FCthr',FC_cutoff,'.png',sep=''))
    print(A1)
    dev.off()
    plot_L= c(plot_L,'A1'=list(A1))
    
    #make network directed by adding reverse of each edge
    #from now on use '_from/to' instead of 'node1/2'
    STRunexpand_directed_f= STRunexpand_new
    STRunexpand_directed_f$type_from= rep('out-term',dim(STRunexpand_directed_f)[1])
    STRunexpand_directed_f$type_to= rep('in-term',dim(STRunexpand_directed_f)[1])
    colnames(STRunexpand_directed_f)[colnames(STRunexpand_directed_f)=='node1']= 'node_from'
    colnames(STRunexpand_directed_f)[colnames(STRunexpand_directed_f)=='node2']= 'node_to'
    colnames(STRunexpand_directed_f)[colnames(STRunexpand_directed_f)=='sigf_node1']= 'sigf_from'
    colnames(STRunexpand_directed_f)[colnames(STRunexpand_directed_f)=='sigf_node2']= 'sigf_to'
    colnames(STRunexpand_directed_f)[colnames(STRunexpand_directed_f)=='TP_node1']= 'TP_from'
    colnames(STRunexpand_directed_f)[colnames(STRunexpand_directed_f)=='TP_node2']= 'TP_to'
    STRunexpand_directed_copy= STRunexpand_directed_f   #copy that remains unaltered
    STRunexpand_directed_r= STRunexpand_directed_f
    STRunexpand_directed_r$node_from= STRunexpand_directed_copy$node_to
    STRunexpand_directed_r$node_to= STRunexpand_directed_copy$node_from
    STRunexpand_directed_r$sigf_from= STRunexpand_directed_copy$sigf_to
    STRunexpand_directed_r$sigf_to= STRunexpand_directed_copy$sigf_from
    STRunexpand_directed_r$TP_from= STRunexpand_directed_copy$TP_to
    STRunexpand_directed_r$TP_to= STRunexpand_directed_copy$TP_from
    STRunexpand_directed_r$type_from= STRunexpand_directed_copy$type_to
    STRunexpand_directed_r$type_to= STRunexpand_directed_copy$type_from
    STRunexpand_directed= rbind.data.frame(STRunexpand_directed_f,STRunexpand_directed_r)
    
    #reduce the DF to edges between neighboring TPs or within the same TP
    STRunexpand_directed_tempOrd= STRunexpand_directed[
      (STRunexpand_directed$TP_to - STRunexpand_directed$TP_from) ==  TP_interval |
        (STRunexpand_directed$TP_to - STRunexpand_directed$TP_from) == 0,]
    
    #extract paths with connections only between adjacent TPs
    paths_unexpand= extract_paths(currEdges_exp=STRunexpand_directed_tempOrd,max_St_betwSliceTerm=max_St_betwSliceTerm,
                                  max_pathLenSt=max_pathLenSt,end_prot=NA)
    (unex_compl_pathsN= paths_unexpand[[1]])
    unex_compl_pathsTPs= paths_unexpand[[2]]
    unex_compl_pathsTypes= paths_unexpand[[3]]
    length(unex_compl_pathsN)
    #proteins on paths may not be in right temporal order; 
    #node order of paths may be duplicated (with different temporal order);
    #path may be sub-path of other
    
    #DO NOT trim paths here, as a node is out-term and in-term at the same TP, from which it is randomly selected
     #you may therefore loose nodes that are selected as out-term instead of in-term at the end of the path
      #but require that the last TP of path is later than first TP (not the same)
    unex_compl_pathsN2= list()
    unex_compl_pathsTPs2= list()
    unex_compl_pathsTypes2= list()
    for (i in seq_along(unex_compl_pathsN)){
      if(unex_compl_pathsTPs[[i]][length(unex_compl_pathsTPs[[i]])] > unex_compl_pathsTPs[[i]][1]) {
        unex_compl_pathsN2[[length(unex_compl_pathsN2)+1]]= unex_compl_pathsN[[i]]
        unex_compl_pathsTPs2[[length(unex_compl_pathsTPs2)+1]]= unex_compl_pathsTPs[[i]]
        unex_compl_pathsTypes2[[length(unex_compl_pathsTypes2)+1]]= unex_compl_pathsTypes[[i]]
      }
    }
    
    #collapse overlapping paths
    reducOut= reduce_paths(paths= unex_compl_pathsN2,TPs= unex_compl_pathsTPs2,Types= unex_compl_pathsTypes2)
    (unex_compl_pathsN3= reducOut[[1]])
    unex_compl_pathsTPs3= reducOut[[2]]
    unex_compl_pathsTypes3= reducOut[[3]]
    
    plotSuffixUnexp= paste(plotSuffix,'_unexp',sep='')
    outChronUnexp= chronoPaths( nodePaths_list=unex_compl_pathsN3, TPsPaths_list= unex_compl_pathsTPs3, 
                 typesPaths_list= unex_compl_pathsTypes3, sigDF=firstSigDF_prot, 
                 maxLen= Inf,plotmode='combPlot',plotSuffixUnexp)
    A2= outChronUnexp[[2]]
    plot_L= c(plot_L,'A2'=list(A2))
    
    #### evaluate results...
    
    
     #find how many of the interactions in each pathway are a protein-kinase followed by changing p-site
      #kinase data; https://thebiogrid.org/project/2/s-cerevisiae-kinome.html
    kinaseDF= read.table(paste(loc,'BIOGRID-PROJECT-kinome_project_sc-LATEST/BIOGRID-PROJECT-kinome_project_sc-GENES-4.4.213.projectindex.txt',sep=''),
                        sep='\t',header=TRUE,stringsAsFactors=FALSE)
    kinaseDF$OFFICIAL.SYMBOL= stringr::str_to_title(kinaseDF$OFFICIAL.SYMBOL)
    kinases= kinaseDF$OFFICIAL.SYMBOL[kinaseDF$CATEGORY.VALUES=="Kinase"]
    
    #kinase - terminal connections on paths
     #first, put edges on path into DF
    EdgesOnPaths_unexp= data.frame(matrix(NA,nrow=0,ncol=4))
    length(unex_compl_pathsN3)
    for(i in seq_along(unex_compl_pathsN3)) {
      for (nIdx in 1:(length(unex_compl_pathsN3[[i]])-1)){
        #there is edge between each adjacent pair
        EdgesOnPaths_unexp= rbind.data.frame(EdgesOnPaths_unexp,
                                       c(unex_compl_pathsN3[[i]][c(nIdx,nIdx+1)],
                                         unex_compl_pathsTypes3[[i]][c(nIdx,nIdx+1)]))
      }
    }
    colnames(EdgesOnPaths_unexp)= c('node1','node2','type_from','type_to')
    
    noEdgesToTerm_inPaths_unexp= 0  #edges to terminal
    noEdgesKinToTerm_inPaths_unexp= 0 #edges from kinase to terminal
    kin_term_paths_unexp= data.frame(matrix(NA,nrow=0,ncol=7))  #dataframe listing paths from kinase to terminal
    for (i in seq_along(EdgesOnPaths_unexp$node1)) {
      if(EdgesOnPaths_unexp$type_to[i]== 'in-term') {
        noEdgesToTerm_inPaths_unexp= noEdgesToTerm_inPaths_unexp+1
        if(EdgesOnPaths_unexp$node1[i] %in% kinases) {
          noEdgesKinToTerm_inPaths_unexp= noEdgesKinToTerm_inPaths_unexp+1
          kin_term_paths_unexp= rbind.data.frame(kin_term_paths_unexp,EdgesOnPaths_unexp[i,])
        }
      }
    }
    noEdgesToTerm_inPaths_unexp
    noEdgesKinToTerm_inPaths_unexp
    kin_term_paths_unexp
    
    #..in underlying network: edges between kinase and significant node
    #STRING network is undirected, therefore you cannot tell if edge is coming from or going to kinase
    STRunexpand[1:3,]
    noEdgesNeighTerm_underl_unexp= 0  #edges neighboring terminal
    noEdgesKinNeighTerm_underl_unexp= 0 #edges between kinase and terminal
    for (i in seq_along(STRunexpand$node1)) {
      #"forward direction"
      if(STRunexpand$sigf_node2[i]== TRUE) {
        noEdgesNeighTerm_underl_unexp= noEdgesNeighTerm_underl_unexp+1
        if(STRunexpand$node1[i] %in% kinases) {
          noEdgesKinNeighTerm_underl_unexp= noEdgesKinNeighTerm_underl_unexp+1
          #print(paste('node1:',STRunexpand$node1[i],'node2:',STRunexpand$node2[i],sep=' '))
          next #you don`t want to count a kinase-kinase edge twice 
        }
      }
      #"backward direction"
      if(STRunexpand$sigf_node1[i]== TRUE) {
        noEdgesNeighTerm_underl_unexp= noEdgesNeighTerm_underl_unexp+1
        if(STRunexpand$node2[i] %in% kinases) {
          noEdgesKinNeighTerm_underl_unexp= noEdgesKinNeighTerm_underl_unexp+1
          #print(paste('node2:',STRunexpand$node2[i],'node1:',STRunexpand$node1[i],sep=' '))
        }
      }
    }
    noEdgesNeighTerm_underl_unexp
    noEdgesKinNeighTerm_underl_unexp
  
    
  ### Expand network ###
    #C: Get PPI network --------------------------------------------------------
    #generate STRING network from all changing proteins
    #write.csv(firstSigDF_prot,paste(loc,today,'_allSignif',FC_cutoff,'Xcutoff.csv',sep=''),row.names = FALSE)
     #expand one level: up to 500 proteins: textmining, experiments, databases; high confidence
      #(a few nodes remain unconnected); minimum required interaction score: 0.7
     #download one-way edges .tsv
    STRING_allSigf_allEvid= STRING_1shell #loaded previously and GENE names renamed to protein names
    #STRING_allSigf_allEvid= read.table(paste(loc,'20221222_allSignif_STRING1exp_',FC_cutoff,'Xcutoff_0.9.tsv',sep=''),
     #                                  header=TRUE,stringsAsFactors = FALSE)
    # STRING_allSigf_allEvid= read.table(paste(loc,'20221217_allSignif_STRING1exp_',FC_cutoff,'Xcutoff.tsv',sep=''),
    #                                    header=TRUE,stringsAsFactors = FALSE)
    #STRING_allSigf_allEvid$node1= stringr::str_to_title(STRING_allSigf_allEvid$node1)   #if read in new
    #STRING_allSigf_allEvid$node2= stringr::str_to_title(STRING_allSigf_allEvid$node2)   #if read in new
    dim(STRING_allSigf_allEvid) #allEvid network size: edges
    STRING_allSigf_allEvid[1:3,]
    length(unique(c(STRING_allSigf_allEvid$node1,STRING_allSigf_allEvid$node2))) #allEvid network size: nodes
    
    # D: Filter by 1, Interaction confidence -----------------------------------
    #build network based on 'experimentally_determined_interaction' (physical and
     #genetic interactions; also high-throughput); test using 'database_annotated'
     #(e.g. KEGG, BioCyc) and 'automated_textmining' (titles+abstracts)
    #filter to 'experimentally_determined_interaction' score !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    STRING_allSigf= STRING_allSigf_allEvid[
      STRING_allSigf_allEvid$experimentally_determined_interaction> expEvid_cut,]   #>0.7: 6092
    dim(STRING_allSigf)
    
    # E, Map phosphoproteomics to PPI network -----------------------------------
    makeNetw_out= make_underlNet(edgeDF=STRING_allSigf, sigfDF=firstSigDF_prot,
                               maxAbslogDF=data2_prot)
    STRING_net= makeNetw_out[[1]]
    STRING_nodes= makeNetw_out[[2]]
    STRING_allSigf= makeNetw_out[[3]]
    vertex.attributes(STRING_net)
    
    # ---------- Analysis of underlying network -------------------------------
    length(V(STRING_net))  #network size (nodes)
    length(E(STRING_net))  #network size (edges)
    par(mar=c(0,0,0,0))
    set.seed(1)
    lay_STRING= layout_with_fr(STRING_net)
    nodCol= rep(NA,dim(STRING_nodes)[1])
    nodCol[STRING_nodes$sigf]='red'
    V(STRING_net)$color= nodCol
    
    #dev.new(width=7,height=7,unit='in')
    par(mar=c(1,1,1,1))
    plot(STRING_net, layout= lay_STRING, vertex.size=3,vertex.label= NA)  #main='Underl. network',axes=TRUE
    B1= recordPlot()
    png(paste(plotDir,today,'_STRING_FCthr',FC_cutoff,'_expEvid',expEvid_cut,'.png',sep=''))
    print(B1)
    dev.off()
    plot_L= c(plot_L,'B1'=list(B1))
    
    if(simulAns == 'n'){
      #divide graph into communities: plot; plot within-cluster-SSQ
      out = mark_communities(STRING_net, plot_L)
    STRING_net= out[[1]]
    plot_L= out[[2]]
    }    
    
    #Summarize underlying network
    par(mfrow=c(2,1),mar=c(4,5,1,1))
    #get degree distribution
    degr_underl= degree(graph=STRING_net,mode='all')  #network is undirected
    hist(degr_underl,main=NULL,breaks=70,xlab='Degree')
    median(degr_underl)
    mean(degr_underl)
    
    #get betweenness values of underlying network (to later compare to values of inferred networks)
    findSubgraphs(STRING_net) #how many subgraphs with how many nodes?
    #maxPaths_underl= (length(V(STRING_net))-30)*(length(V(STRING_net))-30)-39  #maximum number of paths in underlying network; of 20/01/23
    btwns_underl= betweenness(graph=STRING_net,directed=TRUE)
    par()
    hist(btwns_underl,main=NULL,breaks=70,xlab='Betweenness Centrality')
    median(btwns_underl)
    mean(btwns_underl)
    
    B2= recordPlot() #record the plot first
    png(paste(plotDir,today,'_DegrBtwnDists_STRING_expEvid',expEvid_cut,'.png',sep=''))
    print(B2)
    dev.off()
    plot_L= c(plot_L,'B2'=list(B2))
    
    
    # (F) For Classical Steiner approach: Find subgraphs of STRING network --------
    if( PC_ans== 'n') {
      subgraph_list1= findSubgraphs(G=STRING_net) 
       #limit to subgraphs of certain size
      subgraph_list= subgraph_list1[unlist(lapply(subgraph_list1,FUN=function(x){length(x)>10}))]
    } else {  #no division into subgraphs necessary if PRICE COLLECTING tree
      subgraph_list= list()
      subgraph_list[[1]]= seq_along(V(STRING_net))
    }
    #continue only with subgraphs of at least 5 nodes - not active
    #subgraph_lengths= unlist(lapply(subgraph_list,FUN=length))
    #subgraph_list= subgraph_list[[which(subgraph_lengths>4)]]
    
    par(mar=c(0,0,2,0),mfrow=c(1,1))
    for(i in seq_along(subgraph_list)) {
      subgrP= subgraph(STRING_net,subgraph_list[[i]])  #subset graph with node indices
      print(plot(subgrP,vertex.size=3,vertex.label= NA,main= paste('Subgraph',i)))
    }
    length(V(subgraph(STRING_net,subgraph_list[[1]])))
    
    # F, Steiner Trees on time slices of subgraphs ----------------------------
      #first generate Steiner tree (minimum spanning tree) for each pair two adjacent
       #time points (do this individually for each subgraph)
    TPs
    (nodesNotInSTRING= firstSigDF_prot[is.na(match(firstSigDF_prot$Gene,STRING_nodes$name)),]) 
                 #diff. phos proteins not found by STRING
    par(mar=c(0,0,1.5,0))
    allUnalloc= c()  
    all_t_graph_l= list()
    Edges_l= list() #edges in steiner tree: list for each subgraph
    runtimeSteinerMtx= matrix(NA,nrow=length(subgraph_list),ncol=length(TPs))  #in min; rows by sG_idx, cols by t_idx
    runtimeLoopMtx= matrix(NA,nrow=length(subgraph_list),ncol=length(TPs))  #in min; rows by sG_idx, cols by t_idx
    sigfCountInOut_mtxL= list()
    treePlot_listList= list()
    for (sG_idx in seq_along(subgraph_list)) {
      print(paste('sG_idx:',sG_idx))
      all_t_graph_l[[sG_idx]]= make_empty_graph(n=0,directed=TRUE)
      Edges_l[[sG_idx]]= data.frame(matrix(NA,nrow=0,ncol=8)) #!!
      colnames(Edges_l[[sG_idx]])= c('1','2','TP_from','TP_to','sigf_from','sigf_to','type_from','type_to')
      start_systimeLoop= Sys.time()
      sigfCountInOut_mtxL[[sG_idx]]= matrix(NA,nrow=(length(TPs)-1),ncol=2) #list significant nodes in input and in tree at each TP
      treePlot_list= list()
      for (t_idx in seq.int(from=1,to= length(TPs)-1) ) {  #T=0 is not part of TPs; lenght-1, because
                                             #I here take the firs timepoint of the pair
        t1= TPs[t_idx]
        t2= TPs[t_idx+1]
        print(paste('t1:',t1,',t2:',t2,sep=' '))
        nodesAtt1= firstSigDF_prot$Gene[firstSigDF_prot$firstSigTP==t1]
        nodesAtt2= firstSigDF_prot$Gene[firstSigDF_prot$firstSigTP==t2]
        nodesAtt_DF= rbind.data.frame(
         cbind.data.frame(nodes=nodesAtt1,t=rep(t1,length(nodesAtt1))),
         cbind.data.frame(nodes=nodesAtt2,t=rep(t2,length(nodesAtt2))) )
        nodesAtt= nodesAtt_DF$nodes
        curr_nodes_DF= nodesAtt_DF[nodesAtt_DF$nodes %in% STRING_nodes$name,]
        
        #pathLengthDist_currT= calcPathsLen(nodeDF= curr_nodes_DF, graph= STRING_net, t1=t1, t2=t2)
          #disabled, because takes very long to run
  
        #summarize degrees => which path lengths are expected
        degree(graph=STRING_net,v=curr_nodes_DF$nodes[curr_nodes_DF$t==t1],mode='out')
        degree(graph=STRING_net,v=curr_nodes_DF$nodes[curr_nodes_DF$t==t2],mode='in')
        
        start_systimeSteiner= Sys.time()
        if(PC_ans== 'n'){
          gS_out= getSteiner(underlNetw=STRING_net,
                  subgr=subgraph_list[[sG_idx]],t1,t2,subgr_idx= sG_idx,
                  design= 'fullCon_visitOne',FC_thres=FC_cutoff,expEvid_cut, prizeToCost,
                  max_St_betwTerm=max_St_betwSliceTerm,gSol_exist,testmodeAns=testmodeAns)  #'fullCon_visitAll'
        } else if (PC_ans== 'y'){
          gS_out= getSteiner(underlNetw=STRING_net,
                  subgr=subgraph_list[[sG_idx]],t1,t2,subgr_idx= sG_idx,
                  design= 'prize_collecting',FC_thres=FC_cutoff,expEvid_cut, prizeToCost,
                  max_St_betwTerm=max_St_betwSliceTerm,gSol_exist,testmodeAns=testmodeAns)
        }
        
        end_systimeSteiner= Sys.time()
        runtimeSteiner= end_systimeSteiner - start_systimeSteiner
        if(attributes(runtimeSteiner)$units== 'mins'){
          runtimeSteiner_min= as.numeric(runtimeSteiner)
        } else if (attributes(runtimeSteiner)$units== 'secs'){
          runtimeSteiner_min= as.numeric(runtimeSteiner)/60
        } else {
          stop('Unrecognized runtime for Steiner')
        }
        runtimeSteinerMtx[sG_idx,t_idx]= runtimeSteiner_min
        
        sigfCountInOut_mtxL[[sG_idx]][t_idx,]= c(gS_out[[3]],gS_out[[4]])
        
        tree_out= gS_out[[1]]
        if(length(V(tree_out))==0){
          print("NO TREE FOUND!")
        } else {
          V(tree_out)$time.slice= t_idx  #no underscore, because this is used for splitting later
        }
        unallocNodes= gS_out[[2]] #unallocated at current TP and subgraph
        
        allUnalloc= c(allUnalloc,unallocNodes)
    
        redrawTree= tree_out
        
        if(length(V(redrawTree))>1){
          # add plot of tree to a collection to eventually plot them in a panel
          treePlot_list= addToTreePlotList(redrawTree, treePlot_list)
          
          ST_sG_list= findSubgraphs(G=redrawTree)
          
          # Filter:
          #PCST algorithm contains code to avoid loops below certain size (6), but at some
          #point this becomes computationally too expensive
          #remove subgraphs from solution that do not have at least one in- and one out-terminal
          filtGraph_list= list()
          for(i in seq_along(ST_sG_list)){
            ST_sG= subgraph(redrawTree,ST_sG_list[[i]])  #subset graph with node indices
            print(plot(ST_sG,vertex.size=10,arrow.width=0.1,main= paste('Subgraph',i)))
            #find out if in-term and out-term are represented in subgraph
            nodeTypesInsG= unique(V(ST_sG)$type)
            if("in-term" %in% nodeTypesInsG & "out-term" %in% nodeTypesInsG){
              filtGraph_list[[length(filtGraph_list)+1]]= ST_sG
            }
          }
          
          for(subgraphNo in seq_along(filtGraph_list)) {
            sG_ST_sG_filt= filtGraph_list[[subgraphNo]]
            
            ### if needed ###
            #svg(paste(plotDir,today,'_',t1,'_to_',t2,'.svg',sep=''))
            plot_timeOrderedSlice(slice_graph= sG_ST_sG_filt, t1 = t1, t2 = t2)
            ### dev.off() ###
            
          }  #end of loop over subgraphs within slice
          if (length(filtGraph_list)>1){
            if( !is.named(filtGraph_list[[1]])){ stop(paste('filtGraph_list is not named! (l.1389)'))}
            ST_sG_filt_1= graph.union(filtGraph_list) #combine the filtered trees for given subgraph
            ST_sG_filt= clean_union_attributes(tree= ST_sG_filt_1) #restore attributes
          } else {
            ST_sG_filt= redrawTree
          }
        } else { #redrawTree didn`t have any nodes
          
          ST_sG_filt= as.directed(redrawTree)
        }
        
        #collect vertices/edges from different timepoints in dataframe
        newNodes_DF= data.frame(vertex.attributes(ST_sG_filt),stringsAsFactors = FALSE) #!!
        edge_ends= ends(ST_sG_filt,es=E(ST_sG_filt))
        #newEdges_DF= data.frame(edge.attributes(ST_sG_filt),stringsAsFactors = FALSE) #!!
        sigf_from= newNodes_DF$sigf[match(edge_ends[,1],newNodes_DF$name)]
        sigf_to= newNodes_DF$sigf[match(edge_ends[,2],newNodes_DF$name)]
        type_from= newNodes_DF$type[match(edge_ends[,1],newNodes_DF$name)]
        type_to= newNodes_DF$type[match(edge_ends[,2],newNodes_DF$name)]
              #no longer use 'type' of edge, but only of nodes
        TP_from= newNodes_DF$t[match(edge_ends[,1],newNodes_DF$name)]
        TP_to= newNodes_DF$t[match(edge_ends[,2],newNodes_DF$name)]
        Edges_l[[sG_idx]]= rbind.data.frame(
          cbind.data.frame(edge_ends,TP_from,TP_to,sigf_from,sigf_to,type_from,type_to),
          Edges_l[[sG_idx]]) #!!
    
        #create union graph for all timepoints of given subgraph
        if(length(ST_sG_filt)>0) { #if there are nodes in graph
          if(length(V(all_t_graph_l[[sG_idx]]))==0){ #initialize
              all_t_graph_l[[sG_idx]]= ST_sG_filt
          } else { #add
            if( !is.named(all_t_graph_l[[sG_idx]]) ) {stop(paste('all_t_graph_l[[sG_idx]] is not named!'))}
            if(!is.named(ST_sG_filt)){stop(paste('ST_sG_filt is not named!'))}
            all_t_graph_l[[sG_idx]]= graph.union(all_t_graph_l[[sG_idx]],ST_sG_filt)
            all_t_graph_l[[sG_idx]]= clean_union_attributes(tree=all_t_graph_l[[sG_idx]]) #this function keeps the
             #first of two conflicting attribute values; as trees are added in temporal order, the earlier TP is kept
          }
        }
       end_systimeLoop= Sys.time()
       runtimeLoop= end_systimeLoop-start_systimeLoop
        if(attributes(runtimeLoop)$units== 'mins'){
          runtimeLoop= as.numeric(runtimeLoop)
        } else if (attributes(runtimeLoop)$units== 'secs') {
          runtimeLoop= as.numeric(runtimeLoop)/60
        } else {
          stop('Loop runtime unit not recognized!')
        }
      runtimeLoopMtx[sG_idx,t_idx]= runtimeLoop 
      } #end of loop over TPs
    treePlot_listList[[sG_idx]]= treePlot_list
    } #end loop over subgraphs
    allUnalloc= unique(allUnalloc)
    print(sigfCountInOut_mtxL)
    plot_L= c(plot_L,'treePlot_list'=treePlot_list)
    
    dev.new(height=3,width=6,noRStudioGD = TRUE)
    par(mar=c(5,4,1,1))
    #par(mai=c(1,1,0,0))
    barplot(sigfCountInOut_mtxL[[1]][,1],names.arg= TPslices, las=2,
            ylab='# Sign. proteins')
    barplot(sigfCountInOut_mtxL[[1]][,2],col='blue', las=2,add=TRUE)
    title(xlab='time slice', line=4)
    legend(1,30,c('input','in tree'),fill=c('grey','blue'))
    D_bar= recordPlot()
    #box(lwd=1,lty=1)
    dev.off()
    png(file=paste(plotDir,today,'_Tree_InOut_',plotSuffix,'.png',sep=''),width=10,height=7,
        units='in',res=300)
    D_bar
    dev.off()
    plot_L= c(plot_L,'D_bar'=list(D_bar))
    
    avgFracRetained= mean(sigfCountInOut_mtxL[[1]][,2]/sigfCountInOut_mtxL[[1]][,1])
    
    # saveRDS(Edges_l,file=paste(loc,today,'_Edges_l.rds',sep=''))
    # Edges_l= readRDS(file=paste(loc,'20230109_Edges_l.rds',sep='')) #laod packages; define loc and functions
    # saveRDS(firstSigDF_prot,file=paste(loc,today,'_firstSigDF_prot.rds',sep=''))
    # firstSigDF_prot= readRDS(file=paste(loc,'20230109_firstSigDF_prot.rds',sep='')) #laod packages; define loc and functions
    # saveRDS(STRING_allSigf_allEvid,file=paste(loc,today,'_STRING_allSigf_allEvid.rds',sep=''))
    # STRING_allSigf_allEvid= readRDS(file=paste(loc,'20230109_STRING_allSigf_allEvid.rds',sep='')) #laod packages; define loc and functions
    
    # TPs_wBetw= sort(c(TPs,TPs[1:(length(TPs)-1)]+2.5))
    # V(all_t_graph_l[[1]])$lvl= match(V(all_t_graph_l[[1]])$t,TPs_wBetw)
    # #all_t_graph_xy <- layout_as_multilevel(all_t_graph_l[[1]],type = "all", alpha = 25, beta = 45)
    # for(alph in seq(0,90,10)){
    #   for(bet in seq(0,90,10)){
    #     all_t_graph_xy <- layout_as_multilevel(all_t_graph_l[[1]],type = "all",alpha=alph,beta=bet)
    #     pdf(file=paste(plotdir,today,'_treeOverTime_alph',alph,'_bet',bet,'.pdf',sep=''))
    #     print(
    #       ggraph(all_t_graph_l[[1]], "manual", x = all_t_graph_xy[, 1], y = all_t_graph_xy[, 2]) +
    #         geom_edge_link(edge_width = 0.3,arrow= arrow(),end_cap = circle(3, 'mm') ) + #do not end in center of node
    #         geom_node_point(aes(color = as.factor(lvl)), size = 3) +
    #         #geom_node_label(aes(label = name,color = V(all_t_graph_l[[1]])$color),  #color = as.factor(lvl)
    #          #               repel = FALSE,cex=2.5 ) +
    #         coord_cartesian(clip = "off", expand = TRUE) +
    #         ggtitle( "Tree over time")+
    #         theme(legend.position = "none",plot.title = element_text(hjust = 0.5)) )  #family="TT Times New Roman"
    #     dev.off()
    #   }
    # }
    # 
    # all_t_graph_xy <- layout_as_multilevel(all_t_graph_l[[1]],type = "all",
    #                                        project2D=FALSE)
    # all_t_graph_l[[1]]$layout <- all_t_graph_xy
    # graphjs(g=all_t_graph_l[[1]],vertex.size=0.5)
    # 
    # plot(0,0,type='n',bty='n',axes=FALSE)
    ######### END OF DISPLAY of graph union over time ############ 
    
    #an edge/node may have been added twice to Steiner tree -> I keep both
  } else if (testmodeAns== 'y') {
    #################
    test_out= testcase_fun()
    g_test= test_out[[1]]   #the test-graph is both the underlying network..
    
    all_t_graph_l= list(g_test)  #..and the graph inferred over time
    term_test= c(test_out[[2]],test_out[[3]])
    max_St_betwSliceTerm= 10
    PCST_OUT= PCST_DIR(expEvid_cut=expEvid_cut,
                testmodeAns= testmodeAns,testcase_fun=testcase_fun,testString=testString)    #file set to testcases; asks for Gurobi results
       #C:\Users\wrtlb\Desktop\20220909_Kanshin_local\20221024_PCST_ILP\20221209_testing
    testTree= PCST_OUT[[1]]
    if(length(V(testTree))== 0){
      stop("Aborting testcase: No tree found!")   #exit test if no tree found
    }
    print( paste('Runtime',PCST_OUT[[2]],sep=' ') )
    Edges_l= list()
    edge_ends= data.frame(ends(testTree,es=E(testTree)),stringsAsFactors = FALSE)
    colnames(edge_ends)= c('node_from','node_to')
    sigf_from= V(testTree)$sigf[match(edge_ends[,1],V(testTree)$name)]
    sigf_to= V(testTree)$sigf[match(edge_ends[,2],V(testTree)$name)]
    type_from= V(testTree)$type[match(edge_ends[,1],V(testTree)$name)]
    type_to= V(testTree)$type[match(edge_ends[,2],V(testTree)$name)]
    #no longer use 'type' of edge, but only of nodes
    TP_from= V(testTree)$t[match(edge_ends[,1],V(testTree)$name)]
    TP_to= V(testTree)$t[match(edge_ends[,2],V(testTree)$name)] 
    Edges_l[[1]]= cbind.data.frame(edge_ends,TP_from,TP_to,sigf_from,sigf_to,type_from,type_to)
    firstSigDF_prot= cbind.data.frame( gene=term_test, 
                     firstSigTP= V(g_test)$t[match(term_test, V(g_test)$name)]  )
    
    firstSigDF_prot$comm= rep(1, dim(firstSigDF_prot)[1])
    firstSigDF_prot$commCol= rep(NA, dim(firstSigDF_prot)[1])
    STRING_allSigf_allEvid=cbind.data.frame( ends(g_test,es=E(g_test)) )
    colnames(STRING_allSigf_allEvid)= c('node1','node2')
    STRING_allSigf= cbind.data.frame(STRING_allSigf_allEvid,
      experimentally_determined_interaction= rep(NA,dim(STRING_allSigf_allEvid)[1]),
      database_annotated= rep(NA,dim(STRING_allSigf_allEvid)[1]),
      automated_textmining= rep(NA,dim(STRING_allSigf_allEvid)[1]),
      combined_score= rep(NA,dim(STRING_allSigf_allEvid)[1]) )
    STRING_nodes= tibble( name= unique(as.character(ends(g_test,es=E(g_test)))),
                           sigf= rep(FALSE,length(unique(as.character(ends(g_test,es=E(g_test)))))) )
    STRING_nodes$sigf[STRING_nodes$name %in% firstSigDF_prot$Gene]= TRUE
    
    STRING_links= cbind.data.frame(STRING_allSigf[,1:2],
                                   costs= 1-STRING_allSigf$experimentally_determined_interaction)
    STRING_net= graph_from_data_frame(d= STRING_links, vertices= STRING_nodes,
                                      directed= FALSE)
    STRING_net$name= STRING_nodes
    commInfo= cbind.data.frame(names= STRING_net$name, commNo= rep(1,length(STRING_net)),
                               commColor= rep('black',length(STRING_net)))
    
    l=layout_with_fr(STRING_net)
    l_nodes= V(STRING_net)$name
    print( plot(STRING_net,layout=l, vertex.size=20,arrow.width=0.1,main=paste('test',testString,sep='')))
    
    btwns_underl= betweenness(graph=test_out[[1]],directed=TRUE)
    ###############
  }
  
  out_atg= plot_graph_overTime(all_t_graph_l,plotSuffix=plotSuffix) #return new 'all_t_graph_l', because some attributes added
  all_t_graph_l2= out_atg[[1]]
  all_t_graph_l2ayout= out_atg[[2]]
  plot_L= c(plot_L,'E1'= list(out_atg[[3]]))  #recorded all_t_graph_l2
  length( (V(all_t_graph_l2[[1]]))[V(all_t_graph_l2[[1]])$sigf==TRUE])  #number of significant nodes
  length( (V(all_t_graph_l2[[1]])))  #number of vertices
  length( E(all_t_graph_l2[[1]]))   #number of edges
  length(findSubgraphs(all_t_graph_l2[[1]]))  #number of subgraphs
  
  # saveRDS(all_t_graph_l2,file=paste(loc,today,'_all_t_graph_l2.rds',sep=''))
  # all_t_graph_l2= readRDS(file=paste(loc,'20230322_all_t_graph_l2.rds',sep='')) #laod packages; define loc and functions
  
  # TPs_wBetw= sort(c(TPs,TPs[1:(length(TPs)-1)]+2.5))
  # V(all_t_graph_l2[[1]])$lvl= match(V(all_t_graph_l2[[1]])$t,TPs_wBetw)
  # #all_t_graph_xy <- layout_as_multilevel(all_t_graph_l2[[1]],type = "all", alpha = 25, beta = 45)
  # for(alph in seq(0,90,10)){
  #   for(bet in seq(0,90,10)){
  #     all_t_graph_xy <- layout_as_multilevel(all_t_graph_l2[[1]],type = "all",alpha=alph,beta=bet)
  #     pdf(file=paste(plotdir,today,'_treeOverTime_alph',alph,'_bet',bet,'.pdf',sep=''))
  #     print(
  #       ggraph(all_t_graph_l2[[1]], "manual", x = all_t_graph_xy[, 1], y = all_t_graph_xy[, 2]) +
  #         geom_edge_link(edge_width = 0.3,arrow= arrow(),end_cap = circle(3, 'mm') ) + #do not end in center of node
  #         geom_node_point(aes(color = as.factor(lvl)), size = 3) +
  #         #geom_node_label(aes(label = name,color = V(all_t_graph_l2[[1]])$color),  #color = as.factor(lvl)
  #          #               repel = FALSE,cex=2.5 ) +
  #         coord_cartesian(clip = "off", expand = TRUE) +
  #         ggtitle( "Tree over time")+
  #         theme(legend.position = "none",plot.title = element_text(hjust = 0.5)) )  #family="TT Times New Roman"
  #     dev.off()
  #   }
  # }
  # 
  # all_t_graph_xy <- layout_as_multilevel(all_t_graph_l2[[1]],type = "all",
  #                                        project2D=FALSE)
  # all_t_graph_l2[[1]]$layout <- all_t_graph_xy
  # graphjs(g=all_t_graph_l2[[1]],vertex.size=0.5)
  # 
  # plot(0,0,type='n',bty='n',axes=FALSE)
  ######### END OF DISPLAY of graph union over time ############ 
  
  
  ################
  ########!!!!!! only subgraph 1 ###########
  #prepare extraction of paths
  currEdges_exp= Edges_l[[1]]
  #currEdges= currEdges[order(currEdges$t),]
  colnames(currEdges_exp)[1:2]= c('node_from','node_to')
  #subset(currEdges_exp,select=-t)  #remove edge-time; deprecated
  currEdges_exp[1:4,]
  #currEdges[currEdges$node_from=='SNF8',]
  
  #give each node an index (for plotting)
  allNfromE= unique(c(currEdges_exp$node_from, currEdges_exp$node_to)) #all nodes from edge-list
  length(allNfromE)  #nodes in Steiner tree
  node_idx_DF= cbind.data.frame(node=allNfromE, idx=seq_along(allNfromE))
  currEdges_exp$idx_from= node_idx_DF$idx[match(currEdges_exp$node_from,node_idx_DF$node)]
  currEdges_exp$idx_to= node_idx_DF$idx[match(currEdges_exp$node_to,node_idx_DF$node)]
  currEdges_exp[1:4,]
  dim(currEdges_exp)  #edges in Steiner tree
  
  #manual check
  #currEdges_exp[currEdges_exp$node_from=='HOG1' | currEdges_exp$node_to=='HOG1',]
  #currEdges_exp[currEdges_exp$node_from=='RIM101',]
  #currEdges_exp[currEdges_exp$node_to=='SNF8',]
  
  #extract paths
  #(I could also extract paths by creating a graph and using all_simple_paths, however, this would
  #not take into account temporal order)
  exPat= extract_paths(currEdges_exp=currEdges_exp,max_St_betwSliceTerm=max_St_betwSliceTerm,
                       max_pathLenSt=max_pathLenSt,end_prot= "")
  (compl_pathsN= exPat[[1]])
  compl_pathsTPs= exPat[[2]]
  (compl_pathsType= exPat[[3]])
  compl_pathsIdxes= exPat[[4]]
  
  if(length(compl_pathsN)== 0) {
    stop('No paths at N step!')
  }
  
  
  #plot paths (untrimmed)
  #chronOUT_untr = chronoPaths( nodePaths_list=compl_pathsN, TPsPaths_list= compl_pathsTPs, 
  #             typesPaths_list= compl_pathsType, sigDF=firstSigDF_prot, maxLen= max_pathLenSt,
  #             plotmode='combPlot',plotSuffix)
  #node_idx_untrim= chronOUT_untr[[1]]
  
  OUT_trimList= trimPaths(compl_pathsN, compl_pathsTPs, compl_pathsType,max_St_betwSliceTerm)
  (compl_pathsN2= OUT_trimList[[1]])
  compl_pathsTPs2= OUT_trimList[[2]]
  compl_pathsType2= OUT_trimList[[3]]
  length(compl_pathsN2)
  
  if(length(compl_pathsN2)==0) {
    stop('No paths at N2 step!')
  }
  
  #there may now be duplicated paths, because two paths may have ended in different protein that was then
  #trimmed off, but may otherwise be the same
  #or one path may be subpath of other, because it started at later node
  # -> remove shorter or duplicated path
  reducOut= reduce_paths(paths= compl_pathsN2,TPs= compl_pathsTPs2,Types= compl_pathsType2)
  compl_pathsN3= reducOut[[1]]
  compl_pathsTPs3= reducOut[[2]]
  compl_pathsType3= reducOut[[3]]
  
  if(length(compl_pathsN3)==0) {
    stop('No paths at N3 step!')
  }
  
  #check how often nodes are unique in network to decide how many plots you need
  allNodes= unlist(compl_pathsN3)
  length(allNodes) #all nodes
  length(unique(unlist(compl_pathsN3))) #unique nodes
  duplNodes= allNodes[which(duplicated(allNodes))]
  length(duplNodes) #duplicated nodes
  length(allNodes[is.na(match(allNodes,duplNodes))]) #nodes only occurring once
  
  ## COUNT where losses occurred
  if(testmodeAns== 'n') {  #only if not in testmode
    noSigfProt= dim(firstSigDF_protA)[1]  #how many significant proteins originally
    #noInSTRING= dim(firstSigDF_prot)[1] #in STRING DB; this is not useful, because STRING output is already filtered
    firstSigDF_prot$Gene [firstSigDF_prot$Gene %in% 
                            unique(c(STRING_allSigf_allEvid$node1,STRING_allSigf_allEvid$node2)) == FALSE]
    noSTRINGfilt= sum(firstSigDF_prot$Gene %in% STRING_nodes$name) #how many in STRING network (filtered
    #for on interaction evidence scores)
    nodesInTree= unique(c(currEdges_exp$node_from,currEdges_exp$node_to))
    noInTree= sum(firstSigDF_prot$Gene %in% nodesInTree) #how many part of tree
    sigfInTree= firstSigDF_prot$Gene[firstSigDF_prot$Gene %in% nodesInTree]
    length(allUnalloc)  #not allocated to tree
    noPathsUntrim= sum(firstSigDF_prot$Gene %in% unique(unlist(compl_pathsN)))  #on paths before trimming
    sigfOnPaths= firstSigDF_prot$Gene[firstSigDF_prot$Gene %in% unique(unlist(compl_pathsN))]
    sigfInTree[sigfInTree %in% sigfOnPaths == FALSE]  #in tree, but not on path
    noPathsTrim= sum(firstSigDF_prot$Gene %in% unique(unlist(compl_pathsN3 ))) #on paths after trimming
    firstSigDF_prot$Gene[firstSigDF_prot$Gene %in% unique(unlist(compl_pathsN3 ))]
    par(mar=c(4,4,1,1))
    
    #png(paste(plotDir,today,'_losses_',plotSuffix,'.png',sep=''))
    x= barplot(c(noSigfProt,noSTRINGfilt,noInTree,noPathsTrim),ylab='# proteins',
               ylim= c(0,noSigfProt+2)) #,noInSTRINGnoPathsUntrim,
    names.arg =c('experiment','STRING filt.','on trees','on paths')  #'STRING','on trim. paths',
    text(cex=1, x=x, y=-3.2, names.arg, xpd=TRUE, srt=90)
    text(4,38,paste('max. Steiner betw. term.:',max_St_betwSliceTerm))
    text(x,c(noSigfProt,noSTRINGfilt,noInTree,noPathsTrim)+1,   #,noInSTRING
         c(noSigfProt,noSTRINGfilt,noInTree,noPathsTrim))  #,noInSTRING
  
    #dev.off()
    
    #write 'losses' to file to combine CST and PCST
    # write.csv(
    #   cbind.data.frame(
    #     var= c("noSigfProt","noSTRINGfilt","noInTree","noPathsTrim"),
    #     count= c(noSigfProt,noSTRINGfilt,noInTree,noPathsTrim),
    #     PC= rep(PC_ans,4)
    #   ),
    #   file= paste(today,'_losses_PC',PC_ans,'.csv',sep=''),row.names = FALSE
    # )
  }
  
  #plot paths 
  #you could still try: each path gets its own color; sorting of points to minimize intersections
  ###### if needed####
  chronOUT= chronoPaths( nodePaths_list=compl_pathsN3, TPsPaths_list= compl_pathsTPs3, 
                              typesPaths_list= compl_pathsType3, sigDF=firstSigDF_prot, maxLen= max_pathLenSt,
                              plotmode='combPlot',plotSuffix)
  chronoPaths( nodePaths_list=compl_pathsN3, TPsPaths_list= compl_pathsTPs3, 
                              typesPaths_list= compl_pathsType3, sigDF=firstSigDF_prot, maxLen= max_pathLenSt,
                              plotmode='indivPlots',plotSuffix)
  node_idx_trim= chronOUT[[1]]
  plot_L= c(plot_L,'E2'= list(chronOUT[[2]]))
  
  #check number of subgraphs in chronopaths
    #DOESN'T WORK LIKE THIS, BECAUSE NODES OCCURING AT MULTIPLE TPs CAUSE CONNECTION
  # chronoPathsDF= data.frame(matrix(NA,nrow=0,ncol=2))
  # for(i in seq_along(compl_pathsN3)) {
  #   for (j in 1:(length(compl_pathsN3[[i]])-1)){
  #    #print(paste(compl_pathsN3[[i]][j],compl_pathsN3[[i]][j+1]))
  #    chronoPathsDF= rbind.data.frame(chronoPathsDF,c(compl_pathsN3[[i]][j],compl_pathsN3[[i]][j+1])) 
  #   }
  # }
  # colnames(chronoPathsDF)=c('node1','node2')
  # g_fromChronoP= graph_from_data_frame(d=chronoPathsDF,directed=TRUE)
  # plot(g_fromChronoP)
  # findSubgraphs(g_fromChronoP)
  
  #color according to communities
  plotSuffix_comm= paste('comm',plotSuffix,sep='_')
  commColor_l= 
    lapply(compl_pathsN3,function(x){V(STRING_net)$commCol[match(x, V(STRING_net)$name)]})
  chronOUT_comm= chronoPaths( nodePaths_list=compl_pathsN3, TPsPaths_list= compl_pathsTPs3, 
                              typesPaths_list= compl_pathsType3, sigDF=firstSigDF_prot, maxLen= max_pathLenSt,
                              plotmode='combPlot',plotSuffix_comm,altColor_l= commColor_l)
  node_idx_trim= chronOUT_comm[[1]]
  chronoPaths( nodePaths_list=compl_pathsN3, TPsPaths_list= compl_pathsTPs3, 
                              typesPaths_list= compl_pathsType3, sigDF=firstSigDF_prot, maxLen= max_pathLenSt,
                              plotmode='indivPlots',plotSuffix_comm,altColor_l= commColor_l)
  
  #select protein Bem2, Ipl1
  compl_pathsN3[which(unlist(lapply(compl_pathsN3,function(x){any('Bem2' %in% x)})))]
  #compl_pathsN3[which(unlist(lapply(compl_pathsN3,function(x){any('Ipl1' %in% x)})))]
  compl_pathsN2[which(unlist(lapply(compl_pathsN2,function(x){any('Bem2' %in% x)})))]
  compl_pathsN[which(unlist(lapply(compl_pathsN,function(x){any('Ipl1' %in% x)})))]
  compl_pathsType[which(unlist(lapply(compl_pathsN,function(x){any('Ipl1' %in% x)})))]
  
  #store all paths in DF
  allTrimPathsDF= data.frame(matrix(NA,nrow=0,ncol=6))
  for(i in seq_along(compl_pathsN3)){
    allTrimPathsDF= rbind.data.frame(allTrimPathsDF,
          cbind.data.frame(rep(i,length(compl_pathsN3[[i]])),compl_pathsN3[[i]],
              node_idx_trim$idx[match(compl_pathsN3[[i]],node_idx_trim$node)], compl_pathsTPs3[[i]],
              compl_pathsType3[[i]], compl_pathsN3[[i]] %in% firstSigDF_prot$Gene) )
  }
  colnames(allTrimPathsDF)= c('pathNo','node','node_idx','TP','type','sigf')
  
  #summarize path lengths: all nodes
  png(paste(plotDir,today,'_pathLenHistAll_',plotSuffix,'.png',sep=''))
  avg_nodInPath= mean(unlist(lapply(compl_pathsN3,length)))
  hist(unlist(lapply(compl_pathsN3,length)),breaks=seq(-0.5,max(unlist(lapply(compl_pathsN3,length)))+0.5),
       xlab='Path lengths',main='')
  dev.off()
  
  #summarize path lengths: termini
  png(paste(plotDir,today,'_pathLenHistTermini_',plotSuffix,'.png',sep=''))
  hist(unlist(lapply(compl_pathsType3, function(x) {sum(x=='out-term' | x=='in-term')})),
       breaks=seq(-0.5,max(unlist(lapply(compl_pathsType3, function(x) {sum(x=='out-term' | x=='in-term')})))+0.5),
       xlab='# sigf. prot. in paths',main='')
  dev.off()
  
  allOutTerm= unique(allTrimPathsDF$node[allTrimPathsDF$type=='out-term'])
  allInTerm= unique(allTrimPathsDF$node[allTrimPathsDF$type=='in-term'])
  extremOutTerm= allOutTerm[allOutTerm %in% allInTerm == FALSE]
  extremInTerm= allInTerm[allInTerm %in% allOutTerm == FALSE]
  
  #the longest paths:
  longPathsOrd= order(unlist(lapply(compl_pathsN2,length)),decreasing=TRUE)
  (compl_pathsN2[longPathsOrd])[1:10]
  (compl_pathsType2[longPathsOrd])[1:10]
  
  nodesDF_N3= data.frame(matrix(NA,nrow=0,ncol=2))
  for (lst in seq_along(compl_pathsN3)) {
    #cat('path: ',compl_pathsN3[[lst]],'\n')
    for(lc in 1:(length(compl_pathsN3[[lst]])-1)) {
      #cat(' element: ',compl_pathsN3[[lst]][lc],'\n')
      nodesDF_N3= rbind.data.frame(nodesDF_N3,c(compl_pathsN3[[lst]][lc],compl_pathsN3[[lst]][lc+1]))
    }
  }
  nodesDF_N3= unique(nodesDF_N3)
  
  # WOULD IT BE BETTER TO MAKE THE ALL-T NETWORK FROM ALL_T_GRAPH[[1]]? WHY NOT?
  #get betweenness centrality in all-t network   
  all_t_graph_sG1_N3= graph_from_data_frame(d=nodesDF_N3,directed=TRUE)
  idx_in_l2= match(V(all_t_graph_sG1_N3)$name,V(all_t_graph_l2[[1]])$name)
  V(all_t_graph_sG1_N3)$t = V(all_t_graph_l2[[1]])$t[idx_in_l2]
  V(all_t_graph_sG1_N3)$type = V(all_t_graph_l2[[1]])$type[idx_in_l2]
  V(all_t_graph_sG1_N3)$type = V(all_t_graph_l2[[1]])$type[idx_in_l2]
  V(all_t_graph_sG1_N3)$sigf = V(all_t_graph_l2[[1]])$sigf[idx_in_l2]
  V(all_t_graph_sG1_N3)$color = V(all_t_graph_l2[[1]])$color[idx_in_l2]
  V(all_t_graph_sG1_N3)$frameCol = V(all_t_graph_l2[[1]])$frameCol[idx_in_l2]
  V(all_t_graph_sG1_N3)$labCol = V(all_t_graph_l2[[1]])$labCol[idx_in_l2]
  if(length(V(all_t_graph_l2[[1]])$commCol)>0){
    V(all_t_graph_sG1_N3)$commCol= 
      V(all_t_graph_l2[[1]])$commCol[idx_in_l2] 
  }
  if(length(V(all_t_graph_l2[[1]])$commColWSigf)>0){
    V(all_t_graph_sG1_N3)$commColWSigf= 
      V(all_t_graph_l2[[1]])$commColWSigf[match(V(all_t_graph_sG1_N3)$name,V(all_t_graph_l2[[1]])$name)]
  }
  if(length(V(all_t_graph_l2[[1]])$time.slice)>0){
    V(all_t_graph_sG1_N3)$time.slice= 
      V(all_t_graph_l2[[1]])$time.slice[match(V(all_t_graph_sG1_N3)$name,V(all_t_graph_l2[[1]])$name)]
  }
  

  if(testmodeAns== 'y'){
    stop('Stopping: No betweenness and evaluation if in testmode!')
  }
  
  btwns= (betweenness(graph=all_t_graph_sG1_N3,v=V(all_t_graph_sG1_N3),directed=TRUE,weights=NULL)) 
  #equivalent paths are counted as fractions; e.g. if there are two shortest paths between two nodes, each counts 1/2 
  x=hist(btwns,main=NULL,breaks=70)
  V(all_t_graph_sG1_N3)$btwns= btwns
  hist(V(all_t_graph_sG1_N3)$btwns[V(all_t_graph_sG1_N3)$type=='in-term' | 
                                     V(all_t_graph_sG1_N3)$type=='out-term'],breaks=70,col=rgb(0,0,1,0.25),alpha=0.5,add=TRUE)

  #correlate betweenness in current network with betweenness in underlying network for current nodes
  if(testmodeAns== 'n'){
    corresp_btwns_underl= btwns_underl[match(V(all_t_graph_sG1_N3)$name, V(STRING_net)$name)]
  } else if (testmodeAns== 'y') {
    corresp_btwns_underl= btwns_underl[match(V(all_t_graph_sG1_N3)$name, V(g_test)$name)]
  }
  log10corresp_btwns_underl= log(corresp_btwns_underl,10)
  log10_0subs_corresp_btwns_underl= log10corresp_btwns_underl
  log10_0subs_corresp_btwns_underl[log10_0subs_corresp_btwns_underl<0]=0
  log10_btwns= log(btwns,10)
  log10_0subs_btwns= log10_btwns
  log10_0subs_btwns[log10_0subs_btwns<0]=0
  if(all(log10_0subs_corresp_btwns_underl==0)){log10_0subs_corresp_btwns_underl= 
    log10_0subs_corresp_btwns_underl+runif(n=length(log10_0subs_corresp_btwns_underl),min=0,max=0.2)}
  if(all(log10_0subs_btwns==0)){log10_0subs_btwns= 
    log10_0subs_btwns+runif(n=length(log10_0subs_btwns),min=0,max=0.2)}
  ########## if needed ############
  svg(paste(plotDir,today,'_btwns_corr_',plotSuffix,'.svg',sep=''))
  par(pty='s')
  plotmax= max(c(log10_0subs_corresp_btwns_underl,log10_0subs_btwns))
  plot(log10_0subs_corresp_btwns_underl,log10_0subs_btwns,xlim=c(0,plotmax),ylim=c(0,plotmax),
       xlab='log10(btwns) Underlying network',ylab='log10(btwns) Inferred network')
  lines(x=c(-1,plotmax+1),y=c(-1,plotmax+1),col='green',lty=2,lwd=2)
  corUnderlInf_btwns= cor(log10_0subs_corresp_btwns_underl,log10_0subs_btwns)
  
  btwn_DF= cbind.data.frame(underl= log10_0subs_corresp_btwns_underl,inferr=log10_0subs_btwns)
  (btwn_DF= btwn_DF[order(btwn_DF$inferr-btwn_DF$underl),])
  
  #nodes for which betweenness centrality in underlying network was <2 and is >10 in inferred network
  #constitute hubs only in inferred network
  infHub_idx= which(corresp_btwns_underl<=2 & btwns>10)  #inferred hub indices ("hub" is not actually a good word)
  V(all_t_graph_sG1_N3)$name[infHub_idx]
  points(log10_0subs_corresp_btwns_underl[infHub_idx],log10_0subs_btwns[infHub_idx],col='red')
  dev.off()
  ###############################
  
  #mark nodes which have high betweenness in inferred compared to underlying network
  all_t_graph_sG1_N3_wHubs= all_t_graph_sG1_N3
  V(all_t_graph_sG1_N3_wHubs)$frameCol= rep('black',length(V(all_t_graph_sG1_N3_wHubs)))
  V(all_t_graph_sG1_N3_wHubs)$frameCol[infHub_idx]= 'red'
  #sorting just so that 'Hub' nodes appear on top
  permIdxes= 1:length(V(all_t_graph_sG1_N3_wHubs))
  permIdxes[! permIdxes %in% infHub_idx]= 1:length(permIdxes[! permIdxes %in% infHub_idx])  #fill up non-hub idxes from 1
  permIdxes[infHub_idx]= (length(permIdxes)-length(infHub_idx)+1):length(permIdxes) #add higher numbers for hubs
  
  all_t_graph_sG1_N3_wHubs_sort= permute.vertices(all_t_graph_sG1_N3_wHubs,permutation=permIdxes) #permutation is new index of this position
  #permute keeps attributes
  graph.isomorphic(all_t_graph_sG1_N3_wHubs,all_t_graph_sG1_N3_wHubs_sort) #needs to be TRUE
  l_orig= layout.fruchterman.reingold(all_t_graph_l2[[1]])
  rescale_mtx= matrix( c(rep(1/20,dim(l_orig)[1]),rep(1/20,dim(l_orig)[1])),ncol=2,byrow = FALSE ) #rescale manually to fit plot
  l= l_orig*rescale_mtx
  l_nodes= V(all_t_graph_l2[[1]])$name
  l_sG1_N3_wHubs_sort= l[match(V(all_t_graph_sG1_N3_wHubs_sort)$name,l_nodes),] #new layout: remove nodes from layout that are no longer present and sort
  ### if needed
  if (PC_ans== 'y'){
    svg(paste(plotDir,today,'_incrBtwnsNetwOvrTime_',plotSuffix,'.svg',sep=''))  #increased betweenness in inferred compared to underlying graph
  } else {
    svg(paste(plotDir,today,'_incrBtwnsNetwOvrTime_',plotSuffix,'.svg',sep=''))
  }
  ###
  #par(mar=c(2,2,0,0))
  leg_TPs= TPs
  par(mar=c(0,0,0,0))
  plot(all_t_graph_sG1_N3_wHubs_sort,layout=l_sG1_N3_wHubs_sort,vertex.size=10,arrow.width=0.1,
       vertex.label.cex=0.5,vertex.label.color=V(all_t_graph_sG1_N3_wHubs_sort)$labCol,
       vertex.frame.color=V(all_t_graph_sG1_N3_wHubs_sort)$frameCol,
       vertex.frame.width=2,rescale=FALSE)  #, axes=TRUE
  if (PC_ans== 'n'){
    legend(0.9,-0.1,legend=leg_TPs,pt.bg=rgb(0.3,leg_TPs/60,0.6),  title='t/s',
           pch=21,col="#777777", pt.cex=2, cex=.8, bty="n", ncol=1)     #
    legend(0.3,0.9,c('other','increased betw. cent.'),  title='Node type',
           pch=1,col=c('black','red'), pt.cex=2, cex=.8, bty="n", ncol=1) 
  } else {
    legend(0.7,1.1,,legend=leg_TPs,pt.bg=rgb(0.3,leg_TPs/60,0.6),  title='t/s',
           pch=21,col="#777777", pt.cex=2, cex=.8, bty="n", ncol=1)     #
    legend(0.3,1.1,c('other','increased betw. cent.'),  title='Node type',
           pch=1,col=c('black','red'), pt.cex=2, cex=.8, bty="n", ncol=1) }
  dev.off()
  ###
  #tkplot(graph= all_t_graph_sG1_N3_wHubs_sort) #interactive
  
  
  #sorting just so that nodes with highest btwns appear on top
  btwns_plusNoise= V(all_t_graph_sG1_N3)$btwns + 
    sample(seq(0,0.01,0.00001),length(V(all_t_graph_sG1_N3)),replace=FALSE) #random noise to avoid ties
  btwnsIdx= match(btwns_plusNoise, sort(btwns_plusNoise))  #get indices: match(A,sort(A))
  
  all_t_graph_sG1_N3_sort=permute.vertices(graph= all_t_graph_sG1_N3,permutation=btwnsIdx)
  l_sG1_N3_sort= l[match(V(all_t_graph_sG1_N3_sort)$name,l_nodes),] #remove nodes from layout that are no longer present and apply sorting
  V(all_t_graph_sG1_N3_sort)$labCol= rep('white',length(V(all_t_graph_sG1_N3_sort)))
  V(all_t_graph_sG1_N3_sort)$labCol[V(all_t_graph_sG1_N3_sort)$t/60 >0.3 | 
                                      V(all_t_graph_sG1_N3_sort)$btwns/300 <5]='black'
  
  plotSuffixBtwn= paste('btwn_',plotSuffix,sep='')
  #reduce the layout to the nodes in graph
  all_t_graph_l2ayout_btw= all_t_graph_l2ayout[
    match(V(all_t_graph_sG1_N3_sort)$name, row.names(all_t_graph_l2ayout)),]
  
  vertexSize_scaleFac= 150 #300
  out_atg_btw= plot_graph_overTime(all_t_graph_l=list(all_t_graph_sG1_N3_sort),plotSuffix=plotSuffixBtwn,
                      layout= all_t_graph_l2ayout_btw,vertexColors=NULL, 
                      vertexSizes=V(all_t_graph_sG1_N3_sort)$btwns/vertexSize_scaleFac,
                      vertexLabCol='black',vertexLabCex=.7,
                      vertexLabDist=.5, vertexLabDegree=pi/2,vertexSize_scaleFac)
  plot_L= c(plot_L,'F1'= list(out_atg_btw[[3]]))  #recorded plot
  
  ###########################
  #separate plots for each terminal (terminal included in path)
  # (note: this is different than just plotting each path individually)
  ###terminals= node_idx_trim$node[node_idx_trim$type=='terminal']
  terminals= unique(allTrimPathsDF$node[allTrimPathsDF$type=='out-term' | allTrimPathsDF$type=='in-term'])
  
  plotIndivTermPaths(terminals = terminals, nodeType_paths_list= compl_pathsType3,
                               nodeIdxDF= node_idx_trim, leg_TPs= leg_TPs) 
    
  
  
  #manual check
  targ= 'CLC1'
  targ_idxes= which(!is.na(unlist(lapply(compl_pathsN3,function(x){match(targ,x)}))))
  for(Ti in targ_idxes) {
    print(compl_pathsN3[[Ti]])
    print(compl_pathsTPs3[[Ti]]) }
  
  #EVALUATION
  
  #for edges on inferred paths: what are STRING scores?
  EdgesOnPaths= data.frame(matrix(NA,nrow=0,ncol=7))
  #columns: node_from, node_to,type_from,type_to expScore, DBscore, textmScore
  for(i in 1:length(compl_pathsN3)) {
    for (nIdx in 1:(length(compl_pathsN3[[i]])-1)){
      #there is edge between each adjacent pair
      EdgesOnPaths= rbind.data.frame(EdgesOnPaths,
                                     c(compl_pathsN3[[i]][c(nIdx,nIdx+1)],
                                       compl_pathsType3[[i]][c(nIdx,nIdx+1)],
                                       NA,NA,NA))
    }
  }
  dim(EdgesOnPaths)
  EdgesOnPaths= unique(EdgesOnPaths)
  dim(EdgesOnPaths)
  colnames(EdgesOnPaths)= c('node1','node2','type_from','type_to',
                            'expScore','DBscore','textmScore')
  
  #get Scores
  for (i in 1:dim(EdgesOnPaths)[1]) {
    STRINGidx= which( (STRING_allSigf[,1]==EdgesOnPaths[i,1] & 
                         STRING_allSigf[,2]==EdgesOnPaths[i,2]) |
                        (STRING_allSigf[,1]==EdgesOnPaths[i,2] & 
                           STRING_allSigf[,2]==EdgesOnPaths[i,1]) )#STRING edges are undirected
    EdgesOnPaths[i,5:7]= 
      c(STRING_allSigf$experimentally_determined_interaction[STRINGidx],
        STRING_allSigf$database_annotated[STRINGidx],
        STRING_allSigf$automated_textmining[STRINGidx])
  }
  EdgesOnPaths$expScore= as.numeric(EdgesOnPaths$expScore)
  EdgesOnPaths$DBscore= as.numeric(EdgesOnPaths$DBscore)
  EdgesOnPaths$textmScore= as.numeric(EdgesOnPaths$textmScore)
  
  #plot scores for edges on paths
  EdgesOnPaths_melt= melt(EdgesOnPaths,,id.vars = c('node1','node2'),
                          measure.vars= c('expScore','DBscore','textmScore'),variable.name='score_type',
                          value.name='score')
  
  dev.new(width=3,height=1.3)
  par(mai=c(0,0.5,0,0))
  pathPlot= ggplot(data= EdgesOnPaths_melt)+
    geom_histogram(aes(x=score,fill=score_type),
                   position = "identity",bins=30,alpha=0.5)+
    theme_classic()+
    ggtitle('Edges on Paths')
  
  
  #plot scores for all edges in underlying network (thresholded to experimentally
  #determined interaction score)
  STRING_allSigf_copy= STRING_allSigf
  colnames(STRING_allSigf_copy)[match(c('experimentally_determined_interaction',
                                        'database_annotated','automated_textmining'),colnames(STRING_allSigf_copy))]=
    c('expScore','DBscore','textmScore')
  STRING_allSigf_melt= melt(STRING_allSigf_copy,id.vars = c('node1','node2'),
                            measure.vars= c('expScore','DBscore','textmScore'),variable.name='score_type',
                            value.name='score')
  dev.off()
  
  dev.new(width=3,height=1.3)
  par(mai=c(0,0.5,0,0))
  undNet_plot= ggplot(data= STRING_allSigf_melt)+
    geom_histogram(aes(x=score,fill=score_type),
                   position = "identity",bins=30,alpha=0.5)+
    theme_classic()+
    ggtitle('Edges in underlying network')
  dev.off()
  
  #plot scores for Steiner tree
  expScore=numeric()
  DBscore=numeric()
  textmScore=numeric()
  for (i in 1:dim(currEdges_exp)[1]) {
    STRINGidx= which( (STRING_allSigf[,1]==currEdges_exp[i,1] & 
                         STRING_allSigf[,2]==currEdges_exp[i,2]) |
                        (STRING_allSigf[,1]==currEdges_exp[i,2] & 
                           STRING_allSigf[,2]==currEdges_exp[i,1]) )#STRING edges are undirected
    expScore[i]= STRING_allSigf$experimentally_determined_interaction[STRINGidx]
    DBscore[i]= STRING_allSigf$database_annotated[STRINGidx]
    textmScore[i]= STRING_allSigf$automated_textmining[STRINGidx]
  }
  currEdges_exp2= cbind.data.frame(currEdges_exp,expScore,DBscore,textmScore)
  currEdges_exp2_melt= melt(currEdges_exp2,id.vars = c('node_from','node_to'),
                            measure.vars= c('expScore','DBscore','textmScore'),
                            variable.name='score_type',value.name='score')
  
  dev.new(width=3,height=1.3)
  par(mai=c(0,0.5,0,0))
  treePlot= ggplot(data= currEdges_exp2_melt)+
    geom_histogram(aes(x=score,fill=score_type),
                   position = "identity",bins=30,alpha=0.5)+
    theme_classic()+
    ggtitle('Edges in Steiner Tree')
  dev.off()
  
  # plot(x=0,y=0,type='n',axes=FALSE,bty='n',xlab='',ylab='')
  # legend(-0.5,0,levels(currEdges_exp2_melt$score_type),pch=22,pt.bg=c('red','green','blue'),
  #        pt.cex=2, cex=.8, bty="n", ncol=1)
  
  #plot all three into one plot (http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/)
  plot_L= c(plot_L,'G1'=list(undNet_plot),'G2'=list(treePlot),'G3'=list(pathPlot))
  svg(paste(plotDir,today,'_edgeScores_',plotSuffix,'.svg',sep=''))
  Gcomb= ggpubr::ggarrange(undNet_plot,
                           treePlot,
                           pathPlot, 
                           #labels = c("A", "B", "C"),
                           ncol = 1, nrow = 3, widths=c(3),heights=c(1.3,1.3,1.3),
                           common.legend = TRUE,legend= 'bottom')
  #ignoring "partial match" warning
  print(Gcomb)
  dev.off()
  
  #test Steiner Tree vs underlying network curated database/textmining scores
  wilcox.test(x=currEdges_exp2$DBscore,y=STRING_allSigf_copy$DBscore,
              alternative='greater')
  wilcox.test(x=currEdges_exp2$textmScore,y=,STRING_allSigf_copy$textmScore,
              alternative='greater')
  #test inferred paths vs underlying network curated database/textmining scores
  wilcox.test(x=EdgesOnPaths$DBscore,y=STRING_allSigf_copy$DBscore,
              alternative='greater')
  wilcox.test(x=EdgesOnPaths$textmScore,y=STRING_allSigf_copy$textmScore,
              alternative='greater')
  #test inferred paths vs Steiner Tree curated database/textmining scores
  wilcox.test(x=EdgesOnPaths$DBscore,y=currEdges_exp2$DBscore,
              alternative='greater')
  wilcox.test(x=EdgesOnPaths$textmScore,y=currEdges_exp2$textmScore,
              alternative='greater')
  
  median(STRING_allSigf_copy$textmScore)
  median(currEdges_exp2$textmScore)
  median(EdgesOnPaths$textmScore)
  
  
  #################################
  #manual investigation
   #which edges are in tree, but not on paths
  currEdges_exp2$node1_2= apply(currEdges_exp2,1,function(x){paste(x[1],x[2],sep='_')})
  EdgesOnPaths$node1_2= apply(EdgesOnPaths,1,function(x){paste(x[1],x[2],sep='_')})
  length(currEdges_exp2$node1_2)
  length(EdgesOnPaths$node1_2)
  onBoth_idx= match(EdgesOnPaths$node1_2,currEdges_exp2$node1_2)
  length(onBoth_idx)
  onlyOnPaths= !(EdgesOnPaths$node1_2 %in% currEdges_exp2$node1_2)  #should be none
  onlyOnTree= !(currEdges_exp2$node1_2 %in% EdgesOnPaths$node1_2)
  #view(currEdges_exp2[onlyOnTree,])
  currEdges_exp2$onPath_ind= (currEdges_exp2$node1_2 %in% EdgesOnPaths$node1_2) #indicator if on path
  currEdges_exp2$Fun12Rps31_ind= str_detect(currEdges_exp2$node1_2,'Fun12') | str_detect(currEdges_exp2$node1_2,'Rps31')

  # arranging multiple plots: https://bookdown.org/ndphillips/YaRrr/arranging-plots-with-parmfrow-and-layout.html
  #dev.new(width=7,height=5,noRStudioGD=TRUE)
  #png(paste(plotDir,today,'_scoreScatterWBoxes_',plotSuffix,'.png',sep=''),width=7,height=5,units='in',res=500)
  #jpeg(paste(plotDir,today,'_scoreScatterWBoxes_',plotSuffix,'.jpg',sep=''),width=7,height=5,units='in',res=500)
  pdf(paste(plotDir,'supplFig1.pdf',sep=''),width=7,height=5)
  layMtx=matrix(c(1,2,3),ncol=3)
  layout(mat=layMtx,heights= 5,width=c(5,1,1)) #one row, two columns
  #layout.show()
  par(pty='s',mar=c(4,4,1,1))
  plot(currEdges_exp2$expScore,currEdges_exp2$textmScore,xlim=c(0,1),ylim=c(0,1),
       xlab='expScore',ylab='textmScore')
  points(currEdges_exp2$expScore[currEdges_exp2$Fun12Rps31_ind],
         currEdges_exp2$textmScore[currEdges_exp2$Fun12Rps31_ind], col='blue',pch=19)
  par(pty='m',mar=c(4,0,1,1))
  boxplot(currEdges_exp2$textmScore[!currEdges_exp2$Fun12Rps31_ind],ylim=c(0,1),yaxt='n',frame=FALSE)
  boxplot(currEdges_exp2$textmScore[currEdges_exp2$Fun12Rps31_ind],ylim=c(0,1),col='blue',yaxt='n',frame=FALSE)
  dev.off()
  #########################
  
  #evaluate via number of kinases connected to terminal nodes
  #..on paths
  kinases[1:3]
  EdgesOnPaths[1:3,]
  
  noEdgesToTerm_inPaths= 0  #edges to terminal
  noEdgesKinToTerm_inPaths= 0 #edges from kinase to terminal
  kin_term_paths= data.frame(matrix(NA,nrow=0,ncol=7))  #dataframe listing paths from kinase to terminal
  for (i in 1:dim(EdgesOnPaths)[1]) {
    if(EdgesOnPaths$type_to[i]== 'in-term') {
      noEdgesToTerm_inPaths= noEdgesToTerm_inPaths+1
      if(EdgesOnPaths$node1[i] %in% kinases) {
        noEdgesKinToTerm_inPaths= noEdgesKinToTerm_inPaths+1
        kin_term_paths= rbind.data.frame(kin_term_paths,EdgesOnPaths[i,])
      }
    }
  }
  noEdgesToTerm_inPaths
  noEdgesKinToTerm_inPaths
  kin_term_paths
  
  #..in underlying network: edges from kinase to any node (irrespective if significant or not)
  #STRING network is undirected, therefore you cannot tell if edge is coming from or going to kinase
  STRING_allSigf[1:3,]
  noEdges_inUndNet= dim(STRING_allSigf)[1]  #edges
  noEdgesNeighKin_inUndNet= 0 #edges from/to kinase (i.e. neighbors of kinases)
  for (i in 1:dim(STRING_allSigf)[1]) {
    if(STRING_allSigf$node1[i] %in% kinases | STRING_allSigf$node2[i] %in% kinases) {
      noEdgesNeighKin_inUndNet= noEdgesNeighKin_inUndNet+1
    }
  }
  noEdges_inUndNet
  noEdgesNeighKin_inUndNet
  
  #..in underlying network: edges between kinase and significant node
  #STRING network is undirected, therefore you cannot tell if edge is coming from or going to kinase
  STRING_allSigf[1:3,]
  noEdgesNeighTerm_underl= 0  #edges neighboring terminal
  noEdgesKinNeighTerm_underl= 0 #edges between kinase and terminal
  for (i in 1:dim(STRING_allSigf)[1]) {
    #"forward direction"
    if(STRING_allSigf$sigf_node2[i]== TRUE) {
      noEdgesNeighTerm_underl= noEdgesNeighTerm_underl+1
      if(STRING_allSigf$node1[i] %in% kinases) {
        noEdgesKinNeighTerm_underl= noEdgesKinNeighTerm_underl+1
      }
    }
    #"backward direction"
    if(STRING_allSigf$sigf_node1[i]== TRUE) {
      noEdgesNeighTerm_underl= noEdgesNeighTerm_underl+1
      if(STRING_allSigf$node2[i] %in% kinases) {
        noEdgesKinNeighTerm_underl= noEdgesKinNeighTerm_underl+1
      }
    }
  }
  noEdgesNeighTerm_underl
  noEdgesKinNeighTerm_underl
  
  
  #on Steiner tree
  kinases[ kinases %in% V(all_t_graph_l2[[1]])$name ]   #which kinases
  currEdges_exp2[1:3,]
  noEdgesToTerm_inSteiner= 0  #edges to terminal
  noEdgesKinToTerm_inSteiner= 0 #edges from kinase to terminal
  for (i in 1:dim(currEdges_exp2)[1]) {
    if(currEdges_exp2$type_to[i]== 'in-term') {
      noEdgesToTerm_inSteiner= noEdgesToTerm_inSteiner+1
      if(currEdges_exp2$node_from[i] %in% kinases) {
        noEdgesKinToTerm_inSteiner= noEdgesKinToTerm_inSteiner+1
      }
    }
  }
  currEdges_exp2[currEdges_exp2$node_from %in% kinases & currEdges_exp2$type_to== 'in-term',c(1,2)] #inspect edges
  noEdgesToTerm_inSteiner
  noEdgesKinToTerm_inSteiner
  
  #test: Steiner vs. underl. network
  kin_2way= matrix(c(noEdgesToTerm_inSteiner,noEdgesKinToTerm_inSteiner,
                     noEdgesNeighTerm_underl,noEdgesKinNeighTerm_underl),nrow=2,ncol=2)
  fisher.test(kin_2way,alternative="less")
  
  #test: paths vs. underl. network
  kin_2way= matrix(c(noEdgesToTerm_inPaths,noEdgesKinToTerm_inPaths,
                     noEdgesNeighTerm_underl,noEdgesKinNeighTerm_underl),nrow=2,ncol=2)
  fisher.test(kin_2way,alternative="less")

plot_LoL[[PC_ans_idx]]= plot_L
}



######## ARRANGED PLOTS ##########
plot_LoL

A1 = plot_LoL[[1]][['A1']]
A2 = plot_LoL[[1]][['A2']]
#png(paste(loc,'plots/',today,'_FigA.png',sep=''),width=5,height=20,units='in',res=500)
pdf(paste(loc,'plots/','FigA.pdf',sep = ''), width = 5, height = 20)
ggpubr::ggarrange(A1, A2, nrow = 2, ncol = 1, labels = c('A','B'),
                  widths = 3.5, heights = c(3.5,10.5))
dev.off()

# losses figure (largely manual)
losses_CST= read.csv(file=paste(loc,'20230411_losses_PCn.csv',sep=''))
losses_PCST= read.csv(file=paste(loc,'20230411_losses_PCy.csv',sep=''))
losses= rbind.data.frame(losses_CST,losses_PCST[3:4,])
losses$PC[c(1,2)]='both'
losses$var = factor(
  c('experiment','STRING filt.','CST\non trees','CST\non paths','PCST\non trees','PCST\non paths'),
  levels= c('experiment','STRING filt.','CST\non trees','CST\non paths','PCST\non trees','PCST\non paths'),
  ordered= is.ordered(c('experiment','STRING filt.','CST\non trees','CST\non paths','PCST\non trees','PCST\non paths')))

par(mar=c(4,5,1,1))
B3= ggplot(data=losses, aes(x=var,y=count,group=PC,fill=PC))+
  geom_bar(stat='identity',position=position_dodge(),color='black')+
  scale_fill_manual(values= c('#9e9e9d','#595957','#2e2e2d') )+  #grDevices::terrain.colors(3) 
  geom_text(aes(x=1:6,y=count+5,label= losses$count))+
  scale_y_continuous(name='# proteins')+
  theme_classic()+
  theme(axis.title.x= element_blank(), text=element_text(size=15),
        legend.position = "none")
#B3= recordPlot()
png(paste(plotDir,today,'_losses.png',sep=''))
print(B3)
dev.off()

# fig. B
B1= plot_LoL[[1]][['B1']]
B2= plot_LoL[[1]][['B2']]
#png(paste(loc,'plots/',today,'_FigB.png',sep=''),width=21,height=6,units='in',res=500)
pdf(paste(loc,'plots/','FigB.pdf',sep=''),width=21,height=6)
ggarrange(B1,B2,B3,nrow=1,ncol=3,labels=c('A','B','C'),widths=c(3,2.1,2.9),heights=2.5)
dev.off()


# fig. C
C1= plot_LoL[[1]][['C1']]
C2= plot_LoL[[1]][['C2']]
C3= plot_LoL[[1]][['C3']]
#png(paste(loc,'plots/',today,'_FigC.png',sep=''),width=10,height=8,units='in',res=500)
pdf(paste(loc,'plots/','FigC.pdf',sep=''),width=10,height=8)
ggarrange(C1,
          ggarrange(C2,C3,nrow=2,ncol=1,labels=c('B','C'),widths=c(4),heights=c(4,4)),
nrow=1,ncol=2,widths=c(6,4),heights=c(8),labels=c('A') )
dev.off()

# fig. D
D1a= plot_LoL[[1]][['treePlot_list1']]
D2a= plot_LoL[[1]][['treePlot_list2']]
D3a= plot_LoL[[1]][['treePlot_list3']]
D4a= plot_LoL[[1]][['treePlot_list4']]
D5a= plot_LoL[[1]][['treePlot_list5']]
D6a= plot_LoL[[1]][['treePlot_list6']]
D7a= plot_LoL[[1]][['treePlot_list7']]
D8a= plot_LoL[[1]][['treePlot_list8']]
D9a= plot_LoL[[1]][['treePlot_list9']]
D10a= plot_LoL[[1]][['treePlot_list10']]
D_bar_a= plot_LoL[[1]][['D_bar']]
D1b= plot_LoL[[2]][['treePlot_list1']]
D2b= plot_LoL[[2]][['treePlot_list2']]
D3b= plot_LoL[[2]][['treePlot_list3']]
D4b= plot_LoL[[2]][['treePlot_list4']]
D5b= plot_LoL[[2]][['treePlot_list5']]
D6b= plot_LoL[[2]][['treePlot_list6']]
D7b= plot_LoL[[2]][['treePlot_list7']]
D8b= plot_LoL[[2]][['treePlot_list8']]
D9b= plot_LoL[[2]][['treePlot_list9']]
D10b= plot_LoL[[2]][['treePlot_list10']]
D_bar_b= plot_LoL[[2]][['D_bar']]

#png(paste(loc,'plots/',today,'_FigD.png',sep=''),width=12.5,height=20,unit='in',res=500)
pdf(paste(loc,'plots/','FigD.pdf',sep=''),width=12.5,height=20)
ggarrange(
  ggarrange(D1a,D2a,D3a,D4a,
            nrow=1,ncol=4,widths=c(3,3,3,3),heights=c(3),labels=c('A','','')),
  ggarrange(D5a,D6a,D7a,D8a,
            nrow=1,ncol=4,widths=c(3,3,3,3),heights=c(3)),
  ggarrange(
    ggarrange(D9a,D10a,nrow=1,ncol=2,widths = c(3,3),heights=c(3)), 
    D_bar_a,nrow=1,ncol=2,widths=c(6,6),heights=c(3),labels=c('','B') ),
  
  ggarrange(D1b,D2b,D3b,D4b,
            nrow=1,ncol=4,widths=c(3,3,3,3),heights=c(3),labels=c('C','','')),
  ggarrange(D5b,D6b,D7b,D8b,
            nrow=1,ncol=4,widths=c(3,3,3,3),heights=c(3)),
  ggarrange(
    ggarrange(D9b,D10b,nrow=1,ncol=2,widths = c(3,3),heights=c(3)), 
    D_bar_b,nrow=1,ncol=2,widths=c(6,6),heights=c(3),labels=c('','D') ),
  
  nrow=6,ncol=1,widths=c(12),heights=c(3,3,3,3,3,3) )
dev.off()

# fig. E
E1a= plot_LoL[[1]][['E1']]
E2a= plot_LoL[[1]][['E2']]
E1b= plot_LoL[[2]][['E1']]
E2b= plot_LoL[[2]][['E2']]
svg(paste(loc,'plots/',today,'_FigE.svg',sep=''))
#png(paste(loc,'plots/',today,'_FigE.png',sep=''),width=6,height=11.5,units='in',res=500)
#pdf(paste(loc,'plots/','FigE.pdf',sep=''),width=6,height=11.5)
ggarrange(
  ggarrange(E1a,E2a,nrow=2,ncol=1,labels=c('A','C'),widths=c(3),heights=c(3,8.5)),
  ggarrange(E1b,E2b,nrow=2,ncol=1,labels=c('B','D'),widths=c(3),heights=c(3,8.5)),
nrow=1,ncol=2,widths=c(3,3),heights=c(11.5) )
dev.off()

# fig. F
F1a= plot_LoL[[1]][['F1']]
F1b= plot_LoL[[2]][['F1']]
svg(paste(loc,'plots/',today,'_FigF.svg',sep=''))
#png(paste(loc,'plots/',today,'_FigF.png',sep=''),width=8,height=4,units='in',res=500)
#pdf(paste(loc,'plots/','FigF.pdf',sep=''),width=8,height=4)
ggarrange(F1a,F1b,nrow=1,ncol=2,labels=c('A','B'),widths=c(4,4),heights=c(4))
dev.off()

# fig. G
G1a= plot_LoL[[1]][['G1']]
G2a= plot_LoL[[1]][['G2']]
G3a= plot_LoL[[1]][['G3']]
G1b= plot_LoL[[2]][['G1']]
G2b= plot_LoL[[2]][['G2']]
G3b= plot_LoL[[2]][['G3']]
#png(paste(loc,'plots/',today,'_FigG.png',sep=''),width=6,height=6,units='in',res=500)
pdf(paste(loc,'plots/','FigG.pdf',sep=''),width=6,height=6)
#dev.new(width=6,height=6,noRStudioGD = TRUE)
leg= get_legend(G1a,position='bottom')
ggarrange(
  ggarrange(G1a,G2a,G3a,nrow=3,ncol=1,widths=c(3),heights=c(2,2,2),legend="none"),
  ggarrange(G1b,G2b,G3b,nrow=3,ncol=1,widths=c(3),heights=c(2,2,2),legend="none"),
nrow=1,ncol=2,widths=c(3,3),heights=c(6),labels=c('A','B'),
common.legend = TRUE,legend= 'bottom',legend.grob=leg)
dev.off()



########## some manual explorations ##############

##### find literature proteins of HOG-pathway
litProt= c('Hrk1','Msb2','Sln1','Ypd1','Ssk1','Ssk2','Ssk22','Sho1','Cdc42','Opy2','Ste20','Ste50','Ste11',
           'Pbs2','Hog1','Cla4','Crm1','Ptc1', 'Ptp2', 'Ptp3',
           'Nha1', 'Tok1', 'Sic1', 'Msn2', 'Msn4', 'Sko1', 'Hot1', 'Smp1')
litProt[litProt %in% V(all_t_graph_l2[[1]])$name]
litProt[litProt %in% unique(unlist(compl_pathsN3))]

Ptp3_idxes= which(unlist(lapply(compl_pathsN3,function(x){'Ptp3' %in% x})))
compl_pathsN3[Ptp3_idxes]
Hog1_idxes= which(unlist(lapply(compl_pathsN3,function(x){'Hog1' %in% x})))
compl_pathsN3[Hog1_idxes]

litProt[litProt %in% firstSigDF_prot$Gene]

litExpDiffPhos= c('Ssk2','Ssk22','Ste11', 
                  'Pbs2','Hog1',
                  'Nha1', 'Tok1', 'Sic1', 'Msn2', 'Msn4', 'Sko1', 'Hot1', 'Smp1')  #proteins with expected observable phosphorylation difference
litExpDiffPhos[litExpDiffPhos %in% firstSigDF_prot$Gene]  #only sigf. prot. found in STRING
#litExpDiffPhos[litExpDiffPhos %in% firstSigDF_protA$Gene]

compl_pathsN[ which(unlist(lapply(compl_pathsN, function(x){'Ptc1' %in% x}))) ]

#rate of mapping of significant proteins
(mappedSig= firstSigDF_prot$Gene[ firstSigDF_prot$Gene %in% V(all_t_graph_l2[[1]])$name ])
noSig= length(firstSigDF_prot$Gene)
noMapSig= length(mappedSig)
noMapSig/noSig

# proteins that are not expected to be differentially phosphorylated
litNotDiffPhos= litProt[ litProt %in% litExpDiffPhos == FALSE ]
litNotDiffPhos [litNotDiffPhos %in% V(all_t_graph_l2[[1]])$name ]

#proteins with expected, but unobserved differential phosphorylation
expButNotObs= litExpDiffPhos[litExpDiffPhos %in% firstSigDF_prot$Gene == FALSE]
expButNotObs[expButNotObs %in% V(all_t_graph_l2[[1]])$name]

#GO analysis
#done externally: https://yeastgenome.org/goTermFinder
targets= V(all_t_graph_l2[[1]])$name
# write.table(V(all_t_graph_l2[[1]])$name,paste(loc,today,'_protIn_all_t_graph_',plotSuffix,'.txt',sep=''),
#             ,row.names = FALSE)
underlProt_1shell_filt= V(STRING_net)$name
 # write.table(V(STRING_net)$name,paste(loc,today,'_protInUnderl_expEvid',expEvid_cut,'_',
 #                                      plotSuffix,'.txt',sep=''),row.names = FALSE)
posRegCellCyc= str_to_title(c('GLC7', 'EDE1', 'RIM15', 'MIH1', 'CLB2', 'CDC28', 'IGO1', 'ZDS1', 'DBF4', 'IGO2', 
                              'CDC5', 'KSS1', 'CDC55'))
posRegCellCyc[ posRegCellCyc %in% firstSigDF_prot$Gene]

#color chronopath by posRegCellCyc
plotSuffix_CellCyc= paste('CellCyc',plotSuffix,sep='_')
cellCycColor_l= lapply(compl_pathsN3,function(x){
  A= rep('black',length(x))
  A[x %in% posRegCellCyc] = 'brown'
  return(A)  })
chronOUT_cc= chronoPaths( nodePaths_list=compl_pathsN3, TPsPaths_list= compl_pathsTPs3, 
                            typesPaths_list= compl_pathsType3, sigDF=firstSigDF_prot, maxLen= max_pathLenSt,
                            plotmode='combPlot',plotSuffix= plotSuffix_CellCyc,altColor_l= cellCycColor_l)
node_idx_trim= chronOUT_cc[[1]]
   
Rpo21_pathIdxes= which(unlist(lapply(compl_pathsN3, function(x){'Rpo21' %in% x})))
compl_pathsN3[Rpo21_pathIdxes]

##########
which(unlist(lapply(compl_pathsN3, function(x){'Fun12' %in% x})))
compl_pathsN3[ which(unlist(lapply(compl_pathsN3, function(x){'Fun12' %in% x}))) ]
compl_pathsN2[ which(unlist(lapply(compl_pathsN2, function(x){'Fun12' %in% x}))) ]
compl_pathsType2[ which(unlist(lapply(compl_pathsN2, function(x){'Fun12' %in% x}))) ]
compl_pathsN[ which(unlist(lapply(compl_pathsN, function(x){'Fun12' %in% x}))) ]
nodesAroundFun12= unique(unlist(compl_pathsN[ which(unlist(lapply(compl_pathsN, function(x){'Fun12' %in% x}))) ]))
compl_pathsType[ which(unlist(lapply(compl_pathsN, function(x){'Fun12' %in% x}))) ]
V(all_t_graph_l2[[1]])$name[ match(nodesAroundFun12, V(all_t_graph_l2[[1]])$name) ]
V(all_t_graph_l2[[1]])$prizes[ match(nodesAroundFun12, V(all_t_graph_l2[[1]])$name) ]

edxIdxAroundFun12= unique( c(match(nodesAroundFun12, ends(all_t_graph_l2[[1]],es=E(all_t_graph_l2[[1]]))[,1]),
                             match(nodesAroundFun12, ends(all_t_graph_l2[[1]],es=E(all_t_graph_l2[[1]]))[,2]) ) )
E(all_t_graph_l2[[1]])$costs[edxIdxAroundFun12]

which(unlist(lapply(compl_pathsN3,function(x){'Hsl1' %in% x})))
which(unlist(lapply(compl_pathsN3,function(x){'Gpa2' %in% x})))
which(unlist(lapply(compl_pathsN3,function(x){'Sko1' %in% x})) &
      unlist(lapply(compl_pathsN3,function(x){'Hog1' %in% x})) &
      unlist(lapply(compl_pathsN3,function(x){'Ptp3' %in% x})) )

which(unlist(lapply(compl_pathsN3,function(x){'Pbs2' %in% x})) &
        unlist(lapply(compl_pathsN3,function(x){'Hog1' %in% x})) &
        unlist(lapply(compl_pathsN3,function(x){'Ptp3' %in% x})) )

which(unlist(lapply(compl_pathsN3,function(x){'Ssn2' %in% x})) &
        unlist(lapply(compl_pathsN3,function(x){'Rgr1' %in% x})) &
        unlist(lapply(compl_pathsN3,function(x){'Rtf1' %in% x})) )

which(unlist(lapply(compl_pathsN3,function(x){'Ssn2' %in% x})) &
        unlist(lapply(compl_pathsN3,function(x){'Srb7' %in% x})) &
        unlist(lapply(compl_pathsN3,function(x){'Paf1' %in% x})) )

which(  unlist(lapply(compl_pathsN3,function(x){'Rpo21' %in% x})) &
        unlist(lapply(compl_pathsN3,function(x){'Hog1' %in% x})))

which(unlist(lapply(compl_pathsN3,function(x){'Ras2' %in% x})) &
        unlist(lapply(compl_pathsN3,function(x){'Tpk1' %in% x})) &
        unlist(lapply(compl_pathsN3,function(x){'Pfk26' %in% x})) )

c('Ssn2','Rgr1','Spt15','Bdp1','Rpo26','Spt5','Rtf1','Rpo21') %in% firstSigDF_protA$Gene

#enrichment of proteins of HOG-pathway in network vs. enrichment of significant proteins
prot_underl=  length(V(STRING_net))   #proteins in underlying network
sigf_underl= dim(firstSigDF_protA)[1] #number of significant proteins in underlying network
HOG_underl= sum(litProt %in% V(STRING_net)$name)  #HOG-pathway proteins in underlying network

prot_infer=  length(V(all_t_graph_l2[[1]]))   #proteins in inferred network
sigf_infer= sum(firstSigDF_protA$Gene %in% V(all_t_graph_l2[[1]])$name) #number of significant proteins in inferred network
HOG_infer= sum(litProt %in% V(all_t_graph_l2[[1]])$name)  #HOG-pathway proteins in inferred network

prot_paths= length(unique(unlist(compl_pathsN3)))  #proteins on paths
sigf_paths= sum(firstSigDF_protA$Gene %in% unique(unlist(compl_pathsN3))) #number of significant proteins on paths
HOG_paths= sum(litProt %in% unique(unlist(compl_pathsN3)))  #HOG-pathway proteins on paths

(sigf_enrich= (sigf_infer/sigf_underl) / (prot_infer/prot_underl))
(HOG_enrich= (HOG_infer/HOG_underl) / (prot_infer/prot_underl))
(HOG_enrichPaths= (HOG_paths/HOG_underl) / (prot_paths/prot_underl))



################
#for all significant proteins:
 #put node1 and node2 each in both column1 and column2; use dcast to create a symmetric matrix with all
 #proteins in all rows and columns
 #cluster; then remove one triangle
  #fill the other triangle with tree proteins

#1, reverse node1 with node2 column
dim(STRING_allSigf)
STRING_allSigf_rev= STRING_allSigf
STRING_allSigf_rev$node1= STRING_allSigf$node2
STRING_allSigf_rev$node2= STRING_allSigf$node1
STRING_allSigf_origRev= rbind.data.frame(STRING_allSigf,STRING_allSigf_rev)
dim(STRING_allSigf_origRev)

#2, create matrix with dcast
length(unique(c(STRING_allSigf$node1,STRING_allSigf$node2)))  #you should end up with mtx of these dimensions
STRING_allSigf_mtx_orig= dcast(data=STRING_allSigf_origRev,node1~node2,value.var='automated_textmining') #'database_annotated'
dim(STRING_allSigf_mtx_orig)
STRING_allSigf_mtx_orig[1:4,1:4]
STRING_allSigf_mtx= as.matrix(STRING_allSigf_mtx_orig[,2:dim(STRING_allSigf_mtx_orig)[2]])
rownames(STRING_allSigf_mtx)= as.character(STRING_allSigf_mtx_orig$node1)
STRING_allSigf_mtx[is.na(STRING_allSigf_mtx)]= 0  #heatmap.2 has problems with NA
STRING_allSigf_mtx[1:4,1:4]

#3, cluster
library('gplots')
library('RColorBrewer')
clust= hclust( dist(STRING_allSigf_mtx,method='euclidean'), method='ward.D2')
 # https://www.biostars.org/p/398548/
clustNo= 4 #8
clust_STRINGallSigf= cutree(clust,k=clustNo)
clustCol= brewer.pal(8, 'Accent')[clust_STRINGallSigf]  #rainbow(8)
#png('20230421_clustHeatm_DBscore_FC3_PCST_p2c1_maxTree2_maxPath10.png')   #manual
heat_STRING= heatmap.2(STRING_allSigf_mtx,trace='none',col=colorRampPalette(brewer.pal(9,'Blues')),
          dendrogram='none', Rowv=as.dendrogram(clust),Colv=as.dendrogram(clust),
          labRow=FALSE , labCol=FALSE,
          colRow= clustCol,RowSideColors = clustCol, # to add nice colored strips  
          main='textmScore') #'DBscore'
#dev.off()

plot(clust)
head(clust_STRINGallSigf)
head(IDmap)
clustDF= cbind.data.frame(as.data.frame(clust_STRINGallSigf),clustCol)
clustDF$gene= rownames(clustDF)
clustDF2= merge(x=clustDF, y=IDmap, by.x='gene', by.y='Gene.Name',all.x=TRUE)
#write.csv(clustDF2,paste(loc,today,'_clustSTRING_allSigf_DBscore2.csv',sep=''),row.names = FALSE)
#write.csv(clustDF2,paste(loc,today,'_clustSTRING_allSigf_textmScore2.csv',sep=''),row.names = FALSE)

#DB score
manu_clustGO= c('misc',
                'chromatin binding',
                'ribosome',
                'RNA binding',
                'DNA binding',
                'RNA binding',
                'transcription regulator activity',
                'ribosome'
)

#textm score
manu_clustGO= c('misc',
                'ribosome',
                'RNA binding',
                'transcription regulator activity'
)

plot.new()
legend(x=.2,y=.8,fill=brewer.pal(clustNo, 'Accent'),legend=manu_clustGO) #x=-1,y=2
heat_STRING$rowInd
STRING_allSigf_mtx_ord= STRING_allSigf_mtx[heat_STRING$rowInd,heat_STRING$colInd]  #re-order
STRING_allSigf_mtx_ord_lower= STRING_allSigf_mtx_ord
STRING_allSigf_mtx_ord_lower[ upper.tri(STRING_allSigf_mtx_ord_lower) ]= 0  #lower triangle remains


#check heatmap
heatmap.2(STRING_allSigf_mtx_ord_lower,trace='none',col=colorRampPalette(brewer.pal(9,'Blues')),
          dendrogram='none', Rowv=FALSE, Colv=FALSE, labRow=FALSE, labCol=FALSE)   #no re-clustering

#4, fill upper triangle with "tree" values
 #currEdges_exp2 may have several entries for the same edge, because it may exist at several times
  #they will have the same scores -> remove duplicates
dim(currEdges_exp2)
tree_dupRem= currEdges_exp2
tree_dupRem= tree_dupRem[ !duplicated(tree_dupRem[,c(1,2)]),]
dim(tree_dupRem)
length(unique(currEdges_exp2$node1_2))  #double check

 #make full matrix from tree edges (pretend they are undirected)
tree_rev= tree_dupRem
tree_rev$node_from= tree_dupRem$node_to
tree_rev$node_to= tree_dupRem$node_from
tree_rev$node1_2= paste(tree_rev$node_from,tree_rev$node_to,sep='_')
tree_origRev= rbind.data.frame(tree_dupRem,tree_rev)
dim(tree_origRev)
 #there may have been edges in both directions
tree_origRev_dupRem= tree_origRev[ !duplicated(tree_origRev[,c(1,2)]),]
dim(tree_origRev_dupRem)

 #look up edge scores
rNam= rownames(STRING_allSigf_mtx_ord)
cNam= colnames(STRING_allSigf_mtx_ord)
rNam[1:10]
cNam[1:10]
STRING_Tree_mtx= STRING_allSigf_mtx_ord_lower
for (ro in 1:dim(STRING_Tree_mtx)[1]){
  for(co in (ro+1):dim(STRING_Tree_mtx)[2]) {
    idx= which(tree_origRev_dupRem$node_from== rNam[ro] & tree_origRev_dupRem$node_to== cNam[co])
    if(length(idx) == 1){
      STRING_Tree_mtx[ro,co]= tree_origRev_dupRem$DBscore[idx]
    } else if (length(idx) > 1) {
      stop("Multiple entries detected for the same edge in tree!")
    }
  }
}


#5, Show the heatmap with both triangles filled
heatmap.2(STRING_Tree_mtx,trace='none',col=colorRampPalette(brewer.pal(9,'Blues')),
          dendrogram='none', Rowv=FALSE, Colv=FALSE, labRow=FALSE, labCol=FALSE)   #no clustering

# Can we determine what the populations with low and high DB score are
DB_le0.25= STRING_allSigf_copy[STRING_allSigf_copy$DBscore<.25,]
#write.csv(DB_le0.25,paste(loc,today,'DBscore_le0.25.csv',sep=''))
DB_gr0.75= STRING_allSigf_copy[STRING_allSigf_copy$DBscore>0.75,]
#write.csv(DB_gr0.75,paste(loc,today,'DBscore_gr0.75.csv',sep=''))
