PCST_DIR= function(dir_graph,out_term,in_term,t1,t2,FC_thres,expEvid_cut,prizeToCost,max_St_betwTerm,gSol_exist,
                   testmodeAns,testcase_fun,testString) {

  #derived from '20221107_PCST_ILP_DIR.R'
  #plot.new()
  closeAllConnections()
  startTime= Sys.time()
  
  source(paste('./20221024_PCST_ILP',
               '20230310_SUBTOUR_ELIM.R',sep='/'))
  
  ############
  GRAPH_FROM_SIMPPATHS= function(graph,out_terms,in_terms,maxNoEdges) {
    #function to create subgraph from edges between out-terminals and in-terminals
     #using all_simple_paths
     #"simple path": path on which each vertex is visited only once
    print( paste('Graph with',length(V(graph)),'nodes with',
            length(out_terms),'out-terms and',length(in_terms),'in-terms.') )
    allowEdges= data.frame(matrix(NA,ncol=2,nrow=0))  #each entry contains nodes of allowed path
    for(o in out_terms) {
      currAllowPaths=
        all_simple_paths(graph=graph,from=o,to=in_terms,mode="out",maxNoEdges) 
      #cutoff necessary for speed; cutoff is number of edges
      for (p in seq_along(currAllowPaths)) {
        aPn= currAllowPaths[[p]]$name
        #print(aPn)
        allowEdges= unique(rbind.data.frame(allowEdges,
                                            cbind.data.frame(aPn[1:(length(aPn)-1)],aPn[2:length(aPn)])) )
      }
    }
    allowVertices= unique(c(allowEdges[,1],allowEdges[,2]))
    print(paste(length(allowVertices),'allowed vertices;',dim(allowEdges)[1],'allowed edges',sep=' '))
    
    #remove nodes that are not on these paths
    #this also removes termini that are not on path from in-terminus to out-terminus
    subgr_fromPaths= graph_from_data_frame(d=allowEdges,directed=TRUE,vertices=allowVertices)
    if(length(V(subgr_fromPaths))==0)  {
      print('You asked me to produce empty graph!')
      return(make_empty_graph())
    }
    #I didn`t find function that creates 'edge induced subgraph' to keep attributes
    V(subgr_fromPaths)$t= V(graph)$t[match(V(subgr_fromPaths)$name,V(graph)$name)]
    V(subgr_fromPaths)$type= V(graph)$type[match(V(subgr_fromPaths)$name,V(graph)$name)]
    V(subgr_fromPaths)$sigf= V(graph)$sigf[match(V(subgr_fromPaths)$name,V(graph)$name)]
    V(subgr_fromPaths)$prizes= V(graph)$prizes[match(V(subgr_fromPaths)$name,V(graph)$name)]
    V(subgr_fromPaths)$color= V(graph)$color[match(V(subgr_fromPaths)$name,V(graph)$name)]
    
    #the graph may not have 'commCol' or 'commColWSigf' attributes (e.g. from testcase)
    if(!is.null(V(graph)$commCol) & !is.null(V(graph)$commColWSigf)){
      V(subgr_fromPaths)$commCol= V(graph)$commCol[match(V(subgr_fromPaths)$name,V(graph)$name)]
      V(subgr_fromPaths)$commColWSigf= V(graph)$commColWSigf[match(V(subgr_fromPaths)$name,V(graph)$name)]
    }
    E(subgr_fromPaths)$name= paste(ends(subgr_fromPaths,es=E(subgr_fromPaths))[,1],
                                   ends(subgr_fromPaths,es=E(subgr_fromPaths))[,2],sep='_')
    #transferring edge costs from original graph is a bit involved..
    edge_cost_DF_orig= cbind.data.frame(ends(graph,es=E(graph)),costs=edge.attributes(graph)$costs)
    costTrans_FUN= function(x) {
      edge_cost_DF_orig$costs[edge_cost_DF_orig[,1]==x[1] & edge_cost_DF_orig[,2]==x[2]]}
    E(subgr_fromPaths)$costs= apply(ends(subgr_fromPaths,es=E(subgr_fromPaths)),1,costTrans_FUN)
    print( paste('Returning new subgraph with',length(V(subgr_fromPaths)),'nodes.',sep=' ') )
    
    return(subgr_fromPaths)
  }
#############

  locILP= './20221024_PCST_ILP/' 
  
  #prepare ILP to create Prize Collecting Steiner Tree
   #undirected graph
    #this implementation avoids singly-connected nodes that therefore do not
    #contribute to connecting terminal nodes, even if they contribute positive prize
    #it can produce unconnected subgraphs as long as each subgraph has at least two
    #terminals
  
  #node prizes and edge costs both need to be positive!
  
  if(testmodeAns== 'y') {
    ### testcases:
    library(pcSteiner)
    library(stringr)
    library(igraph)
    locILP= './20221024_PCST_ILP/20221209_testing/'
    source(paste(locILP,'testcases.R',sep=''))
    test_out= testcase_fun()
    g_orig= test_out[[1]]
    out_term_orig= test_out[[2]]
    in_term_orig= test_out[[3]]
    t1=1
    t2=2
    max_St_betwTerm= 4
    prizeToCost=1
    gSol_exist='n'
    FC_thres= NA
  } else if (testmodeAns== 'n') {
    ###
    #if not using testcases
    g_orig= dir_graph
    out_term_orig= out_term
    in_term_orig= in_term
    ###
  }
  
  V(g_orig)$color= rep(NA,length(V(g_orig)))
  V(g_orig)[out_term_orig]$color= 'yellow'
  V(g_orig)[in_term_orig]$color= 'orange'
    
  E(g_orig)$name= paste(ends(g_orig,es=E(g_orig))[,1],ends(g_orig,es=E(g_orig))[,2],sep='_')
  
  par(mar=c(0,0,2,0),cex.main=0.8)
  # plot(g_orig,edge.label= attributes(E(g_orig))$vnames,vertex.label.cex=0.5,
  #      edge.label.cex=0.5)
  edgeLabs_orig= unlist(lapply(rbind.data.frame(attributes(E(g_orig))$vnames,E(g_orig)$costs),FUN= 
                            function(x){paste(x,collapse='\n')} ))
  vertexLabs_orig= unlist(lapply(rbind.data.frame(V(g_orig)$name,V(g_orig)$prizes),FUN=
                              function(x){paste(x,collapse='\n')}) )
   # plot(g_orig,edge.label= edgeLabs_orig,vertex.label=vertexLabs_orig,vertex.label.cex=0.6,
   #      edge.label.cex= 0.6,edge.arrow.size=0.5,main='Original Graph') 
   # plot(g_orig,vertex.label.cex=0.5,edge.arrow.size=0.5,vertex.size=10,
   #    main='Original graph')  #w/o edge labels
  # set.seed(1)
  # lay= layout.fruchterman.reingold(g_orig)
  # frame_col= rep('black',length(V(g_orig)))
  # frame_col[V(g_orig)$sigf=='sigf']= 'red'
  # frame_width= rep(1,length(V(g_orig)))
  # frame_width[V(g_orig)$sigf=='sigf']= 2
  # print(
  # plot(g_orig,layout=lay,edge.arrow.size=0.5,vertex.size=4,vertex.label=NA,vertex.frame.color=frame_col,
  #      vertex.frame.width=frame_width, main='Original graph') ) #vertex.label.cex=0.5,   w/o vertex or edge labels
  # #layout = layout.fruchterman.reingold, edge.label= attributes(E(g_orig))$vnames,edge.label.cex=0.5
  # legend("bottomright", c("Source","Sink"), pch=21,
  #        col="#777777", pt.bg=c('yellow','orange'), pt.cex=2, cex=.8, bty="n", ncol=1) #x=-1, y=-0.8
  # legend("bottomleft", c("sigf."), pch=21, col='red', pt.cex=2, cex=.8, bty="n", ncol=1)
  #not yet: https://stackoverflow.com/questions/38999656/increasing-spaces-between-vertices-for-r-igraph
  
  ## INITIAL TRIMMING and possible division into subgraphs
  
  print('  Start trimming..')
  #find all nodes on directed paths between out- and in-termini
   #"simple path": path on which each vertex is visited only once
  g= GRAPH_FROM_SIMPPATHS(graph=g_orig,out_terms=out_term_orig,
                                 in_terms=in_term_orig,maxNoEdges=max_St_betwTerm+1)
  if(length(V(g))==0) {
    print("NO TREE FOUND!")
    return(list(make_empty_graph(),0))
  }
    
  # allowEdges= data.frame(matrix(NA,ncol=2,nrow=0))  #each entry contains nodes of allowed path
  # for(o in out_term_orig) {
  #   #print('Trimming to paths between terminals: Please wait!')
  #   currAllowPaths=
  #     all_simple_paths(graph=g_orig,from=o,to=in_term_orig,mode="out",cutoff=max_St_betwTerm+1) 
  #       #cutoff necessary for speed; cutoff is number of edges
  #   for (p in seq_along(currAllowPaths)) {
  #     aPn= currAllowPaths[[p]]$name
  #     #print(aPn)
  #     allowEdges= unique(rbind.data.frame(allowEdges,
  #                         cbind.data.frame(aPn[1:(length(aPn)-1)],aPn[2:length(aPn)])) )
  #   }
  # }
  # allowVertices= unique(c(allowEdges[,1],allowEdges[,2]))
  
  #remove nodes that are not on these paths
   #this also removes termini that are not on path from in-terminus to out-terminus
  # g= graph_from_data_frame(d=allowEdges,directed=TRUE,vertices=allowVertices)
  # V(g)$sigf= V(g_orig)$sigf[match(V(g)$name,V(g_orig)$name)]
  # V(g)$prizes= V(g_orig)$prizes[match(V(g)$name,V(g_orig)$name)]
  # V(g)$color= V(g_orig)$color[match(V(g)$name,V(g_orig)$name)]
  # #transferring edge costs from original graph is a bit involved..
  # edge_cost_DF_orig= cbind.data.frame(ends(g_orig,es=E(g_orig)),costs=edge.attributes(g_orig)$costs)
  # costTrans_FUN= function(x) {
  #   edge_cost_DF_orig$costs[edge_cost_DF_orig[,1]==x[1] & edge_cost_DF_orig[,2]==x[2]]}
  # E(g)$costs= apply(ends(g,es=E(g)),1,costTrans_FUN)
  print('  Trimming complete!')
  
  #edgeLabs= unlist(lapply(rbind.data.frame(attributes(E(g))$vnames,E(g)$costs),FUN= 
   #                 function(x){paste(x,collapse='\n')} ))
  edgeLabs= sapply(E(g)$costs,round,digits=1)
  # vertexLabs= unlist(lapply(rbind.data.frame(V(g)$name,V(g)$prizes),FUN=
  #                   function(x){paste(x,collapse='\n')}) )
  vertexLabs= unlist(lapply(rbind.data.frame(V(g)$name,round(V(g)$prizes,2)),FUN=
                              function(x){paste(x,collapse='\n')}) )
  # plot(g,edge.label= edgeLabs,vertex.label=vertexLabs,vertex.label.cex=0.6,
  #       edge.label.cex= 0.6,edge.arrow.size=0.5,edge.curved=1,main='Trimmed graph') 
  # plot(g,vertex.label=NA,vertex.size=5,edge.arrow.size=0.5,edge.curved=1,
  #      main='Trimmed graph') #w/o labels
  #lg= layout_as_tree(g)
  # plot(g,edge.label= NULL,vertex.label=NA,vertex.label.cex=0.6,vertex.size=5,
  #     edge.label.cex= 0.6,edge.arrow.size=0.5,edge.curved=1,main='Trimmed graph (no labels)')
  
  out_term= V(g)$name[V(g)$name %in% out_term_orig]
  in_term= V(g)$name[V(g)$name %in% in_term_orig]
  # distances(graph=g,v = out_term,to = in_term,mode='out') #shortest paths!!!
  # #distance_table(graph=g)  #shortest paths!!!
  # singNodeTrim= GRAPH_FROM_SIMPPATHS(graph=g,out_terms=out_term[1],in_terms=in_term[1],
  #                                    maxNoEdges=max_St_betwSliceTerm+1)
  # lay_tree= layout_as_tree(singNodeTrim, root= V(singNodeTrim)[V(singNodeTrim)$name==out_term[1]])
  # lay_tree[V(singNodeTrim)$name==in_term[1],2]=-1  #y-coordinate of in-term
  #some custom visualizations
  # firstNeighborsV= neighbors(graph=singNodeTrim,v='LAS17',mode='out') #first neighbors of out-term
  # V(singNodeTrim)$color[firstNeighborsV]='green'
  # secondNeighborsI= unlist(adjacent_vertices(graph=singNodeTrim,v=firstNeighborsV,mode='out'))
  # V(singNodeTrim)$color[secondNeighborsI]= 'red'
  # V(singNodeTrim)$color[V(singNodeTrim)$name=='LSB3']='blue'
  # V(singNodeTrim)$color[V(singNodeTrim)$name=='UBI4']='violet'
  # V(singNodeTrim)$color[V(singNodeTrim)$name=='RSP5']='pink'
  # plot(singNodeTrim,edge.label= NULL,layout=lay_tree,vertex.label=NA,vertex.label.cex=0.6,vertex.size=5,
  #      edge.label.cex= 0.6,edge.arrow.size=0.5,edge.curved=1,main='Trimmed graph (no labels)')
  
  if(length(out_term)==0 | length(in_term)==0) {
    print('EXITING: No more out- and/or in-terminals after trimming!')
    return(make_empty_graph(n=0,directed=TRUE)) #empty graph (directed!)
  }
  
  potSteinerN= V(g)$name[is.na(match(V(g)$name,c(out_term,in_term)))]
    ###### arrives here in reasonable time
  
  #root= 1
  #pcs.tree(graph=g,terminals=c(out_term,in_term),lambda=1,root=root,depth=5, eps=1e-3, max_iter=10)
  #tr= pcs.tree(graph=g,terminals=c(out_term,in_term),lambda=1,root=root,depth=5, eps=-1, max_iter=10)
  
  #create ILP formulation
  #https://stackoverflow.com/questions/2470248/write-lines-of-text-to-a-file-in-r/2470277#2470277
  if(testmodeAns== 'n'){
    out= paste(locILP,'Steiner1_dir_FC',FC_thres,'_PtoC',prizeToCost,'_evid',expEvid_cut,
             '_maxSt',max_St_betwTerm,'_',t1,'_',t2,'.lp',sep='')
  } else if (testmodeAns== 'y'){
    out= paste(locILP,'test',testString,'_FC',FC_thres,'_PtoC',prizeToCost,
               '_maxSt',max_St_betwTerm,'.lp',sep='')
  }
  sink(out)
  
  #Decision variables: edge-selected: x_ij; node-selected: y_j
  cat('\\This implementation avoids singly-connected nodes that therefore do not
   \\contribute to connecting terminal nodes, even if they contribute positive prize
   \\it can produce unconnected subgraphs as long as each subgraph has at least two terminals\n')
  cat('\\Decision variables: edge-selected: x_ij; node-selected: y_j')
  
  
  edg_vars= paste('x_',unlist(lapply(X=data.frame(
    t(matrix(as_edgelist(graph=g),ncol=2)),stringsAsFactors = FALSE),FUN=
                                 function(x){paste(x,collapse='_')})),sep='')
  C_vec= E(g)$costs
  node_vars= paste('y_',V(g)$name,sep='')
  P_vec= prizeToCost * V(g)$prizes
  cat('\n\nMinimize ',
    paste(paste(C_vec,edg_vars,sep=' '),collapse=' + '),' - ',
    paste(paste(P_vec,node_vars,sep=' '),collapse=' - '),sep='' )
    
  cat('\n\nsubject to')
  
   #out-terminals: need to have at least one outgoing edge
  out_neighIdx_list= adjacent_vertices(graph=g,
                        v= out_term,mode="out") #reference by name
  for(i in seq_along(out_neighIdx_list)){
    if(length(out_neighIdx_list[[i]])>0) {
      cat('\n',  paste(paste('x',out_term[i],(out_neighIdx_list[[i]])$name,sep='_'),
          collapse=' + '),' >= 1 \\edges from terminal',sep='')
    } else {
      stop(paste("ERROR: No potential connection from out-terminal",out_term[i])) #initial trimming should have taken care of this
    }
  }
  
  in_neighIdx_list= adjacent_vertices(graph=g, v= in_term,mode="in")
  for(j in seq_along(in_neighIdx_list)){
    if(length(in_neighIdx_list[[j]])>0) {
      cat('\n',  paste(paste('x',(in_neighIdx_list[[j]])$name,in_term[j],sep='_'),
          collapse=' + '),' >= 1 \\edges to terminal',sep='')
    }else {
      stop(paste("ERROR: No potential connection to in-terminal",out_term[i])) #initial trimming should have taken care of this
    }
  }
  
  #terminal nodes are selected by default
  cat('\n')
  cat(paste(paste('y_',out_term,' >= 1',sep=''),' \\terminal',collapse='\n'))
  cat('\n')
  cat(paste(paste('y_',in_term,' >= 1',sep=''),' \\terminal',collapse='\n'))
  
  #out-edges: y_i <= sum over i (x_ij)
  for (i in potSteinerN){
    out_neigh= (neighbors(graph=g,v=i,mode="out"))$name
    if(length(out_neigh)>0) {
      cat('\ny_',i,' - ',paste(paste('x', i,out_neigh,sep='_'),collapse=' - '),
          ' <= 0 \\Steiner node out',sep='') 
    } else {
      stop(paste('No potential out-connection from node',i)) #initial trimming should have taken care of this
    }
  }
  #in-edges: y_j <= sum over j (x_ij)
  for (j in potSteinerN){
    in_neigh= (neighbors(graph=g,v=j,mode="in"))$name
    if(length(in_neigh)>0){
      cat('\ny_',j,' - ',paste(paste('x',in_neigh,j,sep='_'),collapse=' - '),
          ' <= 0 \\Steiner node in',sep='')
    } else {
      stop(paste('No potential in-connection to node',j)) #initial trimming should have taken care of this
    }
  }
  
  #connected nodes are selected (this is necessary even though node prizes are
   #positive, because allowing edge addition may allow addition of nodes that 
   #would not otherwise be connected)
  #y_i >= x_ij  for all j   => y_i - x_ij >=0 for all j  /outgoing
  #y_j >= x_ij  for all i   => y_j - x_ij >=0 for all i  /incoming
  for (i in V(g)$name){
    out_neigh= neighbors(graph=g,v=i,mode="out")
    for(j in out_neigh$name){
      cat('\ny_',i,' - x_',i,'_',j,' >= 0 \\select node if connected (outgoing)',sep='')
    }
  }
  for (j in V(g)$name){
    in_neigh= neighbors(graph=g,v=j,mode="in")
    for(i in in_neigh$name){
      cat('\ny_',j,' - x_',i,'_',j,' >= 0 \\select node if connected (incoming)',sep='')
    }
  }
  
  endTime= Sys.time()
  runTime= endTime-startTime
  # print(runTime)
  # if(attributes(runTime)$units== 'mins'){
  #   runtime_min= as.numeric(runTime)
  # } else if (attributes(runTime)$units== 'secs'){
  #   runTime_min= as.numeric(runTime)/60
  # } else {
  #   stop('Unrecognized runtime for Steiner')
  # }
  
  #avoid circles of up to 6 nodes (Subtour Elimination)
   #starting from base node, find all out-nodes; if any of the out nodes corresponds
    #to in-node, circle is formed
   #continue with out nodes of out nodes (level2 out-nodes); if any of them corresponds
    #to in-node, circle is formed
   #etc. for higher level out nodes
  #print(' Avoiding circles..')    #gets printed to file
  reduc_graph= g #originally contains all nodes
   for(base in V(reduc_graph)$name){
     out_l1= neighbors(graph=reduc_graph,v=base,mode='out')$name
     in_l1= neighbors(graph=reduc_graph,v=base,mode='in')$name
     for(out_l1_idx in which(out_l1 %in% in_l1)){
       cat('\n',  paste( paste('x',base,out_l1[out_l1_idx],sep='_'),
                         paste('x',out_l1[out_l1_idx],base,sep='_'),sep=' + '),' <= 1 \\avoid 2-node cycles',sep='')
     }
     for(o_l1 in out_l1[out_l1 != base]){
       out_l2= neighbors(graph=reduc_graph,v=o_l1,mode='out')$name
       for(out_l2_idx in which(out_l2 %in% in_l1)){
         cat('\n',  paste( paste('x',base,o_l1,sep='_'),
                            paste('x',o_l1,out_l2[out_l2_idx],sep='_'),
                            paste('x',out_l2[out_l2_idx],base,sep='_'),sep=' + '),' <= 2 \\avoid 3-node cycles',sep='')
       }
       for(o_l2 in out_l2[! out_l2 %in% c(base,o_l1)]){
         out_l3= neighbors(graph=reduc_graph,v=o_l2,mode='out')$name
         for(out_l3_idx in which(out_l3 %in% in_l1)){
           cat('\n',  paste( paste('x',base,o_l1,sep='_'),
                              paste('x',o_l1,o_l2,sep='_'),
                              paste('x',o_l2,out_l3[out_l3_idx],sep='_'),
                              paste('x',out_l3[out_l3_idx],base,sep='_'),sep=' + '),' <= 3 \\avoid 4-node cycles',sep='')
         }
         for(o_l3 in out_l3[! out_l3 %in% c(base,o_l1,o_l2)]){
           out_l4= neighbors(graph=reduc_graph,v=o_l3,mode='out')$name
           for(out_l4_idx in which(out_l4 %in% in_l1)){
             cat('\n',  paste( paste('x',base,o_l1,sep='_'),
                                paste('x',o_l1,o_l2,sep='_'),
                                paste('x',o_l2,o_l3,sep='_'),
                                paste('x',o_l3,out_l4[out_l4_idx],sep='_'),
                                paste('x',out_l4[out_l4_idx],base,sep='_'),sep=' + '),' <= 4 \\avoid 5-node cycles',sep='')
           }
           for(o_l4 in out_l4[! out_l4 %in% c(base,o_l1,o_l2,o_l3)]){
             out_l5= neighbors(graph=reduc_graph,v=o_l4,mode='out')$name
             for(out_l5_idx in which(out_l5 %in% in_l1)){
               cat('\n',  paste( paste('x',base,o_l1,sep='_'),
                                  paste('x',o_l1,o_l2,sep='_'),
                                  paste('x',o_l2,o_l3,sep='_'),
                                  paste('x',o_l3,o_l4,sep='_'),
                                  paste('x',o_l4,out_l5[out_l5_idx],sep='_'),
                                  paste('x',out_l5[out_l5_idx],base,sep='_'),sep=' + '),' <= 5 \\avoid 6-node cycles',sep='')
             }
           }
         }
       }
     }
     reduc_graph= delete.vertices(graph= reduc_graph,v=base)
   }
  #print(' Circle calculations complete!')   #gets printed to file
   
  cat('\n\nbinary')
  cat('\n')
  cat( paste('y_',V(g)$name,sep=''),sep='\n')
  out_neighIdx_list_all= adjacent_vertices(graph=g,v=V(g),mode="out")
  for(i in 1:length(out_neighIdx_list_all)){
    if(length(out_neighIdx_list_all[[i]])>0) {
      cat( paste('x',V(g)$name[i],(out_neighIdx_list_all[[i]])$name,sep='_'), sep='\n')
    }
  }   #you don`t need to do the same for incoming edges, as all edges are already listed
  
  
  cat('\nEnd')
  closeAllConnections()
  
##############
  #display ILP results from Gurobi resultfile  
  if(gSol_exist== 'y') {
    resFileName= paste('Steiner1_dir_FC',FC_thres,'_PtoC',prizeToCost,'_evid',expEvid_cut,
                       '_maxSt',max_St_betwTerm,'_',t1,'_',t2,'.sol',sep='')
  } else {
    resFileName= readline("Give filename of Gurobi resultfile: ")
    while(! str_detect(resFileName,'\\.sol$') ){
      resFileName= readline("Incorrect Gurobi result file! Try again:")
    }
  }
  
  res= read.table(paste(locILP,resFileName,sep=''),stringsAsFactors = FALSE)
  res_chosen= res[res[,2]==1,]
  chos_nod= vapply(res_chosen[,1],function(x){(str_match(x,'y_(.*)'))[2]},
                   FUN.VALUE = character(length= 1))
  (chos_nod= chos_nod[!is.na(chos_nod)])
  chos_edg= vapply(res_chosen[,1],function(x){(str_match(x,'x_(.*_.*)'))[2]},
                   FUN.VALUE = character(length= 1))
  (chos_edg= chos_edg[!is.na(chos_edg)])
  chos_edg_DF= data.frame(matrix(NA,nrow=0,ncol=2))
  for(i in seq_along(chos_edg)){
    chos_edg_DF[i,]= str_split(chos_edg[i],'_')[[1]]
  }
  # chos_edg_DF$sigf_from= chos_edg_DF$sigf[match(chos_edg_DF[,1],V(g)$name)]
  # chos_edg_DF$sigf_to= chos_edg_DF$sigf[match(chos_edg_DF[,2],V(g)$name)]
  # chos_edg_DF$t_from= chos_edg_DF$t[match(chos_edg_DF[,1],V(g)$name)]
  # chos_edg_DF$t_to= chos_edg_DF$t[match(chos_edg_DF[,2],V(g)$name)]
  # chos_edg_DF$type_from= chos_edg_DF$type[match(chos_edg_DF[,1],V(g)$name)]
  # chos_edg_DF$type_to= chos_edg_DF$type[match(chos_edg_DF[,2],V(g)$name)]
  
  V(g)$shape= rep("circle",length(V(g)))
  V(g)$shape[V(g)$name %in% chos_nod]= "square"
  E(g)$color= rep('grey',length(E(g)))
  E(g)$color[ unlist(lapply(data.frame(t(ends(graph=g, es=E(g),names=TRUE)),
              stringsAsFactors = FALSE),
                function(x){paste(x,collapse='_')}))  %in% chos_edg] ='red'
  
  #plot(g,vertex.shape=V(g)$shape)
  plot(g,vertex.shape=V(g)$shape,vertex.label.cex=0.4,edge.arrow.size=0.4)
  legend('bottomleft', c("not chosen","chosen"),col=c('grey','red'),lty=c(1,1))
  legend('topright', c("not chosen","chosen"),pch=c(1,0))
  
  g_chos= graph_from_data_frame(d= chos_edg_DF,directed= TRUE,
            vertices= chos_nod) #new graph only from chosen edges
  #distances(graph= g_chos,v='LAS17',to=c('NSR1','HOG1','SNU66'),mode='out') #t=5-10
  
  V(g_chos)$t= V(g)$t[match(V(g_chos)$name,V(g)$name)]
  V(g_chos)$type= V(g)$type[match(V(g_chos)$name,V(g)$name)]
  V(g_chos)$sigf= V(g)$sigf[match(V(g_chos)$name,V(g)$name)]
  V(g_chos)$prizes= V(g)$prizes[match(V(g_chos)$name,V(g)$name)]
  V(g_chos)$color= NA
  V(g_chos)$color[V(g_chos)$name %in% out_term]= 'yellow'
  V(g_chos)$color[V(g_chos)$name %in% in_term]= 'orange'
  if( !is.null(V(g)$commCol) & !is.null(V(g)$commColWSigf) ){
    V(g_chos)$commCol= V(g)$commCol[match(V(g_chos)$name,V(g)$name)]
    V(g_chos)$commColWSigf= V(g)$commColWSigf[match(V(g_chos)$name,V(g)$name)]
  }
  E(g_chos)$name= paste(ends(g_chos,es=E(g_chos))[,1],ends(g_chos,es=E(g_chos))[,2],sep='_')
  E(g_chos)$costs= E(g)$costs[match(E(g_chos)$name,E(g)$name)]
  # 
  # vLabs= apply(cbind.data.frame(V(g_chos)$name,V(g_chos)$prizes),1,paste,collapse='\n')
  # plot(g_chos,vertex.shape=V(g_chos)$shape,vertex.label.cex=0.4,edge.arrow.size=0.4,
  #      vertex.label=vLabs)
  
  return(list(g_chos,runTime))
}

#some helpful igraph functions
# get.vertex.attribute(g)
# V(g)$name
# incident / incident_edges
# neighbors / adjacent_vertices
#delete_vertices
# vertex.attributes

#https://stackoverflow.com/questions/53971847/tidying-an-igraph-plot

#change log:
#20221028a: #no longer using root node
            #test cases in separate file
            #automatic construction of E_mtx and C_mtx
#20221028b: #further testing
#20221028_PCST_ILP_DIR.R: for directed graph
#20221029:  #more functions from iGraph, instead of manual coding
            #x_12 -> x_1_2
            #y_i <=0 for non-terminal, non-Steiner nodes
#20221029b: test cases with non-numeric vertex names
#20221029c: 'select node if connected' restriction removed
#20221029d: remove early on non-terminal nodes that cannot be Steiner nodes
#20221031:  #use iGraph functions for initial network trimming
#20221101:  #'select node if connected' restriction added back in corrected form
#20221104:  #avoid cycles of up to 5 nodes
#20221105:  #display Gurobi results
#20221107:  #throw error when printing output from nodes that lack potential connections
             #(initial trimming should have taken care of this)
#20221229: #instead of trimming to all nodes that are on paths with the allowed number of Steiner nodes,
            #I now trim to the edges (and nodes) on these paths