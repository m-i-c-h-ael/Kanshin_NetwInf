library(igraph)

####### undirected graphs #########
testcase1= function() {
  g <- graph('Bull')
  E(g)$costs <- c(3, 3, 3, 3, 3)
  V(g)$prizes <- c(10, 2, 2, 2, 2)
  term= c(4,5)
  
  return(list(g,term))
  #expected: obj = -6; all nodes selected; all edges except [2,3] selected
}

testcase2= function() {
  E= rbind.data.frame(c(1,2),c(1,3),c(1,4),c(2,5),c(4,6),c(5,6),c(5,7),c(6,8))
  colnames(E)= c('node1','node2')
  E$costs= c(1,2,2,3,1,3,2,2)
  # E$names= paste(E$node1,E$node2,sep="_")
  E
  g= graph_from_data_frame(d=E,directed=FALSE,vertices=1:8)
  V(g)$prizes= c(1,2,3,4,5,6,7,8)
  term= c(1,8)
  
  return(list(g,term))
  #expected: obj: -14; all nodes except 3 and 7 selected; 
  #all edges, except for those to nodes 3 and 7 selected
}

testcase3= function() {
  set.seed(1)
  E= rbind.data.frame(c(1,2),c(1,3),c(1,4),c(2,5),c(4,6),c(5,6),c(5,7),c(6,8))
  colnames(E)= c('node1','node2')
  E$costs= sample(1:10,8)
  # E$names= paste(E$node1,E$node2,sep="_")
  E
  g= graph_from_data_frame(d=E,directed=FALSE,vertices=1:8)
  V(g)$prizes= sample(1:10,8)
  term= sample(1:8,2)
  
  return(list(g,term))
  #expected: obj: -12; nodes 2,5,7 and edges between them selected
}

testcase4= function() {
  set.seed(2)
  E= rbind.data.frame(c(1,2),c(1,3),c(1,4),c(2,5),c(4,6),c(5,6),c(5,7),c(6,8))
  colnames(E)= c('node1','node2')
  E$costs= sample(1:10,8)
  # E$names= paste(E$node1,E$node2,sep="_")
  E
  g= graph_from_data_frame(d=E,directed=FALSE,vertices=1:8)
  V(g)$prizes= sample(1:10,8)
  numTerm= sample(1:8,1)
  term= sample(1:8,numTerm)
  
  return(list(g,term))
  #expected: obj: -15; all nodes selected; edges [1,3],[1,4],[2,5],[5,7] and
   #[6,8] selected
}

testcase5= function() {
  set.seed(5)
  E= rbind.data.frame(c(1,2),c(1,3),c(1,4),c(2,5),c(4,6),c(5,6),c(5,7),c(6,8))
  colnames(E)= c('node1','node2')
  E$costs= sample(1:10,8)
  # E$names= paste(E$node1,E$node2,sep="_")
  E
  g= graph_from_data_frame(d=E,directed=FALSE,vertices=1:8)
  V(g)$prizes= sample(1:10,8)
  numTerm= sample(1:8,1)
  term= sample(1:8,numTerm)
  
  return(list(g,term))
  #expected: obj: -7; all nodes except 7 selected; edges [1,3],[1,4],[2,5],[4,6]
   # and [6,8] selected
}

##### directed graphs ############
testcase2_dir= function() {
  E= rbind.data.frame(c(1,2),c(1,3),c(1,4),c(2,5),c(4,6),c(5,6),c(5,7),c(8,6))
  colnames(E)= c('node1','node2')
  E$costs= c(1,2,2,3,1,3,2,2)
  # E$names= paste(E$node1,E$node2,sep="_")
  E
  g= graph_from_data_frame(d=E,directed=TRUE,vertices=1:8)
  V(g)$prizes= c(1,2,3,4,5,6,7,8)
  out_term= 1
  in_term= 8  #make sure no terminal node is out-term and in-term at once
  
  return(list(g,out_term,in_term))
  #expected: fails, because node 8 has no incoming edge
  #expected: 'noTerm'
   #penalty of 10: obj: 0: no node selected
   #penalty of 0: obj: -28: all nodes chosen; all except 1 are terminals
   #penalty of 3: obj: -10: nodes 5,6,7,8 chosen; all terminals
}

testcase6_dir= function() {
  set.seed(6)
  E= rbind.data.frame(c(1,2),c(1,3),c(1,4),c(2,5),c(4,6),c(5,6),c(5,7),c(8,6))
  colnames(E)= c('node1','node2')
  E$costs= c(1,2,2,3,1,3,2,2)
  # E$names= paste(E$node1,E$node2,sep="_")
  E
  g= graph_from_data_frame(d=E,directed=TRUE,vertices=1:8)
  V(g)$prizes= c(1,2,3,4,5,6,7,8)
  out_term= sample(E[,1],1)
  in_term= sample(E[,2],1)  
  in_term= in_term[is.na(match(in_term,out_term))]#make sure no terminal node is out-term and in-term at once
  
  return(list(g,out_term,in_term))
  #expected: no valid solution
}

testcase7_dir= function() {
  set.seed(7)
  E= rbind.data.frame(c(1,2),c(1,4),c(2,5),c(4,6),c(5,6),c(5,7),c(8,6))
  colnames(E)= c('node1','node2')
  E$costs= c(1,2,2,1,3,2,2)
  E
  g= graph_from_data_frame(d=E,directed=TRUE,vertices=1:8)
  V(g)$prizes= c(1,2,3,4,5,6,7,8)
  out_term= 2
  in_term= 7 
  in_term= in_term[is.na(match(in_term,out_term))]#make sure no terminal node is out-term and in-term at once
  
  return(list(g,out_term,in_term))
  #node #3 is unconnected
  #expected: obj.val= -10, nodes 2,5,7 and edges [2,5] and [5,7] are selected
}


testcase8_dir= function() {
  set.seed(8)
  vert= c('A','B','C','D','E','F','G','H')
  E_idx= rbind.data.frame(c(1,2),c(1,4),c(2,5),c(4,6),c(5,6),c(5,7),c(8,6))
  E= cbind.data.frame(vert[E_idx[,1]],vert[E_idx[,2]])
  colnames(E)= c('node1','node2')
  E$costs= c(1,2,2,1,3,2,2)
  E
  g= graph_from_data_frame(d=E,directed=TRUE,vertices=vert)
  V(g)$prizes= c(1,2,3,4,5,6,7,8)
  out_term= 'B'
  in_term= 'G' 
  in_term= in_term[is.na(match(in_term,out_term))]#make sure no terminal node is out-term and in-term at once
  
  return(list(g,out_term,in_term))
  #node #C is unconnected; only E can be connected between B and G
  #expected: obj.val= -10, nodes B,E,G and edges [B,E] and [E,G] are selected
}


testcase9_dir= function() {
  set.seed(9)
  n= 50
  g= sample_pa(n,directed=TRUE)  #Barabasi-Albert scale-free
  E(g)$costs= sample(x=1:10,size=length(E(g)),replace=TRUE)
  V(g)$name= paste('a',V(g),sep='')
  V(g)$prizes= sample(x=1:10,size=length(V(g)),replace=TRUE)
  numb_out_term= sample(1:5,1)
  numb_in_term= sample(1:5,1)
  out_term= sample(V(g)$name,numb_out_term)
  in_term= sample(V(g)$name,numb_in_term)
  in_term= in_term[is.na(match(in_term,out_term))]#make sure no terminal node is out-term and in-term at once
  
  return(list(g,out_term,in_term))
  #expected: no solutions: no connections to in-terminals
}

testcase8b_dir= function() {
  #contains long circles that would be preferred due to low cost, but that should be removed in ILP
  set.seed(8)
  vert= c('A','B','C','D','E','F','G','H','I')
  E_idx= rbind.data.frame(c(2,1),c(1,4),c(2,5),c(4,6),c(5,6),c(5,7),c(8,6),
                          c(7,8),c(8,5),c(4,2),c(4,1),c(4,7),c(2,8),c(2,9),c(9,1))
  E= cbind.data.frame(vert[E_idx[,1]],vert[E_idx[,2]])
  colnames(E)= c('node1','node2')
  E$costs= c(1,2,200,1,3,2,2,1,1,1,1,100,100,2,2)
  E
  g= graph_from_data_frame(d=E,directed=TRUE,vertices=vert)
  V(g)$prizes= c(1,2,3,4,5,6,7,8,9)
  out_term= 'B'
  in_term= 'G' 
  in_term= in_term[is.na(match(in_term,out_term))]#make sure no terminal node is out-term and in-term at once
  
  return(list(g,out_term,in_term))
  #expected: obj.val: 81; B - H - E - G
}


testcase9b_dir= function() {
  #as 9a, but all edges in both directions
  set.seed(9)
  n= 50
  g1= sample_pa(n,directed=TRUE)  #Barabasi-Albert scale-free
  g1r= reverse_edges(g1)
  g= union(g1,g1r)
  E(g)$costs= sample(x=1:10,size=length(E(g)),replace=TRUE)
  V(g)$name= paste('a',V(g),sep='')
  V(g)$prizes= sample(x=1:10,size=length(V(g)),replace=TRUE)
  numb_out_term= sample(1:5,1)
  numb_in_term= sample(1:5,1)
  out_term= sample(V(g)$name,numb_out_term)
  in_term= sample(V(g)$name,numb_in_term)
  in_term= c(in_term,'a3')
  in_term= in_term[is.na(match(in_term,out_term))]#make sure no terminal node is out-term and in-term at once
  
  return(list(g,out_term,in_term))
  #expected: obj: -11: all nodes selected; all edges from source to sinks
  #expected 'noTerm':
   #penalty of 0: -144
   #penalty of 10: -32
}

testcase8a_dir= function() {
  #contains long circles that would be preferred due to low cost, but that should be removed in ILP
  set.seed(8)
  vert= c('A','B','C','D','E','F','G','H')
  E_idx= rbind.data.frame(c(2,1),c(1,4),c(2,5),c(4,6),c(5,6),c(5,7),c(8,6),
                          c(7,8),c(8,5),c(4,2),c(4,1),c(4,7),c(2,8))
  E= cbind.data.frame(vert[E_idx[,1]],vert[E_idx[,2]])
  colnames(E)= c('node1','node2')
  E$costs= c(1,2,200,1,3,2,2,1,1,1,1,100,100)
  E
  g= graph_from_data_frame(d=E,directed=TRUE,vertices=vert)
  V(g)$prizes= c(1,2,3,4,5,6,7,8)
  out_term= 'B'
  in_term= 'G' 
  in_term= in_term[is.na(match(in_term,out_term))]#make sure no terminal node is out-term and in-term at once
  
  return(list(g,out_term,in_term))
  #expected for "no term": 
   #penalty of 10: no nodes selected
   #penalty of 0: obj=-23; nodes A,B,D,E,G,H: all except H are terminals
   #penalty of 3: obj=-12; nodes E,G,H: G and E are terminals
}


testcase9c_dir= function() {          #like 9a, but with further attributes
  set.seed(9)
  n= 50
  g= sample_pa(n,directed=TRUE)  #Barabasi-Albert scale-free
  E(g)$costs= sample(x=1:10,size=length(E(g)),replace=TRUE)
  V(g)$name= paste('a',V(g),sep='')
  V(g)$prizes= sample(x=1:10,size=length(V(g)),replace=TRUE)
  numb_out_term= sample(1:5,1)
  numb_in_term= sample(1:5,1)
  V(g)$type= rep('Steiner',length(V(g)))
  V(g)$type[sample(seq_along(V(g)),numb_out_term)]= 'out-term'
  V(g)$type[sample(seq_along(V(g)),numb_in_term)]= 'in-term'
  V(g)$type[V(g)$name=='a3']= 'in-term'
  V(g)$t= rep(15,length(V(g)))
  V(g)$t[V(g)$type=='out-term']=10
  V(g)$t[V(g)$type=='in-term']=20
  V(g)$sigf= rep('n',length(V(g)))
  V(g)$sigf[V(g)$type=='out-term' | V(g)$type=='in-term']= TRUE
  
  vertex.attributes(g)
  
  out_term= V(g)$name[V(g)$type=='out-term']
  in_term= V(g)$name[V(g)$type=='in-term']
  return(list(g,out_term,in_term))
  #expected: obj: -13; all nodes and edges after trimming selected: a28,a45,a19,a23,a3
}

###########
#testcases for path-extraction: it is required that TPs are visited in right order, but you don`t need terminals
 #at each TP
testcase10_dir_multTP= function() {
  set.seed(1)
  #non-random edges
  out_term= c('START1','START2')
  in_term= c('END1','END2')
  S_vec= c('START1','SFP1','SKP1','SNU66','SPE4','END1')
  M_vec= c('START2','MMS1','MRS6','MTR2','END2')
  edgeDF1= cbind.data.frame( node_from= S_vec[1:5],node_to= S_vec[2:6])
  edgeDF2= cbind.data.frame(node_from= M_vec[1:4],node_to= M_vec[2:5] )
  
  nodeDF= read.csv('C:/Users/wrtlb/Desktop/20220909_Kanshin_local/allSignif.csv',stringsAsFactors = FALSE)
  colnames(nodeDF)= c('name','t')
  nodeDF= rbind.data.frame( nodeDF, 
                cbind.data.frame(name=c('START1','START2','END1','END2'),t=c(10,10,60,50)) )
  nodeDF$t[match(S_vec,nodeDF$name)]= seq(10,60,10)
  nodeDF$t[match(M_vec,nodeDF$name)]= seq(10,50,10)
  nodeDF$sigf= rep(NA,dim(nodeDF)[1])
  nodeDF$sigf[match(c('START1','START2','END1','END2'),nodeDF$name)]= TRUE
  nodeDF$type= rep("Steiner",dim(nodeDF)[1])
  nodeDF$type[match(c('START1','START2'),nodeDF$name)]= 'out-term'
  nodeDF$type[match(c('END1','END2'),nodeDF$name)]= 'in-term'
  nodeDF$prizes= sample(x=1:10,size=dim(nodeDF)[1],replace=TRUE)
  
  #random edges
  n_rand= 10
  sigf_freq= 0.1
  edgeDF_rand= cbind.data.frame(node_from=sample(nodeDF$name,n_rand),node_to=sample(nodeDF$name,n_rand) )
  nodeDF$sigf[sample(seq_along(nodeDF$sigf),round(dim(nodeDF)[1]*sigf_freq))]= TRUE  #random significance
  
  edgeDF= rbind.data.frame(edgeDF1,edgeDF2,edgeDF_rand)
  g= graph_from_data_frame(d=edgeDF,directed=TRUE,vertices=nodeDF)
  E(g)$costs= sample(x=1:10,size=dim(edgeDF)[1],replace=TRUE)
  vertex.attributes(g)
  edge.attributes(g)
  return(list(g,out_term,in_term))
  #expected: both S_vec and M_vec are used as paths (prize-cost = 34)
  }

testcase10_dir_2TP= function() {
  set.seed(1)
  #non-random edges
  out_term= c('START1','START2')
  in_term= c('END1','END2')
  S_vec= c('START1','SFP1','SKP1','SNU66','SPE4','END1')
  M_vec= c('START2','MMS1','MRS6','MTR2','END2')
  edgeDF1= cbind.data.frame( node_from= S_vec[1:5],node_to= S_vec[2:6])
  edgeDF2= cbind.data.frame(node_from= M_vec[1:4],node_to= M_vec[2:5] )
  
  nodeDF= read.csv('C:/Users/wrtlb/Desktop/20220909_Kanshin_local/allSignif.csv',stringsAsFactors = FALSE)
  colnames(nodeDF)= c('name','t')
  nodeDF$t= rep(15,dim(nodeDF)[1])
  nodeDF= rbind.data.frame( nodeDF, 
                            cbind.data.frame(name=c('START1','START2','END1','END2'),t=c(10,10,20,20)) )
  nodeDF$sigf= rep(NA,dim(nodeDF)[1])
  nodeDF$sigf[match(c('START1','START2','END1','END2'),nodeDF$name)]= TRUE
  nodeDF$type= rep("Steiner",dim(nodeDF)[1])
  nodeDF$type[match(c('START1','START2'),nodeDF$name)]= 'out-term'
  nodeDF$type[match(c('END1','END2'),nodeDF$name)]= 'in-term'
  nodeDF$prizes= sample(x=1:10,size=dim(nodeDF)[1],replace=TRUE)
  
  #random edges
  n_rand= 10
  sigf_freq= 0.1
  edgeDF_rand= cbind.data.frame(node_from=sample(nodeDF$name,n_rand),node_to=sample(nodeDF$name,n_rand) )
  nodeDF$sigf[sample(seq_along(nodeDF$sigf),round(dim(nodeDF)[1]*sigf_freq))]= TRUE  #random significance
  
  edgeDF= rbind.data.frame(edgeDF1,edgeDF2,edgeDF_rand)
  g= graph_from_data_frame(d=edgeDF,directed=TRUE,vertices=nodeDF)
  E(g)$costs= sample(x=1:10,size=dim(edgeDF)[1],replace=TRUE)
  vertex.attributes(g)
  edge.attributes(g)
  return(list(g,out_term,in_term))
  #expected: both S_vec and M_vec are used as paths (prize-cost= 34)
}


testcase11_dir_multTP= function() {
#join the paths
  set.seed(1)
  #non-random edges
  out_term= c('START1','START2')
  in_term= c('END1','END2')
  S_vec= c('START1','SFP1','SKP1','SNU66','SPE4','END1')
  M_vec= c('START2','MMS1','MRS6','MTR2','END2')
  edgeDF1= cbind.data.frame( node_from= S_vec[1:5],node_to= S_vec[2:6])
  edgeDF2= cbind.data.frame(node_from= c(M_vec[1:4],'SNU66'),node_to= c(M_vec[2:5],'END2') )
  
  nodeDF= read.csv('C:/Users/wrtlb/Desktop/20220909_Kanshin_local/allSignif.csv',stringsAsFactors = FALSE)
  colnames(nodeDF)= c('name','t')
  nodeDF= rbind.data.frame( nodeDF, 
                            cbind.data.frame(name=c('START1','START2','END1','END2'),t=c(10,10,60,50)) )
  nodeDF$t[match(S_vec,nodeDF$name)]= seq(10,60,10)
  nodeDF$t[match(M_vec,nodeDF$name)]= seq(10,50,10)
  nodeDF$sigf= rep(NA,dim(nodeDF)[1])
  nodeDF$sigf[match(c('START1','START2','END1','END2'),nodeDF$name)]= TRUE
  nodeDF$type= rep("Steiner",dim(nodeDF)[1])
  nodeDF$type[match(c('START1','START2'),nodeDF$name)]= 'out-term'
  nodeDF$type[match(c('END1','END2'),nodeDF$name)]= 'in-term'
  nodeDF$prizes= sample(x=1:10,size=dim(nodeDF)[1],replace=TRUE)
  
  #random edges
  n_rand= 10
  sigf_freq= 0.1
  edgeDF_rand= cbind.data.frame(node_from=sample(nodeDF$name,n_rand),node_to=sample(nodeDF$name,n_rand) )
  nodeDF$sigf[sample(seq_along(nodeDF$sigf),round(dim(nodeDF)[1]*sigf_freq))]= TRUE  #random significance
  
  edgeDF= rbind.data.frame(edgeDF1,edgeDF2,edgeDF_rand)
  g= graph_from_data_frame(d=edgeDF,directed=TRUE,vertices=nodeDF)
  E(g)$costs= sample(x=1:10,size=dim(edgeDF)[1],replace=TRUE)
  vertex.attributes(g)
  edge.attributes(g)
  return(list(g,out_term,in_term))
  
  #expected: START1 to END1, START2 to END2 (prize-cost= 34)
}


testcase12_dir_multTP= function() {
  #join the paths; increase cost of MMS1-MRS6 edge
  set.seed(1)
  #non-random edges
  out_term= c('START1','START2')
  in_term= c('END1','END2')
  S_vec= c('START1','SFP1','SKP1','SNU66','SPE4','END1')
  M_vec= c('START2','MMS1','MRS6','MTR2','END2')
  edgeDF1= cbind.data.frame( node_from= S_vec[1:5],node_to= S_vec[2:6])
  edgeDF2= cbind.data.frame(node_from= c(M_vec[1:4],'SNU66'),node_to= c(M_vec[2:5],'END2') )
  
  nodeDF= read.csv('C:/Users/wrtlb/Desktop/20220909_Kanshin_local/allSignif.csv',stringsAsFactors = FALSE)
  colnames(nodeDF)= c('name','t')
  nodeDF= rbind.data.frame( nodeDF, 
                            cbind.data.frame(name=c('START1','START2','END1','END2'),t=c(10,10,60,50)) )
  nodeDF$t[match(S_vec,nodeDF$name)]= seq(10,60,10)
  nodeDF$t[match(M_vec,nodeDF$name)]= seq(10,50,10)
  nodeDF$sigf= rep(NA,dim(nodeDF)[1])
  nodeDF$sigf[match(c('START1','START2','END1','END2'),nodeDF$name)]= TRUE
  nodeDF$type= rep("Steiner",dim(nodeDF)[1])
   #make signif nodes out-term
  nodeDF$type[nodeDF$sigf==TRUE]= 'out-term'
   #then overwrite for START and END nodes
  nodeDF$type[match(c('START1','START2'),nodeDF$name)]= 'out-term'
  nodeDF$type[match(c('END1','END2'),nodeDF$name)]= 'in-term'
  nodeDF$prizes= sample(x=1:10,size=dim(nodeDF)[1],replace=TRUE)
  
  #random edges
  n_rand= 10
  sigf_freq= 0.1
  edgeDF_rand= cbind.data.frame(node_from=sample(nodeDF$name,n_rand),node_to=sample(nodeDF$name,n_rand) )
  nodeDF$sigf[sample(seq_along(nodeDF$sigf),round(dim(nodeDF)[1]*sigf_freq))]= TRUE  #random significance
  
  edgeDF= rbind.data.frame(edgeDF1,edgeDF2,edgeDF_rand)
  g= graph_from_data_frame(d=edgeDF,directed=TRUE,vertices=nodeDF)
  E(g)$costs= sample(x=1:10,size=dim(edgeDF)[1],replace=TRUE)
  E(g)$costs[7]=100
  vertex.attributes(g)
  edge.attributes(g)
  
  plot(g,vertex.size=4)
  
  return(list(g,out_term,in_term))
  
  #expected: START1 to END1, START2 to END2: terminals need to be included (prize-cost= -65)
  #betweenness: START1, START2, END1, END2:0; MMS1, MTR2: 3; MRS6: 4;
                   #SFP1, SPE4: 4; SKP1, SNU66: 6
}


testcase13_dir_multTP= function() {
  set.seed(1)
  #join the paths
  TPs= seq(10,60,10)
  #non-random edges
  out_term= c('START1','START2')
  in_term= c('END1','END2')
  S_vec= c('START1','SFP1','SKP1','SNU66','SPE4','END1')
  M_vec= c('START2','MMS1','MRS6','MTR2','END2')
  edgeDF1= cbind.data.frame( node_from= S_vec[1:5],node_to= S_vec[2:6])
  edgeDF2= cbind.data.frame(node_from= c(M_vec[1:4],'SNU66'),node_to= c(M_vec[2:5],'END2') )
  
  nodeDF= read.csv('C:/Users/wrtlb/Desktop/20220909_Kanshin_local/allSignif.csv',stringsAsFactors = FALSE)
  colnames(nodeDF)= c('name','t')
  nodeDF= rbind.data.frame( nodeDF, 
                            cbind.data.frame(name=c('START1','START2','END1','END2'),t=c(10,10,60,50)) )
  nodeDF$t[match(S_vec,nodeDF$name)]= seq(10,60,10)
  nodeDF$t[match(M_vec,nodeDF$name)]= seq(10,50,10)
  nodeDF$sigf= rep(NA,dim(nodeDF)[1])
  nodeDF$sigf[match(c('START1','START2','END1','END2'),nodeDF$name)]= TRUE
  nodeDF$type= rep("Steiner",dim(nodeDF)[1])
  nodeDF$type[match(c('START1','START2'),nodeDF$name)]= 'out-term'
  nodeDF$type[match(c('END1','END2'),nodeDF$name)]= 'in-term'
  nodeDF$prizes= sample(x=1:10,size=dim(nodeDF)[1],replace=TRUE)
  
  #random edges
  n_rand= 100
  sigf_freq= 0.1
  rand_from= sample(nodeDF$name,n_rand,replace=TRUE)
  rand_to= c() #match 
  for(i in seq_along(rand_from)){
    TP_from= nodeDF$t[match(rand_from[i],nodeDF$name)]
    t_idx= which(TPs==TP_from)+1 #index in order of TPs; to-index is one above from-index
    toCand= nodeDF$name[which(nodeDF$t==TPs[t_idx])]
    if(length(toCand)>0){
      rand_to[i]=sample(toCand,1)
    } else {
      rand_to[i]= NA
    }
  }
  edgeDF_rand= cbind.data.frame(node_from=rand_from,node_to=rand_to)
  edgeDF_rand= edgeDF_rand[!is.na(edgeDF_rand$node_to) & edgeDF_rand$node_to != edgeDF_rand$node_from,]   #no self-edges
  
  nodeDF$sigf[sample(seq_along(nodeDF$sigf),round(dim(nodeDF)[1]*sigf_freq))]= TRUE  #random significance
  
  edgeDF= rbind.data.frame(edgeDF1,edgeDF2,edgeDF_rand)
  edgeDF= unique(edgeDF) #no duplicates
  g= graph_from_data_frame(d=edgeDF,directed=TRUE,vertices=nodeDF)
  E(g)$costs= sample(x=1:10,size=dim(edgeDF)[1],replace=TRUE)
  vertex.attributes(g)
  edge.attributes(g)
  cbind.data.frame(edgeDF,E(g)$costs)
  
  plot(g)
  
  return(list(g,out_term,in_term))
  
  #expected: START1-MMS1-MRS6-VPS53-END2 START2-SFP1-ART5-MTR2-SPE4-END1 =>
   #better prize-cost than 'predefined' paths over S-vec and M-vec (prize: 47)
}


testcase14_dir_multTP= function() {
  #compared to 'testcase13_dir_multTP', enforce use of paths through S-vec and M-vec by increasing prize of
   #nodes on paths -> prizes of 10 were not sufficient
  set.seed(1)
  TPs= seq(10,60,10)
  #non-random edges
  out_term= c('START1','START2')
  in_term= c('END1','END2')
  S_vec= c('START1','SFP1','SKP1','SNU66','SPE4','END1')
  M_vec= c('START2','MMS1','MRS6','MTR2','END2')
  edgeDF1= cbind.data.frame( node_from= S_vec[1:5],node_to= S_vec[2:6])
  edgeDF2= cbind.data.frame(node_from= c(M_vec[1:4],'SNU66'),node_to= c(M_vec[2:5],'END2') )
  
  nodeDF= read.csv('C:/Users/wrtlb/Desktop/20220909_Kanshin_local/allSignif.csv',stringsAsFactors = FALSE)
  colnames(nodeDF)= c('name','t')
  nodeDF= rbind.data.frame( nodeDF, 
                            cbind.data.frame(name=c('START1','START2','END1','END2'),t=c(10,10,60,50)) )
  nodeDF$t[match(S_vec,nodeDF$name)]= seq(10,60,10)
  nodeDF$t[match(M_vec,nodeDF$name)]= seq(10,50,10)
  nodeDF$sigf= rep(NA,dim(nodeDF)[1])
  nodeDF$sigf[match(c('START1','START2','END1','END2'),nodeDF$name)]= TRUE
  nodeDF$type= rep("Steiner",dim(nodeDF)[1])
  nodeDF$type[match(c('START1','START2'),nodeDF$name)]= 'out-term'
  nodeDF$type[match(c('END1','END2'),nodeDF$name)]= 'in-term'
  nodeDF$prizes= sample(x=1:10,size=dim(nodeDF)[1],replace=TRUE)
  nodeDF$prizes[match(c(S_vec,M_vec),nodeDF$name)]= rep(10,length(c(S_vec,M_vec)))
  
  #random edges
  n_rand= 100
  sigf_freq= 0.1
  rand_from= sample(nodeDF$name,n_rand,replace=TRUE)
  rand_to= c() #match 
  for(i in seq_along(rand_from)){
    TP_from= nodeDF$t[match(rand_from[i],nodeDF$name)]
    t_idx= which(TPs==TP_from)+1 #index in order of TPs; to-index is one above from-index
    toCand= nodeDF$name[which(nodeDF$t==TPs[t_idx])]
    if(length(toCand)>0){
      rand_to[i]=sample(toCand,1)
    } else {
      rand_to[i]= NA
    }
  }
  edgeDF_rand= cbind.data.frame(node_from=rand_from,node_to=rand_to)
  edgeDF_rand= edgeDF_rand[!is.na(edgeDF_rand$node_to) & edgeDF_rand$node_to != edgeDF_rand$node_from,]   #no self-edges
  
  nodeDF$sigf[sample(seq_along(nodeDF$sigf),round(dim(nodeDF)[1]*sigf_freq))]= TRUE  #random significance
  
  edgeDF= rbind.data.frame(edgeDF1,edgeDF2,edgeDF_rand)
  edgeDF= unique(edgeDF) #no duplicates
  g= graph_from_data_frame(d=edgeDF,directed=TRUE,vertices=nodeDF)
  E(g)$costs= sample(x=1:10,size=dim(edgeDF)[1],replace=TRUE)
  vertex.attributes(g)
  edge.attributes(g)
  
  plot(g,vertex.size=5)
  tkplot(g)
  
  return(list(g,out_term,in_term))
  
  #expected: START1-MMS1-MRS6-VPS53-END2; START2-SFP1-ART5-MTR2-SPE4-END1 (prize-cost= 67)
}


testcase15_dir= function() {
  #case with one long but cheap and one short, but expensive path; both going from same out-term to same in-term
   #to test application of max_Steiner limit
  
  out_term= c('START1')
  in_term= c('END1')
  S_vec= c('START1','SFP1','SFP2','SFP3','SFP4','SFP5','SFP6','SFP7','SFP8','SFP9','SFP10','END1')
  M_vec= c('START1','MFP1','MFP2','MFP3','END1')
  edgeDF1= cbind.data.frame( node_from= S_vec[1:11],node_to= S_vec[2:12])
  edgeDF2= cbind.data.frame(node_from= M_vec[1:4],node_to= M_vec[2:5] )
  edgeDF= rbind.data.frame(edgeDF1,edgeDF2)
  
  nodes= c('START1','SFP1','SFP2','SFP3','SFP4','SFP5','SFP6','SFP7','SFP8','SFP9','SFP10','END1','MFP1','MFP2','MFP3')
  TPs= c(1,rep(1.5,10),2,rep(1.5,3))
  nodeDF= cbind.data.frame(nodes,TPs)
  colnames(nodeDF)= c('name','t')
  nodeDF$sigf= c('sig',rep(FALSE,10),TRUE,rep(FALSE,3))
  nodeDF$type= c('out-term',rep('Steiner',10),'in-term',rep('Steiner',3))
  nodeDF$color= c('yellow',rep(NA,10),'orange',rep(NA,3))
  nodeDF$prizes= rep(1,dim(nodeDF)[1])
  
  g= graph_from_data_frame(d=edgeDF,directed=TRUE,vertices=nodeDF)
  E(g)$costs= c(rep(1,12),rep(100,3))
  vertex.attributes(g)
  edge.attributes(g)
  
  plot(g)
  
  return(list(g,out_term,in_term))
  
  #expected: START1-MFP1-MFP2-MFP3-END1 and if max. Steiner <10 and >4 (prize-cost= -296)
  #betweenness: START1, END1: 0; MFP1,MFP3: 3; MFP2: 4
}


testcase16_dir= function() {
  #even though every node is on a path with max. three Steiner nodes between 'out' and 'in', longer paths
  #can be chosen
  endsDF= rbind.data.frame(
    c('out','A'),c('A','B'),c('B','C'),c('C','D'),c('C','in'),c('D','in'),c('D','E'),c('E','F'),c('F','in'),
    c('out','D'))
  #(nodes= unique(c(endsDF[,1],endsDF[,2])) )
  g= graph_from_data_frame(d=endsDF,directed=TRUE)
  V(g)
  lay= matrix( c(2,5,
                 3,4,
                 3,3,
                 3,2,
                 1,4,
                 1,3,
                 1,2,
                 2,1),ncol=2,byrow=TRUE )
  #plot(g,layout= lay)
  V(g)$sigf= c(TRUE,rep(FALSE,6),TRUE)
  V(g)$prizes= rep(10,8)
  V(g)$type= c('out-term',rep('Steiner',6),'in-term')
  V(g)$t= rep(15,length(V(g)))
  V(g)$t[V(g)$type=='out-term']=10
  V(g)$t[V(g)$type=='in-term']=20
  V(g)$color= c('yellow',rep(NA,6),'orange')  
  E(g)$costs= c(rep(1,10))

  in_term= V(g)$name[V(g)$type=='in-term']
  out_term= V(g)$name[V(g)$type=='out-term']
  
  plot(g)
  
  return(list(g,out_term,in_term))
  #expected: at p2c=1: out-A-B-C-D-E-F-in
  #betweenness centrality inferred:   out   A   B   C   D   E   F  in -> 0   6  10  12  12  10   6   0 
  #betweenness centrality underlying: out   A   B   C   D   E   F  in -> 0   2   6   8   9   5   1   0 
}


testcase17_dir= function() {
  g1= make_full_graph(n=4,directed= TRUE)
  V(g1)$name= c('a1','b1','c1','d1')
  g2= make_full_graph(n=5,directed= TRUE)
  V(g2)$name= c('a2','b2','c2','d2','e2')
  (edgFram= cbind.data.frame(c('b2','i'),c('i','a1')))
  i_graph= graph_from_data_frame(edgFram)
  g= g1+g2+i_graph
  
  out_term= 'a2'
  in_term='d1'
  V(g)$sigf= rep(FALSE,length(V(g)))
  V(g)$sigf[match(c(out_term,in_term),V(g)$name)]= TRUE
  V(g)$prizes= rep(1,length(V(g)))
  V(g)$type= rep('Steiner',length(V(g)))
  V(g)$type[V(g)$name==in_term]= 'in-term'
  V(g)$type[V(g)$name==out_term]= 'out-term'
  V(g)$t= rep(15,length(V(g)))
  V(g)$t[V(g)$type=='out-term']=10
  V(g)$t[V(g)$type=='in-term']=20
  V(g)$color= rep(NA,length(V(g)))
  V(g)$color[V(g)$type=='out-term']='yellow'
  V(g)$color[V(g)$type=='in-term']='orange'
  E(g)$costs= rep(1,length(E(g)))
  E(g)$costs[c(21,26)]=0.9  #avoid ties between equivalent graphs
  
  plot(g)
  
  return(list(g,out_term,in_term))
  
  #expected:
     # if max. Steiner <5: no solution
     # if >=5 allowed Steiner: at p2c=2: a2 - c2 - b2 - i - a1 - c1 - d1
  #underlying betweenness: i,b2=20;a1=18; others: 0     (~l. 1572)
  #inferred betweenness: i=9;a1,b2=8;c1,c2:5; a2,d1:0   (~l. 1915)
}


testcase18_dir= function() {
  g1= make_full_graph(n=4,directed= TRUE)
  V(g1)$name= c('a1','b1','c1','d1')
  g2= make_full_graph(n=5,directed= TRUE)
  V(g2)$name= c('a2','b2','c2','d2','e2')
  (edgFram= cbind.data.frame(c('b2','i'),c('i','a1')))
  i_graph= graph_from_data_frame(edgFram)
  g= g1+g2+i_graph
  
  out_term= c('a2','b2')
  in_term= c('d1','i')
  V(g)$sigf= rep(FALSE,length(V(g)))
  V(g)$sigf[match(c(out_term,in_term),V(g)$name)]= TRUE
  V(g)$prizes= rep(1,length(V(g)))
  V(g)$type= rep('Steiner',length(V(g)))
  V(g)$type[V(g)$name %in% in_term]= 'in-term'
  V(g)$type[V(g)$name %in% out_term]= 'out-term'
  V(g)$t= rep(15,length(V(g)))
  V(g)$t[V(g)$type=='out-term']=10
  V(g)$t[V(g)$type=='in-term']=20
  V(g)$color= rep(NA,length(V(g)))
  V(g)$color[V(g)$type=='out-term']='yellow'
  V(g)$color[V(g)$type=='in-term']='orange'
  E(g)$costs= rep(1,length(E(g)))
  E(g)$costs[c(21,26)]=0.9  #avoid ties between equivalent graphs
  
  plot(g)
  
  return(list(g,out_term,in_term))
  
  #expected:  ### This testcase demonstrate that I need introduce time when node enters tree, instead of only when it is significant
  # if max. Steiner <2: no solution
  # if >=2 allowed Steiner: at p2c=2: a2 - c2 - b2 - i - a1 - c1 - d1
  #underlying betweenness: i,b2=20;a1=18; others: 0     (~l. 1572)
  #inferred betweenness: i=9;a1,b2=8;c1,c2:5; a2,d1:0   (~l. 1915)
}

testcase19_dir= function() {
  #testcase for a situation where an out-term at t1 needs to connect via another out-term at t1
    # (i.e. is signaling within a time-slice possible)
  d= data.frame( matrix(c('A','B',
                          'B','C'),ncol=2,byrow=TRUE))
  g= graph_from_data_frame(d,directed=TRUE)
  
  out_term= c('A','B')
  in_term= 'C'
  V(g)$sigf= rep(FALSE,length(V(g)))
  V(g)$sigf[match(c(out_term,in_term),V(g)$name)]= TRUE
  V(g)$prizes= rep(1,length(V(g)))
  V(g)$type= rep('Steiner',length(V(g)))
  V(g)$type[V(g)$name %in% in_term]= 'in-term'
  V(g)$type[V(g)$name %in% out_term]= 'out-term'
  V(g)$t= rep(15,length(V(g)))
  V(g)$t[V(g)$type=='out-term']=10
  V(g)$t[V(g)$type=='in-term']=20
  V(g)$color= rep(NA,length(V(g)))
  V(g)$color[V(g)$type=='out-term']='yellow'
  V(g)$color[V(g)$type=='in-term']='orange'
  E(g)$costs= rep(1,length(E(g)))
  
  plot(g)
  
  return(list(g,out_term,in_term))
  
  #expected: A - B - C, even though B is 'out-term'
}