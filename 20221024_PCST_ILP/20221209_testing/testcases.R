#Testcases are designed to test the creation of trees on a single timeslice
 #testing of the complete workflow, including stiching together slices is not
 #implemented (as of 11.11.2023)

library(igraph)

assign_attrib= function(g, in_term, out_term) {
  #assign default attributes to vertices and edges
  
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
  
  return(g)
}

##### directed graphs ############
testcase2_dir= function() {
  E= rbind.data.frame(c(1,2),c(1,3),c(1,4),c(2,5),c(4,6),c(5,6),c(5,7),c(8,6))
  colnames(E)= c('node1','node2')
  E$costs= c(1,2,2,3,1,3,2,2)
  # E$names= paste(E$node1,E$node2,sep="_")
  E
  g= graph_from_data_frame(d=E,directed=TRUE,vertices=1:8)
  
  out_term= 1
  in_term= 8  #make sure no terminal node is out-term and in-term at once
  g= assign_attrib(g, in_term, out_term)
  V(g)$prizes= c(1,2,3,4,5,6,7,8)
  E(g)$costs= c(1,2,2,3,1,3,2,2)
  
  distances(g, v= out_term, to= in_term, weights= E(g)$costs, mode= 'out')
  par(mar=c(0,0,0,0))
  plot(g)
  
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
  
  out_term= sample(E[,1],1)
  in_term= sample(E[,2],1)  
  in_term= in_term[is.na(match(in_term,out_term))]#make sure no terminal node is out-term and in-term at once
  g= assign_attrib(g, in_term, out_term)
  V(g)$prizes= c(1,2,3,4,5,6,7,8)
  E(g)$costs= c(1,2,2,3,1,3,2,2)
  
  distances(g, v= out_term, to= in_term, weights= E(g)$costs, mode= 'out')
  par(mar=c(0,0,0,0))
  plot(g)
  
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
  
  out_term= 2
  in_term= 7 
  in_term= in_term[is.na(match(in_term,out_term))]#make sure no terminal node is out-term and in-term at once
  g= assign_attrib(g, in_term, out_term)
  V(g)$prizes= c(1,2,3,4,5,6,7,8)
  E(g)$costs= c(1,2,2,1,3,2,2)
  
  distances(g, v= out_term, to= in_term, weights= E(g)$costs, mode= 'out')
  shortest_paths(g, from= out_term, to= in_term, weights= E(g)$costs, mode= 'out')
  par(mar=c(0,0,0,0))
  plot(g)
  
  return(list(g,out_term,in_term))
  #node #3 is unconnected
  #expected: obj.val= -10, nodes 2,5,7 and edges [2,5] and [5,7] are selected
}  #passed 15/11/23


testcase8_dir= function() {
  set.seed(8)
  vert= c('A','B','C','D','E','F','G','H')
  E_idx= rbind.data.frame(c(1,2),c(1,4),c(2,5),c(4,6),c(5,6),c(5,7),c(8,6))
  E= cbind.data.frame(vert[E_idx[,1]],vert[E_idx[,2]])
  colnames(E)= c('node1','node2')
  E$costs= c(1,2,2,1,3,2,2)
  E
  g= graph_from_data_frame(d=E,directed=TRUE,vertices=vert)
  
  out_term= 'B'
  in_term= 'G' 
  in_term= in_term[is.na(match(in_term,out_term))]#make sure no terminal node is out-term and in-term at once
  g= assign_attrib(g, in_term, out_term)
  V(g)$prizes= c(1,2,3,4,5,6,7,8)
  E(g)$costs= c(1,2,2,1,3,2,2)
  
  shortest_paths(g, from= out_term, to= in_term, weights= E(g)$costs, mode= 'out')
  par(mar= c(0,0,0,0))
  plot(g)
  
  return(list(g,out_term,in_term))
  #node #C is unconnected; only E can be connected between B and G
  #expected: obj.val= -10, nodes B,E,G and edges [B,E] and [E,G] are selected
} #passed 15/11/23

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
  out_term= 'B'
  in_term= 'G' 
  in_term= in_term[is.na(match(in_term,out_term))]#make sure no terminal node is out-term and in-term at once
  
  g= assign_attrib(g, in_term, out_term)
  V(g)$prizes= c(1,2,3,4,5,6,7,8)
  E(g)$costs= c(1,2,200,1,3,2,2, 1,1,1,1,100,100)
  
  distances(g, v= out_term, to= in_term, weights= E(g)$costs, mode= 'out')
  
  # B - A - D -G
  d_BA= distances(g, v= 'B', to= 'A', weights= E(g)$costs, mode= 'out')
  d_AG= distances(g, v= 'A', to= 'G', weights= E(g)$costs, mode= 'out')
  p_BADG= sum( V(g)$prizes[match(c('B', 'A', 'D', 'G'), V(g)$name)] )
  (pc_BADG= p_BADG - d_BA - d_AG)
  
  d_BH= distances(g, v= 'B', to= 'H', weights= E(g)$costs, mode= 'out')
  d_HG= distances(g, v= 'H', to= 'G', weights= E(g)$costs, mode= 'out')
  p_BHEG= sum( V(g)$prizes[match(c('B', 'H', 'E', 'G'), V(g)$name)] )
  (pc_BHEG= p_BHEG - d_BH - d_HG)
  
  all_shortest_paths(g, from= out_term, to= in_term, weights= E(g)$costs, mode= 'out')
    # does not take into account prizes
  plot(g, layout= layout.fruchterman.reingold(g))
  
  return(list(g,out_term,in_term))
  #expected: B - H - E - G
} #15/11/23: new result - semi-verified

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
  g= assign_attrib(g, in_term, out_term)
  E(g)$costs= c(1,2,200,1,3,2,2,1,1,1,1,100,100,2,2)
  
  par(mar= c(0,0,0,0))
  plot(g)
  
  return(list(g,out_term,in_term))
  #expected: obj.val: 99; B - A - D - G
} # 15/11/23: new result; semi-verified

testcase9_dir= function() {
  set.seed(9)
  n= 50
  g= sample_pa(n,directed=TRUE)  #Barabasi-Albert scale-free
  V(g)$name= paste('a',V(g),sep='')
  
  numb_out_term= sample(1:5,1)
  numb_in_term= sample(1:5,1)
  out_term= sample(V(g)$name,numb_out_term)
  in_term= sample(V(g)$name,numb_in_term)
  in_term= in_term[is.na(match(in_term,out_term))]#make sure no terminal node is out-term and in-term at once
  g= assign_attrib(g, in_term, out_term)
  E(g)$costs= sample(x=1:10,size=length(E(g)),replace=TRUE)
  V(g)$prizes= sample(x=1:10,size=length(V(g)),replace=TRUE)
  
  return(list(g,out_term,in_term))
  #expected: no solutions: no connections to in-terminals
}  #passed 15/11/23

testcase9b_dir= function() {
  #as 9a, but all edges in both directions
  set.seed(9)
  n= 50
  g1= sample_pa(n,directed=TRUE)  #Barabasi-Albert scale-free
  g1r= reverse_edges(g1)
  g= graph.union(g1,g1r)
  V(g)$name= paste('a',V(g),sep='')
  
  numb_out_term= sample(1:5,1)
  numb_in_term= sample(1:5,1)
  out_term= sample(V(g)$name,numb_out_term)
  in_term= sample(V(g)$name,numb_in_term)
  in_term= c(in_term,'a3')
  in_term= in_term[is.na(match(in_term,out_term))]#make sure no terminal node is out-term and in-term at once
  
  g= assign_attrib(g, in_term, out_term)
  E(g)$costs= sample(x=1:10,size=length(E(g)),replace=TRUE)
  V(g)$prizes= sample(x=1:10,size=length(V(g)),replace=TRUE)
  
  par(mar=c(0,0,0,0))
  plot(g, edge.arrow.size=.5)
  
  return(list(g,out_term,in_term))
  #expected: obj: -11: all nodes selected; all edges from source to sinks
  # 9 paths of length 4 vertices; from a17/a11/a12 to a3/a7/a33 (not a42)
} # 15/11/23 changed result; seems OK

testcase9c_dir= function() {          #like 9a, but with further attributes
  set.seed(9)
  n= 50
  g= sample_pa(n,directed=TRUE)  #Barabasi-Albert scale-free
  
  V(g)$name= paste('a',V(g),sep='')
  numb_out_term= sample(1:5,1)
  numb_in_term= sample(1:5,1)
  out_term= sample(V(g)$name,numb_out_term)
  in_term= 'a3'

  g= assign_attrib(g, in_term, out_term)
  vertex.attributes(g)
  V(g)$prizes= sample(x=1:10,size=length(V(g)),replace=TRUE)
  E(g)$costs= sample(x=1:10,size=length(E(g)),replace=TRUE)
  
  distances(g, v= out_term, to= in_term, weights= E(g)$costs, mode= 'out')
  par(mar= c(0,0,0,0))
  plot(g)
  
  return(list(g,out_term,in_term))
  #expected: no path
} # passed 15/11/23


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
  nodeDF$sigf= c(TRUE,rep(FALSE,10),TRUE,rep(FALSE,3))
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
  g = assign_attrib(g, in_term, out_term)
  E(g)$costs[c(21,26)]=0.9  #avoid ties between equivalent graphs
  
  par(mar = c(0, 0, 0, 0))
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
  g= assign_attrib(g, in_term, out_term)
  E(g)$costs[c(21,26)]=0.9  #avoid ties between equivalent graphs
  
  par(mar=c(0,0,0,0))
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
  g= assign_attrib(g, in_term, out_term)
  
  plot(g)
  
  return(list(g,out_term,in_term))
  
  #expected: A - B - C, even though B is 'out-term'
}
