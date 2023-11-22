# Network simulations

#Start from an underlying network; for each TP, randomly select terminals;
#draw a random number of edges between terminals at adjacent and random number
#of edges within the same TP
#Calculate sensitivity and specificity

library('igraph')

cat('\014')
rm(list=ls())
set.seed(19102023)

# only for testing: get graphs from testcases
source('C:/Users/wrtlb/Dropbox/20221006_Kanshin_NetwInf/github/20221024_PCST_ILP/20221209_testing/testcases.R')
g1= testcase14_dir_multTP()[[1]]  #testcase16_dir()[[1]]
g1= as.undirected(g1)  #first make undirected, so that then each edge in both directions

#base simulation on underlying network from '20230313_allSignif_STRING1exp_3Xcutoff_TexExDB0.7.tsv'
# undNet= read.table('C:/Users/wrtlb/Dropbox/20221006_Kanshin_NetwInf/github/20231019_simulation/20230313_allSignif_STRING1exp_3Xcutoff_TexExDB0.7.tsv',
#                    sep='\t',header=TRUE)
# head(undNet)
# dim(undNet)  #9324 edges
# 
# g1= graph_from_data_frame(undNet[,1:2],directed=FALSE)  #undirected
V(g1)  #620 vertices

#make it directed
g= as.directed(g1, mode='mutual')

#simulation
TPs= seq(0,60,5)
meanVmin1_perTP= 2  #poisson lambda-1 of edges significant per TP (I add 1 to each lambda draw, so that at least one vertex is selected)
meanE_betwTP= 5
meanE_withinTP= 1

#sample signif. vertices at TP1 using poisson distribution with selected mean
V_TP1= sample(V(g),min(length(V(g)),rpois(n=1,lambda=meanVmin1_perTP)+1),replace=FALSE)
V_list= list(V_TP1)

# Color vertices by the TP at which they are selected
V(g)$colr= rep(NA,length(V(g)))
pal= c(palette.colors()[2:9],'yellowgreen','green','purple','magenta')
V(g)$colr[V(g) %in% V_TP1]= 'pink'
E(g)$colr= rep('black',length(E(g)))
#E(g)$numName= as.character(E(g)) #name that remains stable when edges are deleted

# sample vertices at the next TP from all neighbors of vertices at current TP
g_copy= g  #make a copy of g, so I can delete reverse edges of selected edges
V_prevTP= V_TP1
par(mar=c(0,0,3,0))
lay= layout.fruchterman.reingold(g_copy)
plot(g_copy, layout=lay,vertex.label=V(g_copy), edge.label= E(g_copy)) #edge.label= E(g_copy)$numName,
for(j in 2:length(TPs)){
  plot(g_copy, vertex.size=14, vertex.color= V(g_copy)$colr, vertex.label=V(g_copy),
       edge.arrow.size=.5, main=paste(j-1),layout=lay, edge.color=E(g_copy)$colr, edge.label= E(g_copy))  #, edge.label= E(g_copy)$numName, edge.arrow.width= 0.5, 
  
  potE_currTP= c()
  for(k in V_prevTP){
    #print(k)
    #List the potential edges attached to k: outgoing
    potE_currTP= c(potE_currTP, incident(graph=g_copy, v=k, mode='out') )  #potential / candidate vertices
  }
  
  E_currTP= sample(potE_currTP, 
                   min(length(potE_currTP),rpois(n=1,lambda=meanVmin1_perTP)+1),replace=FALSE)
  
  
  
  #List the vertices at the current TP
  V_currTP= unique(as.numeric(ends(graph=g_copy, E_currTP, names=FALSE )))
  V_currTP= V_currTP[!V_currTP %in% V_prevTP]  #remove vertices selected at previous TP
  V(g_copy)$colr[V_currTP]= pal[j]
  E(g_copy)$colr[E_currTP]= pal[j]
  
  #Delete the corresponding incoming edges from the graph, so that signal does not immediately flow back
  bothEnds= ends(graph=g_copy, E_currTP, names=FALSE )

  
  #mark edges one at a time, so that not all edges between sets of vertices are deleted
   #mark them first instead of deleting, so that deletion of one edge doesn`t change the number for the next
  #loop over elements of 'bothEnds'; 'bothEnds' may be empty
  for(i in seq_along(bothEnds[,1])){
    if(bothEnds[i,1] %in% V_currTP){
      E(g_copy)[ bothEnds[i,1] %->% bothEnds[i,2] ]$colr= 'red'
    } 
    if (bothEnds[i,2] %in% V_currTP){   #choose if instead of 'else if', so that also edges between two 'V_currTP' vertices are deleted
      E(g_copy)[ bothEnds[i,2] %->% bothEnds[i,1] ]$colr= 'red'
    }
  }
  plot(g_copy, vertex.size=14, vertex.color= V(g_copy)$colr, vertex.label=V(g_copy),
       edge.arrow.size=.5, main=paste(j),layout=lay, edge.color=E(g_copy)$colr, edge.label= E(g_copy))  #, edge.label= E(g_copy)$numName, edge.arrow.width= 0.5, 
  
  g_copy= delete_edges(g_copy, E(g_copy)[E(g_copy)$colr=='red'  ]) #delete all at once, otherwise edge numbers change if deleting one at a time
  
  
  V_prevTP= V_currTP #current becomes previous
}

par(mar=c(0,0,0,0))
plot(g_copy, vertex.size=3,vertex.label= NA, vertex.color= V(g_copy)$colr,vertex.label=V(g_copy),
     edge.arrow.width= 0.5, edge.arrow.size=.5, edge.label= E(g_copy), layout= lay)  #, edge.label= E(g_copy)$numName
legend('bottomleft',legend= TPs, fill= pal[1:length(TPs)])
