# Network simulations

#Start from an underlying network; for each TP, randomly select terminals;
#draw a random number of edges between terminals at adjacent and random number
#of edges within the same TP
#Calculate sensitivity and specificity of network inference algorithm

## TODO:
 # Hht1 duplicated: OK?

library('igraph')

cat('\014')
rm(list=ls())
set.seed(19102023)

loc= './'
simulNo= 1  #simulation number

#base simulation on underlying network from '20230313_allSignif_STRING1exp_3Xcutoff_TexExDB0.7.tsv'
undNet= read.table('C:/Users/wrtlb/Dropbox/20221006_Kanshin_NetwInf/github/20231019_simulation/20230313_allSignif_STRING1exp_3Xcutoff_TexExDB0.7.tsv',
                    sep='\t',header=TRUE)
head(undNet)
dim(undNet)  #9324 edges
# 
g1= graph_from_data_frame(undNet[,1:2],directed=FALSE)  #undirected
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
V_list= list(V_TP1$name)

# Color vertices by the TP at which they are selected
V(g)$colr= rep(NA,length(V(g)))
pal= c(palette.colors()[2:9],'yellowgreen','green','purple','magenta','orange')
V(g)$colr[V(g) %in% V_TP1]= 'pink'
E(g)$colr= rep('black',length(E(g)))
#E(g)$numName= as.character(E(g)) #name that remains stable when edges are deleted

# sample vertices at the next TP from all neighbors of vertices at current TP
g_copy= g  #make a copy of g, so I can delete reverse edges of selected edges
V_prevTP= V_TP1
print("TP: 0")
print( paste( "   Selected", length(V_TP1), "nodes"))

# plotting is slow for large networks
#par(mar=c(0,0,3,0))
lay= layout.fruchterman.reingold(g_copy)
#plot(g_copy, layout=lay,vertex.label=V(g_copy)) #, edge.label= E(g_copy), edge.label= E(g_copy)$numName,

for(j in 2:length(TPs)){
  print(paste("TP:", TPs[j]))
  
  # plotting is slow for large networks
  # plot(g_copy, vertex.size=14, vertex.color= V(g_copy)$colr, vertex.label= NA,
    #   edge.arrow.size=.5, main=paste(j-1),layout=lay, edge.color=E(g_copy)$colr)  #, edge.label= E(g_copy), edge.label= E(g_copy)$numName, edge.arrow.width= 0.5, 
  
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
  V_list[[j]]= V(g)$name[V_currTP]
  
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
  # plotting is slow for large networks
  # plot(g_copy, vertex.size=14, vertex.color= V(g_copy)$colr, vertex.label=V(g_copy),
    #   edge.arrow.size=.5, main=paste(j),layout=lay, edge.color=E(g_copy)$colr, edge.label= E(g_copy))  #, edge.label= E(g_copy)$numName, edge.arrow.width= 0.5, 
  
  g_copy= delete_edges(g_copy, E(g_copy)[E(g_copy)$colr=='red'  ]) #delete all at once, otherwise edge numbers change if deleting one at a time
  
  print(paste("   Vertices:", paste(V_currTP, collapse = ',')))
  V_prevTP= V_currTP #current becomes previous
}

par(mar=c(0,0,0,0))
plot(g_copy, vertex.size=3,vertex.label= NA, vertex.color= V(g_copy)$colr,vertex.label = NA,
     edge.arrow.width= 0.5, edge.arrow.size=.5, layout= lay)  #, edge.label= E(g_copy), edge.label= E(g_copy)$numName
legend('bottomleft',legend= TPs, fill= pal[1:length(TPs)])


#----------- construct the output DF of selected vertices --------------------# 
  # contains gene name, ORF, pSite, log2FC, log10int
V_vec = unlist(V_list)
V_number= length(V_vec)
log2FC_mtx= matrix(0, nrow= V_number, ncol= length(TPs))
log10int_mtx= matrix(1, nrow= V_number, ncol= length(TPs))
for (i in 1:length(V_list)){
  curr_V= V_list[[i]]
  # assign '3', i.e. 8-fold change to signif. vertices
  log2FC_mtx[match(curr_V, V_vec), i] = rep(3, length(curr_V)) 
}
colnames(log2FC_mtx) = paste(rep('log2_T', length(TPs)), 
                             str_pad(TPs, 2, pad = "0"), sep = '')
colnames(log10int_mtx) = paste(rep('log10int_T', length(TPs)), 
                             str_pad(TPs, 2, pad = "0"), sep = '')

# look up ORF name
IDmap= read.csv(paste(loc,'20150727_UniProt_yeast_IDmapping.csv',sep=''), stringsAsFactors = FALSE)
ORF_vec= IDmap$ORF.Name[match(V_vec, IDmap$Gene.Name)]

outDF= cbind.data.frame(
  Gene= V_vec,
  ORF= ORF_vec,
  pSite = rep(1,length(V_vec)),
  log2FC_mtx,
  log10int_mtx
)

write.csv(outDF, 
          paste(loc,'20231019_simulation/simulTC_', 
                str_pad(simulNo, 3, pad = "0"), '.csv', sep=''),
          row.names = FALSE)
