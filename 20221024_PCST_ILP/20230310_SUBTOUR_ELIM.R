SUBTOUR_ELIM= function(g) {
  
  #Eliminate subtours (circles) up to a certain size, n (e.g. 6), via a brute force algorithm:
    #for all combinations of n of nodes, list the potential edges and make sure maximally
    # n-1 are selected

  circSize= 6 #size of circle to avoid
  
  out_idx= which(V(g)$type== 'out-term')
  in_idx= which(V(g)$type== 'in-term')
  
  for(curr_o in out_idx) {  #select one out-term and one in-term and all other nodes
    for(curr_i in in_idx) {
      
      otherNodes= V(g)$name[-c(curr_o,curr_i)]
      
      allCombs_woOthers= combn(otherNodes,(circSize-2) )  #all combinations
      allCombs= rbind.data.frame(
        rep(V(g)$name[curr_o],dim(allCombs_woOthers)[2]),
        allCombs_woOthers,
        rep(V(g)$name[curr_i],dim(allCombs_woOthers)[2])
      )
      dim(allCombs)
      
      for(colNo in 1:dim(allCombs)[2]){
        sG= subgraph(graph=g, vids= allCombs[,colNo])
        endsDF= ends(graph=sG, es=E(sG))
        
        print(
          paste(
            paste( paste('y',apply(endsDF,1,paste,collapse='_'),sep='_'), collapse=' + '),
            ' < ',circSize-1,' \avoid circle of size ',circSize,sep='')   #Gurobi takes '<' as '<='  
        )
      }
    }
  }
}