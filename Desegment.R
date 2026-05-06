#Finding the closest cell to a cell in a cell-feature table

###improvements: 1-remove duplicated distance calculation 
###2- have it parallel
###3- add in the the nth number closest
###4- add radius search
fst.colsest<- function(cellID, data, R){
  #closest<- numeric(n)
  #dist<- numeric(n)
  X<-data$X_centroid
  Y<-data$Y_centroid
  #A<- data$Area
  a<- as.numeric(data[which(data$CellID==cellID), c(which(colnames(data)=="X_centroid"), which(colnames(data)=="Y_centroid"))])
  a.x<- a[1]
  a.y<- a[2]
  r<- 0.00001
  find.table<- data[which(data$Y_centroid < a.y+r & data$Y_centroid > a.y-r & data$X_centroid > a.x-r & data$X_centroid < a.x+r), c(which(colnames(data)=="CellID"),which(colnames(data)=="X_centroid"),which(colnames(data)=="Y_centroid"))]
    i<-1
  while(nrow(find.table)==1){
    r<- (i*R)
    find.table<- data[which(data$Y_centroid < a.y+r & data$Y_centroid > a.y-r & data$X_centroid > a.x-r & data$X_centroid < a.x+r), c(which(colnames(data)=="CellID"),which(colnames(data)=="X_centroid"),which(colnames(data)=="Y_centroid"))]
    i<- 1+i
  }
  dist<- numeric(dim(find.table)[1])
  for (j in 1:dim(find.table)[1]){
    dist[j] <-  sqrt(((find.table[j,2])-a[1])^2+((find.table[j,3])-a[2])^2)    
  }
  
  closest<-find.table[which(dist==sort(dist)[2]), 1]
  distance<- sort(dist)[2]
  
  out<- cbind(closest, distance)
  colnames(out)<- c("CellID", "Distance")
  return(out)
}



##using the function

N_cell<- numeric(dim(G)[1])
N_distance<- numeric(dim(G)[1])


for(c in G$CellID[82290:dim(G)[1]]){
  i<-which(G$CellID==c)
  #if(N_cell[i]==0){
    N_cell[i]<-as.numeric(fst.colsest(c, G, 5)[,1])
    #N_cell[which(G$CellID==N_cell[i])]<- which(N_cell==N_cell[i])
    N_distance[i]<- as.numeric(fst.colsest(c, G, 20)[,2])
    #N_distance[which(G$CellID==N_distance[i])]<- which(N_distance==N_distance[i])
  #}
}
G<-cbind(G, N_cell, N_distance)

