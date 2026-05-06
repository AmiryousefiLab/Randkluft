
setwd("/Users/bwilliams/git-hub/amiryousefi-lab/Creavasse")
df <- read.csv("data/scaled.csv")
column_index = c(1:4, 6:8, 10:12, 14:21, 23:24)
unique_id = unique(df$imageid)
names_columns = names(df[,column_index])
i=1



  par(mfrow=c(5,4), mar = c(5.1, 3, 4.1, 2))
  for (j in column_index){
    target<- df[df$imageid == unique_id[i], j]
    target<- remove_outliers(target)
    marker<- paste0(c(unique_id[i], names_columns[j] ),collapse="; ")
    output <- skew_gate(target, 0.8)
    output["maker"] = marker
    plot_gating(target, output, marker)
  }

