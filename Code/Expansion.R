
library(dplyr)

setwd("C:/Users/glori/Desktop/LBDM/Expansion_files/")

x<-c()
y<-c()
index=1

for (file in list.files()){
  print(file)
  expanded<-read.csv(file,sep =',',header = T)
  present<-expanded %>% dplyr::filter(expanded$X0 == "pomzp3")
  if (!is_empty(present$X)){
  x[index]<-c(present$X)
  y[index]<-c(file)
  index=index+1
  }
}

new_y<-str_remove(y, ".csv")

data<-data.frame(X=x, Y=new_y)

write.csv(data,file='Where_POMZP3.csv',row.names = F)


