
library(dplyr)
library(rlang)
library(stringr)
# Set the working directory to the folder containing the expansion files

# Initialize vectors to store the X values and corresponding file names
x<-c()
y<-c()
index=1

# Loop through each file in the directory
for (file in list.files()){
  print(file)
  # Read the current file
    expanded<-read.csv(file,sep =',',header = T)
  # Filter for rows where X0 is "pomzp3"
  present<-expanded %>% dplyr::filter(expanded$X0 == "pomzp3")
  # Check if the 'present' dataframe is not empty
  if (!is_empty(present$X)){
  # Add the X value and file name to the vectors
  x[index]<-c(present$X)
  y[index]<-c(file)
  index=index+1
  # }
  }
}

# Remove the ".csv" extension from the file names
new_y<-str_remove(y, ".csv")

# Create a new dataframe with the X values and modified file names
data<-data.frame(X=x, Y=new_y)

# Write the new dataframe to a CSV file without row names
write.csv(data,file='Where_POMZP3.csv',row.names = F)

#POMZP3 <- read.csv('Where_POMZP3.csv') no POMZP3 for 99%, also no at 90% and neither at 80% 
