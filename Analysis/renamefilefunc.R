
getwd()

rm(list= ls())
setwd("../..")

folder_name <- list.files("~/Desktop/Lizards")
dir         <- paste0("~","/Desktop/Lizards/",folder_name, "/")

for(i in 1:length(folder_name)){
	liz       <- folder_name[i]
	directory <- dir[i]
	fileList  <- list.files(directory)
	
	setwd(directory)
		for(j in 1:length(fileList)){
			filedat <- read.table(fileList[j])
			write.table(filedat, file = paste0(liz, ".", fileList[j]), col.names = FALSE, row.names = FALSE, sep = "\t")
		}
}


