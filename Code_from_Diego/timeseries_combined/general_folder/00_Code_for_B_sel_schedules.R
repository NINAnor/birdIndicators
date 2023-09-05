#We first generate path to where our data is storaged, and to were the schedule is
#POZOR!!! I'm using the backup folder to overwrite the schedule, be sure this is what you want
wd <- "P:/12580000_tov_fauna/TOV-E/2022_Combining_trends/Norway_combination/Norway_combination/"

setwd(paste0(wd, 'NorwayOld_converted'))

#path.TT.data="D:\\Prace_2022\\RSWAN-Combine\\Zoom\\working_folder\\"
#schedule.path="D:\\Prace_2022\\RSWAN-Combine\\Zoom\\general_folder\\"

path.TT.data = paste0(wd,"working_folder/")
schedule.path= paste0(wd,"general_folder/")


#Rember each species has 3 files of data associated, you just need one of them to avoid repetitions, so select indices_TT files
patron=(".*\\_indices_TT.csv")
vector_files=list.files(path.TT.data,pattern= patron)
vector_files#Now you have all the names of the files

#Generate a data set where you are going to store the each combination of species country
data_sp_country=as.data.frame(matrix(data= NA,ncol = 2,nrow=1))
names(data_sp_country)= c("Sp_code","Country_code")

#This loop selects only thosfiles at national level, and save the species code and the country code in its correct position
for(i in 1:length(vector_files)){
  x=strsplit(vector_files[i],"_")
  if(x[[1]][2]== 1){
    data_sp_country[i,1]=x[[1]][1]
    data_sp_country[i,2]=x[[1]][3]
  }
}

#Only use this code if the schedule has only an example of the levels and does not each species species-level combination
schedule.data=read.table(paste0(schedule.path,"Swan_schedules_comb.csv"),header = T,sep = ";",dec=",")
xx=vector(mode="character")
pp=unique(data_sp_country$Sp_code)
xx=rep(pp,each=nrow(schedule.data))

schedule.skeleton=as.data.frame(matrix(ncol=ncol(schedule.data),nrow=length(xx),NA))
names(schedule.skeleton)=names(schedule.data)
for(i in 1:ncol(schedule.data)){
  schedule.skeleton[,i]=rep(schedule.data[,i],length(pp))
}
schedule.skeleton$Species_nr=xx
#?rep

#Now is time to load the schedules files
#if you have used the previous code use the second lane

schedule.data=read.table(paste0(schedule.path,"Swan_schedules.csv"),header = T,sep = ";",dec=",")
schedule.data=schedule.skeleton
for (i in 1:nrow(schedule.data)) {
  sp.subset=subset(data_sp_country,data_sp_country$Sp_code==schedule.data$Species_nr[i])
  for (j in 1:nrow(sp.subset)) {
    if(schedule.data$Level1[i]==sp.subset$Country_code[j]){
      schedule.data$B_sel[i]=1
    }
    
  }
  
}
write.table(schedule.data, paste0(schedule.path, "Swan_schedules.csv"),row.names=F,col.names=T,sep=";",dec=".")

