## helper_functions_SS3.R
## Anne-Elise Nieblas (adapted from Sylvain Bonhommeau)
## 6/10/2017
## DESCRIPTION: This script is read by app.R for the IOTC_SS3 shiny. It loads the ss3 .nc model outputs, 
## melts the .nc arrays into data frames, adds/calculates new variables and assigns colnames.
## INPUTS: model, run, dir_save (where to find the .nc files)
## OUTPUTS: CPUE data frame (Fleet, Year, Exp, Obs,Dev,StDev);
##          LEN  data frame (Fleet, Year, Bin, Exp, Obs)
##          
## AEN TO DO:
## IMPROVEMENTS: Have to have different functions ?- can put them all together to load data once?
##               Not all dimensions are necessarily accounted for in column naming, e.g., SEASON!
##                make function to replace single dims.
##               ## MUST PARSE UNITS STRING INTO VALUES


# paste(dirData,'data_ss324_',input$Model,'_',input$Run,'.nc',sep='')


for_ncCPUE <- function(x, dir_save){
  nc            <-nc_open(paste(dir_save,x, sep=""))
  exp           <-ncvar_get(nc,'Exp_cpue')
  
  dims          <-list(nc$dim$fleet$vals,nc$dim$year$vals,nc$dim$season$vals)
  coldims       <-c('Fleet','Year','Season')
  dr            <-NULL
  for(d in 1:length(dims)){if(length(dims[[d]])==1){dr=c(dr,d)}}
  dims[dr]      <-NULL
  if(length(dr)>0){coldims       <-coldims[-dr]}
  dimnames(exp) <-dims
  
  CPUE          <-melt(exp)
  obs           <-ncvar_get(nc,'Obs_cpue')
  Obs           <-melt(obs)
  Dev           <-Obs$value-CPUE$value
  stDev         <-Dev/sd(Dev,na.rm=T)
  CPUE          <-cbind(CPUE,Obs$value,Dev,stDev)
  
  colnames(CPUE)<-c(coldims,'Exp','Obs','Dev','stDev')
  CPUE_data     <- CPUE
  return(CPUE_data)
}

read_CPUE_data <- function(model,run, dir_save){
  comb_run     <- expand.grid(model=model,run=run)
  data_runs    <- cbind(File=paste("data_ss324_", comb_run$model,"_", comb_run$run, ".nc", sep=""), comb_run)
  run_name     <- paste(run)
  Model_run    <- paste(model,run)
  CPUE_data    <- ddply(data_runs, .(File), function(x) data.frame(Model_run=paste(x$model, x$run),ldply(lapply(x$File, for_ncCPUE, dir_save))))[,-1]
  return(CPUE_data)
}

for_ncLEN <- function(x, dir_save){
  nc           <-nc_open(paste(dir_save,x, sep=""))
  Lobs         <-ncvar_get(nc,'Obs_len')
  
  dims         <-list(nc$dim$fleet$vals,nc$dim$year$vals,nc$dim$season$vals,nc$dim$lenbin$vals)
  coldims      <-c('Fleet','Year','Season','Bin')
  dr           <-NULL
  for(d in 1:length(dims)){if(length(dims[[d]])==1){dr=c(dr,d)}}
  dims[dr]    <-NULL
  if(length(dr)>0){coldims     <-coldims[-dr]}
  dimnames(Lobs)<-dims
  
  lobs        <-melt(Lobs)
  Lexp        <-ncvar_get(nc,'Exp_len')
  lexp        <-melt(Lexp)
  Len         <-cbind(lexp,lobs$value)
  colnames(Len)<-c(coldims,'Exp','Obs')
  
  LEN_data    <-Len
  return(LEN_data)
}

read_LEN_data <- function(model,run, dir_save){
  comb_run    <- expand.grid(model=model,run=run)
  data_runs   <- cbind(File=paste("data_ss324_", comb_run$model,"_", comb_run$run, ".nc", sep=""), comb_run)
  run_name    <- paste(run)
  Model_run   <- paste(model,run)
  LEN_data    <- ddply(data_runs, .(File), function(x) data.frame(Model_run=paste(x$model, x$run),ldply(lapply(x$File, for_ncLEN, dir_save))))[,-1]
  return(LEN_data)
}


for_ncRecDev  <- function(x, dir_save){
  nc          <-nc_open(paste(dir_save,x, sep=""))
  RD_val      <-ncvar_get(nc,'RecDev_val')
  
  dims        <-list(nc$dim$year$vals)
  dr          <-NULL
  for(d in 1:length(dims)){if(length(dims[[d]])==1){dr=c(dr,d)}}
  dims[dr]    <-NULL
  dimnames(RD_val)<-dims
  
  rd_val      <-melt(RD_val)
  # recruitment deviation lci (recdev$Value - 1.96 * recdev$Parm_StDev)
  RD_lci      <-ncvar_get(nc,'RecDev_lci')
  rd_lci      <-melt(RD_lci)
  ## recruitment deviation uci (recdev$Value + 1.96 * recdev$Parm_StDev)
  RD_uci      <-ncvar_get(nc,'RecDev_uci')
  rd_uci      <-melt(RD_uci)
  
  recdev      <-cbind(rd_val,rd_lci$value,rd_uci$value)
  colnames(recdev)<-c('Year','value','lci','uci')
  
  RecDev_data <- recdev
  return(RecDev_data)
}

read_RecDev_data <- function(model,run, dir_save){
  comb_run    <- expand.grid(model=model,run=run)
  data_runs   <- cbind(File=paste("data_ss324_", comb_run$model,"_", comb_run$run, ".nc", sep=""), comb_run)
  run_name    <- paste(run)
  Model_run   <- paste(model,run)
  RecDev_data <- ddply(data_runs, .(File), function(x) data.frame(Model_run=paste(x$model, x$run),ldply(lapply(x$File, for_ncRecDev, dir_save))))[,-1]
  return(RecDev_data)
}


for_ncR0 <- function(x, dir_save){
  nc          <-nc_open(paste(dir_save,x, sep=""))
  R0<-ncvar_get(nc,'Recruit_0')
  
  ## MUST PARSE UNITS STRING INTO VALUES - ERA
  xx=nc$dim$era$units
  xx.val=strsplit(xx,' ')
  xx.vals=strsplit(xx.val[[1:length(xx.val)]],':')
  era.vals<-NULL
  for(x in 1:length(xx.vals)){era.vals=c(era.vals,xx.vals[[x]][2])}
  
  dims        <-list(nc$dim$area$vals,nc$dim$year$vals,nc$dim$season$vals,era.vals)
  coldims     <-c('Area','Year','Season','Era')
  dr          <-NULL
  for(d in 1:length(dims)){if(length(dims[[d]])==1){dr=c(dr,d)}}
  dims[dr]    <-NULL
  if(length(dr)>0){coldims     <-coldims[-dr]}
  dimnames(R0)<-dims
  
  r0          <-melt(R0)
  colnames(r0)<-c(coldims,'rec0')
  r0$rec0     <-r0$rec0/1000
  R0_data     <-r0
  return(R0_data)
}

read_R0_data <- function(model,run, dir_save){
  comb_run    <- expand.grid(model=model,run=run)
  data_runs   <- cbind(File=paste("data_ss324_", comb_run$model,"_", comb_run$run, ".nc", sep=""), comb_run)
  run_name    <- paste(run)
  Model_run   <- paste(model,run)
  R0_data     <- ddply(data_runs, .(File), function(x) data.frame(Model_run=paste(x$model, x$run),ldply(lapply(x$File, for_ncR0, dir_save))))[,-1]
  return(R0_data)
}

for_ncSSB <- function(x, dir_save){
  nc          <- nc_open(paste(dir_save,x, sep=""))
  SSB         <- ncvar_get(nc,'SpawnBio')
  
  xx=nc$dim$era$units
  xx.val=strsplit(xx,' ')
  xx.vals=strsplit(xx.val[[1:length(xx.val)]],':')
  era.vals<-NULL
  for(x in 1:length(xx.vals)){era.vals=c(era.vals,xx.vals[[x]][2])}
  
  dims        <- list(nc$dim$area$vals,nc$dim$year$vals,nc$dim$season$vals,era.vals)
  coldims     <- c('Area','Year','Season','Era')
  dr          <- NULL
  for(d in 1:length(dims)){if(length(dims[[d]])==1){dr=c(dr,d)}}
  dims[dr]    <- NULL
  if(length(dr)>0){coldims     <- coldims[-dr]}
  dimnames(SSB)<-dims
  
  ssb         <- melt(SSB)
  colnames(ssb)<-c(coldims,'spbio')
  ssb$spbio   <- ssb$spbio/1000
  SSB_data    <- ssb
  return(SSB_data)
}

read_SSB_data <- function(model,run, dir_save){
  comb_run    <- expand.grid(model=model,run=run)
  data_runs   <- cbind(File=paste("data_ss324_", comb_run$model,"_", comb_run$run, ".nc", sep=""), comb_run)
  run_name    <- paste(run)
  Model_run   <- paste(model,run)
  SSB_data    <- ddply(data_runs, .(File), function(x) data.frame(Model_run=paste(x$model, x$run),ldply(lapply(x$File, for_ncSSB, dir_save))))[,-1]
  return(SSB_data)
}

for_ncBio <- function(x, dir_save){
  nc          <- nc_open(paste(dir_save,x, sep=""))
  Bioall      <- ncvar_get(nc,'Bio_all_ts')
  
  ## MUST PARSE UNITS STRING INTO VALUES - ERA
  xx=nc$dim$era$units
  xx.val=strsplit(xx,' ')
  xx.vals=strsplit(xx.val[[1:length(xx.val)]],':')
  era.vals<-NULL
  for(x in 1:length(xx.vals)){era.vals=c(era.vals,xx.vals[[x]][2])}
  
  
  dims        <- list(nc$dim$area$vals,nc$dim$year$vals,nc$dim$season$vals,era.vals)
  coldims     <- c('Area','Year','Season','Era')
  dr          <- NULL
  for(d in 1:length(dims)){if(length(dims[[d]])==1){dr=c(dr,d)}}
  dims[dr]    <- NULL
  if(length(dr)>0){coldims     <- coldims[-dr]}
  dimnames(Bioall)<-dims
  
  bioall      <- melt(Bioall)
  colnames(bioall)<-c(coldims,'bio')
  Bio_data    <- bioall
  return(Bio_data)
}

read_Bio_data <- function(model,run, dir_save){
  comb_run    <- expand.grid(model=model,run=run)
  data_runs   <- cbind(File=paste("data_ss324_", comb_run$model,"_", comb_run$run, ".nc", sep=""), comb_run)
  run_name    <- paste(run)
  Model_run   <- paste(model,run)
  Bio_data    <- ddply(data_runs, .(File), function(x) data.frame(Model_run=paste(x$model, x$run),ldply(lapply(x$File, for_ncBio, dir_save))))[,-1]
  return(Bio_data)
}

for_ncffmsy <- function(x, dir_save){
  nc          <- nc_open(paste(dir_save,x, sep=""))
  FFmsy       <- ncvar_get(nc,'FFmsy')
  
  dims        <- list(nc$dim$year$vals)
  dr          <- NULL
  for(d in 1:length(dims)){if(length(dims[[d]])==1){dr=c(dr,d)}}
  dims[dr]    <- NULL
  dimnames(FFmsy)<-dims
  
  ffmsy       <- melt(FFmsy)
  
  colnames(ffmsy)<-c('Year','FFmsy')
  ffmsy_data  <- ffmsy
  return(ffmsy_data)
}

read_ffmsy_data <- function(model,run, dir_save){
  comb_run    <- expand.grid(model=model,run=run)
  data_runs   <- cbind(File=paste("data_ss324_", comb_run$model,"_", comb_run$run, ".nc", sep=""), comb_run)
  run_name    <- paste(run)
  Model_run   <- paste(model,run)
  ffmsy_data  <- ddply(data_runs, .(File), function(x) data.frame(Model_run=paste(x$model, x$run),ldply(lapply(x$File, for_ncffmsy, dir_save))))[,-1]
  return(ffmsy_data)
}

for_ncSsel <- function(x, dir_save){
  nc          <- nc_open(paste(dir_save,x, sep=""))
  Ssel        <- ncvar_get(nc,'sizeselex')
  end_yr        <- ncatt_get(nc,0)
  end_yr<-end_yr$model_end_year
  
  ## MUST PARSE UNITS STRING INTO VALUES
  xx=nc$dim$size.factor$units
  xx.val=strsplit(xx,' ')
  xx.vals=strsplit(xx.val[[1:length(xx.val)]],':')
  size.vals<-NULL
  for(x in 1:length(xx.vals)){size.vals=c(size.vals,xx.vals[[x]][2])}
  
  dims        <- list(nc$dim$fleet$vals,nc$dim$year$vals,nc$dim$gender$vals,size.vals,nc$dim$sizeselex.size$vals)
  coldims     <- c('Fleet','Year','Gender','Factor','Size')
  dr          <- NULL
  for(d in 1:length(dims)){if(length(dims[[d]])==1){dr=c(dr,d)}}
  dims[dr]    <- NULL
  if(length(dr)>0){coldims     <- coldims[-dr]}
  dimnames(Ssel)<-dims
  
  ssel        <- melt(Ssel)
  ssel$endyr  <- end_yr
  colnames(ssel)<-c(coldims,'sel','endyr')
  Ssel_data   <- ssel
  return(Ssel_data)
}

read_Ssel_data <- function(model,run, dir_save){
  comb_run    <- expand.grid(model=model,run=run)
  data_runs   <- cbind(File=paste("data_ss324_", comb_run$model,"_", comb_run$run, ".nc", sep=""), comb_run)
  run_name    <- paste(run)
  Model_run   <- paste(model,run)
  Ssel_data   <- ddply(data_runs, .(File), function(x) data.frame(Model_run=paste(x$model, x$run),ldply(lapply(x$File, for_ncSsel, dir_save))))[,-1]
  return(Ssel_data)
}



######## sylvain's scripts ###########
for_data <- function(x, dir_save){
  load(paste(dir_save,x, sep=""))
  data_plot <- outputs_VPA$main
  data_plot$variable <- as.character(data_plot$variable)
  iwh2 <- which(is.na(data_plot$variable))
  if (length(iwh2)>0){data_plot$variable[iwh2] <- "Fplusgroup"} ## small patch due to data import for first get_runs_data function
  return(data_plot)
}

read_data <- function(run, seed_nb, dir_save){
  comb_run  <- expand.grid(run=run, seed_nb=seed_nb)
  data_runs <- cbind(File=paste("data_retro_", comb_run$run, "_", comb_run$seed_nb, ".Rdata", sep=""), comb_run)
  run_name  <- paste(run)
  seed_nb   <- paste(seed_nb)
  Run_seed  <- paste(run, seed_nb)
  data_plot <- ddply(data_runs, .(File), function(x) data.frame(Run_seed=paste(x$run, x$seed_nb),ldply(lapply(x$File, for_data, dir_save))))[,-1]
  return(data_plot)
}


plot_output_VPA <- function(data_runs, vari){
  data_sub <- data_runs[data_runs$variable==vari,]
  #paste(runs_id$path)
  h1 <- hPlot(value ~ Year , data = data_sub, type = "line", group = "Run_seed")
  h1$addParams(dom = vari)
  h1$chart(zoomType = "xy")
  #h1$save('essai.html', standalone=TRUE)
  return(h1)
}

# ##### RELOADS UPDATED ss3.24.R to INFRASTRUCTURE ####
# overwrite<-T #SET TO "TRUE" IF THE FILES ALREADY ON THE WORKSPACE SHOULD BE OVERWRITTEN
# outputs_WS <- paste("/Home",username,"Workspace/IOTC_SS3_Shiny/Rscripts",sep="/")
# listWS(outputs_WS) #GET THE LIST OF FILES AND FOLDERS IN ONE SUB-FOLDER
# helper_SS3=paste('/home/anne.elise.nieblas/SS3/IOTC_SS3_Shiny/Rscripts/helper_functions_SS3.R',sep='') # FILE WITH THE FUNCTION TO WRITE OGC 19115 metadata
# uploadWS(outputs_WS,helper_SS3,overwrite)
# #####################################################
