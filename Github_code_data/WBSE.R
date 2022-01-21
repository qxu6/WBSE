# This code is to calculate burn severity, emissions, and daily emissions
# Author: Qingqing Xu
# Contact: qxu6@ucmerced.edu
# Last updated: 01/08/2022
# Citation: 

# ATTENTION! 
# 1) Need to first download CBI regression models from this GitHub repository: CBSE/Github_code_data/git_data
# 2) To generate data based on Monitoring Trends in Burn Severity (MTBS) shapefile and dNBR/NBR images, you should go to https://mtbs.gov/viewer/index.html to download wildfire data and save the data under the directory of this project: Data/MTBS
# 3) To generate data based on CAL FIRE shapefile and composite images generated with Google Earth Engine, you should contact the authors for access and save the dNBR images under the directory of this project: Data/CALFIRE
# 4) To calculate daily emissions, you should contact the authors for access and save the data to this folder: Data/rasterfile/dob_raster
# 5) User needs to update information in Section 2 to run the script based on input fire data 
# 6) Only run once of some sections of code to save processing time. To comment or uncomment multiple line, use the key combination "command/control + shift + c"
options(warn=1) # print out warnings

##---------------------------------------------------------------
##    Load packages, create functions, and define variables    --
##---------------------------------------------------------------
# 1. Load packages
fxns_pkg <- function(p){
  if(!is.element(p,installed.packages()[,1]))
  {install.packages(p, dep = T)}
  require(p,character.only = T)
} # a function to install and load required packages

sapply(
  c("raster","sf", "sp","rgdal","reshape2","foreign","dplyr","tidyverse","foreach",
    "data.table","stringr","ggmap","compare","purrr","ncdf4","parallel","doParallel",
    "maptools", "gridExtra","maps","rasterVis","docxtractr","magrittr"),
  fxns_pkg)

# 2. Create some functions for use to run this script
# (1) a function to get the full directories given a pattern
fxns_lfiles <- function(dir,pat){
  list.files(dir,pattern = pat, recursive = TRUE,full.names = TRUE)
} 

# (2) functions to use CBI models
fxns_para2 <- function(objs) {
  my_list <- list(a = cbi_model[which(cbi_model$level2 == objs),"A"],
                  b = cbi_model[which(cbi_model$level2 == objs),"B"],
                  c = cbi_model[which(cbi_model$level2 == objs),"C"])
  return(my_list) 
} # a function to get the a, b, c parameters of the CBI equation

fxns_sigB <- function(t,h) {
  t$a/(1.0+exp(-1.0*(h-t$b)/t$c))} # a function to apply the "pyeq3.Models_2D.Sigmoidal.SigmoidB" model

fxns_polyL <- function(t,h) {
  t$a+t$b*h} # a function to apply the "pyeq3.Models_2D.Polynomial.Linear" model

# (3) a function to calculate emissions
fxns_em <- function(vegR,sevR,grassE,shrubE,forestE){ #vegR is the vegetation classification raster; sevR is the severity classification raster; grassE, shrubE, and forestE are emission factors for grasslands, shrublands and forest
  overlay(vegR, sevR,fun=function(x, y) { ifelse ((x==1 & y %in% c(2:5)) | y==5, 4*fxns_emu*grassE, #x==1 means grass veg, y %in% c(2:5) means in low,mod, high or grass burn; y==5 means grass burn
                                                  ifelse (x==2 & y %in% c(2:4), 47*fxns_emu*shrubE, #x==2 means shrub, y %in% c(2:4) means in low or mod or high 
                                                          ifelse (x==3 & y ==2, 69*fxns_emu*forestE, #x==3 means Forest<5500', y ==2 means low severity
                                                                  ifelse (x==3 & y ==3, 82*fxns_emu*forestE, #x==3 means Forest<5500', y ==3 means mod severity
                                                                          ifelse (x==3 & y ==4, 89*fxns_emu*forestE, #x==3 means Forest<5500', y ==4 means high severity
                                                                                  ifelse (x==4 & y ==2, 75*fxns_emu*forestE, #x==4 means Forest 5500'-7500', y ==2 means low severity
                                                                                          ifelse (x==4 & y ==3, 92*fxns_emu*forestE, #x==4 means Forest 5500'-7500', y ==3 means mod severity
                                                                                                  ifelse (x==4 & y ==4, 101*fxns_emu*forestE,  #x==4 means Forest 5500'-7500', y ==4 means high severity
                                                                                                          ifelse (x==5 & y ==2, 45*fxns_emu*forestE,  #x==5 means Forest>7500', y ==2 means low severity
                                                                                                                  ifelse (x==5 & y ==3, 56*fxns_emu*forestE,#x==5 means Forest>7500', y ==3 means mod severity
                                                                                                                          ifelse (x==5 & y ==4, 62*fxns_emu*forestE, 0)))))))))))} ) #x==5 means Forest>7500', y ==4 means high severity
  
  
  
  
  
}

# (4) functions to convert units
fxns_mh<- function(x) as.numeric(x)*30*30/10000 # a function to convert 30m*30m area to hectares
# fuel consumed (t/ha)
# Emission Factors (g emitted/kg burned)
fxns_emu <- 30*30*0.0001*1000*0.001  # a function to convert the units into kg for final emission calculation results. 30m*30m*(square meter to hectare index 0.0001)*t to kg(1000)* (g to kg 0.001), so the final *unit will be kg

# 3. Set CBI burn severity threshold, need to ask Josh about vegetative CBI values
reclass_sev <- c(0,0.1, 1, # unchanged, the reference table set (0,0.1] for unchanged
                 0.1,1.25,2,   # low severity
                 1.25,2.25,3,  # moderate severity
                 2.25,3,4,     # high severity
                 3,4, 5)       # grass total burn, I set CBI value to be 3.5 to classify as grass in calculation
#CBI threshold reference: https://www.researchgate.net/publication/261559367_A_Mixed-Effects_Heterogeneous_Negative_Binomial_Model_for_Post-fire_Conifer_Regeneration_in_Northeastern_California_USA
class.matrix <- matrix(reclass_sev, ncol = 3, byrow = TRUE)

# 4. Define the crs that we will use for this analysis
mycrs <- "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs" # the same as MTBS dNBR CRS

##---------------------------------------------------------------
##          Section 0: Download data--just run one time!       --
##---------------------------------------------------------------
# # select all the codes below in this section and use the key combination "command/control + shift + c" to uncomment the codes in this section to to run for the first time

# # This section is to download a few publicly available GIS datasets using url for use in this project
# # 1. LANDFIRE Biophysical Settings (BPS)
# # 2. Monitoring Trends in Burn Severity (MTBS) Burned Areas Boundaries Dataset 
# # 3. Ecoregions of North America (Omernik)
# # 4. California state noundary
# # 5. CAL FIRE perimeters
# # 6. LANDFIRE Existing Vegetation Type
# options(timeout=1000) # The default time 60 seconds and that would cause Error in download.file "Timeout of 60 seconds was reached". So we set it to be 1000 seconds
# 
# fxns_dzip <- function(url,dir){
#   temp=tempfile()
#   download.file(url, temp)
#   unzip(temp,exdir=dir)
# } #a function that I created to download file from url and unzip it
# 
# #Do this just once! Create Directories
# # dir.create("Data")
# # dir.create("Data/rasterfile")
# # dir.create("Data/shapefile")
# 
# #1.Load functions and packages
# 
# # 1. LANDFIRE BPS remap raster data -- https://landfire.gov/version_download.php
# url_bps_remap <- "https://landfire.gov/bulk/downloadfile.php?FNAME=US_200_mosaic-LF2016_BPS_200_CONUS.zip&TYPE=landfire"
# sapply(url_bps_remap,fxns_dzip,dir="Data/rasterfile") 
# 
# # 2. MTBS Burned Areas Boundaries Dataset -- https://mtbs.gov/direct-download
# url_mtbs <- "https://edcintl.cr.usgs.gov/downloads/sciweb1/shared/MTBS_Fire/data/composite_data/burned_area_extent_shapefile/mtbs_perimeter_data.zip"
# sapply(url_mtbs,fxns_dzip,dir="Data/shapefile") 
# 
# # 3. Omernik vector data -- https://www.epa.gov/eco-research/ecoregions-north-america
# omernik3 <- "https://gaftp.epa.gov/EPADataCommons/ORD/Ecoregions/cec_na/NA_CEC_Eco_Level3.zip"
# sapply(url_omernik3,fxns_dzip,dir="Data/shapefile") # Omernick III layer includes Omernick I and Omernick II data
# 
# # 4. California state boundary -- https://data.ca.gov/dataset/ca-geographic-boundaries/resource/3db1e426-fb51-44f5-82d5-a54d7c6e188b
# url_CA_boundary <- "https://data.ca.gov/dataset/e212e397-1277-4df3-8c22-40721b095f33/resource/3db1e426-fb51-44f5-82d5-a54d7c6e188b/download/ca-state-boundary.zip"
# sapply(url_CA_boundary,fxns_dzip,dir="Data/shapefile")
# 
# # 5. CAL FIRE perimeters -- https://frap.fire.ca.gov/frap-projects/fire-perimeters/
# url_frap <- "https://frap.fire.ca.gov/media/50dgwqrb/fire20_1.zip" #the link to download FRAP GIS data. For details: https://frap.fire.ca.gov/frap-projects/fire-perimeters/
# sapply(url_frap,fxns_dzip,dir="Data/shapefile") #download zip file and unzip
# 
# ogrListLayers("Data/shapefile/fire20_1.gdb") #read gdb data
# firep20_1 = st_read("Data/shapefile/fire20_1.gdb","firep20_1") #read interested data layer "firep20_1" which depicts wildfire perimeters from contributing agencies current as 2020
# st_write(firep20_1, "Data/shapefile/firep20_1.shp")
# 
# # 6. LANDFIRE EVT data -- https://landfire.gov/version_download.php
# url_evt_remap<- "https://landfire.gov/bulk/downloadfile.php?FNAME=US_200_mosaic-LF2016_EVT_200_CONUS.zip&TYPE=landfire"
# url_evt_2014 <- "https://landfire.gov/bulk/downloadfile.php?FNAME=US_140_mosaic-US_140EVT_20180618.zip&TYPE=landfire"
# url_evt_2012 <- "https://landfire.gov/bulk/downloadfile.php?FNAME=US_130_Mosaic-US_130_EVT_04232015.zip&TYPE=landfire"
# url_evt_2010 <- "https://landfire.gov/bulk/downloadfile.php?FNAME=US_120_Mosaic-US_120_EVT_06142017.zip&TYPE=landfire"
# url_evt_2001 <- "https://landfire.gov/bulk/downloadfile.php?FNAME=US_105_mosaic_Refresh-US_105evt_09122104.zip&TYPE=landfire"
# 
# sapply(c(url_evt_remap,url_evt_2014,url_evt_2012,url_evt_2010,url_evt_2001),fxns_dzip,dir="Data/rasterfile/LANDFIRE_EVT") 

# # subset LANDFIRE EVT data to CA region, and reclassify EVT to general classification
# # this process will take a while. Just run once and load the data generated for future use

# # (1) read original EVT data
# dir_evt <- "Data/rasterfile/LANDFIRE_EVT"
# list_evt <- fxns_lfiles(dir_evt,"w001001.adf|.tif$") #[1:5]
# evt_raster <- lapply(list_evt, raster)

# # (2) subset to CA region
# # this process will take a while. Just run once and load the data generated for future use
# for (j in 1:5){
#  evt_raster[j][[1]] <-  evt_raster[j][[1]]%>% 
#                           crop(CA) %>%
#                           mask(CA)
#  }
# names(evt_raster) = c("evt200_ca","evt105_ca","evt120_ca","evt130_ca","evt140_ca")
# list2env(evt_raster,envir=.GlobalEnv)
# writeRaster(evt200_ca,"Data/Rdata/evt_ca/evt200_ca.tif")
# writeRaster(evt105_ca,"Data/Rdata/evt_ca/evt105_ca.tif")
# writeRaster(evt120_ca,"Data/Rdata/evt_ca/evt120_ca.tif")
# writeRaster(evt130_ca,"Data/Rdata/evt_ca/evt130_ca.tif")
# writeRaster(evt140_ca,"Data/Rdata/evt_ca/evt140_ca.tif")

## load the data 
## (3) get the raster value frequency
## this process will take a while. Just run once and load the final data generated for future use
# fxns_cb <- function(x){
#   as.data.frame(table(getValues(x)))
# }

# evt_values <- lapply(evt_raster, fxns_cb)
# names(evt_values) = c("evt200_value","evt105_value","evt120_value","evt130_value","evt140_value")
# list2env(evt_values,envir=.GlobalEnv) 

## (4) adding information from EVT CSV files
# list_csv <- fxns_lfiles(dir_evt,".csv$")
# evt_csv = lapply(list_csv, read.csv)
# names(evt_csv) = c("evt200","evt105","evt120","evt130","evt140")

# colnames(evt_csv$evt105)[which(colnames(evt_csv$evt105) == "Value")] <- "VALUE" # rename some columns to join with other dataframes
# evt_csv$evt105$EVT_FUEL <- evt_csv$evt105$VALUE
#
# for (k in 1:5){
#   if (any(colnames(evt_csv[k][[1]]) %in% c("EVT_Name","CLASSNAME"))){
#      colnames(evt_csv[k][[1]])[which(colnames(evt_csv[k][[1]]) %in% c("EVT_Name","CLASSNAME"))] <- "EVT_NAME"
#      colnames(evt_csv[k][[1]])[which(colnames(evt_csv[k][[1]]) == "EVT_Fuel")] <- "EVT_FUEL"
#   }
# }
#
# lapply(evt_csv,colnames)

# # (5) classify EVT to general classification
# load("Github_code_data/git_data/updated crosswalk of LANDFIRE vegetation types to five generic land covers.rda")
# 
# evt_class <- list()
# evt_class_freq <- list()
# evt_reclass <- list() 
# for (k in 1:5){
#   evt_class[k][[1]] <- evt_values[k][[1]] %>% 
#     mutate(Var1 = as.numeric(as.character(Var1))) %>% 
#     left_join(evt_csv[k][[1]],by=c("Var1"="VALUE")) %>%  
#     left_join(EVT_class)
#   evt_class_freq[k][[1]] <- table(evt_class[k][[1]]$veg_code,useNA = "always")
#   reclass_group <- data.matrix(evt_class[k][[1]][,c("Var1","veg_code")])
#   evt_reclass[k][[1]] <- reclassify(evt_raster[k][[1]],reclass_group)
#   writeRaster(evt_reclass[k][[1]],paste0("Data/Rdata/evt_reclass/",names(evt_csv)[k],"_reclass.tif"), overwrite=TRUE)
# } 
# 
# names(evt_reclass) = c("evt200_reclass","evt105_reclass","evt120_reclass","evt130_reclass","evt140_reclass")
# list2env(evt_reclass,envir=.GlobalEnv)

##---------------------------------------------------------------
##                    Section 1: Load data                     --
##---------------------------------------------------------------
# 1. CBI to dNBR/NBR regression models from Picotte et al 2021
load("Github_code_data/git_data/cbi_dnbr_model.rda")
load("Github_code_data/git_data/cbi_nbr_model.rda")

# 2. BPS raster
bps_remap <- raster("Data/rasterfile/LF2016_BPS_200_CONUS/Tif/LC16_BPS_200.tif") 
bps_csv <- read.csv("Data/rasterfile/LF2016_BPS_200_CONUS/CSV_Data/LF16_BPS_200.csv") 
bps_subcol <- bps_csv[,c("VALUE","BPS_CODE","BPS_NAME","GROUPVEG")] %>% # subset to columns of interests
                mutate(BPS_CODE2 = ifelse(nchar(BPS_CODE) == 5 & str_sub(BPS_CODE,-1,) == 0, str_sub(BPS_CODE,1,4), BPS_CODE)) # add a column to convert BPS_CODE to the same format used in CBI models

# 3. Omernik3
omernik3 <- st_read("Data/shapefile/NA_CEC_Eco_Level3.shp") # this also include Omernik 2 information
omernik3_repo <- st_transform(omernik3,mycrs)

# 4. California State boundary
CA <- st_read("Data/shapefile/CA_State_TIGER2016.shp")  %>% st_transform(crs(mycrs))

# 5. LANDFIRE EVT data and reclassified raster in CA region
dir_evt <- "Data/Rdata/evt_ca"  #get the file directories 
list_evt <- fxns_lfiles(dir_evt,".tif$") #[1:5]
evt_raster <- lapply(list_evt, raster)

dir_evt_reclass <- "Data/Rdata/evt_reclass"  #get the file directories 
list_evt_reclass <- fxns_lfiles(dir_evt_reclass,".tif$") #[1:5]
evt_reclass <- lapply(list_evt_reclass, raster)

# (6) emission factor
load("Github_code_data/git_data/emission_factor.rda")

##---------------------------------------------------------------
##            Section 2: User defined input/output            --
##---------------------------------------------------------------

# *****Attention here! Only run the command lines interested 
input_data_source <- "MTBS" # run this if using MTBS data
input_data_source <- "CALFIRE" # run this if using CAL FIRE data

img_pattern <- "_dnbr.tif$" # run this if using dNBR images 
img_pattern <- "_nbr.tif$" # run this if using NBR images (a small number of wildfires in MTBS dataset)
# *****end

if (input_data_source == "MTBS"){
  mtbs_perims <- st_read("Data/shapefile/mtbs_perims_DD.shp") # MTBS images and shapefiles for wildfires during 1984-2019
}

dir_input <- paste0("Data/", input_data_source) # to generate data based on MTBS data

cbi_model <- data.frame()
if(img_pattern == "_dnbr.tif$"){
     cbi_model <- dnbr_model
     }else{cbi_model <- nbr_model} 

list_tif <- fxns_lfiles(dir_input, img_pattern) # burn severity of 1564 wildfire events were classified using dNBR images

if(input_data_source == "CALFIRE"){
     fire2020 = st_read("Github_code_data/git_data/fraps_ver20211206/fraps_ver20211206.shp") %>% 
       filter(YEAR_ == 2020) %>%
       mutate(FID = 0:499) #read interested data layer "firep20_1" which depicts wildfire perimeters from contributing agencies current as 2020
     # split shapefile to get individual shapefile for each fire event in 2020
      dir.create("Data/shapefile/fire2020all")
      for (i in 1:nrow(fire2020)){
             tmp <- fire2020[i, ] 
             st_write(tmp, paste0("Data/shapefile/fire2020all/", fire2020$FID[i],".shp"))
      }  
     
     list_shp  <- fxns_lfiles("Data/shapefile/fire2020all/",".shp$")
     }else{list_shp  <- fxns_lfiles(dir_input,"_burn_bndy.shp$")}

list_dob <- fxns_lfiles("Data/rasterfile/dob_raster",".tif$") 
list_dup <- list_dob[basename(list_dob) %like% "\\(1\\)"] # 637 obs
ealier_list <- list_dob[basename(list_dob) %like% "^CA"] # 40 obs
list_dob <- setdiff(list_dob,c(list_dup,ealier_list)) # 6952 obs

list_ttmean <-list_dob[basename(list_dob) %like% "tt_mean"] # 1075 obs
ttID <- gsub("tt_mean_|_tt_mean|.tif","",basename(list_ttmean))
ttID <- ifelse(str_detect(ttID,"_"),sapply(strsplit(ttID, "_", fixed = TRUE),function(i) paste(tail(head(i,2),1), collapse = "_")),ttID)

list_ttmean[(basename(list_ttmean) %like% "MOSAIC")] # there were 2 fires that have MOSAIC daily progression images
# f_fid == c("247,"300)

## Do this just once! Create Directories
# dir.create("Data/Output_data/burn_severity")
# dir.create("Data/Output_data/burn_severity_CA")
# dir.create("Data/Output_data/emission_CA")
# dir.create("Data/Output_data/daily_emission_CA")

output_bs <- "Data/Output_data/burn_severity"
output_bs_ca <- "Data/Output_data/burn_severity_CA"
output_em_ca <- "Data/Output_data/emission_CA"
output_de_ca <- "Data/Output_data/daily_emission_CA"

##---------------------------------------------------------------
##                    Section 3: Calculation                   --
##---------------------------------------------------------------

#Create some empty lists/data frames to store information generated in the for loop

total_bs <- data.frame() # a data frame of area burned in different burn severity 
total_bs_ca <- data.frame() # a data frame of area burned in different burn severity of wildfire events for areas only in California region

total_em_class <- data.frame() # a data frame of emissions of each fire by vegetation and severity group
total_em <- data.frame() # a data frame of emission of each fire
total_bps_evt <- data.frame()
total_pm_day <- data.frame() # an empty data frame for calculating daily emissions

# Run on all historical wildfires
# you can try to set i = 1 and jump to start from # (1) to test before running the loop
for (i in 1:length(list_tif)){

    ##---------------------------------------------------------------
    # 1. Burn severity calculation
    ##---------------------------------------------------------------
    
    if(input_data_source == "MTBS"){
      # (1) read sev(dNBR/NBR raster)
      sev <- raster(list_tif[i])
      
      # (2) read shp(shapefile)
      sevdir <- str_remove(list_tif[i], img_pattern)
      shp <- st_read(paste0(sevdir,"_burn_bndy.shp")) %>% 
        st_transform(crs(sev))
      
      if(any(colnames(shp) %in% "Fire_ID")){
        colnames(shp)[which(colnames(shp) =="Fire_ID")] <- "Event_ID"
      } # MTBS uses either "Fire_ID" or "Event_ID". We will convert "Fire_ID" to "Event_ID" to be consistent
      
      if (nrow(shp)>0){
        shp <- shp %>% group_by(Event_ID) %>% 
          summarise()
      } # if there are multiple shape files for a wildfire event, we will merge them into one shape file
      
      #add mask
      mask <- st_read(paste0(sevdir,"_mask.shp"))%>% 
        st_transform(crs(sev))
      #if there is a mask area, mask it out from the shape file
      if (nrow(mask)>0){
        st_agr(shp) = "constant"
        st_agr(mask)= "constant"
        shp <- st_difference(st_make_valid(shp), st_make_valid(st_union(mask)))
      }
    }else{
       sev <- raster(list_tif[i]) %>% 
               projectRaster(crs=mycrs,res = c(30,30))
       file_name <- gsub("_dnbr.tif","",basename(list_tif[i]))
       f_fid <- sapply(strsplit(basename(list_tif[i]), "_"), "[[", 2)
       firename <- gsub("2020_|_dnbr.tif","",basename(list_tif[i]))
    
    # (2) read shapefile
    shp <- st_read(paste0("Data/shapefile/fire2020all","/",f_fid,".shp")) %>% 
      st_transform(crs(mycrs))
    }
    
    # (3) crop raster image to fire perimeter
    sev_crop <- sev %>% 
      crop(shp) %>%
      mask(shp) # mask dNBR images to the shapefile region (wildfire perimeter)
    
    if(any(origin(sev) !=c(15,15))){ 
      sev3 <- sev 
      origin(sev3) <- c(15,15) # since BPS raster also has the origin (15,15), we will convert it to be the same
      sev_crop3 <- crop(sev3,shp)
      sev_crop3 <- mask(sev_crop3, shp)
      sev_crop <- resample(sev_crop, sev_crop3,method='ngb') # assign the nearest neighbor raster value
    }
    
    # (4) apply offset value to raster
    if(input_data_source == "MTBS"){
         fireID <- str_sub(list_tif[i],16,36)
         offset <- mtbs_perims %>% 
                filter(Event_ID==toupper(fireID)) %>%
                pull(dNBR_offst)
         sev_off <- sev_crop - ifelse(length(offset)==0,0,as.numeric(offset))
    } else (sev_off <- sev_crop)
    # (5) find bps info for the raster, reclassify bps raster to convert values to be the same as the BPS code in CBI models
    bps_crop <- bps_remap %>%
                  crop(shp) %>%
                  mask(shp) 
    
    bps_value <- data.frame(table(getValues(bps_crop))) # get bps raster value and freq
    bps_value <- bps_value %>%
                   mutate(Var1 = as.numeric(as.character(Var1))) %>%
                   rename(VALUE = Var1) %>%
                   left_join(bps_subcol, by = "VALUE") %>% # add associated information about "VALUE" "BPS_CODE" "BPS_NAME" "GROUPVEG"
                   mutate(BPS_CODE2 = as.numeric(BPS_CODE2))
    
    bps_value$fireID <- fireID
    
    reclass_b <- data.matrix(bps_value[,c(1,6)])
    bps_model <- reclassify(bps_crop,reclass_b) # this is a raster layer with bps codes like those in CBI models, instead of original BPS raster value
    
    veg <- data.frame(table(getValues(bps_model)))
    
    veg_class  <- bps_value %>%
      select(BPS_CODE2,GROUPVEG) %>%
      rename(code=BPS_CODE2,veg=GROUPVEG)  
    
    
    # (6) classify forest/shrub and others 
    veg_class$class <- ifelse(veg_class$veg %in% c("Barren-Rock/Sand/Clay","Open Water","Perennial Ice/Snow"),1, # assign vegetation class to be 1 for unchanged group
                              ifelse(veg_class$veg %in% c("Conifer","Hardwood","Hardwood-Conifer","Shrubland","Riparian"),2, # assign vegetation class to be 2 for forest/shrub group
                                     ifelse(veg_class$veg %in% c("Savanna","Sparse","Grassland"),3,0 ))) # assign vegetation class to be 3 for grass group and the rest ("NoData", "") to be 0
    
    for (j in 1:nrow(veg_class)){
      if(veg_class$class[j]==2 & veg_class$code[j] %in% cbi_model$level) {veg_class$class[j] <- 0.02}
    } # if a pixel belongs to the category 2 (forest/shrub) and the BPS code has a valid equation, assign the vegetation classification code to be "0.02"
    
    reclass_c <- data.matrix(veg_class[,c(1,3)])
    bps_class <- reclassify(bps_model,reclass_c) 
    class <- as.data.frame(table(getValues(bps_class)))
    class$fireID <- fireID
    
    # (7) add Omernik 2 and Omernik 3 info
    st_agr(shp) = "constant"
    st_agr(omernik3_repo) = "constant"
    shp_O3 <- st_intersection(st_make_valid(shp),omernik3_repo) %>%
                mutate(percent = st_area(geometry)/st_area(shp)) %>%  
                mutate(fireID = fireID)
    fire_O2 <- unique(shp_O3$NA_L2CODE)
    fire_O3 <- unique(shp_O3$NA_L3CODE)
    
    # (8) pixel level calculation
    #If a forest/shrub BPS pixel has a valid equation, keep the BPS codes; 
    #Otherwise, go to step 1, check if any area burned in Omernik 3 regions that have a valid equation. If so, replace the BPS codes with Omernik 3 codes
    #Then go to step 2, check if any area burned in Omernik 2 regions that have a valid equation. If so, replace the BPS codes with Omernik 2 codes
    #Next, go to step 3, any forest/shrub BPS pixels that do not have a BPS, Omernik 3 or Omernik 2 valid equations, will apply the CONUS equation
    #Lastly, go to step 4, for any non forest/shrub pixels, we will assign other arbitrary values
    
    level_mix <- bps_model # level_mix will be the final model codes raster
    
    #step 1: crop to valid O3 region(s)
    O3valid_inter <- shp_O3[which(shp_O3$NA_L3CODE %in% cbi_model$level),] 
    
    if (nrow(O3valid_inter) > 0 ){
      for (t in 1:nrow(O3valid_inter)){
        O3_inter <- O3valid_inter[t,]
        O3_fire  <- bps_model %>%
          crop(O3_inter) %>%
          mask(O3_inter)
        O3_class <- bps_class %>%
          crop(O3_inter) %>%
          mask(O3_inter)
        ix <- getValues(O3_class)==2 
        O3_fire[ix] = cbi_model %>% filter(level==O3_inter$NA_L3CODE) %>% pull(level2)  
        level_mix  <- merge(O3_fire,level_mix )
      } 
    }
    
    #Step 2: crop to valid O2 region(s)
    
    shp_O2 <- shp_O3[which(!(shp_O3$NA_L3CODE %in% cbi_model$level)),]
    
    O2valid_inter <- shp_O2[which(shp_O2$NA_L2CODE %in% cbi_model$level),] 
    if (nrow(O2valid_inter) > 0 ){
      for (k in 1:nrow(O2valid_inter)){
        O2_inter <- O2valid_inter[k,]
        O2_fire  <- bps_model %>%
          crop(O2_inter) %>%
          mask(O2_inter)
        O2_class <- bps_class %>%
          crop(O2_inter) %>%
          mask(O2_inter)
        ix <- getValues(O2_class)==2 
        O2_fire[ix] = cbi_model %>% filter(level==O2_inter$NA_L2CODE) %>% pull(level2)  
        level_mix  <- merge(O2_fire,level_mix )
      } 
    }
    
    #Step 3: crop to the rest region
    shp_all <- shp_O2[which(!(shp_O2$NA_L2CODE %in% cbi_model$level)),]
    
    if (nrow(shp_all) > 0 ) {
      all_fire  <- bps_model %>%
        crop(shp_all) %>%
        mask(shp_all)
      
      all_class <- bps_class %>%
        crop(shp_all) %>%
        mask(shp_all)
      
      ix <- getValues(all_class)==2
      all_fire[ix] = cbi_model %>% filter(level=="All") %>% pull(level2)  
      level_mix  <- merge(all_fire,level_mix )
    }
    
    #Step 4: other veg classes
    ix <- getValues(bps_class)==0 
    level_mix[ix] = 0 # assign raster pixel to be 0 if belongs to the "NA" group 
    
    ix <- getValues(bps_class)==1 
    level_mix[ix] = 1 # assign raster pixel to be 1 if belongs to the "unchanged" group
    
    ix <- getValues(bps_class)==3 
    level_mix[ix] = 3 # assign raster pixel to be 3 if belongs to the "grass" group
    
    #Step 5: pixel level calculation. Need some help to reduce redundancy. I will send out some codes along with the email
    if(img_pattern == "_dnbr.tif$"){
    cbi_value <- overlay(sev_off,level_mix, fun=function(x, y) { ifelse(y==1018, fxns_polyL(fxns_para2(1018),x),
                                                                        ifelse(y==1053, fxns_polyL(fxns_para2(1053),x),
                                                                               ifelse(y==11791, fxns_polyL(fxns_para2(11791),x),
                                                                                      ifelse(y==11792, fxns_polyL(fxns_para2(11792),x),
                                                                                             ifelse(y==1315, fxns_polyL(fxns_para2(1315),x),
                                                                                                    ifelse(y==101, fxns_sigB(fxns_para2(101),x),
                                                                                                           ifelse(y==1011, fxns_sigB(fxns_para2(1011),x),
                                                                                                                  ifelse(y==1016, fxns_sigB(fxns_para2(1016),x),
                                                                                                                         ifelse(y==10321, fxns_sigB(fxns_para2(10321),x),
                                                                                                                                ifelse(y==10322, fxns_sigB(fxns_para2(10322),x),
                                                                                                                                       ifelse(y==1045, fxns_sigB(fxns_para2(1045),x),
                                                                                                                                              ifelse(y==10452, fxns_sigB(fxns_para2(10452),x),
                                                                                                                                                     ifelse(y==1046, fxns_sigB(fxns_para2(1046),x),
                                                                                                                                                            ifelse(y==10471, fxns_sigB(fxns_para2(10471),x),
                                                                                                                                                                   ifelse(y==1051, fxns_sigB(fxns_para2(1051),x),
                                                                                                                                                                          ifelse(y==1052, fxns_sigB(fxns_para2(1052),x),
                                                                                                                                                                                 ifelse(y==10531, fxns_sigB(fxns_para2(10531),x),
                                                                                                                                                                                        ifelse(y==1054, fxns_sigB(fxns_para2(1054),x),
                                                                                                                                                                                               ifelse(y==1055, fxns_sigB(fxns_para2(1055),x),
                                                                                                                                                                                                      ifelse(y==1056, fxns_sigB(fxns_para2(1056),x),
                                                                                                                                                                                                             ifelse(y==1061, fxns_sigB(fxns_para2(1061),x),
                                                                                                                                                                                                                    ifelse(y==1126, fxns_sigB(fxns_para2(1126),x),
                                                                                                                                                                                                                           ifelse(y==1159, fxns_sigB(fxns_para2(1159),x),
                                                                                                                                                                                                                                  ifelse(y==1166, fxns_sigB(fxns_para2(1166),x),
                                                                                                                                                                                                                                         ifelse(y==125, fxns_sigB(fxns_para2(125),x),
                                                                                                                                                                                                                                                ifelse(y==126, fxns_sigB(fxns_para2(126),x),
                                                                                                                                                                                                                                                       ifelse(y==127, fxns_sigB(fxns_para2(127),x),
                                                                                                                                                                                                                                                              ifelse(y==128, fxns_sigB(fxns_para2(128),x),
                                                                                                                                                                                                                                                                     ifelse(y==129, fxns_sigB(fxns_para2(129),x),
                                                                                                                                                                                                                                                                            ifelse(y==130, fxns_sigB(fxns_para2(130),x),
                                                                                                                                                                                                                                                                                   ifelse(y==131, fxns_sigB(fxns_para2(131),x),
                                                                                                                                                                                                                                                                                          ifelse(y==132, fxns_sigB(fxns_para2(132),x),
                                                                                                                                                                                                                                                                                                 ifelse(y==133, fxns_sigB(fxns_para2(133),x),
                                                                                                                                                                                                                                                                                                        ifelse(y==134, fxns_sigB(fxns_para2(134),x),
                                                                                                                                                                                                                                                                                                               ifelse (y==3, 3.5, #assign CBI value 5  to grass
                                                                                                                                                                                                                                                                                                                       ifelse (y==1, 0.01, NA) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) )}) #assign CBI value 0.05 (mean unchanged) to Veg class 1 or NA
    }else{cbi_value <- overlay(sev_crop,level_mix, fun=function(x, y) { ifelse(y==1011, fxns_polyL(fxns_para2(1011),x),
                                                                              ifelse(y==1016, fxns_polyL(fxns_para2(1016),x),
                                                                                     ifelse(y==1018, fxns_polyL(fxns_para2(1018),x),
                                                                                            ifelse(y==1027, fxns_polyL(fxns_para2(1027),x),
                                                                                                   ifelse(y==1028, fxns_polyL(fxns_para2(1028),x),
                                                                                                          ifelse(y==1051, fxns_polyL(fxns_para2(1051),x),
                                                                                                                 ifelse(y==1052, fxns_polyL(fxns_para2(1052),x),
                                                                                                                        ifelse(y==1053, fxns_polyL(fxns_para2(1053),x),
                                                                                                                               ifelse(y==1061, fxns_polyL(fxns_para2(1061),x),
                                                                                                                                      ifelse(y==1080, fxns_polyL(fxns_para2(1080),x),
                                                                                                                                             ifelse(y==11791, fxns_polyL(fxns_para2(11791),x),
                                                                                                                                                    ifelse(y==11792, fxns_polyL(fxns_para2(11792),x),
                                                                                                                                                           ifelse(y==1369, fxns_polyL(fxns_para2(1369),x),
                                                                                                                                                                  ifelse(y==135, fxns_polyL(fxns_para2(135),x),
                                                                                                                                                                         ifelse(y==136, fxns_polyL(fxns_para2(136),x),
                                                                                                                                                                                ifelse(y==139, fxns_polyL(fxns_para2(139),x),
                                                                                                                                                                                       ifelse(y==142, fxns_polyL(fxns_para2(142),x),
                                                                                                                                                                                              ifelse(y==101, fxns_sigB(fxns_para2(101),x),
                                                                                                                                                                                                     ifelse(y==10321, fxns_sigB(fxns_para2(10321),x),
                                                                                                                                                                                                            ifelse(y==10322, fxns_sigB(fxns_para2(10322),x),
                                                                                                                                                                                                                   ifelse(y==1045, fxns_sigB(fxns_para2(1045),x),
                                                                                                                                                                                                                          ifelse(y==10451, fxns_sigB(fxns_para2(10451),x),
                                                                                                                                                                                                                                 ifelse(y==10452, fxns_sigB(fxns_para2(10452),x),
                                                                                                                                                                                                                                        ifelse(y==1046, fxns_sigB(fxns_para2(1046),x),
                                                                                                                                                                                                                                               ifelse(y==10471, fxns_sigB(fxns_para2(10471),x),
                                                                                                                                                                                                                                                      ifelse(y==10531, fxns_sigB(fxns_para2(10531),x),
                                                                                                                                                                                                                                                             ifelse(y==1054, fxns_sigB(fxns_para2(1054),x),
                                                                                                                                                                                                                                                                    ifelse(y==1055, fxns_sigB(fxns_para2(1055),x),
                                                                                                                                                                                                                                                                           ifelse(y==1056, fxns_sigB(fxns_para2(1056),x),
                                                                                                                                                                                                                                                                                  ifelse(y==1125, fxns_sigB(fxns_para2(1125),x),
                                                                                                                                                                                                                                                                                         ifelse(y==1126, fxns_sigB(fxns_para2(1126),x),
                                                                                                                                                                                                                                                                                                ifelse(y==1166, fxns_sigB(fxns_para2(1166),x),
                                                                                                                                                                                                                                                                                                       ifelse(y==1315, fxns_sigB(fxns_para2(1315),x),
                                                                                                                                                                                                                                                                                                              ifelse(y==1330, fxns_sigB(fxns_para2(1330),x),
                                                                                                                                                                                                                                                                                                                     ifelse(y==1480, fxns_sigB(fxns_para2(1480),x),
                                                                                                                                                                                                                                                                                                                            ifelse(y==132, fxns_sigB(fxns_para2(132),x),
                                                                                                                                                                                                                                                                                                                                   ifelse(y==133, fxns_sigB(fxns_para2(133),x),
                                                                                                                                                                                                                                                                                                                                          ifelse(y==134, fxns_sigB(fxns_para2(134),x),
                                                                                                                                                                                                                                                                                                                                                 ifelse(y==138, fxns_sigB(fxns_para2(138),x),
                                                                                                                                                                                                                                                                                                                                                        ifelse(y==140, fxns_sigB(fxns_para2(140),x),
                                                                                                                                                                                                                                                                                                                                                               ifelse(y==143, fxns_sigB(fxns_para2(143),x),
                                                                                                                                                                                                                                                                                                                                                                      ifelse(y==148, fxns_sigB(fxns_para2(148),x),
                                                                                                                                                                                                                                                                                                                                                                             ifelse (y==3, 3.5, #assign CBI value 5  to grass
                                                                                                                                                                                                                                                                                                                                                                                     ifelse (y==1, 0.01, NA) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) )  }) #assign CBI value 0.05 (mean unchanged) to Veg class 1 or NA
  }
    
    
    ix <- getValues(cbi_value) <= 0 
    cbi_value[ix] = 0.001
    
    ix <- getValues(cbi_value) > 3 &  getValues(cbi_value)!= 3.5
    cbi_value[ix] = 3
    if(input_data_source == "CALFIRE"){
      cal_fire_record <- fire2020 %>% st_drop_geometry() %>% filter(FID == f_fid)
      out_filename <- paste0(cal_fire_record$STATE,format(cal_fire_record$ALARM_DATE, "%Y%m%d"),"_",format(cal_fire_record$CONT_DATE, "%Y%m%d"),"_UCM",firename,"_CF",cal_fire_record$INC_NUM)
      writeRaster(cbi_value,paste0(output_bs,"/",out_filename,"_cbi_value.tif"), overwrite=TRUE)
      sev_classified <- reclassify(cbi_value,class.matrix) # classify CBI value to burn severity classfication
      writeRaster(sev_classified,paste0(output_bs,"/",out_filename,"_cbi_severity.tif"), overwrite=TRUE)
      fireID <- out_filename
      dob_filename <- paste0("CA",
             gsub("[^0-9.-]", "", sapply(strsplit(fireID, "_", fixed = TRUE),function(i) paste(tail(i,1), collapse = "_"))),
             "_2020",
             sapply(strsplit(firename, "_", fixed = TRUE),function(i) paste(tail(i,-1), collapse = "_"))
               ) 
    }else{
      writeRaster(cbi_value,paste0(output_bs,"/",fireID,"_cbi_value.tif"),overwrite = TRUE)
      sev_classified <- reclassify(cbi_value,class.matrix) # classify CBI value to burn severity classfication
      writeRaster(sev_classified,paste0(output_bs,"/",fireID,"_cbi_severity.tif"),overwrite = TRUE)
      
    }
    
    cat_sev<- data.frame(rbind(table(values(sev_classified)))) #freq of different severity 
    colnames(cat_sev) = gsub("X", "Sev", colnames(cat_sev))
    
    cat_sev <- cat_sev %>%
                 mutate(Event_ID = toupper(fireID),
                        burned_areaHa = fxns_mh(sum(cat_sev)),
                        index = i) 
    total_bs <-  bind_rows(total_bs,cat_sev)
    
    bs_raster <- sev_classified
    
   
      # fire area only in California region 
      bs_raster_ca <- bs_raster %>% crop(CA) %>% mask(CA)
      writeRaster(bs_raster_ca,paste0(output_bs_ca,"/",fireID,"_burn_severity_CA.tif"),overwrite = TRUE)
      
      bs_ca <- data.frame(rbind(table(values(bs_raster_ca)))) #freq of different severity 
      colnames(bs_ca) = gsub("X", "Sev", colnames(bs_ca))
      bs_ca <- bs_ca %>%
                 mutate(Event_ID = toupper(fireID),
                        burned_areaHa = fxns_mh(sum(bs_ca)),
                        index = i)  
      total_bs_ca <-  bind_rows(total_bs_ca,bs_ca)
##---------------------------------------------------------------
              # 2. Emission calculation
##---------------------------------------------------------------
    if(input_data_source == "MTBS"){
      fireYear <- str_sub(list_tif[i],29,32)
    }else{
      fireYear <- str_sub(basename(list_tif[i]),1,4)
    }
    
    
    #(3) read CBI severity
    bs_raster_ca <- raster(paste0(output_bs_ca,"/",fireID,"_burn_severity_CA.tif"))
    
    #(4) choose which LANDFIRE EVT version based on fireYear
    evt_reclass_layer <- ifelse(fireYear <= 2001, evt_reclass[1], #if a fire burned before the end of 2001, we will use LF 1.0.5
                        ifelse(fireYear > 2001 & fireYear <= 2010, evt_reclass[2], #if a fire burned between 2002 and the end of 2010, we will use LF 1.2.0
                               ifelse(fireYear > 2010 & fireYear <= 2012, evt_reclass[3], #if a fire burned between 2011 and the end of 2012, we will use LF 1.3.0
                                      ifelse(fireYear > 2012 & fireYear <= 2014, evt_reclass[4], evt_reclass[5]))))  #if a fire burned between 2013 and the end of 2014, we will use LF 1.4.0; for all the fires burned starting 2015 we will use LF 2.0.0
    
    evt_reclass_crop <- evt_reclass_layer[[1]] %>%
                  crop(shp) %>%
                  mask(shp) 
    
    evt_layer <- ifelse(fireYear <= 2001, evt_raster[1], #if a fire burned before the end of 2001, we will use LF 1.0.5
                        ifelse(fireYear > 2001 & fireYear <= 2010, evt_raster[2], #if a fire burned between 2002 and the end of 2010, we will use LF 1.2.0
                               ifelse(fireYear > 2010 & fireYear <= 2012, evt_raster[3], #if a fire burned between 2011 and the end of 2012, we will use LF 1.3.0
                                      ifelse(fireYear > 2012 & fireYear <= 2014, evt_raster[4], evt_raster[5]))))  #if a fire burned between 2013 and the end of 2014, we will use LF 1.4.0; for all the fires burned starting 2015 we will use LF 2.0.0
    
    evt_crop <- evt_layer[[1]] %>%
                  crop(shp) %>%
                  mask(shp) 
    
    #(5) calculate 10 types of emissions and save as a raster stack
    em_stack=list()
    for(j in 1:10){
      em_stack[[j]] = fxns_em(evt_reclass_crop, bs_raster_ca,emission_factor[3,j+1][[1]],emission_factor[2,j+1][[1]],emission_factor[1,j+1][[1]])
    }
    
    all_stack = stack(em_stack)
    all_stack = stack(evt_reclass_crop, bs_raster_ca,all_stack)
    names(all_stack) <- c("evt_class","cbi_severity",colnames(emission_factor)[2:11])
    writeRaster(all_stack,paste0("Data/Output_data/emission_CA/",fireID,"_stack.grd"), format="raster",overwrite =  TRUE)
    #calculate emissions based on the group conditions of vegetation classification and burn severity (e.g. how much pm2.5 emitted from Forest 5500'-7500' that burned in low severity...)
    all_freq<- as.data.frame(extract(all_stack,c(1:ncell(all_stack)))) # get the values from each cell
    all_freq <-all_freq[complete.cases(all_freq), ] # remove cells with all NA values in each layer
    em_class <- all_freq  %>%
                group_by(evt_class,cbi_severity) %>%
                summarise_all(list(sum)) %>%
                left_join(all_freq  %>%
                            group_by(evt_class,cbi_severity) %>% 
                            summarise(Freq=n()),by = c("evt_class", "cbi_severity"))%>%
                mutate(Event_ID = fireID)  # sum emissions of each group
    total_em_class <- bind_rows(total_em_class,em_class) # store the results from loop
    
    em_ca <- as.data.frame(t(cellStats(all_stack[[3:12]],sum, merge=TRUE))) %>%
                mutate(Event_ID = fireID) # sum emissions of each fire event
    total_em <- bind_rows(total_em,em_ca) # store the results from loop 

    ##---------------------------------------------------------------
    # 3. BPS and EVT classification  per pixel comparison for each fire
    ##---------------------------------------------------------------
    bps_crop_ca <- bps_crop %>% crop(CA) %>% mask(CA)
    bps_evt_stack = stack(bps_crop_ca,evt_crop)
    names(bps_evt_stack) <- c("BPS_value","EVT_value")
    writeRaster(bps_evt_stack,paste0("Data/Output_data/bps_evt_CA/",fireID,"_stack.grd"), format="raster",overwrite =  TRUE)
    
    be_freq<- as.data.frame(extract(bps_evt_stack,c(1:ncell( bps_evt_stack)))) # get the values from each cell
    be_freq <-be_freq[complete.cases(be_freq), ] # remove cells with all NA values in each layer
    colnames(be_freq) <- c("bps","evt")
    be_sum <- be_freq  %>%
      group_by(bps,evt) %>%
      summarise(count=n()) %>%
      mutate(Event_ID = fireID, FireYear = fireYear)  # sum emissions of each group
    
    total_bps_evt <- bind_rows(total_bps_evt,be_sum)
    
    ##---------------------------------------------------------------
    # 4. Daily PM2.5 Emission Calculation based on Time Tuple Mean data
    ##---------------------------------------------------------------    
  if(input_data_source =="MTBS") {
    dobID <- fireID
  } else{dobID <- f_fid}
    
  if(toupper(dobID) %in% ttID) {
    if (fireYear == 2020){
       if (dobID %in% c("247", "300")){
          dob_mean <- do.call(merge, lapply(list_ttmean[ttID == f_fid],raster))
            }else{dob_mean <- raster(list_ttmean[list_ttmean %like% file_name])}
    } else { dob_mean <- raster(list_ttmean[list_ttmean %like% toupper(fireID)])}

    dob_repo <- dob_mean %>%
      projectRaster(crs=mycrs,method = "ngb")
      
    
    dob_floor <- floor(dob_repo)
    #raster to polygon, merge the polygons that on the same date, get the centroids of the date
    DOY <- as(dob_repo,'SpatialPolygonsDataFrame') # convert to SpatialPolygonsDataFrame
    sf_DOY <- st_as_sf(DOY) 
    colnames(sf_DOY)[1] <- "dob"
    sf_DOY$dob <- floor(sf_DOY$dob)
    # sf_DOY <- sf_DOY[which(sf_DOY$dob!=0),]
    
    union <- unionSpatialPolygons(DOY, sf_DOY$dob)
    
    union_repo <- spTransform(union,"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    centroids <- getSpPPolygonsLabptSlots(union_repo)
    
    union_df <- as(union , "SpatialPolygonsDataFrame")
    
    union_df@data$dummy <- as.numeric(rownames(union_df@data))
    colnames(union_df@data)[1] = "doy"
      #(1)pm
      pm <- raster::extract(all_stack$PM2.5, union, fun = sum, na.rm = TRUE ) %>% 
        as.data.frame() %>% 
        mutate(doy = union_df@data$doy ) 
      
      centroids <- centroids %>%
        as.data.frame() %>%
        magrittr::set_colnames(c("Lon","Lat"))
      
      day_emission <- cbind(pm,centroids) %>% 
        dplyr::rename(PM = V1) %>%
        mutate(fireID = fireID,dob_ID = i) 
      
      union_df@data <- union_df@data %>% 
                       left_join(day_emission[,-6])
      writeOGR(union_df, output_de_ca, paste0(fireID,"_daily_emission"), 
                        driver = "ESRI Shapefile",
               overwrite_layer = T)
      total_pm_day <-bind_rows(total_pm_day,day_emission)
    }
    
  }


total_bs$unchanged <- fxns_mh(total_bs$Sev1)
total_bs$low <- fxns_mh(total_bs$Sev2)
total_bs$mod <- fxns_mh(total_bs$Sev3)
total_bs$hig <- fxns_mh(total_bs$Sev4)
total_bs$grass <- fxns_mh(total_bs$Sev5)

save(total_bs,file=paste0("Data/Rdata/total_bs_",input_data_source,str_sub(img_pattern, ,-6),".rda")) #final burn severity data set without using MTBS mask files

total_bs_ca$unchanged <- fxns_mh(total_bs_ca$Sev1)
total_bs_ca$low <- fxns_mh(total_bs_ca$Sev2)
total_bs_ca$mod <- fxns_mh(total_bs_ca$Sev3)
total_bs_ca$hig <- fxns_mh(total_bs_ca$Sev4)
total_bs_ca$grass <- fxns_mh(total_bs_ca$Sev5)
save(total_bs_ca,file=paste0("Data/Rdata/total_bs_ca_",input_data_source,str_sub(img_pattern, ,-6),".rda")) #wildfire events with mask files

save(total_em,file=paste0("Data/Rdata/total_em_",input_data_source,str_sub(img_pattern, ,-6),".rda")) #wildfire events with mask files
save(total_bps_evt,file=paste0("Data/Rdata/total_bps_evt_",input_data_source,str_sub(img_pattern, ,-6),".rda")) #wildfire events with mask files
save(total_em_class,file=paste0("Data/Rdata/total_em_class_",input_data_source,str_sub(img_pattern, ,-6),".rda")) #wildfire events with mask files
save(total_pm_day,file=paste0("Data/Rdata/total_pm_day_",input_data_source,str_sub(img_pattern, ,-6),".rda")) #wildfire events with mask files
