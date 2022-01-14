# This code is to generate burn severity based on CBI values converted from dNBR values
# Author: Qingqing Xu
# Date: 05/13/2021
# Contact: qxu6@ucmerced.edu

# ATTENTION! Prerequisite:
# 1. Download CBI regression models from this GitHub repository: CBSE/Github_code_data/Data
# 2. Run the script "1-download-data.R" to download data for use to run this script
# 3. If you are interested to generate data based on Monitoring Trends in Burn Severity (MTBS) shapefile and dNBR/NBR images, you should go to https://mtbs.gov/viewer/index.html to download wildfire data interested and save the data under the directory of this project: Data/MTBS_CA_19842019
# 4. If you are interested to generate data based on CAL FIRE shapefiles and composite images generated with Google Earth Engine, you should contact the authors for access and save the data under the directory of this project: Data/

options(warn=1) # print out warnings 

##---------------------------------------------------------------
##    Load packages, create functions, and define variables    --
##---------------------------------------------------------------
#(1) Load packages
ndpkg <- function(p){
  if(!is.element(p,installed.packages()[,1]))
  {install.packages(p, dep = T)}
  require(p,character.only = T)
} # a function to install and load required packages

sapply(
  c("raster","sf", "sp","rgdal","reshape2","foreign","dplyr","tidyverse","data.table","stringr","ggmap","compare","purrr"),
  ndpkg)

#(2) Create some functions for use to run this script
fxns_lfiles <- function(dir,pat){
  list.files(dir,pattern = pat, recursive = TRUE,full.names = TRUE)
} # a function to get the full directories given a pattern

fxns_pacre <- function(x) as.numeric(x)*30*30/4047 # a function to convert pixel numbers to acres

fxns_para2 <- function(objs) {
  my_list <- list(a = dnbr_model[which(dnbr_model$level2 == objs),"A"],
                  b = dnbr_model[which(dnbr_model$level2 == objs),"B"],
                  c = dnbr_model[which(dnbr_model$level2 == objs),"C"])
  return(my_list) 
} # a function to get the a, b, c parameters of the CBI equation

fxns_sigB <- function(t,h) {
  t$a/(1.0+exp(-1.0*(h-t$b)/t$c))} # a function to apply the "pyeq3.Models_2D.Sigmoidal.SigmoidB" model

fxns_polyL <- function(t,h) {
  t$a+t$b*h} # a function to apply the "pyeq3.Models_2D.Polynomial.Linear" model

#(3) Set CBI burn severity threshold, need to ask Josh about vegetative CBI values
reclass_sev <- c(0,0.1, 1, # unchanged, the reference table set (0,0.1] for unchanged
                 0.1,1.25,2,   # low severity
                 1.25,2.25,3,  # moderate severity
                 2.25,3,4,     # high severity
                 3,4, 5)       # grass total burn, I set CBI value to be 3.5 to classify as grass in calculation
#CBI threshold reference: https://www.researchgate.net/publication/261559367_A_Mixed-Effects_Heterogeneous_Negative_Binomial_Model_for_Post-fire_Conifer_Regeneration_in_Northeastern_California_USA
class.matrix <- matrix(reclass_sev, ncol = 3, byrow = TRUE)

#(4) Define the crs that we will use for this analysis
mycrs <- "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs" # the same as MTBS dNBR CRS

##---------------------------------------------------------------
##                    Section 1: Load data                     --
##---------------------------------------------------------------
#(1) CBI to dNBR/NBR regression models from Picotte et al 2021
load("Github_code_data/Data/CBSE_cbi_to_dNBR_NBR_regression_models(cbi_eq).rda")

#(2) BPS raster
bps_remap <- raster("Data/rasterfile/LF2016_BPS_200_CONUS/Tif/LC16_BPS_200.tif") 
bps_csv <- read.csv("Data/rasterfile/LF2016_BPS_200_CONUS/CSV_Data/LF16_BPS_200.csv") 
bps_subcol <- bps_csv[,c("VALUE","BPS_CODE","BPS_NAME","GROUPVEG")] %>% # subset to columns of interests
  mutate(BPS_CODE2 = ifelse(nchar(BPS_CODE) == 5 & str_sub(BPS_CODE,-1,) == 0, str_sub(BPS_CODE,1,4), BPS_CODE)) # add a column to convert BPS_CODE to the same format used in CBI models

#(3) Omernik3
omernik3 <- st_read("Data/shapefile/NA_CEC_Eco_Level3.shp") # this also include Omernik 2 information
omernik3_repo <- st_transform(omernik3,mycrs)

##---------------------------------------------------------------
##                    Section 2: Input data                   --
##---------------------------------------------------------------
#(1) MTBS images and shapefiles for wildfires during 1984-2019
mtbs_perims <- st_read("Data/shapefile/mtbs_perims_DD.shp")
#1.Get list of directories of dnbr, shp, mask
mtbs_dir <- "Data/MTBS_CA_19842019"
dnbr_list <- fxns_lfiles(mtbs_dir,"_dnbr.tif$") # burn severity of 1564 wildfire events were classified using dNBR images
nbr_list  <- fxns_lfiles(mtbs_dir, "_nbr.tif$") # burn severity of 59 wildfire events were classified using NBR images
shp_list  <- fxns_lfiles(mtbs_dir,"_burn_bndy.shp$") #1623
dnbr6_list <- fxns_lfiles(mtbs_dir,"_dnbr6.tif$")
nbr6_list <- lfiles_fxns(mtbs_dir,"_nbr6.tif$")

# (2) 
##---------------------------------------------------------------
##                    Section 2: Calculation                   --
##---------------------------------------------------------------

#(1)Create some empty lists/data frames to store information generated in the for loop
myList <- vector(mode = "list") # a list to show the results from each for loop run (warnings and errors)
total_cbi <- data.frame() # a data frame of area burned in different burn severity
total_cbi_mask <- data.frame() # a data frame of area burned in different burn severity of wildfire events with MTBS mask files

mtbs_rat <- data.frame("ID" = 0:6 ) # a data frame we will use to select legend and colors when plotting MTBS severity map
cbi_rat <- data.frame("ID" = 1:5 ) # a data frame we will use to select legend and colors when plotting CBI severity map

#(2)Run on all historical wildfires
for (i in 1:length(dnbr_list)){
  r <- tryCatch({
    #(1)read MTBS sev(dNBR raster), shp(shapefile)
    sev <- raster(dnbr_list[i])
    sevdir <- str_sub(dnbr_list[i],1,-10)
    shp <- st_read(paste0(sevdir,"_burn_bndy.shp")) %>% 
      st_transform(crs(sev))
    
    if(any(colnames(shp) %in% "Fire_ID")){
      colnames(shp)[which(colnames(shp) =="Fire_ID")] <- "Event_ID"
    } # MTBS uses "Fire_ID" or "Event_ID". We will convert "Fire_ID" to "Event_ID" to be consistent
    
    if (nrow(shp)>0){
      shp <- shp %>% group_by(Event_ID) %>% 
        summarise()
    } # if there are multiple shape files for a wildfire event, we will merge them into one shape file
    
    sev_crop <- sev %>% 
      crop(shp) %>%
      mask(shp) # mask dNBR images to the shape file region (wildfire perimeter)
    
    if(any(origin(sev) !=c(15,15))){ #there are two fire events and all fire events in 2019 that the origin is not the same as the rest (15,15)
      sev3 <- sev 
      origin(sev3) <- c(15,15) # since BPS raster also has the origin (15,15), we will convert it to be the same
      sev_crop3 <- crop(sev3,shp)
      sev_crop3 <- mask(sev_crop3, shp)
      sev_crop <- resample(sev_crop, sev_crop3,method='ngb') # assign the nearest neighbor raster value
    }
    
    #(2)apply offset value to dNBR raster
    fireID <- str_sub(dnbr_list[i],28,48)
    offset <- mtbs_perims %>% 
      filter(Event_ID==toupper(fireID)) %>%
      pull(dNBR_offst)
    sev_off <- sev_crop - ifelse(length(offset)==0,0,as.numeric(offset))
    
    #(3)find bps info for the raster, reclassify bps raster to convert values to be the same as the BPS code in CBI models
    bps_crop <- bps_remap %>%
      crop(shp) %>%
      mask(shp) 
    
    bps_value <- data.frame(table(getValues(bps_crop))) #get bps raster value and freq
    bps_value <- bps_value %>%
      mutate(Var1 = as.numeric(as.character(Var1))) %>%
      rename(VALUE = Var1) %>%
      left_join(bps_subcol, by = "VALUE") %>% # add associated information about "VALUE" "BPS_CODE" "BPS_NAME" "GROUPVEG"
      mutate(BPS_CODE2 = as.numeric(BPS_CODE2))
    
    bps_value$fireID <- fireID
    
    reclass_b <- data.matrix(bps_value[,c(1,6)])
    bps_fire <- reclassify(bps_crop,reclass_b) # this is a raster layer with bps codes like those in CBI models, instead of original BPS raster value
    
    veg <- data.frame(table(getValues(bps_fire)))
    
    veg_class  <- bps_value %>%
      select(BPS_CODE2,GROUPVEG) %>%
      rename(code=BPS_CODE2,veg=GROUPVEG) # this is a layer of BPS vegetation group classification 
    
    
    #(4) classify forest/shrub and others 
    veg_class$class <- ifelse(veg_class$veg %in% c("Barren-Rock/Sand/Clay","Open Water","Perennial Ice/Snow"),1, # assign vegetation class to be 1 for unchanged group
                              ifelse(veg_class$veg %in% c("Conifer","Hardwood","Hardwood-Conifer","Shrubland","Riparian"),2, # assign vegetation class to be 2 for forest/shrub group
                                     ifelse(veg_class$veg %in% c("Savanna","Sparse","Grassland"),3,0 ))) # assign vegetation class to be 3 for grass group and the rest ("NoData", "") to be 0
    
    for (j in 1:nrow(veg_class)){
      if(veg_class$class[j]==2 & veg_class$code[j] %in% dnbr_model$level) {veg_class$class[j] <- 0.02}
    } # if a pixel belongs to the category 2 (forest/shrub) and the BPS code has a valid equation, assign the vegetation classification code to be "0.02"
    
    reclass_c <- data.matrix(veg_class[,c(1,3)])
    bps_class <- reclassify(bps_fire,reclass_c) 
    class <- as.data.frame(table(getValues(bps_class)))
    class$fireID <- fireID
    
    #(5)add Omernik 2 and Omernik 3 info
    st_agr(shp) = "constant"
    st_agr(omernik3_repo) = "constant"
    shp_O3 <- st_intersection(st_make_valid(shp),omernik3_repo) %>%
      mutate(percent = st_area(geometry)/st_area(shp)) %>%  
      mutate(fireID = fireID)
    fire_O2 <- unique(shp_O3$NA_L2CODE)
    fire_O3 <- unique(shp_O3$NA_L3CODE)
    
    #(6)pixel level calculation
    #If a forest/shrub BPS pixel has a valid equation, keep the BPS codes; 
    #Otherwise, go to step 1, check if any area burned in Omernik 3 regions that have a valid equation. If so, replace the BPS codes with Omernik 3 codes
    #Then go to step 2, check if any area burned in Omernik 2 regions that have a valid equation. If so, replace the BPS codes with Omernik 2 codes
    #Next, go to step 3, any forest/shrub BPS pixels that do not have a BPS, Omernik 3 or Omernik 2 valid equations, will apply the CONUS equation
    #Lastly, go to step 4, for any non forest/shrub pixels, we will assign other arbitrary values
    
    level_mix <- bps_fire # level_mix will be the final model codes raster
    
    #step 1: crop to valid O3 region(s)
    O3valid_inter <- shp_O3[which(shp_O3$NA_L3CODE %in% dnbr_model$level),] 
    
    if (nrow(O3valid_inter) > 0 ){
      for (t in 1:nrow(O3valid_inter)){
        O3_inter <- O3valid_inter[t,]
        O3_fire  <- bps_fire %>%
          crop(O3_inter) %>%
          mask(O3_inter)
        O3_class <- bps_class %>%
          crop(O3_inter) %>%
          mask(O3_inter)
        ix <- getValues(O3_class)==2 
        O3_fire[ix] = dnbr_model %>% filter(level==O3_inter$NA_L3CODE) %>% pull(level2)  
        level_mix  <- merge(O3_fire,level_mix )
      } 
    }
    
    #Step 2: crop to valid O2 region(s)
    
    shp_O2 <- shp_O3[which(!(shp_O3$NA_L3CODE %in% dnbr_model$level)),]
    
    O2valid_inter <- shp_O2[which(shp_O2$NA_L2CODE %in% dnbr_model$level),] 
    if (nrow(O2valid_inter) > 0 ){
      for (k in 1:nrow(O2valid_inter)){
        O2_inter <- O2valid_inter[k,]
        O2_fire  <- bps_fire %>%
          crop(O2_inter) %>%
          mask(O2_inter)
        O2_class <- bps_class %>%
          crop(O2_inter) %>%
          mask(O2_inter)
        ix <- getValues(O2_class)==2 
        O2_fire[ix] = dnbr_model %>% filter(level==O2_inter$NA_L2CODE) %>% pull(level2)  
        level_mix  <- merge(O2_fire,level_mix )
      } 
    }
    
    #Step 3: crop to the rest region
    shp_all <- shp_O2[which(!(shp_O2$NA_L2CODE %in% dnbr_model$level)),]
    
    if (nrow(shp_all) > 0 ) {
      all_fire  <- bps_fire %>%
        crop(shp_all) %>%
        mask(shp_all)
      
      all_class <- bps_class %>%
        crop(shp_all) %>%
        mask(shp_all)
      
      ix <- getValues(all_class)==2
      all_fire[ix] = dnbr_model %>% filter(level=="All") %>% pull(level2)  
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
    ix <- getValues(cbi_value) <= 0 
    cbi_value[ix] = 0.001
    
    ix <- getValues(cbi_value) > 3 &  getValues(cbi_value)!= 3.5
    cbi_value[ix] = 3
    writeRaster(cbi_value,paste0("Data/CBI_burn_severity/",fireID,"_cbi_value.tif"))
    
    sev_classified <- reclassify(cbi_value,class.matrix) # classify CBI value to burn severity classfication
    writeRaster(sev_classified,paste0("Data/CBI_burn_severity/",fireID,"_cbi_severity.tif"))
    
    cat_sev<- data.frame(rbind(table(values(sev_classified)))) #freq of different severity 
    colnames(cat_sev) = gsub("X", "Sev", colnames(cat_sev))
    cat_sev <- cat_sev %>%
      mutate(Event_ID = toupper(fireID),Total_cbi <- sum(cat_sev)*30*30/4047,
             burned_area = sum(cat_sev)*30*30/4047,
             index = i) 
    total_cbi <-  bind_rows(total_cbi,cat_sev)
    
    #add mask
    mask <- st_read(paste0(sevdir,"_mask.shp"))%>% 
      st_transform(crs(sev))
    #if there is a mask area, mask it out from the shape file
    if (nrow(mask)>0){
      st_agr(shp) = "constant"
      st_agr(mask)= "constant"
      final_shp <- st_difference(st_make_valid(shp), st_make_valid(st_union(mask)))
      
      sev_classified_crop_mask <- sev_classified%>% 
        crop(final_shp) %>%
        mask(final_shp) 
      writeRaster(sev_classified_crop_mask,paste0("Data/CBI_burn_severity/",fireID,"_cbi_severity_crop_mask.tif"))
      
      mask_sev<- data.frame(rbind(table(values(sev_classified_crop_mask)))) #freq of different severity of each fire
      
      colnames(mask_sev) = gsub("X", "Sev", colnames(mask_sev))
      
      mask_sev <- mask_sev %>%
        mutate(Event_ID = toupper(fireID),Total_cbi <- sum(mask_sev)*30*30/4047,
               burned_area = sum(mask_sev)*30*30/4047,
               index = i)
      
      total_cbi_mask <-  bind_rows(total_cbi_mask,mask_sev)
    }
    
    #plot CBI burn severity image and MTBS burn severity image for each wildfire event
    #prepare to plot CBI burn severity image
    cbi_reclass <- ratify(sev_classified) 
    rat <- levels(cbi_reclass )[[1]]
    rat$legend  <- c("Unchanged", "Low_severity","Medium_severity","High_severity","Grass burn")[cbi_rat$ID %in% rat$ID] 
    levels(sev_classified) <- rat
    my_col <- c("Green","cyan","yellow","red","purple")[cbi_rat$ID %in% rat$ID] 
    
    #prepare to plot MTBS burn severity image
    dnbr6 <- raster(paste0(sevdir,"_dnbr6.tif"))
    dnbr6_reclass <- ratify(dnbr6)
    rat6 <- levels(dnbr6_reclass)[[1]]
    rat6$legend  <- c("NA","Unburned to Low", "Low","Moderate","High","Increased Greenness","Non-processing Area Mask")[mtbs_rat$ID %in% rat6$ID] 
    levels(dnbr6) <- rat6
    mtbs_col <- c("grey","dark green","cyan","yellow","red","light green","white")[mtbs_rat$ID %in% rat6$ID]
    
    fire_name <- mtbs_perims$Incid_Name[which(mtbs_perims$Event_ID==toupper(fireID))]
    
    png(paste0("Data/Rplot/CBI_Severity_map_dnbr_19842019/",fireID,".png"), units="px", width=1600, height=1200, res=300)
    par(mfrow=c(1,2))
    plot(sev_classified,col=my_col,legend=F,box=F,axes=F) 
    legend("bottomright", legend =rat$legend, fill = my_col,xpd=TRUE,cex=0.4,bty='n',inset=c(0,0))
    plot(dnbr6,col=mtbs_col,legend=F,box=F,axes=F)
    legend("bottomright", legend =rat6$legend, fill = mtbs_col,xpd=TRUE,cex=0.4,bty='n',inset=c(0,0))
    mtext(fire_name, outer=TRUE,  cex=1, line=-3)
    dev.off()
    
  }, warning = function(w) {
    list(warningmessage = w$message)
  }, error = function(e) {
    list(errormessage = e$message)
  })
  myList[[i]] <- r
}

#convert area units from pixel numbers to Acres
total_cbi$unchanged <- fxns_pacre(total_cbi$Sev1)
total_cbi$low <- fxns_pacre(total_cbi$Sev2)
total_cbi$mod <- fxns_pacre(total_cbi$Sev3)
total_cbi$hig <- fxns_pacre(total_cbi$Sev4)
total_cbi$grass <- fxns_pacre(total_cbi$Sev5)

save(total_cbi,file="Data/Rdata/total_cbi.rda") #final burn severity data set without using MTBS mask files

total_cbi_mask$unchanged <- fxns_pacre(total_cbi_mask$Sev1)
total_cbi_mask$low <- fxns_pacre(total_cbi_mask$Sev2)
total_cbi_mask$mod <- fxns_pacre(total_cbi_mask$Sev3)
total_cbi_mask$hig <- fxns_pacre(total_cbi_mask$Sev4)
total_cbi_mask$grass <- fxns_pacre(total_cbi_mask$Sev5)
save(total_cbi_mask,file="Data/Rdata/total_cbi_mask.rda") #wildfire events with mask files
