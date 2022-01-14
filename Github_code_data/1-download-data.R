# This code is to download a few publicly available gis datasets using url for use in this project
# (1) LANDFIRE Biophysical Settings (BPS)
# (2) Monitoring Trends in Burn Severity (MTBS) Burned Areas Boundaries Dataset 
# (3) Ecoregions of North America (Omernik)
# (4) California state noundary
# (5) US boundary 

# Author: Qingqing Xu
# Conntact: qxu6@ucmerced.edu
# Date: 05/13/2021


options(timeout=1000) # The default time 60 seconds and that would cause Error in download.file "Timeout of 60 seconds was reached". So we set it to be 1000 seconds

fxns_dzip <- function(url,dir){
  temp=tempfile()
  download.file(url, temp)
  unzip(temp,exdir=dir)
} #a function that I created to download file from url and unzip it

#Do this once! Create Directories
dir.create("Data")
dir.create("Data/rasterfile")
dir.create("Data/shapefile")
dir.create("Data/shapefile/FRAP_Data")
dir.create("Data/shapefile/FRAP_Data/fire2019")
dir.create("Data/shapefile/FRAP_Data/fire2020")
#1.Load functions and packages

#(1) BPS remap raster data -- https://landfire.gov/version_download.php
link_bps_remap <- "https://landfire.gov/bulk/downloadfile.php?FNAME=US_200_mosaic-LF2016_BPS_200_CONUS.zip&TYPE=landfire"
sapply(link_bps_remap,fxns_dzip,dir="Data/rasterfile") 

#(2) MTBS Burned Areas Boundaries Dataset -- https://mtbs.gov/direct-download
link_mtbs <- "https://edcintl.cr.usgs.gov/downloads/sciweb1/shared/MTBS_Fire/data/composite_data/burned_area_extent_shapefile/mtbs_perimeter_data.zip"
sapply(link_mtbs,fxns_dzip,dir="Data/shapefile") 

#(3) Omernik vector data -- https://www.epa.gov/eco-research/ecoregions-north-america
omernik3 <- "https://gaftp.epa.gov/EPADataCommons/ORD/Ecoregions/cec_na/NA_CEC_Eco_Level3.zip"
sapply(link_omernik3,fxns_dzip,dir="Data/shapefile") # Omernick III layer includes Omernick I and Omernick II data

#(4) California state boundary -- https://data.ca.gov/dataset/ca-geographic-boundaries/resource/3db1e426-fb51-44f5-82d5-a54d7c6e188b
link_CA_boundary <- "https://data.ca.gov/dataset/e212e397-1277-4df3-8c22-40721b095f33/resource/3db1e426-fb51-44f5-82d5-a54d7c6e188b/download/ca-state-boundary.zip"
sapply(link_CA_boundary,fxns_dzip,dir="Data/shapefile")

#(5) US state boundary
link_US_state <- "https://www2.census.gov/geo/tiger/TIGER2017//STATE/tl_2017_us_state.zip"
sapply(link_US_state,fxns_dzip,dir="Data/shapefile")

# (6) CAL FIRE perimeters -- https://frap.fire.ca.gov/frap-projects/fire-perimeters/
link_frap <- "https://frap.fire.ca.gov/media/50dgwqrb/fire20_1.zip"
  
  "https://frap.fire.ca.gov/media/3nrpp42r/fire20_1.zip" #the link to download FRAP GIS data. For details: https://frap.fire.ca.gov/frap-projects/fire-perimeters/
sapply(link_frap,fxns_dzip,dir="Data/shapefile") #download zip file and unzip

ogrListLayers("Data/shapefile/fire20_1.gdb") #read gdb data
firep20_1 = st_read("Data/shapefile/fire20_1.gdb","firep20_1") #read interested data layer "firep20_1" which depicts wildfire perimeters from contributing agencies current as 2020
st_write(firep20_1, "Data/shapefile/firep20_1.shp")
