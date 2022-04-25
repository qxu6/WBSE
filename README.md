# WBSE
Wildfire Burn Severity and Emission Inventory 
Contact: Qingqing Xu
Email: qxu6@ucmerced.edu
Citation: Xu Q, Westerling A L, Notohamiprodjo A, Wiedinmyer C, Picotte J J, Parks S A, Hurteau M D, Marlier M E, Kolden C A, Sam J A, Baldwin W J and Ade C (In Review). 
Wildfire Burn Severity and Emission Inventory: An example implementation over California. Environmental Research Letters.

Description: 

The Wildfire Burn Severity and Emission Inventory (WBSE) provides estimates of 30-m resolution burn severity, and emissions of CO2, CO,  CH4, non-methane organic compounds (NMOC), SO2, NH3, NO,  NO2, 
nitrogen oxides (NOx = NO + NO2), PM2.5, OC, and BC. Day of burning and daily emissions information is provided for each fire during 2002-2020 at 500-m resolution. 
We implemented WBSE for large wildfires (> 404ha) during 1984-2020 in the state of California, U.S., on a per-fire event scale since 1984 and a daily scale since 2002.  Fire records  for California from 1984-2019 
were retrieved from MTBS (https://mtbs.gov/viewer/index.html) via interactive viewer on May 8th, 2021, resulting in a dataset with a total of 1623 wildfires. We also acquired fire 
perimeters for wildfires (> 404 ha) in 2020 from CAL FIRE (https://frap.fire.ca.gov/frap-projects/fire-perimeters/) and calculated dNBR for each 2020 fire using the dNBR calculation tool with R and Google Earth Engine.

The final output datasets provided here are:

(1) Burn severity data

           Folder name: “burn_severity_CA”
                     
          
(2) Emission data
            
           Folder name: “emission_CA”
     

(3) Day of burning and daily emissions

      
           Folder name: “daily_emission_CA”
     

General information:

File names start with MTBS ID (format: US State Abbreviations (e.g. ca)  numbers) are calculated based on MTBS shapefiles and dNBR images. Total of 1623 large wildfires in the California region during 1984-2019. 
File names not in  MTBS ID format are the 2020 fires based on CAL FIRE shapefiles and our GEE dNBR images. Total of 72 large wildfires in the California region in 2020. 
Day of burning calculated for fires during 2002-2020 (a total of 841 large wildfires)  
Burn  severity classes: 
Raster value 0 -- non-processing area
Raster value 1 -- unburned
Raster value 2 -- low severity
Raster value 3 -- moderate severity
Raster value 4 -- high severity 
Raster value 5 -- grass burn
