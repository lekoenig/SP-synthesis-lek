## StreamPULSE synthesis: Estimate site width
## Last updated June 2021
## LE Koenig

## The objective of this script is to identify potentially erroneous lat/lon in the streampulse synthesis dataset

# Load packages:
library(tidyr)         # data manipulation
library(dplyr)         # data manipulation
library(ggplot2)       # create plots
library(sf)            # work with geospatial data
library(mapview)       # plot geospatial data
library(tigris)        # interface with U.S. census data to retrieve CONUS shapefile


####################################################
##                    Load data                   ##
####################################################

# Load sites with estimated metabolism:
sites_dat_full <- readRDS("./data/lotic_site_info.rds") %>% rename("sitecode"="Site_ID")

# scan unique epsg codes:
print(unique(sites_dat_full$epsg_crs))

# split data by epsg code and recombine as spatial layer:
sites_full_nad83 <- sites_dat_full %>% filter(epsg_crs == 4269) %>% 
                    st_as_sf(.,coords=c("Lon","Lat"),crs=4269) %>% rename("epsg_crs_orig" = "epsg_crs")
sites_full_wgs84 <- sites_dat_full %>% filter(epsg_crs == 4326) %>% 
                    st_as_sf(.,coords=c("Lon","Lat"),crs=4326) %>% rename("epsg_crs_orig" = "epsg_crs") %>% st_transform(.,4269)
sites_full_sp <- bind_rows(sites_full_nad83,sites_full_wgs84) %>% mutate(epsg_crs_full = 4269)

# Load shapefile of U.S. and filter out states/territories outside of continental U.S. (CONUS):
states_omit <- c("American Samoa","Alaska","Hawaii","Commonwealth of the Northern Mariana Islands","United States Virgin Islands",
                "Guam","Puerto Rico")
conus_states <- states(cb = TRUE) %>% filter(.,!NAME %in% states_omit) %>%
                # transform to Albers Equal Area Conic and combine state geometries:              
                st_transform(.,5070) %>% st_union()

####################################################
##              Inspect site locations            ##
####################################################

# Which sites lie outside of CONUS (within 2 km)?
# Expectation: none, with possible exception of sites in Puerto Rico
sites_inspect <- sites_full_sp[which(st_is_within_distance(st_transform(sites_full_sp,5070),conus_states,sparse=FALSE,dist = 2000)=="FALSE"),]
print(sites_inspect[,c("sitecode","Name","Source")])

# Export sites to be inspected
sites_inspect_export <- sites_inspect %>% st_drop_geometry() %>% select(-epsg_crs_full)
write.csv(sites_inspect_export,paste("./output/intermediate/",format(Sys.Date(),"%Y%m%d"),"_check_site_locations.csv",sep=""),row.names = FALSE)

# Manually adjust erroneous CO sites:
sites_full_sp_fix <- sites_full_sp %>% filter(.,grepl("CO",sitecode)) %>%
                     mutate(lon = st_coordinates(.)[,1],
                            lat = st_coordinates(.)[,2]) %>%
                     st_drop_geometry() %>% 
                     # for now, manually adjust CO sites where a "2" was typed in for lat instead of a "3"
                     filter(lat < 30) %>% mutate(lat = lat+10) %>% 
                     st_as_sf(.,coords=c("lon","lat"),crs=.$epsg_crs_full[1])
                     
# Add flag to sites outside of CONUS and export:
sites_full_sp_export <- sites_full_sp %>% filter(!sitecode %in% sites_full_sp_fix$sitecode) %>%
                        bind_rows(.,sites_full_sp_fix) %>%
                        # site CO_Wfresdown1 seems to have two unique siteNames:
                        filter(!duplicated(sitecode)) %>% 
                        mutate(flag = ifelse(sitecode %in% sites_full_sp_fix$sitecode,"site location not where expected - lat/lon manually adjusted",NA))

saveRDS(sites_full_sp_export,"./output/intermediate/lotic_site_info_filtered.rds")




