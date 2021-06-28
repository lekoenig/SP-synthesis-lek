## StreamPULSE synthesis: Estimate site width
## Last updated June 2021
## LE Koenig

## The objective of this script is to join metabolism site locations to national hydrography dataset (nhdplusv2) flowlines

# Load packages:
library(tidyr)         # data manipulation
library(dplyr)         # data manipulation
library(ggplot2)       # create plots
library(sf)            # work with geospatial data
library(mapview)       # plot geospatial data
library(nhdplusTools)  # interface with national hydrography dataset


####################################################
##                    Load data                   ##
####################################################

# Load sites with estimated metabolism:
sites_full_sp <- readRDS("./output/intermediate/lotic_site_info_full_sitelistqaqc.rds") 
unique(sites_full_sp$flag)
dim(sites_full_sp)

# Load site data from Appling et al. 2018:
appling_info <- data.table::fread("./data/Appling2018/site_data.tsv",data.table=FALSE)


####################################################
##       Join NHD: nhdplusTools helper fxns       ##
####################################################

## 1. Find nearest COMID (method used in prior synthesis compilation):

  # find nhdplus flowlines comid using nhdplusTools discover_nhdplus_id:  
  for(i in seq_along(sites_full_sp$sitecode)){
    sites_full_sp$comid[i] <- discover_nhdplus_id(sites_full_sp[i,])
    print(i)
  } 
  sites_full_sp <- sites_full_sp %>% rename("comid_orig" = "comid") %>% select(sitecode, comid_orig,VPU)

  
## 2. Find nearest flowline comid by adjusting discover_nhdplus_id function to look for nwis sites:
  
  sites_full_sp_nwis <- sites_full_sp %>% filter(grepl("nwis",sitecode)) %>% st_drop_geometry() %>% select(sitecode)
  
  for(i in seq_along(sites_full_sp_nwis$sitecode)){
    sitecode_reformat <- paste("USGS-",substr(sites_full_sp_nwis$sitecode[i],start=6,stop=30),sep="")
    nldi_nwis <- list(featureSource = "nwissite", featureID = sitecode_reformat)
    comid <- suppressMessages(discover_nhdplus_id(nldi_feature = nldi_nwis))
    sites_full_sp_nwis$comid_checknwis[i] <- ifelse(length(comid)>0,comid,NA)
    print(i)
  }
  
  
# Access seamless nhdplusv2 (download national seamless data layer if no local copy is present) and check nearest comid:
#dir.create("./output/intermediate/nhd")
#nhdplusTools::download_nhdplusv2("./output/intermediate/nhd")
#st_layers(dsn="./output/intermediate/nhd/NHDPlusNationalData/NHDPlusV21_National_Seamless_Flattened_Lower48.gdb")
flowline <- read_sf(dsn = "./output/intermediate/nhd/NHDPlusNationalData/NHDPlusV21_National_Seamless_Flattened_Lower48.gdb",
                      layer = "NHDFlowline_Network")
  
## 3. Find nearest flowline using get_flowline_index (need local copy of nhdplus flowlines):
  
  nwis_name_strings <- function(sitecode){ # returns a vector of strings to use in pattern matching
    
    omit_strings <- c("river","creek","crk","at","in","near","below","nr","bl","@","r","c","r.","fk","fk.","above","f","run","abv",
                      "s","e","e.","n","w","l","blw","br","rv","ab","a","bk","branch","east","west","south","north","wb","run,","site",
                      "swamp","n.f.","s.f.","eb","creek,","fork","cr","m","st.","ck","of","site-at","creek-upstream","(upper")
    site_nm <- substr(sitecode,start=6,stop=30)
    nwis_nm <- dataRetrieval::readNWISsite(site_nm) %>% select(station_nm) %>% as.character() %>% tolower()
    nwis_strings <- tolower(unlist(strsplit(nwis_nm,split=" ")))[c(1:3)] %>% 
                    #stringr::str_remove(.,stringr::fixed(paste(omit_strings, collapse = "|"),ignore_case = TRUE)) %>%
                    .[. != ""] %>% .[!.%in% omit_strings] 
    return(nwis_strings)
  }

  comid_fline_ls <- list()
  
  for(i in seq_along(sites_full_sp$sitecode)){
    vpu <- sites_full_sp$VPU[i]  
    comid_check <- get_flowline_index(st_transform(filter(flowline,VPUID==vpu),5070),
                                st_transform(sites_full_sp[i,],5070),
                                search_radius = 2000,max_matches=3)

    if(grepl("nwis",sites_full_sp$sitecode[i],ignore.case=TRUE)){
      comid_name_match_strings <- grepl(paste(nwis_name_strings(sites_full_sp$sitecode[i]),collapse="|"),
                                flowline$GNIS_NAME[which(flowline$COMID %in% comid_check$COMID)],ignore.case=TRUE)
      if("TRUE" %in% comid_name_match_strings){
        comid_name_match <- flowline$COMID[which(flowline$COMID %in% comid_check$COMID)][which(comid_name_match_strings=="TRUE")]
        comid_check_filter <- comid_check[which(comid_check$COMID %in% comid_name_match),]
        comid <- comid_check_filter[1,]
      } else {
        comid <- comid_check[1,]
      }
    } else {
      comid <- comid_check[1,]
    }
    }
      comid$sitecode <- sites_full_sp$sitecode[i] 
      comid_fline_ls[[i]] <- comid

    print(i)
  }
  comid_check_fline <- do.call("rbind",comid_fline_ls) %>% select(sitecode,COMID,REACHCODE,offset) %>% 
                       rename("comid_flineindex"="COMID","reachcode_flineindex"="REACHCODE") #%>%
                       # replace comid with NA if near_dist > 200:
                       #mutate_at(c("comid_flineindex","reachcode_flineindex","offset"),~ifelse(offset > 00, NA, .))
    
  
####################################################
##     Join NHD - set max distance to flowline    ##
####################################################  
  
# Find nearest flowline comid by setting max distance (need local copy of NHDPlusV2, download below):
  
# start function to link point to comid (and define distance to flowline):
link_comid <- function(points,flines,vpu,purge_non_dendritic){ 
  
  # Prep flowline:
  names(flines)[1:(length(names(flines))-1)] <- tolower(names(flines)[1:(length(names(flines))-1)])
  flowline.sub <- flines %>% filter(vpuid == as.character(vpu)) %>% st_zm(.,drop=TRUE,what="ZM")
  
  if(purge_non_dendritic == "TRUE"){
    flowline.sub <- filter(flowline.sub, streamorde == streamcalc)
  }
  
  if(grepl("nwis",points$sitecode)){
  
    # Find flowline segment distances from point (automatically uses great circle distances for geographic coords):
    dist <- st_distance(x=st_transform(points,st_crs(flowline.sub)),y=flowline.sub,by_element=FALSE)
    candidate_reaches <- order(dist,decreasing=FALSE)[c(1:3)]
    candidate_reaches_dist <- dist[order(dist,decreasing=FALSE)[c(1:3)]]
    candidate_name_match <- grepl(paste(nwis_name_strings(points$sitecode),collapse="|"),
                                  flowline.sub$gnis_name[candidate_reaches],ignore.case = TRUE)
    
    cv_candidate_dist <- as.numeric(sd(candidate_reaches_dist)/mean(candidate_reaches_dist))
    
    # if there is big discrepancy between top-3 nearest flowlines, choose closest one:
    if(cv_candidate_dist < 0.51){
      
    # check if one of the nearby flowlines matches the gnis name:
    if("TRUE" %in% candidate_name_match){
      reach_matches <- candidate_reaches[which(candidate_name_match=="TRUE")]
      chosen_reach <- ifelse(length(reach_matches)>1,reach_matches[which.min(dist[reach_matches])],reach_matches)
      chosen_fline <- flowline.sub[chosen_reach,]
      near <- dist[chosen_reach]
    } else {
      chosen_fline <- flowline.sub[which.min(dist),]
      near <- dist[which.min(dist)]
    }
    } else {
      chosen_fline <- flowline.sub[which.min(dist),]
      near <- dist[which.min(dist)] 
    }
  } else {

  # Find flowline segment distances from point (automatically uses great circle distances for geographic coords):
  dist <- st_distance(x=st_transform(points,st_crs(flowline.sub)),y=flowline.sub,by_element=FALSE)
  chosen_fline <- flowline.sub[which.min(dist),]
  near <- dist[which.min(dist)]
  }

  # Return nearest flowline segment:
  join.dat <- bind_cols(points %>% select(sitecode),
                        chosen_fline) %>% st_drop_geometry() %>% 
              select(sitecode,comid,streamcalc,streamorde,totdasqkm,qe_ma,vpuid) %>% 
              mutate(near_dist_m = as.numeric(near))
  return(join.dat)
} # end function

# Find nearest comid for each point:
comid_ls <- list()
for(i in seq_along(sites_full_sp$sitecode)){
  comid_ls[[i]] <- link_comid(points = select(sites_full_sp,sitecode)[i,],
                              vpu = sites_full_sp$VPU[i],
                              flines=flowline, 
                              purge_non_dendritic = "FALSE")
  print(i)
}
comid_check <- do.call("rbind",comid_ls) %>% select(sitecode,comid,near_dist_m) %>% rename("comid_checkdist"="comid") #%>%
               # replace comid with NA if near_dist > 200:
               #mutate_at(c("comid_checkdist","near_dist_m"),~ifelse(near_dist_m > 200, NA, .))



####################################################
##              Compare snapped comids            ##
####################################################

check_divergence <- function(comid){
  div <- flowline$Divergence[which(flowline$COMID==comid)]
}

check_dendritic <- function(comid){
  dendr <- ifelse(flowline$StreamCalc[which(flowline$COMID==comid)]==flowline$StreamOrde[which(flowline$COMID==comid)],"yes","no")
}

# compare against results from discover_nhdplus_id above:
compare_comid <- left_join(st_drop_geometry(sites_full_sp),comid_check,by="sitecode") %>%
                 left_join(.,sites_full_sp_nwis,by="sitecode") %>%
                 left_join(.,appling_info[,c("site_name","nhdplus_id")],by=c("sitecode"="site_name")) %>%
                 left_join(.,comid_check_fline,by="sitecode")

# check whether comid's are within the dendritic network, or are divergences:
compare_comid$orig_div <- mapply(check_divergence,compare_comid$comid_orig)
compare_comid$orig_dendr <- mapply(check_dendritic,compare_comid$comid_orig)

# compare orig drain area to nwis reported drain area as a qc check:
compare_comid$comid_orig_areakm2 <- NA
compare_comid$nwis_areakm2 <- NA

for(i in seq_along(compare_comid$sitecode)){
compare_comid$comid_orig_areakm2[i] <- flowline$TotDASqKM[which(flowline$COMID==compare_comid$comid_orig[i])]

if(grepl("nwis",compare_comid$sitecode[i])){
  name <- substr(compare_comid$sitecode[i],start=6,stop=35)
  area_mi2 <- dataRetrieval::readNWISsite(name) %>% select(drain_area_va) %>% as.numeric()
  area_km2 <- area_mi2 * 2.58999
  compare_comid$nwis_areakm2[i] <- area_km2
} else {
  compare_comid$nwis_areakm2[i] <- NA
}
print(i)
}

# filter based on pct diff between two areas:
compare_comid$area_pctdiff <- (compare_comid$nwis_areakm2-compare_comid$comid_orig_areakm2)/compare_comid$nwis_areakm2 * 100
compare_areas <- compare_comid[which(abs(compare_comid$area_pctdiff)>90),]

# export:
saveRDS(compare_comid[,-15],"./output/intermediate/compare_comids.rds")

# (manually) examine any differences:
comid_diff <- compare_comid %>% filter(comid_orig!=comid_checkdist |
                                       comid_orig!=comid_flineindex |
                                       is.na(comid_checkdist) | is.na(comid_flineindex) |
                                       near_dist_m > 200 | offset > 200)

i <- 1
point <- sites_full_sp[which(sites_full_sp$sitecode == comid_diff$sitecode[i]),]
diff <- comid_diff[i,] %>% rename("comid_nwis"="comid_checknwis")
print(diff)

mapview(point,col.region="red") + 
mapview(flowline[which(flowline$COMID==diff$comid_orig),],color="blue")+
mapview(flowline[which(flowline$COMID==diff$nhdplus_id),],color="orange")+
mapview(flowline[which(flowline$COMID==diff$comid_flineindex),],color="purple")+
mapview(flowline[which(flowline$COMID==diff$comid_checkdist),],color="green")


####################################################
##             Export flagged site codes          ##
####################################################

flagged_sites <- c("nwis_01540500","nwis_03298200","nwis_06805500","nwis_07381600",
           "nwis_11126000","nwis_14202650","nwis_12101100","nwis_01184000",
           "nwis_13206400","MD_POBR","NH_BDC","NH_MCQ","NH_HBF")

# Gather watershed areas for sites with manually-adjusted comid (and no nwis area):
area_nwis12101100 <- read_sf("./data/Appling2018/catchment_shapefile/catchment_shapefile.shp") %>% filter(site_nm == "nwis_12101100") %>%
                     st_transform(.,5070) %>% mutate(area_km2 = st_area(.)/(1000*1000)) %>% select(area_km2) %>% st_drop_geometry() %>% 
                     as.numeric() %>% round(.,digits=2)
area_nwis13206400 <- read_sf("./data/Appling2018/catchment_shapefile/catchment_shapefile.shp") %>% filter(site_nm == "nwis_13206400") %>%
                     st_transform(.,5070) %>% mutate(area_km2 = st_area(.)/(1000*1000)) %>% select(area_km2) %>% st_drop_geometry() %>% 
                     as.numeric() %>% round(.,digits=2)
area_nwis07381600 <- read_sf("./data/Appling2018/catchment_shapefile/catchment_shapefile.shp") %>% filter(site_nm == "nwis_07381600") %>%
                     st_transform(.,5070) %>% mutate(area_km2 = st_area(.)/(1000*1000)) %>% select(area_km2) %>% st_drop_geometry() %>% 
                     as.numeric() %>% round(.,digits=2)
area_NH_BDC <- read_sf("./data/NH_site_data/NH_BDC_StreamStats/globalwatershed.shp") %>% st_transform(.,5070) %>%
               mutate(area_km2 = st_area(.)/(1000*1000)) %>% select(area_km2) %>% st_drop_geometry() %>% as.numeric() %>% round(.,digits=2)
area_NH_MCQ <- read_sf("./data/NH_site_data/NH_MCQ_MShattuck/MCQ.shp") %>% st_transform(.,5070) %>%
               mutate(area_km2 = st_area(.)/(1000*1000)) %>% select(area_km2) %>% st_drop_geometry() %>% as.numeric() %>% round(.,digits=2)
area_NH_HBF <- read_sf("./data/NH_site_data/NH_HBF_LTER/hbef_wsheds/hbef_wsheds.shp") %>% filter(WS=="WS3") %>% st_transform(.,5070) %>%
               mutate(area_km2 = st_area(.)/(1000*1000)) %>% select(area_km2) %>% st_drop_geometry() %>% as.numeric() %>% round(.,digits=2)
area_MD_POBR <- dataRetrieval::readNWISsite("01583570") %>% mutate(area_km2 = drain_area_va*2.58999) %>% select(area_km2) %>% as.numeric() %>% round(.,digits=2)

# Prep data frame containing info. for manually-adjusted sites:
compare_comid_subflags <- compare_comid %>% filter(sitecode %in% flagged_sites) %>%
                          mutate(
                          comid_manual = case_when(
                            sitecode=="nwis_01540500" ~ as.character(2604977),
                            sitecode=="nwis_03298200" ~ as.character(10264164),
                            sitecode=="nwis_06805500" ~ as.character(17416032),
                            sitecode=="nwis_11126000" ~ as.character(17610189),
                            sitecode=="nwis_14202650" ~ as.character(23805136),
                            sitecode=="nwis_12101100" ~ as.character(23982435),
                            sitecode=="nwis_13206400" ~ as.character(23400535),
                            sitecode=="nwis_01184000" ~ as.character(7700898),
                            sitecode=="nwis_07381600" ~ as.character(15182073),
                            sitecode=="MD_POBR" ~ "NA",
                            sitecode=="NH_BDC" ~ "NA",
                            sitecode=="NH_MCQ" ~ "NA",
                            sitecode=="NH_HBF" ~ "NA"),
                          manual_areakm2 = case_when(
                              sitecode=="nwis_01540500" ~ flowline$TotDASqKM[which(flowline$COMID==2604977)],
                              sitecode=="nwis_03298200" ~ flowline$TotDASqKM[which(flowline$COMID==10264164)],
                              sitecode=="nwis_06805500" ~ flowline$TotDASqKM[which(flowline$COMID==17416032)],
                              sitecode=="nwis_11126000" ~ flowline$TotDASqKM[which(flowline$COMID==17610189)],
                              sitecode=="nwis_14202650" ~ flowline$TotDASqKM[which(flowline$COMID==23805136)],
                              sitecode=="nwis_12101100" ~ area_nwis12101100,
                              sitecode=="nwis_13206400" ~ area_nwis13206400,
                              sitecode=="nwis_01184000" ~ flowline$TotDASqKM[which(flowline$COMID==7700898)],
                              sitecode=="nwis_07381600" ~ area_nwis07381600,
                              sitecode=="MD_POBR" ~ area_MD_POBR,
                              sitecode=="NH_BDC" ~ area_NH_BDC,
                              sitecode=="NH_MCQ" ~ area_NH_MCQ,
                              sitecode=="NH_HBF" ~ area_NH_HBF),
                          manual_area_src = case_when(
                            sitecode=="nwis_12101100" ~ "Appling2018_StreamStats",
                            sitecode=="nwis_13206400" ~ "Appling2018_StreamStats",
                            sitecode=="nwis_07381600" ~ "Appling2018_USGS2013",
                            sitecode=="nwis_01540500" ~ "NHDV2_vaa",
                            sitecode=="nwis_03298200" ~ "NHDV2_vaa",
                            sitecode=="nwis_06805500" ~ "NHDV2_vaa",
                            sitecode=="nwis_11126000" ~ "NHDV2_vaa",
                            sitecode=="nwis_14202650" ~ "NHDV2_vaa",
                            sitecode=="nwis_01184000" ~ "NHDV2_vaa",
                            sitecode=="MD_POBR" ~ "nwis_site_description",
                            sitecode=="NH_BDC" ~ "StreamStats",
                            sitecode=="NH_MCQ" ~ "localuser_UNHWQAL",
                            sitecode=="NH_HBF" ~ "localuser_HBFLTER")
                            ) %>%
                          select(-c(reachcode_flineindex,orig_div,orig_dendr))

# Export data frame containing info. for manually-adjusted sites:
write.csv(compare_comid_subflags,paste("./output/",format(Sys.Date(),"%Y%m%d"),"_flagged_sites_checkcomid.csv",sep=""),row.names = FALSE)


#########################################################
##     Export new COMIDs based on manual adjustments   ##
#########################################################

# Assign all comids (either keep original comid assignment or, in the case of flagged sites, replace with comid_manual):
for(i in seq_along(compare_comid$sitecode)){
compare_comid$comid_define[i] <- ifelse(compare_comid$sitecode[i] %in% compare_comid_subflags$sitecode,
                                     compare_comid_subflags$comid_manual[which(compare_comid_subflags$sitecode==compare_comid$sitecode[i])],
                                     compare_comid$comid_orig[i])
}
compare_comid$comid_define <- as.integer(compare_comid$comid_define)

# Look up vpu for adjusted comid's:
nhd_subset <- subset_nhdplus(comids = compare_comid$comid_define[which(!is.na(compare_comid$comid_define))],
                             output_file = tempfile(fileext = ".gpkg"),
                             nhdplus_data = "download", 
                             flowline_only = TRUE,
                             return_data = TRUE)
nhd_dat <- st_drop_geometry(nhd_subset$NHDFlowline_Network)

# Prep data frame for export:
comid_fixes_export <- left_join(compare_comid,nhd_dat[,c("comid","totdasqkm","vpuid")],by=c("comid_define"="comid")) %>%
                        mutate(vpuid = ifelse(is.na(vpuid),VPU,vpuid)) %>%
                        select(sitecode,comid_orig,comid_define,vpuid,totdasqkm) %>%
                        rename("comid_old" = "comid_orig","comid_new"="comid_define","totdasqkm_comidnew" = "totdasqkm")
# Add ws areas:
for(i in seq_along(comid_fixes_export$sitecode)){
  comid_fixes_export$ws_area_km2[i] = ifelse(comid_fixes_export$sitecode[i] %in% compare_comid_subflags$sitecode,
                                      compare_comid_subflags$manual_areakm2[which(compare_comid_subflags$sitecode==comid_fixes_export$sitecode[i])],
                                      comid_fixes_export$totdasqkm_comid[i])
  comid_fixes_export$ws_area_src[i] = ifelse(comid_fixes_export$sitecode[i] %in% compare_comid_subflags$sitecode,
                                      compare_comid_subflags$manual_area_src[which(compare_comid_subflags$sitecode==comid_fixes_export$sitecode[i])],
                                      "NHDV2_vaa")
}
                        
# Export data frame containing comid fixes:
write.csv(comid_fixes_export,paste("./output/intermediate/",format(Sys.Date(),"%Y%m%d"),"_synthesis_fullset_comidfixes.csv",sep=""),row.names = FALSE)

