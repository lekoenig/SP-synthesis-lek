## StreamPULSE synthesis: Estimate site width
## Last updated June 2021
## LE Koenig

# Load packages:
library(tidyr)         # data manipulation
library(dplyr)         # data manipulation
library(ggplot2)       # create plots
library(ggrepel)       # format plots
library(patchwork)     # format plots
library(waffle)        # create plots to visualize data sources
library(sf)            # work with geospatial data
library(nhdplusTools)  # interface with national hydrography dataset
library(dataRetrieval) # interface with NWIS
library(dataMeta)      # helper functions to build metadata/data dictionary files

options(scipen=999)

####################################################
##                    Load data                   ##
####################################################

# Load sites within StreamPULSE data portal (as of 25 May 2021):
sites_portal <- read.csv("./data/SP_portal_site_data/all_basic_site_data-2.csv",header=TRUE) %>% 
                select(regionID,siteID,siteName,dataSource,latitude,longitude,USGSgageID) %>% 
                mutate(sitecode = paste(regionID,siteID,sep="_"))
  # Assign gage number for MD site that is co-located with USGS gage (indicated in an embargoed data set and confirmed by lat/lon):
  sites_portal$USGSgageID[which(sites_portal$sitecode=="MD_DRKR")] <- "01589330"

  # Restore leading zeros for usgs gage names shorter than 8-character minimum:
  sites_portal$USGSgageID_fix <- ifelse(nchar(sites_portal$USGSgageID)==7,paste(0,sites_portal$USGSgageID,sep=""),sites_portal$USGSgageID)

# Load StreamPULSE synthesis sites (stat subset):
sites_dat <- read.csv("./data/streampulse_synthesis_statset.csv",header=TRUE,stringsAsFactors = FALSE) %>% rename(.,"width_mcmanamay"="width")
sites_sp <- st_as_sf(sites_dat,coords=c("lon","lat"),crs = 4326) 

# Load sites with estimated metabolism:
lotic_standardized_metabolism <- readRDS("./data/lotic_standardized_metabolism.rds")
sites_dat_full <- readRDS("./output/intermediate/lotic_site_info_filtered.rds") 
unique(sites_dat_full$flag)
# Create spatial df for full synthesis dataset:
sites_full_sp <- sites_dat_full

# Load site data from Appling et al. 2018:
appling_info <- data.table::fread("./data/Appling2018/site_data.tsv",data.table=FALSE)

# Find associated VPUID for each site (needed for regional hydraulic geometry relationships below):
    # 1. find nhdplus flowlines comid:  
    for(i in seq_along(sites_full_sp$sitecode)){
      sites_full_sp$comid[i] <- discover_nhdplus_id(sites_full_sp[i,])
      print(i)
    } 

    # 2. download flowlines for each comid identified above and extract vpu from VAA table:
    subset <- subset_nhdplus(comids = sites_full_sp$comid,
                            output_file = tempfile(fileext = ".gpkg"),
                            nhdplus_data = "download", 
                            flowline_only = TRUE,
                            return_data = TRUE)
    nhd_dat <- st_drop_geometry(subset$NHDFlowline_Network)
    
    # 3. Join sites to NHD information, flags indicating influence of dams, and portal information indicating co-located usgs gage:
    sites_nhd <- left_join(st_drop_geometry(sites_full_sp),nhd_dat[,c("comid","vpuid","totdasqkm")],by="comid") %>% 
                 left_join(.,appling_info[,c("site_name","struct.dam_flag")],by=c("sitecode"="site_name")) %>%
                 left_join(.,sites_portal[,c("sitecode","USGSgageID_fix")],by="sitecode") %>%
                 left_join(.,sites_dat[,c("sitecode","width_mcmanamay")],by="sitecode") %>%
                 # Add gage name in addition to sitecode, which indicates the name used in the metabolism synthesis:
                 mutate(USGSgage_name = ifelse(grepl("nwis",sitecode),substr(sitecode,start=6,stop=30),USGSgageID_fix)) %>%
                 rename("ws_area_km2" = "totdasqkm")
    
# Load manually-estimated widths and join with synthesis site data:
    # Google Earth:
    widths_manual_GE <- read.csv("./data/manual_width_by_imagery/sites_manual_width_ge.csv",header=TRUE)
    widths_GE <- widths_manual_GE %>% select(!c(lat,lon,comments)) %>% 
                  pivot_longer(!sitecode, names_to = "width_xsection", values_to = "value") %>% 
                  group_by(sitecode) %>% summarize(width_m_GE = median(value,na.rm=T))
    
    # ArcMap:
      # Load initial subset (n=50):
      widths_manual_arc <- read.csv("./data/manual_width_by_imagery/sites_manual_width_arc.csv",header=TRUE)
      widths_arc <- widths_manual_arc %>% select(!OID) %>% rename("width_xsection"="xsection","value"="width_m") %>%
                    group_by(sitecode) %>% summarize(width_m_arc = median(value,na.rm=T))
    
    # USGS in situ field measurements:
      # Define function to retrieve usgs-measured widths:
      get_usgs_widths <- function(metab_sitename,gage_sitename){
        
        # find relevant dates:
        metab_dat <- lotic_standardized_metabolism[[which(names(lotic_standardized_metabolism)==metab_sitename)]]
        metab_rng <- range(metab_dat$Date)
        
        # Convert m3/s to cfs and find discharge range: 
        # *Note that Q data for MD_DRKR appears to be in L/s rather than m3/s:
        if(metab_sitename == "MD_DRKR"){
          q_rng <- range(metab_dat$discharge/1000*35.3147,na.rm=T)
        } else {
          q_rng <- range(metab_dat$discharge*35.3147,na.rm=T)
        }
        # Find field measurements within metabolism date/discharge ranges:
        width_df <- readNWISmeas(siteNumbers=gage_sitename,expanded = TRUE) %>%
                    filter(measurement_dt > metab_rng[1] & measurement_dt < metab_rng[2]) %>%
                    #readNWISmeas(siteNumbers=gage_sitename,expanded = TRUE,startDate = metab_rng[1],endDate = metab_rng[2]) %>%
                    # filter measurements that fall within range of metabolism discharge:
                    filter(discharge_va > q_rng[1] & discharge_va < q_rng[2]) %>%
                    # summarize in situ widths:
                    select(measurement_dt,chan_width) %>% mutate(width_m = chan_width*0.3048) %>%
                    summarize(width_m = median(width_m,na.rm=T),n = n())

        return(width_df)
      }
      
      # Retrieve usgs-measured widths for each site (if available):
      usgs_dat <- data.frame(sitecode = sites_nhd$sitecode,
                             gage_sitename = sites_nhd$USGSgage_name,
                             width_m_usgs = rep(NA,length(sites_nhd$sitecode)),
                             n_dates_meas_usgs = rep(NA,length(sites_nhd$sitecode)))
        
      for(i in 1:length(usgs_dat$gage_sitename)){

        width_df <- tryCatch(
          {get_usgs_widths(metab_sitename = usgs_dat$sitecode[i],
                           gage_sitename = usgs_dat$gage_sitename[i])},
          error=function(cond){return(NA)},
          warning=function(cond){return(NA)})
      
      usgs_dat$width_m_usgs[i] <- ifelse(length(width_df)>1,width_df$width_m,NA)
      usgs_dat$n_dates_meas_usgs[i] <- ifelse(length(width_df)>1,width_df$n,NA)
      print(i)
      }

    # StreamPULSE data portal:
    files <- list.files("./data/SP_portal_site_data/",pattern="geomorph")
    
    SP_dat_ls <- list()
    
    for(i in 1:length(files)){
      site <- substr(files[i], start = 1, stop = 2)
      dat <- read.csv(paste("./data/SP_portal_site_data/",files[i],sep=""),header=TRUE,skip = 1) %>% rename("sitecode"="siteID")
      dat_summary <- dat[,c("sitecode","date","wetted_width_m")] %>% group_by(sitecode,date) %>% 
                    summarize(width_m = median(wetted_width_m,na.rm=T)) %>% group_by(sitecode) %>%
                    summarize(width_m_SP = median(width_m,na.rm=T),n_dates_meas_SP=n())
      for(j in 1:length(dat_summary$sitecode)){
        dat_summary$sitecode[j] <- ifelse(grepl(site, dat_summary$sitecode[j], fixed = TRUE),
                                          dat_summary$sitecode[j],
                                          paste(site,"_",dat_summary$sitecode[j],sep=""))
      }
      SP_dat_ls[[i]] <- dat_summary
    }
    
    SP_dat <- do.call(rbind,SP_dat_ls)
    SP_dat$width_m_SP[which(is.nan(SP_dat$width_m_SP))] <- NA
    
    # NH EPSCoR site not on the data portal (NH_DCF):  
    epscor_dat <- read.csv("./data/NH_site_data/EPSCoR_Stream_Dims.csv",header=TRUE) %>%
                  group_by(Site,Date) %>% summarize(MedWidth_m = median(Width_cm/100,na.rm=TRUE)) %>% 
                  summarize(width_m_SP = median(MedWidth_m,na.rm=TRUE),n_dates_meas_SP = n()) %>% 
                  mutate(sitecode = paste("NH_",Site,sep=""))
    
    # NC site data not on the data portal (NC_NHC and NC_UNHC):
    NC_dat <- read.csv("./data/NC_site_data/NC_SP_field_widths.csv",header=TRUE) %>% group_by(sitecode) %>% 
              summarize(width_m_SP = median(Wetted_width_m,na.rm=TRUE)) %>% mutate(n_dates_meas_SP = 1)
    
    # FL site data not on the data portal (FL_SF2800, and FL_WS1500):
    FL_dat <- read.csv("./data/FL_site_data/FL_SP_field_widths.csv",header=TRUE) %>% group_by(sitecode) %>% 
              summarize(width_m_SP = median(Baseflow_width_m,na.rm=TRUE)) %>% mutate(n_dates_meas_SP = 1)
      
    # Combine SP_dat with new EPSCoR data above:
    SP_dat2 <- bind_rows(SP_dat,epscor_dat,NC_dat,FL_dat)
    
    # Combine manual width data:
    widths_manual <- full_join(widths_arc,widths_GE,by="sitecode") %>%
                     full_join(.,SP_dat2,by="sitecode") %>%
                     full_join(.,usgs_dat,by="sitecode")

# Compare manual width estimates (since GE estimates of distance may not be very accurate):
plot(width_m_GE~width_m_arc,data = widths_manual,xlim=c(0,160))
abline(a=0,b=1)
summary(lm(width_m_GE~width_m_arc,data=widths_manual))

# Compare manual width estimates to in situ widths reported by USGS:
plot(width_m_usgs ~ width_m_arc,data=widths_manual,ylim=c(0,max(widths_manual$width_m_arc,na.rm=T)))
abline(a=0,b=1)
summary(lm(width_m_usgs~width_m_arc,data=widths_manual)) 

# Create sites data frame that includes field-measured widths:
sites <- left_join(sites_nhd,widths_manual,by="sitecode")
head(sites)

# Join SP (measured) widths to co-located USGS gage locations (WI_BRW and WI_BEC):
sites$width_m_SP[which(sites$sitecode=="nwis_05406469")] <- sites$width_m_SP[which(sites$sitecode=="WI_BRW")]
sites$width_m_SP[which(sites$sitecode=="nwis_05406457")] <- sites$width_m_SP[which(sites$sitecode=="WI_BEC")]


####################################################
##             Estimate NWIS widths               ##
####################################################

# Define function to calculate width from median discharge using Jud coefficients in Appling dataset:
calc_width_nwis <- function(site){
  
  # Get the site of interest:
  SOI <- lotic_standardized_metabolism[[site]]
  
  # Calculate median discharge:
  medq <- median(SOI[,"discharge"],na.rm=TRUE)
  
  # Get the empirical hydraulic geometry parameters for the site:
  site_params <- appling_info[which(appling_info$site_name==site),]
  a_coeff <- site_params[,"dvqcoefs.a"]
  b_coeff <- site_params[,"dvqcoefs.b"]
  
  # width = a * (Q^b)
  if(length(a_coeff) + length(b_coeff) == 0){
    width <- NA
  }
  
  if(length(a_coeff) + length(b_coeff) != 0){
    width <- width <- a_coeff * (medq^b_coeff)
  }
  
  # Return width:
  return(width)
  
} # End

sites$width_nwis <- apply(sites["sitecode"],1,calc_width_nwis)
length(which(!is.na(sites$width_nwis))) 


####################################################
##                 Widths: Raymond                ##
####################################################

a_coeff_raymond <- 12.88
b_coeff_raymond <- 0.42

calc_width_raymond <- function(site){
  
  # Get the site of interest:
  SOI <- lotic_standardized_metabolism[[site]]
  
  # Calculate median discharge:
  medq <- median(SOI[,"discharge"],na.rm=TRUE)

  # width = a * (Q^b)
  if(is.na(a_coeff_raymond) + is.na(b_coeff_raymond) == 0){
    width <- a_coeff_raymond * (medq^b_coeff_raymond)
  }
  
  if(is.na(a_coeff_raymond) + is.na(b_coeff_raymond) != 0){
    width <- NA
  }
  
  # Return width:
  return(width)
  
} # End

sites$width_raymond <- apply(sites["sitecode"],1,calc_width_raymond)


####################################################
##             Widths: EPA National               ##
####################################################

# Coefficients estimated based on EPA NRSA and WSA data, see: https://github.com/lekoenig/US-hydraulic-geometry.git
sites$width_natlEPA <- 10^((log10(sites$ws_area_km2) * 0.455)-0.06)

# For site in Florida with 0 upstream area, set width to NA:
sites$width_natlEPA[which(sites$ws_area_km2==0)] <- NA


####################################################
##             Widths: EPA Regional               ##
####################################################

# Load regional sma relationships between wetted width and upstream area
# Coefficients estimated based on EPA NRSA and WSA data, see: https://github.com/lekoenig/US-hydraulic-geometry.git
reg_coef <- read.csv("./data/RegionalHydraulicCoef.csv",header=TRUE,stringsAsFactors = FALSE) %>% rename(.,"coef_slope" = "slope","coef_int"="intercept")

# Join sites with regional coefficients and estimate width:
sites <- left_join(sites,reg_coef[,c("HUC02","coef_slope","coef_int")],by=c("vpuid"="HUC02"))
sites$width_reglEPA <- 10^((log10(sites$ws_area_km2) * sites$coef_slope)+sites$coef_int)

# For site in Florida with 0 upstream area, set width to NA:
sites$width_reglEPA[which(sites$ws_area_km2==0)] <- NA


####################################################
##      Compare widths against validation data    ##
####################################################

# define the manual data we're using:
sites$width_manual <- ifelse(is.na(sites$width_m_usgs),sites$width_m_SP,sites$width_m_usgs)
  
# add colors for sites with dam flag > 50:
sites$dam_color <- ifelse(!is.na(sites$struct.dam_flag) & sites$struct.dam_flag>=50,"flag","ok")

# plot estimated widths against manual ~validation data:
plot.mcman <- ggplot() + geom_point(data=sites,aes(x=width_manual,y=width_mcmanamay,color=dam_color),alpha=.75) + geom_abline(slope=1,intercept=0,lty=2)+
      theme_classic() +  scale_color_manual(values=c("#009F73","black"))+
      ggtitle("McManamay") + labs(x=expression(width~manual~(m)),y=expression(width~predicted~(m))) + theme(legend.position="none",title = element_text(size=10)) +
      scale_x_log10() + scale_y_log10()+
      coord_cartesian(xlim=c(0.1,1000),ylim=c(0.1,1000)) +
      #geom_text_repel(size=3,data=subset(sites, width_manual > 100),aes(x = width_manual,y=width_mcmanamay,label=substr(sitecode, start = 6,stop = 30))) + 
      NULL
plot.nwis <- sites %>% ggplot() + geom_point(aes(x=width_manual,y=width_nwis,color=dam_color),alpha=.75) + geom_abline(slope=1,intercept=0,lty=2) + 
      theme_classic()+ scale_color_manual(values=c("#009F73","black"))+ 
      coord_cartesian(xlim=c(0.1,1000),ylim=c(0.1,1000)) +
      ggtitle("NWIS coefficients") + labs(x=expression(width~manual~(m)),y=expression(width~predicted~(m))) + 
      theme(legend.position="none",title = element_text(size=10)) + 
      scale_x_log10() + scale_y_log10()+
      #geom_text_repel(size=3,data=subset(sites, width_manual > 100),aes(x = width_manual,y=width_nwis,label=substr(sitecode, start = 6,stop = 30))) + 
      NULL
plot.raymond <- sites %>% ggplot() + geom_point(aes(x=width_manual,y=width_raymond,color=dam_color),alpha=.75) + geom_abline(slope=1,intercept=0,lty=2) +
      theme_classic() +  scale_color_manual(values=c("#009F73","black"))+
      ggtitle("Raymond coefficients") + labs(x=expression(width~manual~(m)),y=expression(width~predicted~(m))) + theme(legend.position="none",title = element_text(size=10)) +
      scale_x_log10() + scale_y_log10()+
      coord_cartesian(xlim=c(0.1,1000),ylim=c(0.1,1000)) +
      #geom_text_repel(size=3,data=subset(sites, width_manual > 100),aes(x = width_manual,y=width_raymond,label=substr(sitecode, start = 6,stop = 30))) + 
      NULL
plot.natlEPA <- sites %>% ggplot() + geom_point(aes(x=width_manual,y=width_natlEPA,color=dam_color),alpha=.75) + geom_abline(slope=1,intercept=0,lty=2)+
      theme_classic() + scale_color_manual(values=c("#009F73","black"))+
      ggtitle("Natl EPA coefficients") + labs(x=expression(width~manual~(m)),y=expression(width~predicted~(m))) + theme(legend.position="none",title = element_text(size=10)) +
      scale_x_log10() + scale_y_log10() +
      coord_cartesian(xlim=c(0.1,1000),ylim=c(0.1,1000)) + 
      #geom_text_repel(size=3,data=subset(sites, width_manual > 100),aes(x = width_manual,y=width_natlEPA,label=substr(sitecode, start = 6,stop = 30))) + 
      NULL
plot.reglEPA <- sites %>% ggplot() + geom_point(aes(x=width_manual,y=width_reglEPA,color=dam_color),alpha=.75) + geom_abline(slope=1,intercept=0,lty=2) +
      theme_classic() + scale_color_manual(values=c("#009F73","black"))+
      ggtitle("Regional EPA coefficients") + labs(x=expression(width~manual~(m)),y=expression(width~predicted~(m))) + theme(legend.position="none",title = element_text(size=10)) +
      scale_x_log10() + scale_y_log10()+
      coord_cartesian(xlim=c(0.1,1000),ylim=c(0.1,1000)) +
      #geom_text_repel(size=3,data=subset(sites, width_manual > 100),aes(x = width_manual,y=width_reglEPA,label=substr(sitecode, start = 6,stop = 30))) + 
      NULL
plots <- plot.mcman + plot.nwis + plot.raymond + plot.natlEPA + plot.reglEPA + plot_layout(ncol=3)
print(plots)

# Compare linear regression and RMSE:
  fits.df <- data.frame(Width_est_method = c("McManamay","NWIS coefficients","Raymond coefficients","Natl EPA coefficients","Regional EPA coefficients"),
                      slope = NA,
                      r2 = NA,
                      RMSE = NA,MAE = NA,Mean_bias = NA)

  for(i in 1:length(fits.df$Width_est_method)){
    width_man <- sites$width_manual
    width_pred <- case_when(fits.df$Width_est_method[i] == "McManamay" ~ sites$width_mcmanamay,
                            fits.df$Width_est_method[i] == "NWIS coefficients" ~ sites$width_nwis,
                            fits.df$Width_est_method[i] == "Raymond coefficients" ~ sites$width_raymond,
                            fits.df$Width_est_method[i] == "Natl EPA coefficients" ~ sites$width_natlEPA,
                            fits.df$Width_est_method[i] == "Regional EPA coefficients" ~ sites$width_reglEPA)
    fit <- lm(width_pred ~ width_man)
    fits.df$slope[i] <- round(as.numeric(fit$coefficients[2]),2)
    fits.df$r2[i] <- round(as.numeric(summary(fit)$r.squared),2)
    width_diff <- width_pred - width_man
    fits.df$RMSE[i] <- round(sqrt(mean((width_diff)^2,na.rm=T)),1)
    fits.df$MAE[i] <- round(mean(abs(width_diff),na.rm=T),1)
    fits.df$Mean_bias[i] <- round(mean(width_diff,na.rm=T),1)
  }
  print(fits.df[order(fits.df$RMSE,decreasing=FALSE),])

  
####################################################
##     Export width estimates for SP synthesis    ##
####################################################  
  
# If we estimate width using usgs measurements, how many sites are we missing widths for?
n.SP <- length(which(!is.na(sites$width_m_SP)))
if(length(sites$sitecode[which(is.na(sites$width_m_SP) & is.na(sites$width_m_usgs))])>0){
sites$source <- ifelse(sites$sitecode %in% usgs_dat$sitecode[which(!is.na(usgs_dat$width_m_usgs))],"NWIS field measurements",ifelse(sites$sitecode %in% SP_dat2$sitecode[which(!is.na(SP_dat2$width_m_SP))],"StreamPULSE estimates","Need other estimate"))
} else {
sites$source <- ifelse(sites$sitecode %in% usgs_dat$sitecode[which(!is.na(usgs_dat$width_m_usgs))],"NWIS field measurements","StreamPULSE estimates")
}
counts <- sites %>% group_by(source) %>% summarise(n = n()) %>%  mutate(percent = round(n/sum(n)*100))
case_counts <- counts$n
names(case_counts) <- counts$source
print(case_counts)

waffle::waffle(parts=case_counts,rows=10,colors = c("#66c2a5", "#fc8d62","#8da0cb"))

# Which sites coming from the StreamPULSE data portal don't have corresponding widths?  
sites %>% filter(!grepl('nwis', sitecode)) %>% select(sitecode) %>% filter(! sitecode %in% SP_dat2$sitecode)

# Which sites will we need to estimate by satellite imagery or hydraulic geometry?
need_est <- sites[which(is.na(sites$width_m_SP)&is.na(sites$width_m_usgs)),c("sitecode")]
length(need_est) # out of 421 sites in the full dataset
print(need_est)

# Assign estimated widths and export table:
sites_export <- sites %>% select(sitecode,gage_sitename,width_manual,width_reglEPA,source) %>%
                mutate(usgs_gagename = ifelse(!is.na(gage_sitename),paste("nwis_",gage_sitename,sep=""),NA)) %>% 
                mutate(width_export = ifelse(is.na(width_manual),width_reglEPA,width_manual)) %>%
                select(-c(gage_sitename,width_manual,width_reglEPA)) %>%
                mutate(source = gsub("Need other estimate","Regional geomorphic scaling coeff", source)) %>%
                rename("width_src" = "source","width_m" = "width_export")

# Examine distribution:
quantile(sites_export$width_m,na.rm=T)

sites %>% select(sitecode,width_manual,width_mcmanamay) %>% 
          rename("New"="width_manual","McManamay"="width_mcmanamay") %>% 
          pivot_longer(!sitecode,names_to = "width_est", values_to = "width_m") %>%
          ggplot() + geom_histogram(aes(x=width_m,fill=width_est),color="black") + scale_x_log10() + 
          ggtitle(label="Distribution of estimated widths for synthesis dataset",subtitle =  "note log scale") + 
          scale_fill_manual(values=c("#5ab4ac","#d8b365"),name="Width estimate") + theme_classic() 
                    
## ============================================================ ##
##                     Create a metadata file                   ##
## ============================================================ ##

# Define variable descriptions:
var_names <- c("Site name used in metabolism synthesis (Bernhardt et al.)","Associated USGS gage name","Estimated site width in meters","Width data source")
             
# Define variable types (1=include all options for given variable):              
var_type <- c(0,0,1,0)

# Define categorical variables (variables where type == 1):
opt_descrip2 <- c("Site name used in metabolism synthesis (Bernhardt et al.)","Associated USGS gage name","Estimated site width in meters","Width data source","Width data source")

# Create data dictionary:
linker <- build_linker(sites_export, variable_description = var_names, variable_type = var_type)
data_dictionary <- build_dict(my.data = sites_export, linker = linker, prompt_varopts = FALSE) 

# Add data set attributes and join attributes with data set:
data_description = "This data set contains width estimates for the river locations within the StreamPULSE synthesis statset (Bernhardt et al.)"
sites_export2 <- incorporate_attr(my.data = sites_export, data.dictionary = data_dictionary, main_string = data_description)

sites_export_join_ls <- list(names=attributes(sites_export2)$names,
                     class = attributes(sites_export2)$class,
                     main = attributes(sites_export2)$main,
                     dictionary = attributes(sites_export2)$dictionary,
                     last_edit_date = attributes(sites_export2)$last_edit_date,
                     author = attributes(sites_export2)$author,
                     data = sites_export2)

# Export data:
saveRDS(sites_export_join_ls,paste("./output/",format(Sys.Date(),"%Y%m%d"),"_synthesis_fullset_width.rds",sep=""))
write.csv(sites_export_join_ls$data,paste("./output/",format(Sys.Date(),"%Y%m%d"),"_synthesis_fullset_width.csv",sep=""),row.names=FALSE)

# Start writing to an output file
sink(paste("./output/",format(Sys.Date(),"%Y%m%d"),"_synthesis_fullset_width_metadata.txt",sep=""))

cat("=============================\n")
cat("Data set title\n")
cat("=============================\n")
sites_export_join_ls$main

cat("\n")
cat("\n")

cat("=============================\n")
cat("Author\n")
cat("=============================\n")
sites_export_join_ls$author

cat("\n")
cat("\n")

cat("=============================\n")
cat("Last edit date\n")
cat("=============================\n")
sites_export_join_ls$last_edit_date

cat("\n")
cat("\n")

cat("=============================\n")
cat("Column names\n")
cat("=============================\n")
sites_export_join_ls$names

# Stop writing to the file
sink()                             



