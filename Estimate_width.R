## StreamPULSE synthesis: Estimate site width
## Last updated 18 May 2021
## LE Koenig

library(tidyr)        # data manipulation
library(dplyr)        # data manipulation
library(ggplot2)      # create plots
library(patchwork)    # format plots
library(waffle)       # create plots to visualize data sources
library(sf)           # work with geospatial data
library(nhdplusTools) # interface with national hydrography dataset

####################################################
##                    Load data                   ##
####################################################

# Load StreamPULSE synthesis sites:
sites_dat <- read.csv("./data/streampulse_synthesis_statset.csv",header=TRUE,stringsAsFactors = FALSE) %>% rename(.,"width_mcmanamay"="width")
sites_sp <- st_as_sf(sites_dat,coords=c("lon","lat"),crs = 4326)

# Load sites with estimated metabolism:
lotic_standardized_metabolism <- readRDS("./data/lotic_standardized_metabolism.rds")

# Load site data from Appling et al. 2018:
appling_info <- data.table::fread("./data/Appling2018/site_data.tsv",data.table=FALSE)

# Find sites that overlap synthesis and appling datasets:
site_list <- intersect(names(lotic_standardized_metabolism),appling_info[,"site_name"])

# Load regional sma relationships between wetted width and upstream area:
reg_coef <- read.csv("./data/RegionalHydraulicCoef.csv",header=TRUE,stringsAsFactors = FALSE) %>% rename(.,"coef_slope" = "slope","coef_int"="intercept")

# Find associated VPUID for each site:
  # Use hydroregions (note that a local copy of the NHDPlus VPU map is needed to run the lines below):
  #vpu_dat <- read_sf("/Volumes/T7/SpatialData/NHDPlus_HydroRegions/USA_HydroRegions_VPU02.shp") %>% st_transform(.,5070)
  #sites_sp$vpu <- vpu_dat$VPUID[unlist(st_is_within_distance(st_transform(sites_sp,5070),vpu_dat,sparse=TRUE,dist=1000))]
  #sites <- left_join(sites_dat,st_drop_geometry(sites_sp)[,c("sitecode","vpu")],by="sitecode")

  # OR use nhdplus flowlines (advantage: don't have to have a local copy of nhdplus dataset):
    # 1. find comid:  
    for(i in 1:length(sites_sp$sitecode)){
      sites_sp$comid[i] <- discover_nhdplus_id(sites_sp[i,])
    } 

    # 2. download flowlines for each comid identified above and extract vpu from VAA table:
    subset <- subset_nhdplus(comids = sites_sp$comid,
                            output_file = tempfile(fileext = ".gpkg"),
                            nhdplus_data = "download", 
                            flowline_only = TRUE,
                            return_data = TRUE)
    nhd_dat <- st_drop_geometry(subset$NHDFlowline_Network)
    sites <- left_join(st_drop_geometry(sites_sp),nhd_dat[,c("comid","vpuid")],by="comid")
    sites <- left_join(sites,appling_info[,c("site_name","struct.dam_flag")],by=c("sitecode"="site_name"))
    
# Load manually-estimated widths and join with synthesis site data:
    # Google Earth:
    widths_manual_GE <- read.csv("./data/manual_subset/sites_manual_width.csv",header=TRUE)
    widths_GE <- widths_manual_GE %>% select(!c(lat,lon,comments)) %>% 
                  pivot_longer(!sitecode, names_to = "width_xsection", values_to = "value") %>% 
                  group_by(sitecode) %>% summarize(width_m_GE = median(value,na.rm=T))
    
    # ArcMap:
    widths_manual_arc <- read.csv("./data/manual_subset/sites_manual_width_arc.csv",header=TRUE)
    widths_arc <- widths_manual_arc %>% select(!OID) %>% rename("width_xsection"="xsection","value"="width_m") %>%
              group_by(sitecode) %>% summarize(width_m_arc = median(value,na.rm=T))
    
    # StreamPULSE data portal:
    files <- list.files("./data/SP_portal_site_data/",pattern="geomorph")
    
    SP_dat_ls <- list()
    
    for(i in 1:length(files)){
      site <- substr(files[i], start = 1, stop = 2)
      dat <- read.csv(paste("./data/SP_portal_site_data/",files[i],sep=""),header=TRUE,skip = 1) %>% rename("sitecode"="siteID")
      dat_summary <- dat[,c("sitecode","date","wetted_width_m")] %>% group_by(sitecode,date) %>% 
        summarize(width_m = median(wetted_width_m,na.rm=T)) %>% group_by(sitecode) %>%
        summarize(width_m_SP = median(width_m,na.rm=T))
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
    epscor_dat <- read.csv("./data/NH_EPSCoR_site_data/EPSCoR_Stream_Dims.csv",header=TRUE) %>% filter(Site=="DCF") %>%
      group_by(Date) %>% summarize(MedWidth_m = median(Width_cm,na.rm=TRUE)) %>% summarize(width_m_SP = median(MedWidth_m,na.rm=TRUE)/100) %>% mutate(sitecode = "NH_DCF") %>% select(sitecode,width_m_SP)
    
    # NC site data not on the data portal (NC_NHC and NC_UNHC):
    NC_dat <- read.csv("./data/NC_site_data/NC_SP_field_widths.csv",header=TRUE) %>% group_by(sitecode) %>% 
              summarize(width_m_SP = median(Wetted_width_m,na.rm=TRUE))
    
    # Combine SP_dat with new EPSCoR data above:
    SP_dat2 <- bind_rows(SP_dat,epscor_dat,NC_dat)
    
    # Combine manual width data:
    widths_manual <- left_join(widths_GE,widths_arc,by="sitecode") %>% left_join(.,SP_dat2) %>% bind_rows(.,SP_dat2[! SP_dat2$sitecode %in% .$sitecode,])
    
# Compare manual width estimates (since GE estimates of distance may not be very accurate):
plot(width_m_GE~width_m_arc,data = widths_manual)
abline(a=0,b=1)

sites <- left_join(sites,widths_manual,by="sitecode")
head(sites)

####################################################
##             Calculate NWIS widths              ##
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
length(which(!is.na(sites$width_nwis))) # this is 201 sites/215


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

sites <- left_join(sites,reg_coef[,c("HUC02","coef_slope","coef_int")],by=c("vpuid"="HUC02"))

# Coefficients estimated based on EPA NRSA and WSA data, see: https://github.com/lekoenig/US-hydraulic-geometry.git
sites$width_reglEPA <- 10^((log10(sites$ws_area_km2) * sites$coef_slope)+sites$coef_int)

# For site in Florida with 0 upstream area, set width to NA:
sites$width_reglEPA[which(sites$ws_area_km2==0)] <- NA

####################################################
##   Randomly select sites for manual validation  ##
####################################################

#sites.sub <- sites[sample(c(1:length(sites$sitecode)),50,replace = FALSE),c("sitecode","lat","lon")]
#write.csv(sites.sub,"./data/manual_subset/sites_manual_width.csv",row.names = FALSE)


####################################################
##      Compare widths against validation data    ##
####################################################

# define the manual data we're using:
sites$width_manual <- ifelse(!is.na(sites$width_m_SP),sites$width_m_SP,sites$width_m_GE)

# add colors for sites with dam flag > 50:
sites$dam_color <- ifelse(!is.na(sites$struct.dam_flag) & sites$struct.dam_flag>=50,"flag","ok")

# plot estimated widths against manual ~validation data:
plot.mcman <- sites %>% ggplot() + geom_point(aes(x=width_manual,y=width_mcmanamay,color=dam_color),alpha=.75) + geom_abline(slope=1,intercept=0,lty=2)+
      theme_classic() + coord_cartesian(xlim=c(0,300),ylim=c(0,300)) + scale_color_manual(values=c("#009F73","black"))+
      ggtitle("McManamay") + labs(x=expression(width~manual~(m)),y=expression(width~predicted~(m))) + theme(legend.position="none",title = element_text(size=10))
plot.nwis <- sites %>% ggplot() + geom_point(aes(x=width_manual,y=width_nwis,color=dam_color),alpha=.75) + geom_abline(slope=1,intercept=0,lty=2) + 
      theme_classic()+ coord_cartesian(xlim=c(0,300),ylim=c(0,300)) + scale_color_manual(values=c("#009F73","black"))+
      ggtitle("NWIS coefficients") + labs(x=expression(width~manual~(m)),y=expression(width~predicted~(m))) + theme(legend.position="none",title = element_text(size=10))
plot.raymond <- sites %>% ggplot() + geom_point(aes(x=width_manual,y=width_raymond,color=dam_color),alpha=.75) + geom_abline(slope=1,intercept=0,lty=2) +
      theme_classic() + coord_cartesian(xlim=c(0,300),ylim=c(0,300)) + scale_color_manual(values=c("#009F73","black"))+
      ggtitle("Raymond coefficients") + labs(x=expression(width~manual~(m)),y=expression(width~predicted~(m))) + theme(legend.position="none",title = element_text(size=10))
plot.natlEPA <- sites %>% ggplot() + geom_point(aes(x=width_manual,y=width_natlEPA,color=dam_color),alpha=.75) + geom_abline(slope=1,intercept=0,lty=2)+
      theme_classic() + coord_cartesian(xlim=c(0,300),ylim=c(0,300)) + scale_color_manual(values=c("#009F73","black"))+
      ggtitle("Natl EPA coefficients") + labs(x=expression(width~manual~(m)),y=expression(width~predicted~(m))) + theme(legend.position="none",title = element_text(size=10)) 
plot.reglEPA <- sites %>% ggplot() + geom_point(aes(x=width_manual,y=width_reglEPA,color=dam_color),alpha=.75) + geom_abline(slope=1,intercept=0,lty=2) +
      theme_classic() + coord_cartesian(xlim=c(0,300),ylim=c(0,300)) + scale_color_manual(values=c("#009F73","black"))+
      ggtitle("Regional EPA coefficients") + labs(x=expression(width~manual~(m)),y=expression(width~predicted~(m))) + theme(legend.position="none",title = element_text(size=10))
plots <- plot.mcman + plot.nwis + plot.raymond + plot.natlEPA + plot.reglEPA + plot_layout(ncol=3)
print(plots)

# Compare linear regression and RMSE:
fits.df <- data.frame(Width_est_method = c("McManamay","NWIS coefficients","Raymond coefficients","Natl EPA coefficients","Regional EPA coefficients"),
                      slope = NA,
                      r2 = NA,
                      RMSE = NA)

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
  fits.df$RMSE[i] <- round(sqrt(mean(fit$residuals^2)),1)
}
print(fits.df)

# If we estimate width using NWIS at-a-site relationships, how many sites are we missing widths for?
n.SP <- length(which(SP_dat2$sitecode %in% sites$sitecode))
sites$source <- ifelse(sites$sitecode %in% SP_dat2$sitecode,"StreamPULSE estimates",ifelse(sites$sitecode %in% appling_info$site_name,"NWIS estimates","No estimate"))
counts <- sites %>% group_by(source) %>% summarise(n = n()) %>%  mutate(percent = round(n/sum(n)*100))
case_counts <- counts$n
names(case_counts) <- counts$source
print(case_counts)

waffle::waffle(case_counts,rows=10)

# Which sites coming from the StreamPULSE data portal don't have corresponding widths?  
sites %>% filter(!grepl('nwis', sitecode)) %>% select(sitecode) %>% filter(! sitecode %in% SP_dat2$sitecode)



