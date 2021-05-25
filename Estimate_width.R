## StreamPULSE synthesis: Estimate site width
## Last updated 25 May 2021
## LE Koenig

library(tidyr)         # data manipulation
library(dplyr)         # data manipulation
library(ggplot2)       # create plots
library(ggrepel)       # format plots
library(patchwork)     # format plots
library(waffle)        # create plots to visualize data sources
library(sf)            # work with geospatial data
library(nhdplusTools)  # interface with national hydrography dataset
library(dataRetrieval) # interface with NWIS

####################################################
##                    Load data                   ##
####################################################

# Sites within StreamPULSE data portal (as of 25 May 2021):
sites_portal <- read.csv("./data/SP_portal_site_data/all_basic_site_data-2.csv",header=TRUE) %>% 
                select(regionID,siteID,siteName,dataSource,latitude,longitude,USGSgageID) %>% 
                mutate(sitecode = paste(regionID,siteID,sep="_"))

# Excel deleted leading zeros for most usgs gages:
sites_portal$USGSgageID_fix <- ifelse(nchar(sites_portal$USGSgageID)==7,paste(0,sites_portal$USGSgageID,sep=""),sites_portal$USGSgageID)

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
  # Use nhdplus flowlines (advantage: don't have to have a local copy of nhdplus dataset):
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
    sites_nhd <- left_join(st_drop_geometry(sites_sp),nhd_dat[,c("comid","vpuid")],by="comid") %>% 
                 left_join(.,appling_info[,c("site_name","struct.dam_flag")],by=c("sitecode"="site_name")) %>%
                 left_join(.,sites_portal[,c("sitecode","USGSgageID_fix")],by="sitecode") %>%
                 # Add gage name in addition to sitecode, which indicates the name used in the metabolism synthesis:
                 mutate(USGSgage_name = ifelse(grepl("nwis",sitecode),substr(sitecode,start=6,stop=30),USGSgageID_fix))
    
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
    
      # Load widths for statset:
      widths_manual_arc_all <- read.csv("./data/manual_width_by_imagery/widths_manual_arc_statset.csv",header=TRUE)
      
      # Add widths for cross-sections with multiple widths (e.g. split channels) and then summarize:
      widths_arc_all <- widths_manual_arc_all %>% select(!OID) %>% rename("width_xsection"="xsection","value"="width_m") %>%
                        group_by(sitecode,width_xsection) %>% summarize(width_m_xsection = sum(value,na.rm=T)) %>% 
                        group_by(sitecode) %>% summarize(width_m_arc = median(width_m_xsection,na.rm=T),n=n())
      # Remove St. Johns River (arc est. not reliable):
      widths_arc_all$width_m_arc[which(widths_arc_all$sitecode=="nwis_02234000")] <- NA
      
    # USGS:
      # Define function to retrieve usgs-measured widths:
      get_usgs_widths <- function(metab_sitename,gage_sitename){
        
        # find relevant dates:
        metab_dat <- lotic_standardized_metabolism[[which(names(lotic_standardized_metabolism)==metab_sitename)]]
        metab_rng <- range(metab_dat$Date)
        
        width_df <- readNWISmeas(siteNumbers=gage_sitename,expanded = TRUE,startDate = metab_rng[1],endDate = metab_rng[2]) %>%
                    select(measurement_dt,chan_width) %>% mutate(width_m = chan_width*0.3048) %>%
                    summarize(width_m = median(width_m,na.rm=T))
                    
        # find width measurements from usgs:
        #width_df <- readNWISqw(siteNumbers = site2,parameterCd = "00004", 
        #                       startDate = metab_rng[1],endDate = metab_rng[2]) %>%
        #            filter(parm_cd=="00004") %>% mutate(sitecode = paste("nwis_",site_no,sep="")) %>%
        #            mutate(width_m_all = result_va * 0.3048) %>%
        #            summarize(width_m = median(width_m_all,na.rm=T)) 
                 
        return(as.numeric(width_df$width_m))
    
      }
      
      # Retrieve usgs-measured widths for each site (if available):
      
      usgs_dat2 <- data.frame(sitecode=sites_nhd$sitecode,
                              metabname = sites_nhd$sitecode,
                             width_m_usgs = rep(NA,length(sites_nhd$sitecode)))
      usgs_dat2$metabname[which(usgs_dat2$sitecode=="MD_DRKR")] <- "nwis_01589330"
      
      for(i in 1:length(usgs_dat$sitecode)){
        
        usgs_dat2$width_m_usgs[i] <- tryCatch(
                                     {get_usgs_widths(site = usgs_dat$sitecode[i])},
                                      error=function(cond){return(NA)},
                                      warning=function(cond) {return(NA)})
        print(i)
      }
      # remove width from nwis_05599490 (width measured during flood stage):
      usgs_dat2 <- filter(usgs_dat2,sitecode != "nwis_05599490")
       
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
    epscor_dat <- read.csv("./data/NH_site_data/EPSCoR_Stream_Dims.csv",header=TRUE) %>% filter(Site=="DCF") %>%
                  group_by(Date) %>% summarize(MedWidth_m = median(Width_cm,na.rm=TRUE)) %>% 
                  summarize(width_m_SP = median(MedWidth_m,na.rm=TRUE)/100) %>% 
                  mutate(sitecode = "NH_DCF") %>% select(sitecode,width_m_SP)
    
    # NC site data not on the data portal (NC_NHC and NC_UNHC):
    NC_dat <- read.csv("./data/NC_site_data/NC_SP_field_widths.csv",header=TRUE) %>% group_by(sitecode) %>% 
              summarize(width_m_SP = median(Wetted_width_m,na.rm=TRUE))
    
    # FL site data not on the data portal (FL_SF2500, FL_SF2800, and FL_WS1500):
      # FL_SF2800 and SF_WS1500:
      FL_dat1 <- read.csv("./data/FL_site_data/FL_SP_field_widths.csv",header=TRUE) %>% group_by(sitecode) %>% 
              summarize(width_m_SP = median(Baseflow_width_m,na.rm=TRUE))
      # FL_SF2500:
      FL_dat_usgs_headers <- read.csv("./data/FL_site_data/SF2500_channel_USGS.csv",nrows=1,header=FALSE,skip=14)
      FL_dat_usgs <- read.csv("./data/FL_site_data/SF2500_channel_USGS.csv",header=TRUE,skip=16)
      names(FL_dat_usgs) <- FL_dat_usgs_headers
      FL_dat_usgs$date <- as.Date(FL_dat_usgs$measurement_dt,format="%m/%d/%Y")
      FL_metab_daterg <- range(lotic_standardized_metabolism$FL_SF2500$Date)
      FL_metab_qrg <- range(lotic_standardized_metabolism$FL_SF2500$discharge[which(!is.na(lotic_standardized_metabolism$FL_SF2500$GPP))])*35.3147 # present m3/s Q in cfs
      FL_dat2 <- filter(FL_dat_usgs,date > FL_metab_daterg[1] & date < FL_metab_daterg[2] & chan_discharge > FL_metab_qrg[1] & chan_discharge < FL_metab_qrg[2]) %>% 
               summarize(width_ft = median(chan_width,na.rm=TRUE)) %>% 
               mutate(sitecode = "FL_SF2500",width_m_SP = width_ft*0.3048) # convert channel width from ft to meters
    FL_dat <- bind_rows(FL_dat1,FL_dat2[,c("sitecode","width_m_SP")])
      
    # Combine SP_dat with new EPSCoR data above:
    SP_dat2 <- bind_rows(SP_dat,epscor_dat,NC_dat,FL_dat)
    
    # Combine manual width data:
    widths_manual <- left_join(widths_GE,widths_arc,by="sitecode") %>% 
                     left_join(.,SP_dat2) %>% 
                     bind_rows(.,SP_dat2[! SP_dat2$sitecode %in% .$sitecode,]) %>%
                     left_join(.,usgs_dat,by="sitecode")
    
# Compare manual width estimates (since GE estimates of distance may not be very accurate):
plot(width_m_GE~width_m_arc,data = widths_manual,xlim=c(0,160))
abline(a=0,b=1)
summary(lm(width_m_GE~width_m_arc,data=widths_manual))

# Compare manual width estimates to in situ widths reported by USGS:
plot(width_m_usgs ~ width_m_arc,data=widths_manual,xlim=c(0,100))
abline(a=0,b=1)

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
sites$width_manual <- ifelse(!is.na(sites$width_m_SP),sites$width_m_SP,sites$width_m_arc)

# add colors for sites with dam flag > 50:
sites$dam_color <- ifelse(!is.na(sites$struct.dam_flag) & sites$struct.dam_flag>=50,"flag","ok")

# plot estimated widths against manual ~validation data:
plot.mcman <- ggplot() + geom_point(data=sites,aes(x=width_manual,y=width_mcmanamay,color=dam_color),alpha=.75) + geom_abline(slope=1,intercept=0,lty=2)+
      theme_classic() + coord_cartesian(xlim=c(0.1,300),ylim=c(0.1,300)) + scale_color_manual(values=c("#009F73","black"))+
      ggtitle("McManamay") + labs(x=expression(width~manual~(m)),y=expression(width~predicted~(m))) + theme(legend.position="none",title = element_text(size=10)) +
      #scale_x_log10() + scale_y_log10()+
      geom_text_repel(size=3,data=subset(sites, width_manual > 100),aes(x = width_manual,y=width_mcmanamay,label=substr(sitecode, start = 6,stop = 30))) + 
      NULL
plot.nwis <- sites %>% ggplot() + geom_point(aes(x=width_manual,y=width_nwis,color=dam_color),alpha=.75) + geom_abline(slope=1,intercept=0,lty=2) + 
      theme_classic()+ scale_color_manual(values=c("#009F73","black"))+ coord_cartesian(xlim=c(0.1,300),ylim=c(0.1,300)) +
      ggtitle("NWIS coefficients") + labs(x=expression(width~manual~(m)),y=expression(width~predicted~(m))) + 
      theme(legend.position="none",title = element_text(size=10)) + 
      #scale_x_log10() + scale_y_log10()+
      geom_text_repel(size=3,data=subset(sites, width_manual > 100),aes(x = width_manual,y=width_nwis,label=substr(sitecode, start = 6,stop = 30))) + 
      NULL
plot.raymond <- sites %>% ggplot() + geom_point(aes(x=width_manual,y=width_raymond,color=dam_color),alpha=.75) + geom_abline(slope=1,intercept=0,lty=2) +
      theme_classic() + coord_cartesian(xlim=c(0.1,300),ylim=c(0.1,300)) + scale_color_manual(values=c("#009F73","black"))+
      ggtitle("Raymond coefficients") + labs(x=expression(width~manual~(m)),y=expression(width~predicted~(m))) + theme(legend.position="none",title = element_text(size=10)) +
      #scale_x_log10() + scale_y_log10()+
      geom_text_repel(size=3,data=subset(sites, width_manual > 100),aes(x = width_manual,y=width_raymond,label=substr(sitecode, start = 6,stop = 30))) + 
      NULL
plot.natlEPA <- sites %>% ggplot() + geom_point(aes(x=width_manual,y=width_natlEPA,color=dam_color),alpha=.75) + geom_abline(slope=1,intercept=0,lty=2)+
      theme_classic() + coord_cartesian(xlim=c(0.1,300),ylim=c(0.1,300)) + scale_color_manual(values=c("#009F73","black"))+
      ggtitle("Natl EPA coefficients") + labs(x=expression(width~manual~(m)),y=expression(width~predicted~(m))) + theme(legend.position="none",title = element_text(size=10)) +
      #scale_x_log10() + scale_y_log10() +
      geom_text_repel(size=3,data=subset(sites, width_manual > 100),aes(x = width_manual,y=width_natlEPA,label=substr(sitecode, start = 6,stop = 30))) + 
      NULL
plot.reglEPA <- sites %>% ggplot() + geom_point(aes(x=width_manual,y=width_reglEPA,color=dam_color),alpha=.75) + geom_abline(slope=1,intercept=0,lty=2) +
      theme_classic() + coord_cartesian(xlim=c(0.1,300),ylim=c(0.1,300)) + scale_color_manual(values=c("#009F73","black"))+
      ggtitle("Regional EPA coefficients") + labs(x=expression(width~manual~(m)),y=expression(width~predicted~(m))) + theme(legend.position="none",title = element_text(size=10)) +
      #scale_x_log10() + scale_y_log10()+
      geom_text_repel(size=3,data=subset(sites, width_manual > 100),aes(x = width_manual,y=width_reglEPA,label=substr(sitecode, start = 6,stop = 30))) + 
      NULL
plots <- plot.mcman + plot.nwis + plot.raymond + plot.natlEPA + plot.reglEPA + plot_layout(ncol=3)
print(plots)

# Compare linear regression and RMSE:
compare.fits <- function(df){
  fits.df <- data.frame(Width_est_method = c("McManamay","NWIS coefficients","Raymond coefficients","Natl EPA coefficients","Regional EPA coefficients"),
                      slope = NA,
                      r2 = NA,
                      RMSE = NA,MAE = NA,Mean_bias = NA)

  for(i in 1:length(fits.df$Width_est_method)){
    width_man <- df$width_manual
    width_pred <- case_when(fits.df$Width_est_method[i] == "McManamay" ~ df$width_mcmanamay,
                            fits.df$Width_est_method[i] == "NWIS coefficients" ~ df$width_nwis,
                            fits.df$Width_est_method[i] == "Raymond coefficients" ~ df$width_raymond,
                            fits.df$Width_est_method[i] == "Natl EPA coefficients" ~ df$width_natlEPA,
                            fits.df$Width_est_method[i] == "Regional EPA coefficients" ~ df$width_reglEPA)
    fit <- lm(width_pred ~ width_man)
    fits.df$slope[i] <- round(as.numeric(fit$coefficients[2]),2)
    fits.df$r2[i] <- round(as.numeric(summary(fit)$r.squared),2)
    width_diff <- width_pred - width_man
    fits.df$RMSE[i] <- round(sqrt(mean((width_diff)^2,na.rm=T)),1)
    fits.df$MAE[i] <- round(mean(abs(width_diff),na.rm=T),1)
    fits.df$Mean_bias[i] <- round(mean(width_diff,na.rm=T),1)
  }
return(fits.df)
} # End

  # All sites:
  fits_all <- compare.fits(sites)
  print(fits_all[order(fits_all$RMSE,decreasing=FALSE),])

  # What if St. Johns,FL sites are removed?
  sites_omitStJohn <- sites[-which(sites$sitecode %in% c("nwis_02236000","nwis_02234500","nwis_02234000")),]
  fits_omitStJohn <- compare.fits(sites_omitStJohn)
  print(fits_omitStJohn[order(fits_omitStJohn$RMSE,decreasing=FALSE),])
  
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


