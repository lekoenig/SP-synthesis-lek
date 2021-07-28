README file

Last updated by LE Koenig 3 July 2019
Email: Lauren.Koenig@uconn.edu

This folder consists of one file titled "SP_drainagearea" that contains information about the contributing catchments upstream of respective sites in the StreamPULSE database (as of 3 July 2019).

SP_drainage area contains the following variables:

Region			StreamPULSE region or site-grouping
Site               	Abbreviated site code for each StreamPULSE sensor site
DataSource            	StreamPULSE site type within the project database (e.g. core, leveraged, embargoed)
SiteName              	Full name of StreamPULSE sensor site
latitude             	Geographic location of sensor site (EPSG 4269, CRS NAD83)
longitude           	Geographic location of sensor site (EPSG 4269, CRS NAD83)
USGS_STATIONID          NWIS station ID for sensor sites that are co-located with USGS streamflow gages. 
NWIS_DRAINAREA_SQKM     NWIS drainage area in square kilometers for reference
DRAINAREA_SQKM          Upstream drainage area in square kilometers
BOUNDARIES_SOURCE     	Data source for catchment boundaries
NHD_ORDER               Stream order of linked NHD flowline
NHD_NHDPlusID  		Unique NHDPlusID of linked NHD flowline (applicable for NHD_SCALE HR, NHD high-resolution beta product 1:24,000 scale)
NHD_ComID		Unique ComID of linked NHD flowline (applicable for NHD_SCALE MR, NHDV2 medium-resolution 1:100,000 scale)
NHD_SCALE		Scale of best-available NHD dataset for the region
NOTES			Additional notes regarding site naming conventions and/or site location


References:
Falcone, J. A. GAGES-II: Geospatial Attributes of Gages for Evaluating Streamflow. U.S. Geological Survey, http://water.usgs.gov/GIS/metadata/usgswrd/XML/gagesII_Sept2011.xml (2011).

Falcone, J. A., Baker, N. T. & Price, C. V. Watershed boundaries for study sites of the U.S. Geological Survey Surface Water Trends project. U.S. Geological Survey, http://dx.doi.org/10.5066/F78S4N29 (2017).

