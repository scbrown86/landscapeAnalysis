#
# kyle.taylor@pljv.org
# GIS Programmer/Analyst, Playa Lakes Joint Venture
#
# Here is my swiss army knife -- a bunch of functions that I routinely use in processing spatial data.  
#

#
# getPythonPath()
#

.getPythonPath <- function(){
  PYTHON <- unlist(lapply(as.list(unlist(strsplit(Sys.getenv("PATH"),split=":"))), FUN=list.files, pattern="python", full.names=T))
  if(length(PYTHON)>1){ 
    warning("multiple python binaries found in PATH.  Making a guess as to the default to use.")
    PYTHON <- PYTHON[grep(PYTHON, pattern="\\/python$")] # should correspond to /usr/bin/python on nix platforms.  This probably breaks win32 compat.
  } else if(length(PYTHON) == 0){
    stop("couldn't find a python executable in the user's PATH.")
  }
  return(PYTHON)
}

#
# getGDALtoolByName()
#
.getGDALtoolByName <- function(x=NULL){
  if(grepl(tolower(x),pattern="gdal_polygonize")){
    GDAL_POLYGONIZE <- unlist(lapply(as.list(unlist(strsplit(Sys.getenv("PATH"),split=":"))), FUN=list.files, pattern="gdal_polygonize.py", full.names=T))
      if(length(GDAL_POLYGONIZE)<1) stop("couldn't find gdal_polygonize.py tool in PATH")
        return(GDAL_POLYGONIZE)
  }
  return(NULL)
}

#
# rasterToPolygons()
# convert a raster to polygons using either 'R' or 'GDAL'.  GDAL is the (faster) default selection
#

rasterToPolygons <- function(r=NULL, method='gdal'){
 if(grepl(method,pattern='gdal')){ 
   r_name=deparse(substitute(r))
   writeRaster(r,paste(r_name,"tif",sep="."),overwrite=T);
   unlink(paste(r_name,c("shp","xml","shx","prj","dbf"),sep="."))
   system(paste(.getPythonPath(),.getGDALtoolByName("gdal_polygonize"),"-8",paste(r_name,"tif",sep="."),"-f \"ESRI Shapefile\"",paste(r_name,"shp",sep="."),sep=" "))
   return(rgdal::readOGR(".",r_name,verbose=F));
  }
}

#
# gaussianSmoothing()
# Implements a gaussian smoothing window as specified and implemented by Jeff Evans [2014] (see: http://evansmurphy.wix.com/evansspatial#!spatial-smoothing/ch1)
#

gaussianSmoothing <- function(x, s=1, d=5, filename=FALSE, ...) {
    if (!require(sp)) stop("sp PACKAGE MISSING")
    if (!require(raster)) stop("raster PACKAGE MISSING")
    if (!require(rgdal)) stop("rgdal PACKAGE MISSING")
    if (!inherits(x, "RasterLayer")) stop("MUST BE RasterLayer OBJECT") 
       GaussianKernel <- function(sigma=s, n=d) {
          m <- matrix(nc=n, nr=n)
            col <- rep(1:n, n)
            row <- rep(1:n, each=n)
          x <- col - ceiling(n/2)
          y <- row - ceiling(n/2)
         m[cbind(row, col)] <- 1/(2*pi*sigma^2) * exp(-(x^2+y^2)/(2*sigma^2))
        m / sum(m)
       }
   if (filename != FALSE) {
      focal(x, w=GaussianKernel(sigma=s, n=d), filename=filename, ...)
      print(paste("RASTER WRITTEN TO", filename, sep=": "))
        } else {
      return(focal(x, w=GaussianKernel(sigma=s, n=d), ...))
  }
}

#
# cropRasterByPolygons()
# Accepts a raster and SpatialPolygonDataFrame object, iterates over each polygon feature, creating rasters 
# for each step.  A raster list is returned to the user. Useful for parsing out climate/elevation data, county-by-county,
# for an entire state and then processing with the parallel package.
#

cropRasterByPolygons <- function(r=NULL, s=NULL, field=NULL, write=F){

  require(raster)
  require(rgdal)

  rS <- list()

  # sanity checks
  if(is.null(r) || is.null(s)){
    cat(" -- error: r= and s= are required.\n"); stop();
  } 

  # check our polygon data and projections
  if(class(s) == "SpatialPolygonsDataFrame"){
    if(!is.null(field)){
      if(sum(field %in% colnames(s@data)) == 0){ 
        cat(" -- error: field", field, "not found in shapefile.\n",sep=" "); stop();
      }
    }
  }
  s <- spTransform(s, CRS(projection(r))) # are we working with matching projections
  
  # pre-crop our raster to the extent of our shapefile
  r <- raster::crop(r,s)

  for(i in 1:nrow(s)) { 
    focal <- crop(r,s[i,]); 
      focal <- mask(focal,s[i,]);  
    if(write){ 
      writeRaster(focal,as.character(s[i,]@data[,field]),format="GTiff"); 
    } else {
      rS[[length(rS)+1]] <- focal
    }
    cat(".");
  } 
  cat("\n");

  # return the raster stack to the user, if asked
  if(!write) { return(rS) }
}

#
# spatialPointsToPPP()
# accepts a spatial points data frame, converts to PPP data that can be used by spatstat 
#

spatialPointsToPPP <- function(x,extentMultiplier=1.1){
  # default includes
  require(rgdal)
  require(raster)
  require(spatstat)
  
  e <- extent(x)
  
  if(!is.null(extentMultiplier)) { 
    e@xmin <- e@xmin*extentMultiplier
    e@xmax <- e@xmax+abs(e@xmax*(extentMultiplier-1))
    e@ymin <- e@ymin-abs(e@ymin*(extentMultiplier-1)) 
    e@ymax <- e@ymax*extentMultiplier
  }
  
  if(class(x) == "SpatialPointsDataFrame"){
    x <- x@coords
    x <- ppp(x=x[,1], y=x[,2], window=owin(xrange=c(e@xmin,e@xmax), yrange=c(e@ymin,e@ymax)))
  } 
  
  return(x)
}

csvToSpatialDataFrame <- function(path=NULL, proj4string="+init=epsg:4326"){
  require(sp)
  require(raster)
  # attempt to read-in csv data and convert to SpatialPointsDataFrame
  if(file.exists(path)){ t<-read.csv(path); 
  } else { stop(" -- failed to open input csv.\n") }

  coordinates(t) <- ~longitude+latitude
    projection(t) <- projection(proj4string)
      return(t)
}

#
# clusterReclassify()
# stub function to reclassify a raster using snow.  Note, even on machines that have -gt 3 cores,
# subs doesn't usually benefit from including them in the cluster.
#

clusterReclassify <- function(r,t=NULL, n=3){
  if(!(require(snow))) stop("snow package not installed.")
  if(!(require(raster))) stop("raster package not installed.")
  
  # sanity checks
  if(is.null(t)){
  	t <- try(read.csv("~/PLJV/CDL/NASS_ReMapTable_AD.txt",sep=":"))
  	  stopifnot(class(t) != "try-error")
  }
  
  endCluster(); beginCluster(n=3)
  
  r <- raster(r);
    # r <- clusterR(r, fun=subs, args=list(y=t, by=1, which=2, subWithNA=T)) # this doesn't work.  Consider re-writing using reclassify for cluster support
    cat(" -- warning: clustering support disabled.\n")
    r <- subs(r, y=t, by=1, which=2, subsWithNA=T)

  return(r)    
}

#
# findMinExtent()
# Find minimum extent from list of raster or extent objects.  Will return an extent object by default, 
# or the extent object as a SpatialPolygon if the ret='SpatialPolygons' argument is used.
#

findMinExtent <- function(x, ret=NULL){
	require(raster)
	# sanity-check
	if(!is.list(x)) { cat(" -- error: x= parameter should be a list.\n"); stop(); }
  # pre-process: if this is a raster list, let's solve for individual raster extents
  if(class(x[[1]]) == "RasterLayer" || class(x[[1]]) == "SpatialPolygons") {
    if(length(unique(unlist(lapply(x, FUN=projection)))) > 1) { 
    	  cat(" -- error: spatial list contains rasters with different projections.\n"); stop(); 
    }
  	x<-lapply(x,FUN=extent)
  }
  # basic (slow, but functional) step algorithm
  min <- x[[1]]
  for(i in 2:length(x)){
  	for(j in i:length(x)){
      if(x[[j]] < min) min <- x[[j]]
  	}
  } 
  if(is.null(ret)) { # by default, simply return the extent object
  	return(min)
  } else if(ret == "SpatialPolygons"){ # should we return minimum extent as a spatial polygons object?
      min <- as(min, 'SpatialPolygons')
      if(is.na(projection(min))){ projection(min) <- projection("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") }
      return(min)
  } else { stop(" -- unknown ret= parameter specification.") }
    
}

#
# Mode()
# Find the mode of a raster stack. Bogarted from Stack Exchange.
#

Mode <- function(x) {
  ux <- unique(x)
    ux <- ux[!is.na(ux)]
  return(ux[which.max(tabulate(match(x, ux)))])
}

#
# findMaxResolution()
# quick and dirty.  roll this up with above into a findMinExtent function.
#

findMaxResolution <- function(x) {
  require(raster)
  # sanity-check
  if(!is.list(x)) { cat(" -- error: x= parameter should be a list of spatial rasters.\n"); stop(); }
  if(length(unique(unlist(lapply(x, FUN=projection)))) > 1) { 
        cat(" -- error: spatial object list contains rasters with different projections.\n"); stop(); 
  }
  # pre-process: if this is a raster list, let's solve for individual raster extents
  if(class(x[[1]]) == "RasterLayer") {
    x<-lapply(x,FUN=res)
      x<-lapply(x,FUN=prod)
  }
  # basic (slow, but functional) step algorithm
  max <- x[[1]]
  for(i in 2:length(x)){
    for(j in i:length(x)){
      if(x[[j]] > max) max <- x[[j]]
    }
  } 
  return(sqrt(max))
}

#
# mergeRasters()
# quickly merge rasters passed as either list or vector of filenames. Useful for merging SSURGO RASTER data across counties
# implemented by (ex):
# files <- list.files(pattern="tif"); files <- files[!grepl(files,pattern=".tif.")]
# for(f in files){ focal <- list(raster(f), raster(paste("../../TX615/RASTER",f,sep="/"))); mergeRasters(x=focal, output="../../TX615/RASTER"); }
#

mergeRasters <- function(x, output=NULL){
  # convert a vector of filenames to a list of rasters, if the user didn't already create a list of rasters
  if(!is.list(x)) x <- lapply(as.list(x), FUN=raster)
  
  master <- x[[1]]
  if(length(x) > 1){
    cat(" -- processing: ")
    for(i in 2:length(x)){
      master <- try(raster::merge(master,x[[i]])); 
        if(class(master) == "try-error") { cat(" -- error:", master, "\n",sep=" "); stop(); }
      cat(".")
    }
    cat("\n")
  }      
  cat(" -- done\n")
  if(!is.null(output)){
    #writeRaster(master,paste(output,master@data@names,sep="/"),overwrite=T)
    writeRaster(master,filename=output,overwrite=T)
  } else {
    return(master)
  }
}

#
# clusterResample()
# use a snow cluster to reclassify components of a raster series to the minimum extent and maximum resolution of rasters within the series
#

clusterResample <- function(x, extent=NULL, resolution=NULL, n=4){
  require(raster)
  require(snow)

  # sanity checks
  if(!is.list(x)) x <- lapply(as.list(x), FUN=raster)
  if(length(unique(unlist(lapply(x, FUN=projection)))) > 1) { cat(" -- error: raster list contains rasters with different projections.\n"); stop(); }
  res <- lapply(x,FUN=res)
    res<- unlist(lapply(res,FUN=prod))
      ncell <- unlist(lapply(x, FUN=ncell))
        integrated <- unique(res*ncell)
  if(length(x) >1){
    if(length(integrated) == 1){
      cat(" -- warning: the resolution and number of cells in this raster series appear to be in agreement.  No need to run.\n"); 
      return(NULL)
    }
  }

  # build a snow cluster
  endCluster(); 
    beginCluster(n=n);
  # find our minimum extent and maximum resolution
  if(is.null(extent) || is.null(resolution)){
    extent <- findMinExtent(x)
       res <- findMaxResolution(x)  
  } 

  # build a target raster surface to house our reclassified raster data
  focal <- raster(extent, res=res)
    projection(focal) <- projection(x[[1]])

  x<-lapply(x, FUN=resample, y=focal, method="ngb")
    return(x)
}

#
# clusterProject()
# use a snow cluster to reclassify components of a raster series to the minimum extent and maximum resolution of rasters within the series
#

clusterProjectRaster <- function(x, crs=NULL, n=4){
  require(raster)
  require(rgdal)
  require(snow)

  # sanity checks
  if(!is.list(x)) x <- lapply(as.list(x), FUN=raster)
  if(length(unique(unlist(lapply(x, FUN=projection)))) == 1) { # do the projections of rasters in our list actually differ?
    cat(" -- warning: the projections in this raster series appear to be in agreement.  No need to run.\n"); 
    stop();
  } else if(is.null(crs)){ # did the user fail to provide a target CRS?  
    cat(" -- warning: target crs= not specified by user; using the projection of the first raster in the series as a target.\n")
    crs <- CRS(projection(x[[1]]))

  }
  if(length(x) >1){
    if(length(integrated) == 1){
      cat(" -- warning: the resolution and number of cells in this raster series appear to be in agreement.  No need to run.\n"); 
      return(NULL)
    }
  }

  # build a snow cluster
  endCluster(); 
    beginCluster(n=n);
  # find our minimum extent and maximum resolution
  if(is.null(extent) || is.null(resolution)){
    extent <- findMinExtent(x)
       res <- findMaxResolution(x)  
  } 

  # build a target raster surface to house our reclassified raster data
  focal <- raster(extent, res=res)
    projection(focal) <- projection(x[[1]])

  x<-lapply(x, FUN=resample, y=focal, method="ngb")
    return(x)
}

##
# MERGE RASTERS PROGRAM EXAMPLE
#
# require(raster)
# require(rgdal)
#
# argv <- commandArgs(trailingOnly=T);
# r<-c(raster(argv[1]), raster(paste(argv[2],argv[1],sep="/")));
# mergeRasters(r, output=argv[2]);
# quit();
#
## IN BASH
# cd TX135/RASTER; for r in `ls -1 *.tif`; do R --no-save --vanilla --slave --args $r /tmp/output.geomorphed.upsampled/TX606/RASTER/$r < /tmp/process.R; done; cd ../..;
# cd TX495/RASTER; for r in `ls -1 *.tif`; do R --no-save --vanilla --slave --args $r /tmp/output.geomorphed.upsampled/TX615/RASTER/$r < /tmp/process.R; done; cd ../..;
# cd TX461/RASTER; for r in `ls -1 *.tif`; do R --no-save --vanilla --slave --args $r /tmp/output.geomorphed.upsampled/TX618/RASTER/$r < /tmp/process.R; done; cd ../..;
#
# rm -rf TX135 TX495 TX461 
##

##
# FIND MIN EXTENT PROGRAM EXAMPLE
#
# source("~/Cloud/Code/ktaylor_essentialSpatialAddons.R")
# d <- list.dirs();
# d <- paste(d, "RASTER", sep="/")
# for(rs in d){ 
#   setwd("rs"); 
#   rasters<-list.files(pattern="tif$"); 
#   rasters<-lapply(as.list(rasters), FUN=raster); 
#   m<-findMinExtent(rasters, ret="SpatialPolygons"); 
#   setwd("../.."); 
# }
##

##
# RECLASSIFY A LIST OF RASTERS PROGRAM EXAMPLE
#
# argv <- commandArgs(trailingOnly=T)
# setwd(argv[1])
# source("~/Cloud/Code/ktaylor_essentialSpatialAddons.R")
# require(raster)
# cat(" -- constructing a raster list\n")
# files <- list.files(pattern="tif$")
#   files <- lapply(as.list(files), FUN=raster)
# cat(" -- resampling\n")
# o <- clusterResample(files,n=10)
# names <- paste(unlist(lapply(files, "names")), "tif", sep=".")
# cat(" -- writing rasters to disk\n")
# mapply(o,writeRaster,filename=names, overwrite=T)
# cat(" -- done\n")
##
