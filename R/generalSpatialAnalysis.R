#
# kyle.taylor@pljv.org
# GIS Programmer/Analyst, Playa Lakes Joint Venture
#
# Here is my swiss army knife -- a bunch of functions that I routinely use in processing spatial data.  When possible, I
# try to use calls to GDAL or GRASS to handle raster operations.  A full GIS will usually do these things more quickly than
# can be done with pure R (for now).
#

#' built-in (hidden) wrapper function for require that will rudely attempt to install missing packages.
#'
#' @param x specifies the package name
#' @param from specifies whether we are fetching from CRAN or GitHub
#' @param repo specifies the repository URL to use.  Will pick a few good ones by default.
include <- function(x,from="cran",repo=NULL){
  if(from == "cran"){
    if(!do.call(require,as.list(x))) install.packages(x, repos=c("http://cran.revolutionanalytics.com","http://cran.us.r-project.org"));
    if(!do.call(require,as.list(x))) stop("auto installation of package ",x," failed.\n")
  } else if(from == "github"){
    if(!do.call(require,as.list(x))){
      if(!do.call(require,as.list('devtools'))) install.packages('devtools', repos=c("http://cran.revolutionanalytics.com","http://cran.us.r-project.org"));
      require('devtools');
      install_github(paste(repo,x,sep="/"));
    }
  } else{
    stop(paste("could find package:",x))
  }
}

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
  x <- tolower(x)
    x <- unlist(lapply(as.list(unlist(strsplit(Sys.getenv("PATH"),split=":"))), 
          FUN=list.files, pattern=x, full.names=T))
          
  if(length(x)<1){
    warning(paste("couldn't find",x,"tool in PATH",sep=""))
    return(NULL)
  } else if(length(x)>1){
    warning(paste("multiple references to ",x," tool in PATH -- honoring first occurrence.",sep=""))
    return(x[1])
  } 
  
  return(x)
}

#
# parseLayerDsn()
# 
.parseLayerDsn <- function(x=NULL){
  path <- unlist(strsplit(x, split="/"))
    layer <- gsub(path[length(path)],pattern=".shp",replacement="")
      dsn <- paste(path[1:(length(path)-1)],collapse="/")
  return(c(layer,dsn))
}

#
# readOGRfromPath()
# 
.readOGRfromPath <- function(path=NULL){
  landscapeAnalysis:::include('rgdal')
  path <- .parseLayerDsn(path)
   
  layer <- path[1]
    dsn <- path[2]

  return(readOGR(dsn,layer,verbose=F))
}

#
# rasterToPolygons()
# convert a raster to polygons using either 'R' or 'GDAL'.  GDAL is the (faster) default selection
#
# Author: Kyle Taylor
#

rasterToPolygons <- function(r=NULL, method='gdal'){

 cleanUp <- function(n,rmAll=F){
   if(rmAll) {
     unlink(paste(n,c("shp","xml","shx","prj","dbf","tif"),sep="."),force=T,recursive=T)
   } else {
     unlink(paste(n,c("shp","xml","shx","prj","dbf"), sep="."),force=T,recursive=T)
   }
 }

 if(grepl(method,pattern='gdal')){
   r_name=deparse(substitute(r))
     r_name=paste(r_name,sprintf("%.0f", round(as.numeric(runif(n=1,min=0,max=9999999)))),sep="_")
   raster::writeRaster(r,paste(r_name,"tif",sep="."),overwrite=T);
   cleanUp(r_name)
   if(try(system(paste(.getPythonPath(),.getGDALtoolByName("gdal_polygonize"),"-8",paste(r_name,"tif",sep="."),"-f \"ESRI Shapefile\"",paste(r_name,"shp",sep="."),sep=" ")))==0){
     if(class(try(s<-rgdal::readOGR(".",r_name,verbose=F))) != "try-error"){
       cleanUp(r_name,rmAll=T)
       return(s);
     } else {
       warning("gdal_polygonize error : it's possible we tried to polygonize a raster with only NA values")
       return(NA)
     }

   } else {
     warning("gdal_polygonize error")
     return(NA)
   }
  } else {
    return(raster::rasterToPolygons(r, dissolve=T,progress='text'))
  }
}

#
# extractDensities()
# extract quantiles from a continuous raster surface as spatial polygons using GDAL
#
# Author: Kyle Taylor (kyle.taylor@pljv.org)
#

extractDensities <- function(x,s=5,d=15, p=c(0.5,0.9)){
  p <- sprintf("%.2f", round(as.numeric(p),2)) # force a trailing N
  # extract  range contours for raster surface x
  q <- sprintf("%.2f",as.numeric(seq(0,1,0.05)))
    if(sum(as.character(p) %in% as.character(q)) != length(p)) stop("quantiles are typically extracted in 0.05 interval steps")
      q <- as.numeric(q)[!as.numeric(q) %in% c(0,1)]
  # smooth
  smoothed <- gaussianSmoothing(x,s=s)
    h <- hist(smoothed, plot=F)
      smoothed[smoothed<=h$mids[2]] <- NA
        quantiles <- as.vector(quantile(smoothed,probs=q))
  out <- list()
  for(focal in p){
    smoothed_focal <- smoothed>=quantiles[which(as.numeric(q) == as.numeric(focal))]
      smoothed_focal <- raster::match(smoothed_focal,1,nomatch=NA)
        out[[length(out)+1]] <- landscapeAnalysis::rasterToPolygons(smoothed_focal);
  }
  return(out)
}

#
# gaussianSmoothing()
# Implements a gaussian smoothing window as specified and implemented by Jeff Evans [2014] (see: http://evansmurphy.wix.com/evansspatial#!spatial-smoothing/ch1)
# Author: Kyle Taylor (kyle.taylor@pljv.org) [originally J. Evans]
#

gaussianSmoothing <- function(x, s=1, d=5, filename=FALSE, ...) {
  landscapeAnalysis::include('sp')
  landscapeAnalysis::include('raster')
  landscapeAnalysis::include('rgdal')
  if (!inherits(x, "RasterLayer")) stop("x= argument expects a raster* object")
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
# splitExtent()
#
# When working with SDB queries for large areas, we consistently lose a lot of data... particularly for large
# counties.  This gets around that by splitting an extent object into adjacent quarters so that we can download and
# merge our raster segments later.  For small counties, this is inefficient.  But its better than just flatly attempting downloads
#
# Author: Kyle Taylor (kyle.taylor@pljv.org) [2016]
#
splitExtent <- function(e=NULL,multiple=2){
  landscapeAnalysis::include('raster')
  # define our x/y vector ranges
  x <- rep(NA,multiple+1)
  y <- rep(NA,multiple+1)
  # define the x/y range for calculating the size of our extents
  xStep <- diff(c(e@xmin,e@xmax))/multiple
  yStep <- diff(c(e@ymin,e@ymax))/multiple
  # assign vertices to our product vectors
  for(i in 1:(multiple+1)){
    x[i] <- ifelse(i==1,
                   min(e@xmin),
                   x[i-1]+xStep)
    y[i] <- ifelse(i==1,
                   min(e@ymin),
                   y[i-1]+yStep)
  }
  # assign our vertices to extent objects
  extents <- as.list(rep(NA,multiple*multiple))
  # iterate over our extents, assigning as we go
  yStart <- i <- 1;
  while(i <= length(extents)){
    for(j in 1:multiple){ # stagger our y-values
      extents[i] <- extent(c(x[j],x[j+1],y[yStart],y[yStart+1]))
      i <- i+1;
    }
    yStart <- yStart+1;
  }
  return(extents)
}

#
# cropRasterByPolygons()
# Accepts a raster and SpatialPolygonDataFrame object, iterates over each polygon feature, creating rasters
# for each step.  A raster list is returned to the user. Useful for parsing out climate/elevation data, county-by-county,
# for an entire state and then processing with the parallel package.
#

cropRasterByPolygons <- function(r=NULL, s=NULL, field=NULL, write=F, parallel=F){

  landscapeAnalysis::include('raster')
  landscapeAnalysis::include('rgdal')

  rS <- list()

  # sanity checks
  if(is.null(r) || is.null(s)){
    cat(" -- error: r= and s= are required.\n"); stop();
  }

  # check our polygon data and projections
  if(grepl(class(s),pattern="SpatialPolygonsDataFrame")){
    if(!is.null(field)){
      if(sum(field %in% colnames(s@data)) == 0){
        cat(" -- error: field", field, "not found in shapefile.\n",sep=" "); stop();
      }
    }
  }
  s <- sp::spTransform(s, CRS(projection(r))) # are we working with matching projections?
    r <- raster::crop(r,s)
      s <- sp::split(s, f=1:nrow(s)) # split to list by row for apply operations
  # don't try to parallelize with a large number of polygons unless you are on a system with a whole lot of RAM
  if(parallel){
    landscapeAnalysis::include('parallel');
    cl <- makeCluster(getOption("cl.cores", parallel::detectCores()-1),outfile='outfile.log');
     r <- parLapply(cl=cl,Y=s,fun=raster::crop,X=rep(list(r),length(s)))
      rS <- parLapply(cl=cl,Y=s,fun=raster::mask,x=focal)
  } else {
    for(i in 1:length(s)) {
      focal <- crop(r,s[[i]]);
        focal <- mask(focal,s[[i]]);
      if(write){
        writeRaster(focal,as.character(s[[i]]@data[,field]),format="GTiff");
      } else {
        rS[[length(rS)+1]] <- focal
      }
      cat(".");
    }
    cat("\n");
  }

  # return the raster stack to the user, if asked
  if(!write) { return(rS) }
}
#
# multiplyExtent()
#
multiplyExtent <- function(x,extentMultiplier=1.1){
  e <- extent(x)

  if(!is.null(extentMultiplier)) {
    e@xmin <- e@xmin*extentMultiplier
    e@xmax <- e@xmax+abs(e@xmax*(extentMultiplier-1))
    e@ymin <- e@ymin-abs(e@ymin*(extentMultiplier-1))
    e@ymax <- e@ymax*extentMultiplier
  }

  return(e)
}
#
# as.owin()
#
as.owin <- function(x,extentMultiplier=1.1){
    e <- multiplyExtent(x,extentMultiplier=extentMultiplier)
    return(owin(xrange=c(e@xmin,e@xmax), yrange=c(e@ymin,e@ymax)))
}
#
# spatialPointsToPPP()
# accepts a spatial points data frame, converts to PPP data that can be used by spatstat
#
spatialPointsToPPP <- function(x,extentMultiplier=1.1,field=NULL){
  # attribute 'data' to 'marks' for our PPP
  if(grepl(class(x),pattern="SpatialPoints")){
    d <- x@data
    c <- x@coords
    if(grepl(class(x),pattern="SpatialPointsDataFrame")){
      if(!is.null(field)){
        x <- ppp(x=c[,1], y=c[,2], window=as.owin(x,extentMultiplier=extentMultiplier), marks=d[,field])
      } else {
        x <- ppp(x=c[,1], y=c[,2], window=as.owin(x,extentMultiplier=extentMultiplier), marks=d)
      }
    } else {
      x <- ppp(x=c[,1], y=c[,2], window=as.owin(x,extentMultiplier=extentMultiplier))
    }
  }
  return(x)
}

csvToSpatialDataFrame <- function(path=NULL, proj4string="+init=epsg:4326"){
  landscapeAnalysis::include('sp')
  landscapeAnalysis::include('raster')
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
# subs doesn't usually benefit from including them in the cluster.  Might try reimplementing in GDAL
# or splitting a raster into pieces and reclassifying using 'parallel'.
#

clusterReclassify <- function(r,t=NULL, n=3){
  landscapeAnalysis::include('snow')
  landscapeAnalysis::include('raster')

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
# snapTo()
# Ensure absolute consistency between raster objects by cropping,projecting,snapping,and (if asked) resampling
# a raster object using a template
#

snapTo <- function(x,to=NULL,names=NULL,method='bilinear'){
  require(parallel)
  # set-up a cluster for parallelization
  cl <- makeCluster((parallel::detectCores()-2))
  # crop, reproject, and snap our raster to a resolution and projection consistent with the rest our explanatory data
  if(grepl(tolower(class(x)),pattern="character")){ lapply(x,FUN=raster) }
  e <- as(extent(to[[1]]),'SpatialPolygons')
    projection(e) <- CRS(projection(to[[1]]))
  if(class(x) == "list") {
    x <- parLapply(cl,x,fun=raster::crop,extent(spTransform(e,CRS(projection(x[[1]])))))
      x <- parLapply(cl,x,fun=raster::projectRaster,crs=CRS(projection(to[[1]])))
    extents <- lapply(x,alignExtent,to[[1]])
      for(i in 1:length(x)){ extent(x[[i]]) <- extents[[i]] }
    if(!is.null(method)){
      x <- parLapply(cl,x,fun=resample,y=to[[1]],method=method)
    }
  } else {
    x <- raster::crop(x,extent(spTransform(e,CRS(projection(x)))))
      x <- raster::projectRaster(x,crs=CRS(projection(to[[1]])))
    extent <- alignExtent(x,to[[1]])
      extent(x) <- extent
    if(!is.null(method)){
      x <- raster::resample(x,y=to[[1]],method=method)
    }
  }
  endCluster()
  return(x)
}

#
# findMinExtent()
# Find minimum extent from list of raster or extent objects.  Will return an extent object by default,
# or the extent object as a SpatialPolygon if the ret='SpatialPolygons' argument is used.
#

findMinExtent <- function(x, ret=NULL){
  landscapeAnalysis::include('raster')
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
      if(is.na(projection(min))){ projection(min) <- projection(x[[1]]) }
      return(min)
  } else { stop(" -- unknown ret= parameter specification.") }

}

#
# Mode()
# Find the mode of a raster stack. Use the raster::calc() function to apply across a raster stack.
# Bogarded from Stack Exchange.
#

 Mode <- function(x,na.rm=T){
    ux <- unique(x)
    ux=ux[!is.na(ux)]
    ux[which.max(tabulate(match(x, ux)))]
 }

#
# spatialLinesGridToSpatialPolygons()
# accepts an x=SpatialLines* object and converts to polygons.  Useful for generating samples within national grid units.
#

spatialLinesGridToSpatialPolygons <- function(x, res=1000,method="raster"){
  landscapeAnalysis::include('maptools')
  landscapeAnalysis::include('raster')
  if(!inherits(x,'SpatialLines')) stop("x= argument is not of type SpatialLines*")
  # convert to an arbitrary CRS with metric units and decent consistency across North America for maptools::SpatialLinesMidPoints()
  originalCRS <- raster::CRS(raster::projection(x))
  x <- spTransform(x,CRS(projection("+init=epsg:2163")))
    x <- maptools::SpatialLinesMidPoints(x)
  if(grepl(method,pattern="raster")){
    x <- rasterize(x,raster(ext=extent(x),res=res,crs=CRS(projection(x))))
      x[!is.na(x)] <- 1:sum(!is.na(x)) # give unique values to all valid cells to inform rasterToPolygons
        x <- rasterToPolygons(x,na.rm=T,dissolve=FALSE)
  } else if(grepl(method,pattern="rgeos")){
    x<-rgeos::gBuffer(x,width=res/2,capStyle="SQUARE")
  } else {
    stop("unknown method= argument passed to spatialLinesToSpatialPolygons")
  }
  return(x)
}

#
# findMaxResolution()
# quick and dirty.  roll this up with above into a findMinExtent function.
#

findMaxResolution <- function(x) {
  landscapeAnalysis::include('raster')
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
# lMerge()
# quickly merge rasters passed as either list or vector of filenames. Useful for merging SSURGO RASTER data across counties
# implemented by (ex):
# files <- list.files(pattern="tif"); files <- files[!grepl(files,pattern=".tif.")]
# for(f in files){ focal <- list(raster(f), raster(paste("../../TX615/RASTER",f,sep="/"))); mergeRasters(x=focal, output="../../TX615/RASTER"); }
#

lMerge <- function(x, output=NULL, method="R"){
  # sanity checks
  if(length(x)<2){
    stop("x= argument should be a list or vector with a length greater than 1 for a merge operation.")
  }
  # Default "R" method for merging 
  if(grepl(tolower(method),pattern="r")){
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
  # snappier GDAL method of merging, for environments that support it
  } else if(grepl(tolower(method),pattern="gdal")){
    # sanity-checks
    if(class(x) == "list"){
      # assume x= list is a list of rasters.  Get the filenames from that list
      stop("list operation not yet implemented. Gdal method currently requires x= vector()")
    }
    # attempt a merge
    if(is.null(output)) output <- "gdal_merged.tif"
    if(try(system(paste(.getPythonPath(),.getGDALtoolByName("gdal_merge"),"-tap -o",output,"-of GTiff",x,sep=" ")))==0){
      # success
      return(raster(output))
    } else {
      stop("gdal_merge.py failed (non-zero status)")
    }
  }
}

#
# clusterResample()
# use a snow cluster to reclassify components of a raster series to the minimum extent and maximum resolution of rasters within the series.
# Useful for making a list of wildly different raster objects stack()able.
#

clusterResample <- function(x, extent=NULL, resolution=NULL, n=4){
  landscapeAnalysis::include('raster')
  landscapeAnalysis::include('snow')

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
  # snap rasters in our list for consistency
  e <- lapply(x,FUN=alignExtent,x[[1]])
    for(i in 1:length(x)){ extent(x[[i]]) <- e[[i]] }
  # build a snow cluster
  endCluster();
    beginCluster(n=n);
  # find our minimum extent and maximum resolution
  if(is.null(extent) || is.null(resolution)){
    extent <- findMinExtent(x)
       res <- findMaxResolution(x)
  } else {
    res <- resolution
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
  landscapeAnalysis::include('raster')
  landscapeAnalysis::include('rgdal')
  landscapeAnalysis::include('snow')
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

  endCluster();
    beginCluster(n=n);

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

# polygonPathDistance()
#
# calculates path distances in a verroni tesselation.  Function accepts a polygon object
# with an id= argument specify the focal (center) polygon from which path distances are calculated.
# returns attributed polygons with a $class field indicating the path distance as 1,2,3, etc...
# Incorporates a width= argument for buffering to account for slivers and islands.
#
polygonPathDistance <- function(x=null,id=0, width=NULL, quietly=F){
  landscapeAnalysis::include("rgdal")
  landscapeAnalysis::include("rgeos")
  # sanity checks
  if(!inherits(x,"SpatialPolygons")){
    stop("x= argument should specify a SpatialPolygons* object")
  } else if(is.null(x$id)){
    x$id <- 1:length(x)
  }
  if(is.null(width)){
    if(!grepl(projection(x),pattern="+units=m")){
      warning("units for x= argument are not in meters. Assuming 50-meter degrees-equivalent for buffering. Re-run with width=0 to disable buffering.")
      width <- 0.00044915599
    } else {
      warning("null width= argument. Assuming a 50 meter buffer distance for the step algorithm.  Re-run with width=0 to disable buffering.")
      width <- 50
    }
  }
  # iterate over our polygons until every polygon has been attributed with a class
  class <- 0
    x$class <- -1
      x@data[x$id == id,'class'] <- class

  if(!quietly) cat(" -- stepping through polygon adjacencies: ")
  while(sum(x$class == -1)>0){
    # unclassified polygons
    focal <- x[x$class == -1,];
    # parse-out the classified polygons from previous iteration
    if(width>0){
      lastClass <- gBuffer(x[x$class == class,],width=width)
    } else {
      lastClass <- x[x$class == class,]
    }
    # identify overlap between unclassified polygons and last classified polygons (buffered by some distance)
    suppressMessages(rows <- try(which(as.vector(colSums(gOverlaps(focal,lastClass,byid=T)) > 0))))
    # if gOverlaps() was successful, keep on chugging -- otherwise, return what we have to the user for debugging
    if(class(rows)!="try-error"){
      if(length(rows)==0){
        # break on failure to identify an adjacency for next iteration
        warning(paste("no overlapping polygons found at step ",class," -- quiting."))
        return(x);
      }
      focal@data[rows,'class'] <- rep(class+1, length(focal@data[rows,'class']));
      x@data[x$class == -1,] <- focal@data;
      class <- class + 1;
    } else {
      # break on failure to identify an adjacency for next iteration
      return(x);
    }
    if(!quietly) cat(paste("[",class,"]",sep=""));
  };
  if(!quietly) cat("\n");
  return(x)
}
