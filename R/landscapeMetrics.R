require(raster)
require(rgdal)

# sampling addons and other goodies
# source("~/Cloud/Code/ktaylor_essentialSpatialAddons.R")


#
# spatialPolygonsToHexagonalGrid()
# Quick wrapper funtion for HexPoints2SpatialPolygons(), because I can never remember the name of
# the HexPoints2SpatialPolygons function when I need it.
#
# Author: Kyle Taylor (kyle.taylor@pljv.org)
#

spatialPolygonsToHexagonalGrid <- function(s=NULL,n=NULL){
  # default includes
  require(raster)
  require(rgdal)
  # sanity checks
  if(is.null(s) | is.null(n)) { cat(" -- error: no values specified for s= and n=\n"); stop(); }
  sLocations <- spsample(s, n=n, type="hexagonal")
    return(HexPoints2SpatialPolygons(sLocations))
}

#
# subsampleSurface()
# accepts a RasterLayer (and optional input pts=SpatialPoints object) and either generates a
# buffered point pattern process across its surface or buffers the pts= data and performs an extraction
# across the raster surface using those data.  Returns a list containing the extractions.
#
# Author: Kyle Taylor (kyle.taylor@pljv.org)
#

subsampleSurface <- function(x=NULL, pts=NULL, n=100, type='random', width=NULL, DEBUG=F){
  if(DEBUG){ t1 <- Sys.time(); }
  # default includes
  require(rgeos)
  # sanity checks
  if(is.null(x)) {
    cat(" -- error: x= parameter is undefined.\n");
    stop()
  } else if(is.null(width)){
    cat(" -- error: width= parameter is undefined. Buffer distance should be in meters.\n")
    stop()
  }
  # create a series of points covering the raster surface, x
  if(is.null(pts)){
    s <- as(extent(x), 'SpatialPolygons')
      s <- spsample(s, n=n, type=type)
        s <- split(s,f=rep(1:length(s))) # convert to a list
  } else if(class(pts) == "SpatialPoints") {
    # if the user specified spatial points at runtime, let's use them
    if(!is.na(crs(pts))){
      require(rgdal)
      pts <- spTransform(pts, CRS(raster::projection(x)))
    } else {
      cat(" -- error: NA coordinate reference system for pts= data.\n")
      stop()
    }
    s <- split(pts,f=rep(1:length(pts)))
  } else {
    cat(" -- error: unknown input. pts= argument should be SpatialPoints\n");
    stop();
  }
  # convert points to buffers and mask as needed across each buffered point
  s <- lapply(s, FUN=rgeos::gBuffer, width=width)
    r_s <- lapply(s, FUN=raster::crop, x=x)
        s <- mapply(FUN=raster::mask, x=r_s, mask=s);

  rm(r_s); gc(verbose=F)

  if(DEBUG){ cat(" -- debug (elapsed time): ", Sys.time()-t1, "\n"); }
  return(s)
}

#
# lReclass()
# accepts a list of rasters and reclasses each raster to a binary surface where 1=focal cover value and
# 0=everything else
#
# Author: Kyle Taylor (kyle.taylor@pljv.org)
#

lReclass <- function(x=NULL, inValues=NULL, nomatch=NA) lapply(lapply(x,FUN=raster::match,table=inValues,nomatch=nomatch), FUN=calc, fun=function(x,na.rm=F){x>=1})

#
# calcPatchIsolation()
# Calculate the KNN distance between patches in a binary raster (or raster(s); a vector or list) specified by r=RasterLayer.
# You can either return the KNN index or some statistic (e.g., mean) specified by fun=.  The default action is to return the full index.
# I have added multi-core support and use the gdal backend by default for this analysis. It's much faster than doing it in pure 'R'.
#
# Author: Kyle Taylor (kyle.taylor@pljv.org) [2015]
#

calcPatchIsolation <- function(r, fun=NA, k=1, method='gdal', parallel=FALSE){
  .include('raster');
  .include('sp');
  .include('rgeos');
  .include('FNN');
  .include('parallel');
  # set-up our cluster
  if(parallel) cl <- makeCluster(getOption("cl.cores", parallel::detectCores()-1))
  # convert our raster to polygons
  if(parallel){
    r <- parLapply(cl=cl,as.list(r),fun=landscapeAnalysis::rasterToPolygons,method=method)
  } else {
    r <- lapply(as.list(r),FUN=landscapeAnalysis::rasterToPolygons,method=method)
  }
  # check for NA values
  na_values <- as.vector(unlist(lapply(X=as.list(r),FUN=is.na)))
  if(sum(na_values)>0){
    r[na_values] <- r[which(!na_values)[1]] # overwrite our NA values with something valid
      r <- lapply(X=as.list(r),FUN=getSpPPolygonsLabptSlots)
  } else {
    r <- lapply(X=as.list(r),FUN=getSpPPolygonsLabptSlots)
  }
  # do our NN assessment
  d <- function(x,na.rm=F){ o<-try(FNN::knn.dist(x,k=k)); if(class(o) != "try-error") { x <- o; } else { x <- NA }; return(x)}
  r <- lapply(as.list(r),FUN=d);rm(d);
    r <- lapply(as.list(r), FUN=ifelse(is.na(fun),mean,match.fun(fun)))
      r[na_values] <- NA # restore our NA values
  # clean-up
  if(parallel) stopCluster(cl);
  # return the issolation metric for each raster as a vector
  return(as.vector(unlist(r)))
}

#
# calculateLandscapeMetric()
# quick wrapper function that accepts a raster list and applies landscape metrics across the list.  Code for landscape
# metrics is implemented using an approach coded by Jeremy VanDerWal (jjvanderwal@gmail.com), using the SDMTools package.
#
# Author: Kyle Taylor (kyle.taylor@pljv.org)
#

lCalculateLandscapeMetrics <- function(x=NULL, metric=NULL, DEBUG=F){
  if(DEBUG){ t1 <- Sys.time(); }
  # default includes
  require(SDMTools)
  # sanity checks
  if(class(x) != "list") { cat(" -- error: x= object should be a list() of rasters, as returned by subsampleSurface().\n"); stop(); }

  # sanity check our buffer raster surfaces for valid data before passing to PatchStats and ClassStats
  naCheck <- function(x){ sum(values(is.na(x))) == ncell(x) }
    allNa <- unlist(lapply(x, naCheck))

  if(sum(allNa)/length(allNa) == 1){
      if(DEBUG) warning("100% of buffered samples are NA.  Returning zero for all patch stats.")
      vClass <- vPatch <- rep(0, length(x))
      return(list(vClass, vPatch))
  } else if(sum(allNa)>0){
      if(DEBUG) warning((sum(allNa)/length(allNa))*100, "% of the sample buffers are completely NA.  Removing offenders, but it will lower your sample size.")
      x <- x[!allNa]
  }

  # calculate patch and class statistics
  t <- lapply(x, ClassStat)
    t <- list(t,lapply(x, PatchStat)) # make a list of lists containing class stats and patch stats

  # parse out the requested metrics as vectors and return to user, if requested
  if(!is.null(metric)) return(lapply(as.list(metric), FUN=metricsListToVector, x=t))
  # return the full table of stats offered, by default
  if(DEBUG){ cat(" -- debug (elapsed time): ", Sys.time()-t1, "\n"); }
   return(t)
}

#
# lExtract()
# takes spatial features (e.g., SpatialPoints) from y= and extracts across a list of raster scenes specified by X=
# handy for extracting data across the spatial extent of scenes without having to mosaic the rasters into a single surface.
#

lExtract <- function(X=NULL, y=NULL, fun=mean){
  output <- NULL    # we don't know how many points will actually overlap our surface -- will assign output values dynamically
  y$id <- 1:nrow(y) # order our initial point sample by id, so that we can sort effectively later
  cat(" -- extracting across raster scenes: ")
  for(x in X){
    if(nrow(y)>0){ # sanity-check : prevent an attempted extraction with no points
      y <- spTransform(y, CRS(projection(x)))
      # extract overlapping features and add to our merged output
      if(gIntersects(as(extent(x), 'SpatialPolygons'), as(extent(y), 'SpatialPolygons'))){
         vals <- raster::extract(x, y, sp=F, fun=fun)
        # keep non-NA values and add them to output spatial feature -- place NA values back in y for the next iteration of the raster list
        focal <- y[!is.na(vals),]
          focal@data <- cbind(focal@data,data.frame(val=vals[!is.na(vals)]))
        y <- y[is.na(vals),]
        if(is.null(output)){
          output <- focal
        } else {
          output <- rbind(output, spTransform(focal, CRS(projection(output))))
        }
      }
    }; cat(".")
  }; cat("\n")
  # merge in the remaining points (with NA extracted values)
  if(nrow(y)>1){
    y$val <- rep(NA,nrow(y))
    output <- rbind(output,y)
  }
  # return out output points, sorted by initial $id
  return(output[order(output$id),])
}

#
# metricsListToVector()
# parses the list output from lCalculateLandscapeMetrics() by metric name, returning the results as a vector.
#
metricsListToVector <- function(x,metric="total.area") { 
  ret <- as.vector(unlist(x)[grepl(names(unlist(x)),pattern=metric)])
    ret[is.na(ret)] <- 0;
      return(ret);
}
