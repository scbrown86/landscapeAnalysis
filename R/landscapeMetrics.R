require(raster)
require(rgdal)

#' A quick wrapper funtion for HexPoints2SpatialPolygons(), because I can never remember the name of
#' the HexPoints2SpatialPolygons function when I need it.
#' @export
spatialPolygonsToGrid <- function(s=NULL,n=NULL,area=NULL,type="hexagonal"){
  # default includes
  landscapeAnalysis:::include('raster')
  landscapeAnalysis:::include('rgdal')
  landscapeAnalysis:::include('rgeos')
  # sanity checks
  if(is.null(s) && is.null(n)) { cat(" -- error: no values specified for s= and n=\n"); stop(); }
  if(!is.null(area)){ # did the user define an area (in meters) that we can use to define the number of cells in our mesh grid?
    n <- ceiling(rgeos::gArea(spTransform(s,CRS(projection("+init=epsg:2163"))))/area)
  }
  sLocations <- spsample(s, n=n, type=type)
    s <- HexPoints2SpatialPolygons(sLocations)
      s <- SpatialPolygonsDataFrame(s,data=data.frame(id=1:length(s),row.names=paste(sep="","ID",1:length(s))))
  return(s)
}

#' Accepts a RasterLayer (and optional input pts=SpatialPoints object) and either generates a
#' buffered point pattern process across its surface or buffers the pts= data and performs an extraction
#' across the raster surface using those data.  Returns a list containing the extractions.
#' @export
subsampleSurface <- function(x=NULL, pts=NULL, n=100, type='random', width=NULL, DEBUG=F){
  if(DEBUG){ t1 <- Sys.time(); }
  # default includes
  landscapeAnalysis:::include('rgeos')
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
  } else if(class(pts) == "SpatialPoints") {
    # if the user specified spatial points at runtime, let's use them
    if(!is.na(crs(pts))){
      landscapeAnalysis:::include('rgdal')
      s <- spTransform(pts, CRS(raster::projection(x)))
    } else {
      cat(" -- error: NA coordinate reference system for pts= data.\n")
      stop()
    }
  } else {
    cat(" -- error: unknown input. pts= argument should be SpatialPoints\n");
    stop();
  }
  # convert to a list
  out <- list();
    for(i in 1:length(s)){ out[[length(out)+1]] <- s[i,] }
      s <- out; rm(out);
  # convert points to buffers and mask as needed across each buffered point
  s <- lapply(s, FUN=rgeos::gBuffer, width=width)
    r_s <- lapply(s, FUN=raster::crop, x=x)
        s <- mapply(FUN=raster::mask, x=r_s, mask=s);

  rm(r_s); gc(verbose=F)

  if(DEBUG){ cat(" -- debug (elapsed time): ", Sys.time()-t1, "\n"); }
  return(s)
}

#' accepts a list of rasters and reclasses each raster to a binary surface where 1=focal cover value and
#' 0=everything else
#' @export
lReclass <- function(x=NULL, inValues=NULL, nomatch=NA) lapply(lapply(x,FUN=raster::match,table=inValues,nomatch=nomatch), FUN=calc, fun=function(x,na.rm=F){x>=1})

#' Calculate the KNN distance between patches in a binary raster (or raster(s); a vector or list) specified by r=RasterLayer.
#' You can either return the KNN index or some statistic (e.g., mean) specified by fun=.  The default action is to return the full index.
#' I have added multi-core support and use the gdal backend by default for this analysis. It's much faster than doing it in pure 'R'.
#' @param fun a function that is applied to all KNN distances (i.e., for calculating mean/median/mode distances).
#' @param k the K (e.g., 1,2,3) value used for our K nearest-neighbor calculation.
#' @param method string specifying whether 'gdal' or 'R' is used to vectorize our raster patches for the KNN calculation.
#' @param from string specifying whether to calculate patch issolation distances from the 'edges' (i.e., minimum distance) or patch 'centroid'.
#'        default is 'edges'.
#' @param parallel TRUE/FALSE specifying whether to parallelize our vectorization / KNN operations.
#' @export
calcPatchIsolation <- function(r, fun=NULL, k=1, method='gdal',from='edges',parallel=FALSE){
  landscapeAnalysis:::include('raster');
  landscapeAnalysis:::include('sp');
  landscapeAnalysis:::include('rgeos');
  landscapeAnalysis:::include('FNN');
  landscapeAnalysis:::include('parallel');
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
      r <- lapply(X=as.list(r),FUN=sp::coordinates)
  } else {
    r <- lapply(X=as.list(r),FUN=sp::coordinates)
  }
  # do our NN assessment, pre-checking for empty (or single) points in matrices
  d <- function(x,na.rm=F){
    if(nrow(x)>k){ # are there at-least k observations to inform our KNN?
      o <- suppressMessages(try(FNN::knn.dist(x,k=k)));
      if(class(o) != "try-error") {
        x <- o;
      } else {
        x <- NA
      }
    } else {
      x <- NA
    }
    return(x)
  }
  r <- lapply(as.list(r),FUN=d);rm(d);
    if(is.null(fun)){
      r <- lapply(r,FUN=mean)
    } else {
      r <- lapply(r,FUN=match.fun(fun))
    }
  r[na_values] <- NA # restore our NA values
  # clean-up
  if(parallel) stopCluster(cl);
  # return the issolation metric for each raster as a vector
  return(as.vector(unlist(r)))
}

#' a wrapper function that accepts a raster list and applies landscape metrics across the list.  Code for landscape
#' metrics is implemented using an approach coded by Jeremy VanDerWal (jjvanderwal@gmail.com), using the 'SDMTools' package.
#'
#' @param x list object containing raster(s) that will be passed to 'SDMTools' to calculate class/patch statistics.
#' @param metric a vector containing the target class/patch metrics calculated for each raster object in 'x'.
#' @param class raster cell values (usually binary) indicating the focal class.  Default value is 1.
#' @export
lCalculateLandscapeMetrics <- function(x=NULL, metric=NULL, DEBUG=F, class=1){
  if(DEBUG){ t1 <- Sys.time(); }
  # default includes
  landscapeAnalysis:::include('SDMTools')
  # sanity checks
  if(class(x) != "list") { cat(" -- error: x= object should be a list() of rasters, as returned by subsampleSurface().\n"); stop(); }

  # sanity check our buffer raster surfaces for valid data before passing to PatchStats and ClassStats
  # we will not attempt calculations on a raster surface containing only NA values -- this causes
  # problems for 'SDMTools'
  naCheck <- function(x){ sum(values(is.na(x))) == ncell(x) }
    allNa <- unlist(lapply(x, naCheck))

  if(sum(allNa)/length(allNa) == 1){
      if(DEBUG) warning("100% of buffered samples are NA.  Returning zero for all patch stats.")
      vClass <- vPatch <- rep(0, length(x))
        t <- data.frame(do.call(cbind,list(vClass,vPatch)))
          names(t) <- metric
           return(cbind(id=1:length(x),t))
  } else if(sum(allNa)>0){
      if(DEBUG) warning((sum(allNa)/length(allNa))*100, "% of the sample buffers are completely NA.  Removing offenders, but it will lower your sample size.")
      x <- x[!allNa]
  }

  # calculate patch and class statistics
  t <- lapply(x, ClassStat)
    t <- list(t,lapply(x, PatchStat)) # make a list of lists containing class stats and patch stats

  # parse out the requested metrics as vectors and return to user, if requested
  if(!is.null(metric)) {
    t <- lapply(as.list(metric), FUN=landscapeAnalysis:::metricsListToVector, x=t, class=class)
      t <- data.frame(do.call(cbind,t))
        names(t) <- metric
          # back-fill NA rasters with zero values for our focal metric(s)
          t <- cbind(id=which(as.logical(!allNa)),t)
          t_na <- data.frame(matrix(ncol=(length(metric)+1),rep(0,sum(allNa)*(length(metric)+1))))
            names(t_na) <- names(t)
              t_na$id <- which(as.logical(allNa))
                t <- rbind(t,t_na)
                  t <- t[order(t$id),]
  }
  # return the full table of stats offered, by default
  if(DEBUG){ cat(" -- debug (elapsed time): ", Sys.time()-t1, "\n"); }
   return(t)
}

#' takes spatial features (e.g., SpatialPoints) from y= and extracts across a list of raster scenes specified by X=
#' handy for extracting data across the spatial extent of rasters without having to mosaic and create a 'stack' object.
#' @param X list of rasters
#' @param y Spatial* object to extract across the 'X' raster series.
#' @export
lExtract <- function(X=NULL, y=NULL, fun=mean){
  output <- NULL    # we don't know how many points will actually overlap our surface -- will assign output values dynamically
  y$id <- 1:nrow(y) # order our initial point sample by id, so that we can sort effectively later
  cat(" -- extracting across raster(s): ")
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
#' Internal (hidden) function that parses the list output from lCalculateLandscapeMetrics() by metric name, returning the results as a vector (or matrix, if multiple metrics are selected).
#'
metricsListToVector <- function(x,metric="total.area",class=1) {
  # parse the results of a ClassStat/PatchStat list-of-lists
  parseListByMetric <- function(y){
    res <- y[y$class == class,grepl(names(y),pattern=paste(metric,collapse="|",sep=""))];
      if(length(res)==0){
        return(NA);
      } else {
        return(res);
      }
  }
  # parse class stats
  o <-lapply(x[[1]],FUN=parseListByMetric)
  # parse patch stats
  p <- lapply(x[[2]],FUN=parseListByMetric)
  # bind all of our patch/class metrics into a matrix and return only those metrics
  # that have at least one valid value
  out <- cbind(unlist(o),unlist(p))
    return(out[,apply(out,2,sum,na.rm=T)>0])
}
