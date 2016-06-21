# landscapeAnalysis
This is an MVP for a collection of algorithms for processing and analysing raster and vector data in the 'R' language.  I use it for almost every bit of landscape modeling I do to fill-in-the blanks for things that the 'raster' and 'rgeos' packages cannot do -- and plan on keeping it as a labor-of-love for some time to come. 

The code is dynamic (i.e., immature).  It needs documentation, optimization, and some double-checking if you are using it in a production environment or for published research. It's intended for deployment in high-performance UNIX environments, but should work on any computer that supports 'R'. Please file bug reports and issues as they arise.  

# installation
```
library("devtools")
install_github("ktaylora/landscape-analysis")
```

# usage
```
library("landscapeAnalysis")
```

