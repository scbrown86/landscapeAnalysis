# landscapeAnalysis
This is an MVP for a collection of algorithms for processing and analysing raster and vector data in the 'R' language.  I use it for almost every bit of landscape modeling I do (to fill-in-the blanks for things that the 'raster' and 'rgeos' packages cannot do) -- and plan on keeping it as a labor-of-love for some time to come. 

The code is dynamic (i.e., immature).  It needs documentation, optimization, and some double-checking if you are using it in a production environment or for published research. It's intended for deployment in high-performance unix environments, but should work on any computer that supports 'R'. Please commit bug-fixes, new features, or report issues as they arise.  I welcome collaborators.

# installation
```
library("devtools")
install_github("scbrown86/landscapeAnalysis")
```

# usage
```
library("landscapeAnalysis")
```

