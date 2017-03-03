# hopachAdditions

- *Provides additional distances and an option for a user-specified distance;*
- *removes repeated calculation of a distance matrix.*

**Distances:** 
* general binary;
* special case - Jaccard distance;
* another binary that does not fit into a general framework in (1);
* S-function distance that replaces a step cut-off function with a curve (S-shaped); additional parameters provided to vary the shape of the curve. Takes *p*-values as input.
* An option for a user-provided distance function with arbitrary tuning parameters.

### To use in R:
```
library(hopach)

source("hopachAdditions.R")

# example for user-provided distance:

user.distance <- function(x, y, par1, par2, par3) { ... }

hopach(data, d = "user", par1 = mypar1, par2 = mypar2, par3 = mypar3)
```
