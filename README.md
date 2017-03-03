# hopachAdditions

- **Provides additional distances and an option for a user-specified distance;**
- **removes repeated calculation of a distance matrix.**

**Distances:** 
1. general binary;
2. special case - Jaccard distance;
3. another binary that does not fit into a general framework in (1);
4. S-function distance that replaces a step cut-off function with a curve (S-shaped); additional parameters provided to vary the shape of the curve. Takes *p*-values as input.
5. An option for a user-provided distance function with an aritrary number of parameters.

### To use in R:
```
library(hopach)

source("hopachAdditions.R")

# example for user-provided distance:

user.distance <- function(x, y, par1, par2, par3) { ... }

hopach(data, d = "user", par1 = mypar1, par2 = mypar2, par3 = mypar3)
```
