require(rgdal)
library(spdep)
#------------------------------------------------------------------------------#
# Set up -- Shapefiles
#------------------------------------------------------------------------------#

# # Load a shapefile.
# # You can do this easily with the rgdal or sf packages, and read the shape in an object. 
# # For both packages you need to provide dsn - the data source, 
# # which in the case of a shapefile is the directory, 
# # and layer - which is the shapefile name, minus extension:

Neighmap <- readOGR(dsn = 'D:/re-emergent RSV timing/tl_2017_us_state', layer = "tl_2017_us_state")

# leave states except for Florida
Neighmap <-  Neighmap[Neighmap@data$STUSPS%in%states[-c(which(states=="HI"),which(states=="AK"),which(states=="FL"))],]

# Sort the shapefile by state ID.
# This is the most important step!
Neighmap.sort <- Neighmap[order(Neighmap$GEOID),]

#------------------------------------------------------------------------------#
# Create a neighbor matrix
#------------------------------------------------------------------------------#

# Create neighbors from the sorted shapefile 

# The snap dist is governing the NA in neighbor matrix
# Can remove SNAP statement; 

# Neighbors can either be Queen (any zip that touches another zip - even at 
# corners) or Rook neighbors (share full sides -- not corners)

# if Queen = F then single point touch not allowed and Rook distance used.

neighb <- poly2nb(Neighmap.sort, queen = T, snap = sqrt(0.001))

# Check average number of links -- increase snap distance if it says there is 
# a zip with zero neighbors
neighb # Avg Number of links = 4.408163
#two region with no links, rm AK and HI

# Make neighbor matrix 
# if zero.policy FALSE stop with error for any empty neighbour sets, 
# if TRUE permit the weights list to be formed with zero-length weights vectors
# B is the basic binary coding, W is row standardised (both are commonly used)
neighbors_mat <- nb2mat(neighb, zero.policy = T, style = 'B')