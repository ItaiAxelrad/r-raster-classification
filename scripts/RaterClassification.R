#install the packages 
if (!require(raster)) { install.packages("raster") }
if (!require(rgdal)) { install.packages("rgdal") }
if (!require(rasterVis)) { install.packages("rasterVis") }
if (!require(geosphere)) { install.packages("geosphere") }
if (!require(randomForest)) { install.packages("randomForest") }
# load packages
library(rgdal)
library(raster)
library(rasterVis)
library(geosphere)
library(randomForest)
# Load data 
load("data/GewataB2.rda")
load("data/GewataB3.rda")
load("data/GewataB4.rda")

# Check out the attributes
GewataB2
# Some basic statistics using cellStats()
cellStats(GewataB2, stat = max)
cellStats(GewataB2, stat = mean)
# This is equivalent to:
maxValue(GewataB2)
# What is the maximum value of all three bands? 
max(c(maxValue(GewataB2), maxValue(GewataB3), maxValue(GewataB4)))
# summary() is useful function for a quick overview 
summary(GewataB2)
# Put the 3 bands into a RasterBrick object to summarize together 
gewata <- brick(GewataB2, GewataB3, GewataB4)
# 3 histograms in one window (automatic, if a RasterBrick is supplied) 
hist(gewata)
par(mfrow = c(1, 1)) # reset plotting window
hist(gewata, xlim = c(0, 5000), ylim = c(0, 750000), breaks = seq(0, 5000, by = 100))
pairs(gewata)
ndvi <- overlay(GewataB4, GewataB3, fun = function(x, y) {(x - y) / (x + y) })
plot(ndvi)

# Load the data and check it out 
load("data/vcfGewata.rda")
vcfGewata
plot(vcfGewata)
summary(vcfGewata)
hist(vcfGewata)
vcfGewata[vcfGewata > 100] <- NA
plot(vcfGewata)
summary(vcfGewata)
hist(vcfGewata)
gewata <- calc(gewata, fun = function(x) x / 10000)
# Make a new RasterBrick of covariates by adding NDVI and VCF layers 
covs <- addLayer(gewata, ndvi, vcfGewata)
plot(covs)
names(covs) <- c("band2", "band3", "band4", "NDVI", "VCF")
plot(covs)
# Load the training polygons 
load("data/trainingPoly.rda")
# Superimpose training polygons onto NDVI plot 
plot(ndvi)
plot(trainingPoly, add = TRUE)
# Inspect the data slot of the trainingPoly object 
trainingPoly@data
# The 'Class' column is actually an ordered factor type
trainingPoly@data$Class
str(trainingPoly@data$Class)
# We can convert to integer by using the as.numeric() function, 
# which takes the factor levels
trainingPoly@data$Code <- as.numeric(trainingPoly@data$Class)
trainingPoly@data
# Assign 'Code' values to raster cells (where they overlap) 
classes <- rasterize(trainingPoly, ndvi, field = 'Code')
# Plotting
# Define a colour scale for the classes (as above) 
# corresponding to: cropland, forest, wetland 
cols <- c("orange", "dark green", "light blue")
# Plot without a legend
plot(classes, col = cols, legend = FALSE)
# Add a customized legend
legend("topright", legend = c("cropland", "forest", "wetland"), fill = cols, bg = "white")
covmasked <- mask(covs, classes)
plot(covmasked)
# Combine this new brick with the classes layer to make our input training dataset
names(classes) <- "class"
trainingbrick <- addLayer(covmasked, classes)
plot(trainingbrick)
# Extract all values into a matrix 
valuetable <- getValues(trainingbrick)
valuetable <- na.omit(valuetable)
valuetable <- as.data.frame(valuetable)
head(valuetable, n = 10)
tail(valuetable, n = 10)
valuetable$class <- factor(valuetable$class, levels = c(1:3))
val_crop <- subset(valuetable, class == 1)
val_forest <- subset(valuetable, class == 2)
val_wetland <- subset(valuetable, class == 3)
# 1. NDVI
par(mfrow = c(3, 1))
hist(val_crop$NDVI, main = "cropland", xlab = "NDVI", xlim = c(0, 1), ylim = c(0, 4000), col = "orange")
hist(val_forest$NDVI, main = "forest", xlab = "NDVI", xlim = c(0, 1), ylim = c(0, 4000), col = "dark green")
hist(val_wetland$NDVI, main = "wetland", xlab = "NDVI", xlim = c(0, 1), ylim = c(0, 4000), col = "light blue")
par(mfrow = c(1, 1))
# 2. VCF
par(mfrow = c(3, 1))
hist(val_crop$VCF, main = "cropland", xlab = "% tree cover", xlim = c(0, 100), ylim = c(0, 7500), col = "orange")
hist(val_forest$VCF, main = "forest", xlab = "% tree cover", xlim = c(0, 100), ylim = c(0, 7500), col = "dark green")
hist(val_wetland$VCF, main = "wetland", xlab = "% tree cover", xlim = c(0, 100), ylim = c(0, 7500), col = "light blue")
par(mfrow = c(1, 1))
# 3. Bands 3 and 4 (scatterplots)
plot(band4 ~ band3, data = val_crop, pch = ".", col = "orange", xlim = c(0, 0.2), ylim = c(0, 0.5))
points(band4 ~ band3, data = val_forest, pch = ".", col = "dark green")
points(band4 ~ band3, data = val_wetland, pch = ".", col = "light blue")
legend("topright", legend = c("cropland", "forest", "wetland"), fill = c("orange", "dark green", "light blue"), bg = "white")

# Construct a random forest model
# Covariates (x) are found in columns 1 to 5 of valuetable
# Training classes (y) are found in the 'class' column of valuetable 
# Caution: this step takes fairly long!
# but can be shortened by setting importance=FALSE 
modelRF <- randomForest(x = valuetable[, c(1:5)], y = valuetable$class, importance = TRUE)
# Inspect the structure and element names of the resulting model modelRF
class(modelRF)
str(modelRF)
names(modelRF)
# Inspect the confusion matrix of the OOB error assessment
modelRF$confusion
# to make the confusion matrix more readable
colnames(modelRF$confusion) <- c("cropland", "forest", "wetland", "class.error")
rownames(modelRF$confusion) <- c("cropland", "forest", "wetland")
modelRF$confusion
varImpPlot(modelRF)
# Double-check layer and column names to make sure they match 
names(covs)
# [1] "band2" "band3" "band4" "NDVI" "VCF"
names(valuetable)
# [1] "band2" "band3" "band4" "NDVI" "VCF" "class"
# Predict land cover using the RF model
predLC <- predict(covs, model = modelRF, na.rm = TRUE)
# Plot the results
# recall: 1 = cropland, 2 = forest, 3 = wetland
cols <- c("orange", "dark green", "light blue")
plot(predLC, col = cols, legend = FALSE)
legend("bottomright", legend = c("cropland", "forest", "wetland"), fill = cols, bg = "white")

# --- Unsupervised classification: k-means --- #
valuetable <- getValues(covs)
head(valuetable)
km <- kmeans(na.omit(valuetable), centers = 3, iter.max = 100, nstart = 10)
# km contains the clusters (classes) assigned to the cells
head(km$cluster)
unique(km$cluster)
# displays unique values
# Create a blank raster with default values of 0 
rNA <- setValues(raster(covs), 0)
# Loop through layers of covs
# Assign a 1 to rNA wherever an NA is enountered in covs 
for (i in 1:nlayers(covs)) { rNA[is.na(covs[[i]])] <- 1 }
# Convert rNA to an integer vector
rNA <- getValues(rNA)
# Convert valuetable to a data.frame
valuetable <- as.data.frame(valuetable)
# If rNA is a 0, assign the cluster value at that position 
valuetable$class[rNA == 0] <- km$cluster
# If rNA is a 1, assign an NA at that position
valuetable$class[rNA == 1] <- NA
# Create a blank raster
classes <- raster(covs)
# Assign values from the 'class' column of valuetable 
classes <- setValues(classes, valuetable$class)
plot(classes, legend = FALSE, col = c("dark green", "orange", "light blue"))
# Make an NA-value raster based on the LC raster attributes
formask <- setValues(raster(predLC), NA)
# Assign 1 to formask to all cells corresponding to the forest class 
formask[predLC == 2] <- 1
plot(formask, col = "dark green", legend = FALSE)
# Group raster cells into clumps based on the Queen's Case 
if (!file.exists(fn <- "data/clumformask.grd")) { forestclumps <- clump(formask, directions = 8, filename = fn) } else { forestclumps <- raster(fn) }
plot(forestclumps)
# Assign freqency table to a matrix 
clumpFreq <- freq(forestclumps)
head(clumpFreq)
tail(clumpFreq)
# Coerce freq table to data.frame
clumpFreq <- as.data.frame(clumpFreq)
# which rows of the data.frame are only represented by one cell? 
str(which(clumpFreq$count == 1))
# which values do these correspond to? 
str(clumpFreq$value[which(clumpFreq$count == 1)])
# Put these into a vector of clump ID's to be removed
excludeID <- clumpFreq$value[which(clumpFreq$count == 1)]
# Make a new forest mask to be sieved
formaskSieve <- formask
# Assign NA to all clumps whose IDs are found in excludeID 
formaskSieve[forestclumps %in% excludeID] <- NA
# Zoom in to a small extent to check the results
# Note: you can define your own zoom by using 
e <- drawExtent()
e <- extent(c(811744.8, 812764.3, 849997.8, 850920.3))
opar <- par(mfrow = c(1, 2))
# allow 2 plots side-by-side
plot(formask, ext = e, col = "dark green", legend = FALSE)
plot(formaskSieve, ext = e, col = "dark green", legend = FALSE)
par(opar)
# reset plotting window
load("data/lulcGewata.rda")
# Check out the distribution of the values 
freq(lulcGewata)
hist(lulcGewata)
load("data/LUTGewata.rda")
LUTGewata
lulc <- as.factor(lulcGewata)
# assign a raster attribute table (RAT) 
levels(lulc) <- LUTGewata
lulc
classes <- layerize(lulc)
# Layer names follow the order of classes in the LUT 
names(classes) <- LUTGewata$Class
plot(classes, legend = FALSE)
forest <- raster(classes, 5)
# is equivalent to:
forest <- classes[[5]]
# or (since the layers are named): 
forest <- classes$forest
# Replace 0's (non-forest) with NA's
forest[forest == 0] <- NA
plot(forest, col = "dark green", legend = FALSE)
