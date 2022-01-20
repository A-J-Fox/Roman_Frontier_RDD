library(tmap)
library(maptools)
library(raster)
library(SpatialRDD)
library(ggplot2)
library(sf)
library(stars)
library(rgdal)
library(rdrobust)
library(xtable)
library(stargazer)
library(lwgeom)
library(margins)


#loading tif files into stata
hungary <-raster("~/Desktop/DFP/tif files/Hungary_lights_no0.tif")
hungary
plot(hungary)
col = heat.colors(63)
cellStats(hungary, range)
image(hungary, col=col)

# creating plot of hungary lights
hungary2 <- as.data.frame(hungary)
hist(hungary2$Hungary_lights_no0)

# plot of the brightest regions
image(hungary, zlim=c(50,63), col=col)

# loading as a stars object to convert to sf object
hungary_stars <- read_stars("~/Desktop/DFP/tif files/Hungary_lights_no0.tif")
plot(hungary_stars)
# this converts all pixels to point geometries
hungary.sf <- st_as_sf(hungary_stars, as_points = TRUE, merge = FALSE, long = TRUE, crs = 4326, coords = c("x", "y"))
plot(hungary.sf)

#importing discontinuity/danube. For some reason "~" doesn't work here
danube <- readOGR(dsn="/Users/andrewfox/Desktop/DFP/danube", layer="discont")
plot(danube)

danube.sf <- st_read("/Users/andrewfox/Desktop/DFP/danube", layer = "discont")
ggplot()+geom_sf(data=danube.sf)

# loading elevation data and merging with luminosity
hungary_elev_stars <- read_stars("~/Desktop/DFP/tif files/Hungary_elevation.tif")
plot(hungary_elev_stars)
hungary.elev.sf <- st_as_sf(hungary_elev_stars, as_points = TRUE, merge = FALSE, long = TRUE, crs = 4326, coords = c("x", "y"))
hungary.sf <- st_join(hungary.sf, hungary.elev.sf, join = st_nearest_feature)

# loading and merging population density data
hungary_pop_stars <- read_stars("~/Desktop/DFP/hungary.pop.tif")
hungary.pop.sf <- st_as_sf(hungary_pop_stars, as_points = TRUE, merge = FALSE, long = TRUE, crs = 4326, coords = c("x", "y"))
hungary.sf <- st_join(hungary.sf, hungary.pop.sf, join = st_nearest_feature)

# loading and merging soil data
hungary_soil_stars <- read_stars("~/Desktop/DFP/soil.tif")
hungary.soil.sf <- st_as_sf(hungary_soil_stars, as_points = TRUE, merge = FALSE, long = TRUE, crs = 4326, coords = c("x", "y"))
hungary.sf <- st_join(hungary.sf, hungary.soil.sf, join = st_nearest_feature)
# coding as factor variable
hungary.sf$soil.tf.f <- factor(hungary.sf$soil.tif)
is.factor(hungary.sf$soil.tf.f)

# a plot with the boundary
plot(hungary)
plot(danube, add = TRUE)

#importing the treatment polygon cut in QGIS
treat <- readOGR(dsn="/Users/andrewfox/Desktop/DFP/treatment", layer="treat")
plot(hungary)
plot(treat, add = TRUE)

treat.sf <- st_read("/Users/andrewfox/Desktop/DFP/treatment", layer = "treat")
ggplot()+geom_sf(data=treat.sf)
# setting the coordinate system, so they are consistent
st_crs(treat.sf) = 4326
st_crs(hungary.sf) = 4326
st_crs(danube.sf) = 4326

# creating the treatment area
hungary.sf$treated <- assign_treated(hungary.sf, treat.sf, id = "geometry")
# a warning appears, but for this small area, it shouldn't be a problem
summary(hungary.sf$treated)

# creating ln of luminosity
hungary.sf$log_lights <- log(hungary.sf$Hungary_lights_no0.tif)

# a first simple pooled binary variable regression
lm1 <- lm(log_lights ~ treated, data = hungary.sf)
lm2 <- lm(log_lights ~ treated + Hungary_elevation.tif + hungary.pop.tif + soil.tf.f , data = hungary.sf)

# calculating distance between each point and the danube
hungary.sf$dist2cutoff <- as.numeric(sf::st_distance(hungary.sf, danube.sf))
summary(hungary.sf$dist2cutoff)

# simple pooled regression within 3,5,10 km distance
lm4 <- lm(log_lights ~ treated + Hungary_elevation.tif + hungary.pop.tif, data = hungary.sf[hungary.sf$dist2cutoff < 5000, ])
lm5 <- lm(log_lights ~ treated + Hungary_elevation.tif + hungary.pop.tif + soil.tf.f, data = hungary.sf[hungary.sf$dist2cutoff < 10000, ])
lm6 <- lm(log_lights ~ treated + Hungary_elevation.tif + hungary.pop.tif + soil.tf.f, data = hungary.sf[hungary.sf$dist2cutoff < 3000, ]) 
lm7 <- lm(log_lights ~ treated + Hungary_elevation.tif + hungary.pop.tif + soil.tf.f, data = hungary.sf[hungary.sf$dist2cutoff < 5000, ])
lm8 <- lm(log_lights ~ treated + Hungary_elevation.tif + hungary.pop.tif + soil.tf.f, data = hungary.sf[hungary.sf$dist2cutoff < 10000, ])

stargazer(lm1, lm2, lm6, lm7, lm8, title="Pooled Linear Regressions", align=TRUE)

# note that around the border, the effect is positive

lm3 <- lm(Hungary_lights_no0.tif ~ treated, data = hungary.sf[hungary.sf$dist2cutoff < 15000, ])
coef(summary(lm3))[, "Std. Error"]


# re-centering data at zero
hungary.sf$distrunning <- hungary.sf$dist2cutoff
hungary.sf$distrunning[hungary.sf$treated == 0] <- -1 * hungary.sf$distrunning[hungary.sf$treated == 0]

# plot of lights around Danube
ggplot(data = hungary.sf, aes(x = distrunning, y = Hungary_lights_no0.tif)) + geom_point() + geom_vline(xintercept = 0, col = "red")

# "naive" rdd with local linear estimation
covs1 <- cbind(hungary.sf$Hungary_elevation.tif, hungary.sf$hungary.pop.tif)
covs2 <- cbind(hungary.sf$Hungary_elevation.tif, hungary.sf$hungary.pop.tif, hungary.sf$soil.tf.f)
rd1 <- rdrobust(hungary.sf$log_lights, hungary.sf$distrunning, c = 0)
rd1.1 <- rdrobust(hungary.sf$log_lights, hungary.sf$distrunning, c = 0, p = 2)
rd2 <- rdrobust(hungary.sf$log_lights, hungary.sf$distrunning, covs = hungary.sf$Hungary_elevation.tif, c = 0)
rd3 <- rdrobust(hungary.sf$log_lights, hungary.sf$distrunning, covs = covs1, c = 0)
rd4 <- rdrobust(hungary.sf$log_lights, hungary.sf$distrunning, covs = covs2, c = 0)
rd5 <- rdrobust(hungary.sf$log_lights, hungary.sf$distrunning, covs = covs1, c = 0, p=2)
rd6 <- rdrobust(hungary.sf$log_lights, hungary.sf$distrunning, covs = covs2, c = 0, p =2)

rdplot(hungary.sf$log_lights, hungary.sf$distrunning, covs = covs2, c = 0, p =2,
       y.label = "Luminosity", x.label = "distance to border")

rdplot(hungary.sf$log_lights, hungary.sf$distrunning, covs = covs2, c = 0, p =1,
       y.label = "Luminosity", x.label = "distance to border")

rdplot(hungary.sf$log_lights, hungary.sf$distrunning, covs = covs2, c = 0, p =6, h = 60,
       y.label = "Luminosity", x.label = "distance to border")

rdplot(hungary.sf$Hungary_lights_no0.tif, hungary.sf$distrunning, c = 0, ci = 95, 
       kernel = "triangular", y.label = "Luminosity", x.label = "distance to border")

# for a nice 4-part plot
plot(hungary.sf)

# discretize border (approx every 6 kilometers)
borderpoints.sf <- discretise_border(cutoff = danube.sf, n = 50)
borderpoints.sf$id <- 1:nrow(borderpoints.sf)

# approx every 1 kilometer
borderpoints.sf2 <- discretise_border(cutoff = danube.sf, n = 10)
borderpoints.sf2$id <- 1:nrow(borderpoints.sf2)

# looping over boundary points with spatialrd(), with 50 sections. non-parametric estimate
bdd1 <- spatialrd(y = "log_lights", data = hungary.sf, cutoff.points = borderpoints.sf, 
treated = "treated", minobs = 10, spatial.object = F)
knitr::kable(bdd1)
bdd2 <- spatialrd(y = "log_lights", data = hungary.sf, cutoff.points = borderpoints.sf, 
                  treated = "treated", minobs = 10, covs = covs1, spatial.object = F)
bdd3 <- spatialrd(y = "log_lights", data = hungary.sf, cutoff.points = borderpoints.sf, 
                  treated = "treated", minobs = 10, covs = covs2, spatial.object = F)



plotspatialrd(bdd2, map = T)

#calculating average treatment effect
mean(bdd1$Estimate)
mean(bdd1$pvalR)

# poisson regressions
m1 <- glm(Hungary_lights_no0.tif ~ treated + Hungary_elevation.tif + hungary.pop.tif + soil.tf.f, family="poisson", data=hungary.sf)

m2 <- glm(Hungary_lights_no0.tif ~ treated + Hungary_elevation.tif 
          + hungary.pop.tif + soil.tf.f, family="poisson", data=hungary.sf[hungary.sf$dist2cutoff < 3000, ])

m3 <- glm(Hungary_lights_no0.tif ~ treated + Hungary_elevation.tif 
          + hungary.pop.tif + soil.tf.f, family="poisson", data=hungary.sf[hungary.sf$dist2cutoff < 5000, ])

m4 <- glm(Hungary_lights_no0.tif ~ treated + Hungary_elevation.tif 
          + hungary.pop.tif + soil.tf.f, family="poisson", data=hungary.sf[hungary.sf$dist2cutoff < 10000, ])
# seeing marginal effects
margins(m1)
margins(m2)
margins(m3)
margins(m4)

# compiling to LaTeX table!
stargazer(m1,m2,m3,m4)

# woohoo!

