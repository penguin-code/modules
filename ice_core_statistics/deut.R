#############################
# temperature and deuterium: a first look
#
# by James Bernhard
# 2021-02-14
#
# edited by Penny Rowe
# 2021-06-02
#
# data source:
#
# 1975 data:
# Lorius C, Merlivat L. Distribution of mean surface stable isotope values
# in East Antarctica: Observed changes with depth in the coastal area.
# Inter Assoc Sci Hydro Pub, 1977, 118: 127-137
#
# this script requires installation of the lattice package (for the graphics)
# via: 
# > library(package = "lattice")
#############################
tempDeutData <- read.csv("/Users/prowe/GitHub/PENGUIN/ice_core_statistics/deuteriumAndTemperature_1975.csv")

# first do it the way they do it in the paper
lattice::xyplot(deuterium~temperature, data=tempDeutData, type=c("p","r"))
deutModel <- lm(deuterium~temperature, data=tempDeutData)
coef(deutModel)
# this agrees (up to one decimal place) with the result in abstract of the paper
#    (and on page 4), which is: deut = 6.04 temp - 51
# however, in the abstract, they seem to imply this only holds for certain
#    elevations, which I haven't yet been able to derive.


# now look at it the other way around, which is what is needed for predicting temp
lattice::xyplot(temperature~deuterium, data=tempDeutData, type=c("p","r"))
tempModel <- lm(temperature~deuterium, data=tempDeutData)
coef(tempModel)
# So this is temp = 0.16 deut + 7.69 


# out of curiosity, how about just including a linear adjustment for elevation?
elevModel <- lm(temperature~deuterium + elevation, data=tempDeutData)
summary(elevModel)
# interesting: looks like that elevation term perhaps should be included

# check for inclusion of interaction term by comparing AIC with and without it
intModel <- lm(temperature~deuterium*elevation, data=tempDeutData)
summary(intModel)
AIC(elevModel, intModel)
# the interaction term model has the lower AIC so is perhaps to be preferred

# find confidence interval predicted temperatures using the intModel
predict(intModel, newdata=data.frame(deuterium=c(-265, -250), elevation=c(500, 1500)), interval="confidence")


# should we perhaps also look at terms to adjust for distance?
# possibly.




# 1975 data:
# Dahe et al
# This data includes both deltaD and delta18O
# Temperature is mean annual temperature in degrees C, or  (°C)
# Acc rate is in  (g/cm-2 a-1)
#
tempIsotopeData <- read.csv("/Users/prowe/GitHub/PENGUIN/ice_core_statistics/isotopesAndTemperature_1994b.csv")

# scatter plot temperature vs deltaD
lattice::xyplot(deltaD~temperature, data=tempIsotopeData, type=c("p","r"))
deutModel <- lm(deltaD~temperature, data=tempIsotopeData)
coef(deutModel)
# (Intercept) temperature 
# -80.275069    6.364778 
# => 6.36 T - 80.3
# Compare to paper
#    7.00 T - 30.2

# scatter plot temperature vs delta18O
lattice::xyplot(delta18O~temperature, data=tempIsotopeData, type=c("p","r"))
O18Model <- lm(delta18O~temperature, data=tempIsotopeData)
coef(O18Model)  
# (Intercept) temperature 
# -9.8818896   0.8230352  
# => 0.82 T - 9.88
#
# Compare to paper:
#    0.94 T - 4.99

# scatter plot delta18O vs temperature
lattice::xyplot(temperature~delta18O, data=tempIsotopeData, type=c("p","r"))
O18Model <- lm(temperature~delta18O, data=tempIsotopeData)
coef(O18Model)  
# (Intercept) delta18O 
# 10.447903    1.175302 
# => 1.175302 * delta18O + 10.447903
#


# .. Read in the ice core temperature data from the NOAA web site   
#    the result is a tibble file, which is like a fancy dataframe
#    If you have any trouble, reload the library and move on
website_ice_core_T = "https://www1.ncdc.noaa.gov/pub/data/paleo/icecore/antarctica/epica_domec/edc3deuttemp2007.txt"
ice_core <- readr::read_table(website_ice_core_T,
                              skip = 113, 
                              col_names = c('Bag', 'ztop', 'Age', 'Deuterium', 'Temperature'),
                              col_types = "ddddd")

