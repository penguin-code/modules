#############################
# temperature and deuterium: a first look
#
# by James Bernhard
# 2021-02-14
#
# data source:
#
# Lorius C, Merlivat L. Distribution of mean surface stable isotope values
# in East Antarctica: Observed changes with depth in the coastal area.
# Inter Assoc Sci Hydro Pub, 1977, 118: 127-137
#
# this script requires installation of the lattice package (for the graphics)
# via: 
# > library(package = "lattice")
#############################
tempDeutData <- read.csv("/Users/prowe/GitHub/PENGUIN/ice_core_statistics/deuteriumAndTemperature.csv")
#tempDeutData <- read.csv("deuteriumAndTemperature.csv")

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
