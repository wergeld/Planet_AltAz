library(tidyverse)
library(data.table)

## Get Date Difference From Ephemera
DaysSinceEphemera <- function(dateObs) {
  # Date of planetary constants "2000-01-01 12:00:00 UTC"
  dele <- as.POSIXct("2013-08-16 12:00:00", origin = "1970-01-01", tz="UTC")
  
  d <- as.numeric(difftime(dateObs, dele, units="days"))
  
  return(d)
}


Planet_Eph <- function(planetName) {
  PlanetEphData <- Osculating_Elements_AUG_16_2013 %>%
    dplyr::filter(name_planet == planetName)
  
  return(PlanetEphData)
}


# M = n * d + L - p 
# 
# n is daily motion
# d is the number of days since the date of the elements
# L is the mean longitude
# p is the longitude of perihelion 
# 
# M should be in range 0 to 360 degrees, add or subtract
#   multiples of 360 to bring M into this range.

MeanAnomaly_planet <- function(planetName, dayOffset) {
  PlanetEphData <- Planet_Eph(planetName)
  
  x <- PlanetEphData$n_planet * dayOffset + PlanetEphData$L_planet - PlanetEphData$p_planet
  
  return(x)
}


TrueAnomaly_planet <- function(planetName, dayOffset) {
  PlanetEphData <- Planet_Eph(planetName)
  
  MAp_Rad <- MeanAnomaly_planet(planetName, dayOffset) * (pi/180)
  e <- PlanetEphData$e_planet
  
  x <- MAp_Rad + (2 * e - 0.25 * e^3 + 5/96 * e^5) * sin(MAp_Rad) +
    (1.25 * e^2 - 11/24 * e^4) * sin(2*MAp_Rad) + 
    (13/12 * e^3 - 43/64 * e^5) * sin(3*MAp_Rad) + 
    103/96 * e^4 * sin(4*MAp_Rad) + 1097/960 * e^5 * sin(5*MAp_Rad)
  
  return(x) 
}


RadiusVector_planet <- function(planetName, dayOffset) {
  PlanetEphData <- Planet_Eph(planetName)
  
  a <- PlanetEphData$a_planet
  e <- PlanetEphData$e_planet
  v <- TrueAnomaly_planet(planetName, dayOffset)
  
  x <- a * (1 - e^2) / (1 + e * cos(v))
  
  return(x) 
}


HelioCoords_planet <- function(planetName, dayOffset) {
  PlanetEphData <- Planet_Eph(planetName)
  
  i <- PlanetEphData$i_planet * (pi/180)
  o <- PlanetEphData$o_planet * (pi/180)
  p <- PlanetEphData$p_planet * (pi/180)
  r <- RadiusVector_planet(planetName, dayOffset)
  v <- TrueAnomaly_planet(planetName, dayOffset)
  
  Xp <- r * (cos(o) * cos(v + p - o) - sin(o) * sin(v + p - o) * cos(i))
  Yp <- r * (sin(o) * cos(v + p - o) + cos(o) * sin(v + p - o) * cos(i))
  Zp <- r * (sin(v + p - o) * sin(i))
  
  vec <- c(Xp, Yp, Zp)
  return(vec) 
}


HelioCoords_Earth <- function() {
  PlanetEphData <- Planet_Eph("Earth")
  
  p <- PlanetEphData$p_planet * (pi/180)
  r <- RadiusVector_planet("Earth")
  v <- TrueAnomaly_planet("Earth")
  
  Xp <- r * cos(v + p)
  Yp <- r * sin(v + p)
  Zp <- 0
  
  vec <- c(Xp, Yp, Zp)
  return(vec) 
}


GeoEclipticCoords_Planet <- function(planetName, dayOffset) {
  PlanetEphData <- Planet_Eph(planetName)
  
  i <- PlanetEphData$i_planet * (pi/180)
  o <- PlanetEphData$o_planet * (pi/180)
  p <- PlanetEphData$p_planet * (pi/180)
  r <- RadiusVector_planet(planetName, dayOffset)
  v <- TrueAnomaly_planet(planetName, dayOffset)
  
  Xp <- r * (cos(o) * cos(v + p - o) - sin(o) * sin(v + p - o) * cos(i))
  Yp <- r * (sin(o) * cos(v + p - o) + cos(o) * sin(v + p - o) * cos(i))
  Zp <- r * (sin(v + p - o) * sin(i))
  
  PlanetEphDataEarth <- Planet_Eph("Earth")
  
  pe <- PlanetEphDataEarth$p_planet * (pi/180)
  re <- RadiusVector_planet("Earth", dayOffset)
  ve <- TrueAnomaly_planet("Earth", dayOffset)
  
  Xe <- re * cos(ve + pe)
  Ye <- re * sin(ve + pe)
  Ze <- 0
  
  Xgeo <- Xp - Xe
  Ygeo <- Yp - Ye
  Zgeo <- Zp - Ze
  
  vec <- c(Xgeo, Ygeo, Zgeo)
  return(vec)
}


GeoEquitorialCoords_Planet <- function(planetName, dayOffset) {
  vecGeo <- GeoEclipticCoords_Planet(planetName, dayOffset)
  
  Xq <- vecGeo[1]
  Yq <- vecGeo[2] * cos(ec_radians) - vecGeo[3] * sin(ec_radians)
  Zq <- vecGeo[2] * sin(ec_radians) + vecGeo[3] * cos(ec_radians)
  
  vec <- c(Xq, Yq, Zq)
  return(vec)
}


RA_DEC_Planet <- function(planetName, dayOffset) {
  vecGeo <- GeoEquitorialCoords_Planet(planetName, dayOffset)
  Xq <- vecGeo[1]
  Yq <- vecGeo[2]
  Zq <- vecGeo[3]
  
  Xq_rad <- Xq * (pi/180)
  Yq_rad <- Yq * (pi/180)
  Zq_rad <- Zq * (pi/180)
  
  alpha <- atan(Yq_rad/Xq_rad) * (180/pi)
  
  if (Xq < 0) {
    alpha <- alpha + 180
  } else {
    if (Yq < 0) {
      alpha <- alpha + 360
    }
  }
  
  delta <- atan( (Zq_rad) / (sqrt(Xq_rad*Xq_rad + Yq_rad*Yq_rad))) * (180/pi)
  
  distance = sqrt( Xq^2 + Yq^2 + Zq^2)
  
  alpha <- alpha / 15
  
  vec <- c(alpha, delta, distance)
  return(vec)
}


AltAz_Planet <- function(Planet, Lat, Long, Date) {
  
  # Taken from https://gist.github.com/matshofman/4145718
  # Day offset and Local Siderial Time
  d0 <- as.POSIXct("2000-01-01 12:00:00", tz="UTC")
  
  hour = hour(Date)
  minute = minute(Date)
  dayOffset <- as.numeric(difftime(Date, d0, units="days"))
  
  Days_SE <- DaysSinceEphemera(Date)
  RADEC <- RA_DEC_Planet(Planet, Days_SE)
  
  RA <- RADEC[1] * 15
  Dec <- RADEC[2]
  
  #LST <- mod((100.46 + 0.985647 * dayOffset + Long + 15 * (hour + minute / 60) + 360), 360)
  LST <- (100.46 + 0.985647 * dayOffset + Long + 15 * (hour + minute / 60) + 360) %% 360
  
  # Hour Angle
  #HA <- mod((LST - RA + 360), 360)
  HA <- (LST - RA + 360) %% 360
  
  # HA, DEC, Lat to Alt, AZ
  x <- cos(HA * (pi/180)) * cos(Dec * (pi/180))
  y <- sin(HA * (pi/180)) * cos(Dec * (pi/180))
  z <- sin(Dec * (pi/180))
  
  xhor <- x * cos((90 - Lat) * (pi/180)) - z * sin((90 - Lat) * (pi/180))
  yhor <- y
  zhor <- x * sin((90 - Lat) * (pi/180)) + z * cos((90 - Lat) * (pi/180))
  
  az <- atan2(yhor, xhor) * (180/pi) + 180
  alt <- asin(zhor) * (180/pi)
  
  vec <- data.frame(Date, Planet, alt, az)
  colnames(vec) = c("date","Planet","alt","az")
  return(vec)
}

