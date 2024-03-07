#2023 12 06 
#Sophia Gosselin
#Software takes a table of lattitude and longitude values
#Plots those points on a map of the corresponding area.


## Libraries - Necessary
library(tidyverse)
library(sf)
library(mapview)
library(ggplot2)
library(maps)

## Input file of latitude and longitude values
input_name <- "all_metadata.table"
#input_name <- readLines(con = "stdin", n = 1) #should readin input from cmd. Example "cluster_cluster_72_latlong_log.txt"
lat_long_table <- read.table(input_name)


## Data Prep
# Edit column names
new_headers <- c("Phage","LAT","LON")
colnames(lat_long_table) <- new_headers

#for ggplot2 ploting
lat_long_df <- as.data.frame(lat_long_table)
my_sf_object <- st_as_sf(lat_long_df, coords = c('LON','LAT')) #THIS NEEDS TO BE THE ORDERING
my_sf_object <- st_set_crs(my_sf_object, 4326) #4326 is the -180 - 180 system of GPS choords.

#create a map object to plot points on
world1 <- sf::st_as_sf(map('world', plot = FALSE, fill = TRUE))

## BEGIN PLOTTING
# Create base map
mapview(lat_long_table, xcol = "Longitude", ycol = "Latitude", crs = 4269, grid = FALSE)

#with ggplot2
ggplot()+
  geom_sf(data = ) + 
  geom_sf(data = my_sf_object, aes(color = cluster)) #the aes lets you color based on a collumn of your df


