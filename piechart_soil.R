if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(ggplot2, ggtree, dplyr, tidyr, sp,pacman, maps, pipeR, grid, XML, gtable)


library(geojsonio)
library(leaflet.minicharts)
colnames(weather)[c(30,31,32)]=c('% Sand','% Silt','% Clay')

#Separate by years
for (i in 2014:2018) {
  data1=unique(weather[weather$Year==i,c(1,3,4,5,34,30:32)])
  library(leaflet)
  
  ##Option 1
  # transfrom .json file into a spatial polygons data frame
  states <- 
    geojson_read( 
      x = "https://raw.githubusercontent.com/PublicaMundi/MappingAPI/master/data/geojson/us-states.json"
      , what = "sp"
    )
  v <- leaflet(states) %>%
    setView(-96, 37.8, 4) %>%
    addProviderTiles("MapBox", options = providerTileOptions(
      id = "mapbox.light"))
  s<-v %>% addPolygons(color = "#444444", weight = 1, smoothFactor = 0.5,
                       opacity = 1.3, fillOpacity = 0.3)
  s<-s %>%
    addMinicharts(
      data1$long, data1$lat,
      type = "pie",
      chartdata = data1[, c('% Sand','% Silt','% Clay')], 
      colorPalette = colors
    )
  s=s%>%addLegend('bottomright',opacity = 1,title = paste('Soil composition of field trials - ',i,sep = ''),colors = colors,labels =c('% Sand','% Silt','% Clay') )
  
  #Save the leaflet maps
  library(htmlwidgets)
  saveWidget(s, file=paste("soil",i,".html",sep = ''))
  
}


########################
########################
##Option2

tilesURL <- "http://server.arcgisonline.com/ArcGIS/rest/services/Canvas/World_Light_Gray_Base/MapServer/tile/{z}/{y}/{x}"

basemap <- leaflet(width = "100%", height = "400px") %>%
  addTiles(tilesURL)
colors <- c('#c2b280', "#8b4513",'#bc8f8f')

basemap %>%
  addMinicharts(
    data1$long, data1$lat,
    title(main='Soil composition of field trials - 2014'),
    type = "pie",
    chartdata = data1[, c('% Sand','% Silt','% Clay')], 
    colorPalette = colors
  )
