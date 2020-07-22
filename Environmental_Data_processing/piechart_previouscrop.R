if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(ggplot2, ggtree, dplyr, tidyr, sp,pacman, maps, pipeR, grid, XML, gtable)


library(geojsonio)
library(leaflet.minicharts)
weatherdata = read.table(
  '~/Final_datasets_G2F/ALL_WEATHER/environmental_data_processing_1/Weather_soil_processing_1/weather_semihourly.txt',
  header = T,
  sep = '\t',
  na.strings = NA,
  quote = '',
  comment.char = "#",
  stringsAsFactors = F
)

##Option 1

#Separate by years

data1 = unique(weatherdata[, c(1, 2, 4, 5, 36)])
data1$Previous.crop = as.factor(data1$Previous.crop)
data1[which(data1$Previous.crop == ''), 'Previous.crop'] <- NA
data1.df <- split(data1, data1$Year)


library(leaflet)

##Option 1
# transfrom .json file into a spatial polygons data frame
states <-
  geojson_read(x = "https://raw.githubusercontent.com/PublicaMundi/MappingAPI/master/data/geojson/us-states.json"
               , what = "sp")
v <- leaflet(states) %>%
  setView(-96, 37.8, 4) %>%
  addProviderTiles("MapBox", options = providerTileOptions(id = "mapbox.light"))
s <-
  v %>% addPolygons(
    color = "#444444",
    weight = 1,
    smoothFactor = 0.5,
    opacity = 1.3,
    fillOpacity = 0.3
  )

library(viridis)
wardpal <- colorFactor(viridis(11), data1$Previous.crop)

names(data1.df) %>%
  purrr::walk(function(df) {
    s <<- s %>%
      addCircleMarkers(
        data = data1.df[[df]],
        lat = ~ lat,
        lng = ~ long,
        popup = ~ Previous.crop,
        radius = 10,
        fillColor  = ~ wardpal(Previous.crop),
        stroke = F,
        group = df,
        fillOpacity = 0.8,
        labelOptions = labelOptions(noHide = F,
                                    direction = 'auto')
      )
    
  })

s<-s %>%
  addLayersControl(
    baseGroups = names(data1.df),
    options = layersControlOptions(collapsed = FALSE)
  )

s<-s%>% addLegend(
  'bottomright',
  opacity = 1,
  title = c('Previous crop grown in G2F field trials'),
  pal = wardpal,
  values = ~ data1$Previous.crop
)

s



#Save the leaflet maps
library(htmlwidgets)
saveWidget(s, file = paste("previouscrop_layer_years", ".html", sep = ''))











############
############

##Option 2

#Separate by years
for (i in 2014:2018) {
  data1=unique(weatherdata[weatherdata$Year==i,c(1,2,4,5,36)])
  data1$Previous.crop=as.factor(data1$Previous.crop)
  data1[which(data1$Previous.crop==''),'Previous.crop']<-NA
  
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
  
  library(viridis)
  wardpal<-colorFactor(viridis(9)[2:8],data1$Previous.crop)
  
  s<-s %>%
    addCircleMarkers(
      data1$long, data1$lat,radius = 10,
      fillColor = ~wardpal(data1$Previous.crop),
      stroke = F,
      fillOpacity = 0.8,
    )
  
  s=s%>%addLegend('bottomright',opacity = 1,title = paste('Previous crop grown in G2F field trials - ',i,sep = ''),pal = wardpal,values = ~data1$Previous.crop )
 
  #Save the leaflet maps
  library(htmlwidgets)
  saveWidget(s, file=paste("previouscrop",i,".html",sep = ''))
  
}


########################
########################

##Option3

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
