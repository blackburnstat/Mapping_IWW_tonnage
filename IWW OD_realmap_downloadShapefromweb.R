library(sf)
library(tidygraph)
library(igraph)
library(units)
library(osmdata)
library(rgrass7)
library(nabor)
library(tidyverse)
library(rgdal)
library(rgeos)
library(R.utils)
library(rnaturalearth)
library(eurostat)
library(lubridate)
library(stplanr)
library(leaflet)
#this file maps eurostat origin-destination tonnage quantities onto the real shape of the UNECE inland waterway network for 2020

####import Modified IWW network. A messy "download file" appraoch via the UNECE wiki (updated 2022-11-24)####
download.file("https://wiki.unece.org/download/attachments/115540014/Shape.zip?version=5&modificationDate=1668682363048&api=v2","C:/temp/alex3.zip",mode="wb")
unzip("C:/temp/alex3.zip",exdir = "C:/temp/alex4")
IWW_network<-st_read("C:/temp/alex4")
####function for complete conversion of SHP to network####
#Explained in the webpage https://r-spatial.org/r/2019/09/26/spatial-networks.html #
sf_to_tidygraph = function(x, directed = TRUE) {
  edges <- x %>%
    mutate(edgeID = c(1:n()))
  nodes <- edges %>%
    st_coordinates() %>%
    as_tibble() %>%
    rename(edgeID = L1) %>%
    group_by(edgeID) %>%
    slice(c(1, n())) %>%
    ungroup() %>%
    mutate(start_end = rep(c('start', 'end'), times = n()/2)) %>%
    mutate(xy = paste(.$X, .$Y)) %>% 
    mutate(nodeID = group_indices(., factor(xy, levels = unique(xy)))) %>%
    select(-xy)
  
  source_nodes <- nodes %>%
    filter(start_end == 'start') %>%
    pull(nodeID)
  
  target_nodes <- nodes %>%
    filter(start_end == 'end') %>%
    pull(nodeID)
  
  edges = edges %>%
    mutate(from = source_nodes, to = target_nodes)
  
  nodes <- nodes %>%
    distinct(nodeID, .keep_all = TRUE) %>%
    select(-c(edgeID, start_end)) %>%
    st_as_sf(coords = c('X', 'Y')) %>%
    st_set_crs(st_crs(edges))
  
  tbl_graph(nodes = nodes, edges = as_tibble(edges), directed = directed)
  
}
#next line turns the Shapefile into a NETWORK
graph<-sf_to_tidygraph(IWW_network,directed=FALSE)
#next bit defines the length of each edge
graph <- graph %>%
  activate(edges) %>%
  mutate(length = st_length(geometry))
#####Now lets add IWW data to visualise
#first the NUTS2 shapefiles.
map_nuts2<-get_eurostat_geospatial(output_class="spdf",resolution="60",nuts_level="2",year="2016")
centroids<-gCentroid(map_nuts2,byid=TRUE)
centroids<-spTransform(centroids,CRS=CRS("+init=epsg:4326"))
centroids$ID<-map_nuts2$NUTS_ID
rm(map_nuts2) #tidyup
rm(IWW_network)
####now the IWW activity data####
IWWData<-get_eurostat("iww_go_atygofl")
colnames(IWWData)<-c('Good', 'UnloadRegion','LoadRegion','unit','Country','Time',"Values")
IWWData$Time<-year(IWWData$Time)
#only choose proper nuts regions, and total goods, one year
IWWTest<-IWWData %>% 
  filter(unit=="THS_T",str_length(UnloadRegion)==4,str_length(LoadRegion)==4,Time==2020,Good=="TOTAL")
IWWTest$unit<-NULL
rm(IWWData) #tidyup
#next bit takes the maximum value of trade reported A>B (Austria+Slovakia both report AT31>SK01)
IWWTest<-IWWTest%>%
  group_by(LoadRegion,UnloadRegion) %>%
  summarise(bigcount=max(Values))
IntNuts<-merge(IWWTest,centroids,by.x="LoadRegion",by.y="ID")
IntNuts<-merge(IntNuts,centroids,by.x="UnloadRegion",by.y="ID")
colnames(IntNuts)<-c("UnloadXX","LoadXX","tonnage","o_long","o_lat","d_long","d_lat")
#Arrange by size 
IntNuts<-IntNuts %>%
  arrange(desc(tonnage)) 
#06/10/2022 TO include region X to region X values, the following is used.
IntNuts[IntNuts$UnloadXX==IntNuts$LoadXX,]$o_long<-IntNuts[IntNuts$UnloadXX==IntNuts$LoadXX,]$o_long+5
IntNuts[IntNuts$UnloadXX==IntNuts$LoadXX,]$o_lat<-IntNuts[IntNuts$UnloadXX==IntNuts$LoadXX,]$o_lat+5
IntNuts[IntNuts$UnloadXX==IntNuts$LoadXX,]$d_long<-IntNuts[IntNuts$UnloadXX==IntNuts$LoadXX,]$d_long+5
IntNuts[IntNuts$UnloadXX==IntNuts$LoadXX,]$d_lat<-IntNuts[IntNuts$UnloadXX==IntNuts$LoadXX,]$d_lat+5
length_NUTS<-nrow(IntNuts)
#### now translate the OD points onto the real network####
#choose the nodes 
coords <- graph %>%
  activate(nodes) %>%
  as_tibble() %>%
  st_as_sf() %>%
  st_transform(4326) %>%
  st_coordinates()
# coords <- nodes %>%
#   st_coordinates()
#note the 4326 above. it doesn't work with 3857
#make a big table of all the coordinates
node_index_master<-data.frame(rep(1,length_NUTS),rep(2,length_NUTS))
for (i in 1:length_NUTS){
  node_index_master[i,1]<- knn(data = coords, query = IntNuts[i,4:5]%>%
                         as.matrix(nrow=1), k = 1)$nn.idx
} #start point nodes
for (i in 1:length_NUTS){
  node_index_master[i,2]<- knn(data = coords, query = IntNuts[i,6:7]%>%
                                 as.matrix(nrow=1), k = 1)$nn.idx
} #end point nodes
my_big_epath_list<-list()
#next lines can take a few minutes
for (i in 1:length_NUTS){
  my_big_epath_list[i]<-
          shortest_paths(graph = graph,from = node_index_master[i,1], to = node_index_master[i,2], 
          output = 'epath',weights = graph %>%
            activate(edges) %>% 
            pull(SECTION_LE)
  )$epath
}
#final step: the below works but defines thousands of separate objects which is terrible. Can we make this into a list too?
for (i in 1:length_NUTS){
  assign(paste0("path_graph",i),   graph %>%
        subgraph.edges(eids = my_big_epath_list[i] %>% 
                         unlist()) %>%
                              as_tbl_graph() %>%
                       activate(edges) %>%
                         mutate(tonnage=IntNuts$tonnage[i]) %>%
                         as_tibble() %>%
                         st_as_sf())
        }
#now we just need all ~1900 things together as a list to feed into overline
list_paths<-mget(ls(pattern='^path\\_graph\\d'))
#next lines takes 10 minutes+ for total goods (to test the code, use a single good eg "GT01")
combined_graphs<-do.call(rbind,list_paths)
combined_graphs2<-st_set_crs(combined_graphs,3857)
#overline discussed here https://github.com/ropensci/stplanr/blob/HEAD/R/overline.R
overlined_segments<-overline(combined_graphs,attrib="tonnage")
####Map the results####
overlined_segments<-st_transform(overlined_segments, "+init=epsg:4326")
leaflet() %>% 
  addPolylines(data = overlined_segments, weight = ~(tonnage/300)^0.54, opacity = 0.98,fill =FALSE,smoothFactor = 2,label = paste0(overlined_segments$tonnage," thousand tonnes passed in",Year)) %>%
  addProviderTiles('CartoDB.Positron')

