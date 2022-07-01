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
#this file maps eurostat origin-destination tonnage quantities onto the real shape of the inland waterway network

####import Modified IWW network. A messy "download file" appraoch via the UNECE wiki####
download.file("https://wiki.unece.org/download/attachments/115540014/Shape.zip?version=1&modificationDate=1655840838175&api=v2&download=true","C:/temp/alex3.zip",mode="wb")
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
  filter(unit=="THS_T",str_length(UnloadRegion)==4,str_length(LoadRegion)==4,Time==2019,Good=="TOTAL")
IWWTest$unit<-NULL
rm(IWWData) #tidyup
#next bit takes the maximum value of trade reported A>B (Austria+Slovakia both report AT31>SK01)
IWWTest<-IWWTest%>%
  group_by(LoadRegion,UnloadRegion) %>%
  summarise(bigcount=max(Values))
IntNuts<-merge(IWWTest,centroids,by.x="LoadRegion",by.y="ID")
IntNuts<-merge(IntNuts,centroids,by.x="UnloadRegion",by.y="ID")
colnames(IntNuts)<-c("UnloadXX","LoadXX","tonnage","o_long","o_lat","d_long","d_lat")
#Arrange by size and remove A>A volumes
#if we don't filter them out, nothing happens as the map is based on edges and not vertices (not sure)
IntNuts<-IntNuts %>%
  arrange(desc(tonnage)) %>%
  filter(UnloadXX!=LoadXX)
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
#final step: the below works but defines 1909 objects which is terrible. Can we make this into a list too?
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
#now we just need all 1900 things together as a list to feed into overline
list_paths<-mget(ls(pattern='^path\\_graph\\d'))
#next lines takes 15-30 minutes for total goods (to test the code, use a single good eg "GT01")
combined_graphs<-do.call(rbind,list_paths)
combined_graphs2<-st_set_crs(combined_graphs,3857)
# test1<-cc %>%
#   st_as_sf()
 #setting the crs at different points is sometimes necessary, don't know why
#overline discussed here https://github.com/ropensci/stplanr/blob/HEAD/R/overline.R
overlined_segments<-overline(combined_graphs,attrib="tonnage")
####Map the results####
Euro2<-map_data(map="world") #define basemap
#coordinates for a reasonable Euro map
x1<--2
x2<-32
y1<-43
y2<-64
####finally, the graph####
ggplot() + #the ggplot version
  geom_polygon(data=Euro2,mapping=aes(x=long,y=lat,group=group),color="grey",fill="beige")+
  geom_sf(data = graph %>% activate(edges) %>% as_tibble() %>% st_as_sf(), col = 'darkgrey') +
  geom_sf(data = overlined_segments, lwd = (overlined_segments$tonnage^0.41)/30, col = 'red')+
  coord_sf(xlim=c(x1,x2),ylim=c(y1,y2),expand=TRUE,crs = st_crs(4326))+
  theme_classic()
####can we map it interactively in Leaflet??####
#result is partial success, but shows start and end-point connection for some reason.

#necessary to show in leaflet for some reason
overlined_segments<-st_transform(overlined_segments, "+init=epsg:4326")

leaflet(overlined_segments) %>% 
  addProviderTiles("CartoDB.Positron") %>% 
   addPolylines(weight = overlined_segments$tonnage/10000) %>% 
  addPolygons(weight = overlined_segments$tonnage/10000)
#polylines seems to work correctly, but misses out many country data (eg france)
#polygons seems to work correctly but connects start and end points

