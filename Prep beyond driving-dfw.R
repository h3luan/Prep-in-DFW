path <- paste0(here(),"/Data/GTFS")
r5r.core <- setup_r5(data_path = path,verbose = F)

##The zip codes in Dallas and Fort Worth
dfw.zip <- read.csv(paste0(here(),"/Data/DFW-2017-21-new HIV.csv"))
all.zip <- st_read(paste0(here(),"/Data/zcta17-centroids-popweighted.gpkg"))
all.zip$Zip.Code <- as.numeric(all.zip$ZCTA5)

zip.join <- left_join(all.zip,dfw.zip,"Zip.Code")
dfw.join <- filter(zip.join,!is.na(Total.cases))

origins <- data.frame(id=dfw.join$Zip.Code,lon=dfw.join$xcnt_pp,lat=dfw.join$ycnt_pp)

dfw.prep <- st_read(paste0(here(),"/Data/PrEP-DFW-May24.shp"))

destinations <- data.frame(id=1:dim(dfw.prep)[1],lon=dfw.prep$X,lat=dfw.prep$Y)

mode <- "TRANSIT"
max.trip.duration <- 720

departure.time <- as.POSIXct("22-07-2024 9:00:00",format = "%d-%m-%Y %H:%M:%S")

travel.perc <- 50
travel.window <- 420 ##9am-4pm

transit.time <- travel_time_matrix(r5r_core = r5r.core,origins = origins,destinations = destinations,mode = mode,departure_datetime = departure.time,max_trip_duration = max.trip.duration,time_window = travel.window,percentiles = travel.perc)

walking.time <- travel_time_matrix(r5r_core = r5r.core,origins = origins,destinations = destinations,mode = "walk",time_window = travel.window, percentiles = travel.perc)
# 
cycling.time <- travel_time_matrix(r5r_core = r5r.core,origins = origins,destinations = destinations,mode = "bicycle",max_bike_time = 1440,time_window = travel.window, percentiles = travel.perc)

driving.time <- travel_time_matrix(r5r_core = r5r.core,origins = origins,destinations = destinations,mode = "car",time_window = travel.window, percentiles = travel.perc)

driving.nearest <- group_by(driving.time,from_id) %>%
  summarise(.,car.dis = min(travel_time_p50)) %>%
  mutate(.,GEOID=as.numeric(from_id))

transit.nearest <- group_by(transit.time,from_id) %>%
  summarise(.,transit.dis = min(travel_time_p50)) %>%
  mutate(.,GEOID=as.numeric(from_id))

cycle.nearest <- group_by(cycling.time,from_id) %>%
  summarise(.,cycle.dis = min(travel_time_p50)) %>%
  mutate(.,GEOID=as.numeric(from_id))

walk.nearest <- group_by(walking.time,from_id) %>%
  summarise(.,walk.dis = min(travel_time_p50)) %>%
  mutate(.,GEOID=as.numeric(from_id))

## The correlation btw driving and transit time
rcorr(driving.nearest$car.dis,transit.nearest$transit.dis,type = "pearson") ##rho = 0.81, p-value = 0 

quantile(driving.nearest$car.dis,probs = seq(0,1,1/3))
quantile(transit.nearest$transit.dis,probs = seq(0,1,1/3))


# Plot maps of multimodal PrEP accessibility ------------------------------

dfw.zcta21 <- st_read(paste0(here(),"/Data/DFW-zcta21.shp"))
dfw.zcta.join <- left_join(dfw.zcta21,transit.nearest,"GEOID") %>%
  left_join(.,select(driving.nearest,c("GEOID","car.dis"))) %>%
  left_join(.,select(walk.nearest,c("GEOID","walk.dis"))) %>%
  left_join(.,select(cycle.nearest,c("GEOID","cycle.dis")))

dfw.zcta.join$walk.dis[which(is.na(dfw.zcta.join$walk.dis))] <- 9999
dfw.zcta.join$cycle.dis[which(is.na(dfw.zcta.join$cycle.dis))] <- 9999

## Function for plotting maps
plot_map <- function(layer,field,title,classnum,breaks,pal,labels){
  plot <- tm_shape(shp = layer) +
    tm_fill(field,title=title,n=classnum,breaks = breaks,palette = pal,labels = labels,alpha = 0.8,legend.is.portrait = T)+
    tm_borders(lwd = 0.3)+
    tm_layout(frame = F,
              panel.label.size = 1.5)+ ##,legend.position = c(0.25,0.88)
    tm_compass(position = c(0.01,0.8))+
    tm_scale_bar(position = c(0.4,0))
  return(plot)
}

class.num <- 3
class.all <- 5

breaks.transit <- c(0,30,60,720)
breaks.car <- c(0,10,15,30)

breaks.all <- c(0,10.1,20.1,30.1,60.1,1000000)

pal.transit <- rev(brewer.pal(class.all,"Greens"))
pal.car <- rev(brewer.pal(class.all,"Blues"))
pal.walk <- rev(brewer.pal(class.all,"Oranges"))
pal.cycle <- rev(brewer.pal(class.all,"PuRd"))

labels.transit <- c("<=30","30-60",">60")
labels.car <- c("<=10","10-15","15-30")

labels.all <- c("[0,10]","(10,20]","(20,30]","(30,60]",">60")

##Map of PrEP accessibility based on public transit
plot.transit <- plot_map(dfw.zcta.join,"transit.dis","Public transit",class.all,breaks.all,pal.transit,labels.all)

##Map of PrEP accessibility based on driving
plot.car <- plot_map(dfw.zcta.join,"car.dis","Driving",class.all,breaks.all,pal.car,labels.all)

plot.walk <- plot_map(dfw.zcta.join,"walk.dis","Walking",class.all,breaks.all,pal.walk,labels.all)

plot.cycle <- plot_map(dfw.zcta.join,"cycle.dis","Cycling",class.all,breaks.all,pal.cycle,labels.all)

tmap_save(tm = tmap_arrange(plot.transit,plot.car,nrow = 1),filename = "PrEP accessibility-transit-car-Aug 31.png",width = 10,height= 4,dpi = 600)

tmap_save(tm = tmap_arrange(plot.walk,plot.cycle,nrow = 1),filename = "PrEP accessibility-walk-cycle-Aug 31.png",width = 10,height= 4,dpi = 600)

##Bivariate of transit- and driving-based PrEP accessibiilty
##Prepare the data
dfw.zcta.join$transit.class <- cut(dfw.zcta.join$transit.dis,breaks = breaks.transit)
dfw.zcta.join$car.class <- cut(dfw.zcta.join$car.dis,breaks = breaks.car)

dfw.bivariate.labels <- biscale::bi_class_breaks(dfw.zcta.join, x = transit.class, y = car.class,dim = 3, dig_lab = 3, split = FALSE)
dfw.bivariate.labels$bi_x[3] <- ">60"

dfw.bivariate.data <- bi_class(dfw.zcta.join,x="transit.class",y="car.class",style = "quantile",dim = 3)
##Prepare the map
dfw.bivariate.map <- ggplot()+
  geom_sf(data = dfw.bivariate.data,mapping = aes(fill=bi_class),color="white",size=0.1,show.legend = F)+
  bi_scale_fill(pal = "GrPink",dim=3)+
  # labs(title = "New HIV diagnosis VS distance to nearest PrEP clinic")+
  labs(title = "")+
  bi_theme()

dfw.bivariate.legend <- bi_legend(pal = "GrPink",dim = 3,xlab = "Longer pubic transit time",ylab = "Longer driving time",breaks=dfw.bivariate.labels,size = 6)

dfw.bivariate.plot <- ggdraw()+
  draw_plot(dfw.bivariate.map,0,0,1,1)+
  draw_plot(dfw.bivariate.legend,0,0.67,0.27,0.27)

ggsave("DFW_Bivariate_PrEP-Aug 13.png",plot = dfw.bivariate.plot,width=5,height=4,dpi = 600) ##,width = 7,height = 9


# Plot the socioeconomic deprivation index --------------------------------

SDOH20 <- read_excel(paste0(here(),"/Data/SDOH_2020_ZIPCODE_1_0.xlsx"),2)
SDOH.sel <- select(SDOH20, c("ACS_PCT_PERSON_INC_BELOW99_ZC","ACS_PCT_LT_HS_ZC","ACS_MEDIAN_HH_INC_ZC","ACS_PCT_UNEMPLOY_ZC","ACS_GINI_INDEX_ZC")) %>%
  mutate(.,GEOID=as.numeric(SDOH20$ZIPCODE))

SDOH.sel2 <- select(SDOH20,c("ACS_PCT_BLACK_ZC","ACS_PCT_HISPANIC_ZC","ACS_PCT_UNINSURED_ZC","ACS_PCT_ASIAN_ZC","CEN_POPDENSITY_ZC","ACS_PCT_WHITE_ZC")) %>%
  mutate(.,GEOID=as.numeric(SDOH20$ZIPCODE))

dfw.zcta.sdoh <- left_join(dfw.zcta.join,SDOH.sel,"GEOID")

## Impute missing values
which(is.na(dfw.zcta.sdoh$ACS_PCT_PERSON_INC_BELOW99_ZC))
which(is.na(dfw.zcta.sdoh$ACS_MEDIAN_HH_INC_ZC))
which(is.na(dfw.zcta.sdoh$ACS_PCT_LT_HS_ZC))
which(is.na(dfw.zcta.sdoh$ACS_PCT_UNEMPLOY_ZC))
which(is.na(dfw.zcta.sdoh$ACS_GINI_INDEX_ZC))

dfw.zcta.sdoh$ACS_PCT_PERSON_INC_BELOW99_ZC[c(4,8)] <- mean(dfw.zcta.sdoh$ACS_PCT_PERSON_INC_BELOW99_ZC,na.rm=T)
dfw.zcta.sdoh$ACS_MEDIAN_HH_INC_ZC[c(4,8)] <- mean(dfw.zcta.sdoh$ACS_MEDIAN_HH_INC_ZC,na.rm=T)
dfw.zcta.sdoh$ACS_PCT_LT_HS_ZC[8] <- mean(dfw.zcta.sdoh$ACS_PCT_LT_HS_ZC,na.rm=T)
dfw.zcta.sdoh$ACS_PCT_UNEMPLOY_ZC[8] <- mean(dfw.zcta.sdoh$ACS_PCT_UNEMPLOY_ZC,na.rm=T)

dfw.zcta.sdoh$ACS_GINI_INDEX_ZC[c(4,8)] <- mean(dfw.zcta.sdoh$ACS_GINI_INDEX_ZC,na.rm=T)

##PCA to derive socioeconomic deprivation
dfw.deprive <- -prcomp(x = cbind(dfw.zcta.sdoh$ACS_PCT_PERSON_INC_BELOW99_ZC,dfw.zcta.sdoh$ACS_PCT_LT_HS_ZC, dfw.zcta.sdoh$ACS_MEDIAN_HH_INC_ZC,dfw.zcta.sdoh$ACS_PCT_UNEMPLOY_ZC),center = T, scale=T)$x[,1]

rcorr(dfw.deprive, dfw.zcta.sdoh$ACS_MEDIAN_HH_INC_ZC,"spearman")

dfw.zcta.sdoh$deprive <- dfw.deprive

# ## Quartile
# class.num.dep <- 4
# breaks.dep <- c(-4.3,-1.1,-0.3,1.2,9.0)
# pal.dep <- brewer.pal(class.num.dep,"Reds")
# quantile(dfw.zcta.sdoh$deprive,probs = seq(0,1,1/4))
# labels.dep <- c("[-4.3,-1.1]","(-1.1, -0.3]","(-0.3,1,2]","(1.2,9.0]")

## Tertile 
class.num.dep <- 3
breaks.dep <- c(-4.3,-0.8,0.8,9.0)
pal.dep <- brewer.pal(class.num.dep,"Reds")
quantile(dfw.zcta.sdoh$deprive,probs = seq(0,1,1/3))
# labels.dep <- c("[-4.3,-0.8]","(-0.8, 0.8]","(0.8,9.0]")
labels.dep <- c("Low","Moderate","High")

quantile(dfw.zcta.sdoh$ACS_GINI_INDEX_ZC,probs = seq(0,1,1/3))
breaks.gini <- c(0.18,0.39,0.44,0.6)
pal.gini <- brewer.pal(class.num,"BuPu")

plot.dep <- plot_map(dfw.zcta.sdoh,"deprive","Deprivation",class.num.dep,breaks.dep,pal.dep,labels.dep)

plot.gini <- plot_map(dfw.zcta.sdoh,"ACS_GINI_INDEX_ZC","Gini index",class.num,breaks.gini,pal.gini,labels.dep)

tmap_save(tm = tmap_arrange(plot.dep,plot.gini,nrow = 1),filename = "DFW-sdoh-2016-20-tertile.tiff",width = 10,height= 5,units = "in", dpi = 600,compression="lzw")


# Plot the bivariate map: public transit vs. gini -------------------------

##Prepare the data
dfw.zcta.sdoh$transit.class <- cut(dfw.zcta.sdoh$transit.dis,breaks = breaks.transit)
dfw.zcta.sdoh$gini.class <- cut(dfw.zcta.sdoh$ACS_GINI_INDEX_ZC,breaks = breaks.gini)

dfw.bivariate.gini.labels <- biscale::bi_class_breaks(dfw.zcta.sdoh, x = transit.class, y = gini.class,dim = 3, dig_lab = 3, split = FALSE)
dfw.bivariate.gini.labels$bi_x[3] <- ">60"

dfw.bivariate.gini.data <- bi_class(dfw.zcta.sdoh,x="transit.class",y="gini.class",style = "quantile",dim = 3)
##Prepare the map
dfw.bivariate.gini.map <- ggplot()+
  geom_sf(data = dfw.bivariate.gini.data,mapping = aes(fill=bi_class),color="white",size=0.1,show.legend = F)+
  bi_scale_fill(pal = "DkViolet",dim=3)+
  # labs(title = "New HIV diagnosis VS distance to nearest PrEP clinic")+
  labs(title = "")+
  bi_theme()

dfw.bivariate.gini.legend <- bi_legend(pal = "DkViolet",dim = 3,xlab = "Longer pubic transit time",ylab = "Higher income inequailty",breaks=dfw.bivariate.gini.labels,size = 6)

dfw.bivariate.gini.plot <- ggdraw()+
  draw_plot(dfw.bivariate.gini.map,0,0,1,1)+
  draw_plot(dfw.bivariate.gini.legend,0,0.7,0.27,0.27)

ggsave("DFW_Bivariate_PrEP vs Gini-Aug 13.tif",plot = dfw.bivariate.gini.plot,device = "tiff",dpi = 600,width = 7,height = 9)


# Plot the bivariate map: public transit vs. socioeconomic deprivation -------------------------

##Prepare the data
# dfw.zcta.sdoh$transit.class <- cut(dfw.zcta.sdoh$transit.dis,breaks = breaks.transit)
dfw.zcta.sdoh$dep.class <- cut(dfw.zcta.sdoh$deprive,breaks = breaks.dep)

dfw.bivariate.dep.labels <- biscale::bi_class_breaks(dfw.zcta.sdoh, x = transit.class, y = dep.class,dim = 3, dig_lab = 3, split = FALSE)
dfw.bivariate.dep.labels$bi_x[3] <- ">60"
dfw.bivariate.dep.labels$bi_y <- c("Low","Moderate","High")

dfw.bivariate.dep.data <- bi_class(dfw.zcta.sdoh,x="transit.class",y="dep.class",style = "quantile",dim = 3)
##Prepare the map
dfw.bivariate.dep.map <- ggplot()+
  geom_sf(data = dfw.bivariate.dep.data,mapping = aes(fill=bi_class),color="white",size=0.1,show.legend = F)+
  bi_scale_fill(pal = "DkBlue",dim=3)+
  # labs(title = "New HIV diagnosis VS distance to nearest PrEP clinic")+
  labs(title = "")+
  bi_theme()

dfw.bivariate.dep.legend <- bi_legend(pal = "DkBlue",dim = 3,xlab = "Longer pubic transit time",ylab = "Higher deprivation",breaks=dfw.bivariate.dep.labels,size = 6)

dfw.bivariate.dep.plot <- ggdraw()+
  draw_plot(dfw.bivariate.dep.map,0,0,1,1)+
  draw_plot(dfw.bivariate.dep.legend,0,0.67,0.27,0.27)

ggsave("DFW_Bivariate_PrEP vs dep-Aug 13.png",plot = dfw.bivariate.dep.plot,width = 5,height = 4,dpi = 600) ##,width = 7,height = 9

tmap_save(tm = tmap_arrange(plot.transit,plot.car,dfw.bivariate.plot,nrow = 1),filename = "PrEP accessibility-prep-by mode.png",width = 15,height= 4,units = "in", dpi = 600,compression="lzw")

# Map the new HIV diagnosis rate ------------------------------------------

dallas.hiv21 <- read.xlsx(paste0(here(),"/Data/Dallas_AIDSVu_DownloadableDataset_2021.xlsx"),sheetIndex = 1,startRow = 4)
dallas.hiv21 <- dallas.hiv21[-dim(dallas.hiv21)[1],]

fw.hiv21 <- read.xlsx(paste0(here(),"/Data/FortWorth_AIDSVu_DownloadableDataset_2021.xlsx"),sheetIndex = 1,startRow = 4)
fw.hiv21 <- fw.hiv21[-dim(fw.hiv21)[1],]

dfw.hiv21 <- rbind(dallas.hiv21,fw.hiv21)
dfw.hiv21$X5.Year.Cumulative.ZIP.Code.Risk[87:88] <- 0 

quantile(dfw.hiv21$X5.Year.Cumulative.ZIP.Code.Risk,probs = seq(0,1,1/4))
quantile(dfw.hiv21$X5.Year.Cumulative.ZIP.Code.Risk,probs = seq(0,1,1/3))

dfw.zcta.hiv <- left_join(dfw.zcta.sdoh,dfw.hiv21,by=c("GEOID"="Zip.Code"))
dfw.zcta.full <- left_join(dfw.zcta.sdoh,dfw.hiv21,by=c("GEOID"="Zip.Code"))

# dfw.zcta.join$X5.Year.Cumulative.ZIP.Code.Risk[which(dfw.zcta.join$X5.Year.Cumulative.ZIP.Code.Risk<0)] <- NA

##zcta with suppressed HIV data
zcta.missing <- which(dfw.zcta.hiv$X5.Year.Cumulative.ZIP.Code.Risk<0)
length(zcta.missing) ## 12
dfw.zcta.hiv <- dfw.zcta.hiv[-zcta.missing,]

quantile(dfw.zcta.hiv$X5.Year.Cumulative.ZIP.Code.Risk,probs = seq(0,1,1/3))

breaks.hiv <- c(-1,104,182,1011)
labels.hiv <- c("[0-104]","(104-181]","(181,1010]")

# Bivariate map: PrEP vs. HIV ---------------------------------------------

##Prepare the data
# dfw.zcta.sdoh$transit.class <- cut(dfw.zcta.sdoh$transit.dis,breaks = breaks.transit)
dfw.zcta.hiv$hiv.class <- cut(dfw.zcta.hiv$X5.Year.Cumulative.ZIP.Code.Risk,breaks = breaks.hiv)

dfw.bivariate.hiv.labels <- biscale::bi_class_breaks(dfw.zcta.hiv, x = transit.class, y = hiv.class,dim = 3, dig_lab = 3, split = FALSE)
dfw.bivariate.hiv.labels$bi_x[3] <- ">60"
dfw.bivariate.hiv.labels$bi_y <- c("0-104","105-181","182-1010")

dfw.bivariate.hiv.data <- bi_class(dfw.zcta.hiv,x="transit.class",y="hiv.class",style = "quantile",dim = 3)
##Prepare the map
dfw.bivariate.hiv.map <- ggplot()+
  geom_sf(data = dfw.bivariate.hiv.data,mapping = aes(fill=bi_class),color="white",size=0.1,show.legend = F)+
  bi_scale_fill(pal = "Bluegill",dim=3)+
  # labs(title = "New HIV diagnosis VS distance to nearest PrEP clinic")+
  labs(title = "")+
  bi_theme()

dfw.bivariate.hiv.legend <- bi_legend(pal = "Bluegill",dim = 3,xlab = "Longer pubic transit time",ylab = "Higher new HIV diagnosis",breaks=dfw.bivariate.hiv.labels,size = 6)

dfw.bivariate.hiv.plot <- ggdraw()+
  draw_plot(dfw.bivariate.hiv.map,0,0,1,1)+
  draw_plot(dfw.bivariate.hiv.legend,0,0.67,0.27,0.27)

ggsave("DFW_Bivariate_PrEP vs HIV-Aug 13.png",plot = dfw.bivariate.hiv.plot,width=5,height=4,dpi = 600) ##,width = 7,height = 9


# Merge bivariate maps -------------------------------

## HIV and socioeconomic deprivation
dfw.merge <- plot_grid(dfw.bivariate.hiv.plot, dfw.bivariate.dep.plot, ncol = 2, nrow = 1)

ggsave("DFW_Bivariate_PrEP vs HIV & dep-Aug 14.png",plot = dfw.merge,width = 10,height = 4,dpi = 600)

## Public transit- and driving-based PrEP accessibility
dfw.prep.all <- tmap_grob(tmap_arrange(plot.transit,plot.car,nrow = 1))
dfw.prep.merge <- plot_grid(dfw.prep.all,dfw.bivariate.plot, nrow = 1)
ggsave("DFW_prep-all-Aug 14.png",plot = dfw.prep.merge,width = 15,height = 4,dpi = 600)

which(dfw.zcta.hiv$X5.Year.Cumulative.ZIP.Code.Risk>181 & dfw.zcta.hiv$transit.dis>60 & dfw.zcta.hiv$deprive>0.8 & dfw.zcta.hiv$car.dis<15) ## 55  97 114 119 142

which(dfw.zcta.hiv$X5.Year.Cumulative.ZIP.Code.Risk>181 & dfw.zcta.hiv$transit.dis>60 & dfw.zcta.hiv$deprive>0.8) ## 55 97 104 114 119 142

dfw.zcta.hiv$car.dis[c(5,97,104,114,119,142)]

dfw.zcta.hiv$GEOID[c(55,97,104,114,119,142)]  ## zip 75241 75236 75180 76112 75240 75217

length(which(dfw.zcta.hiv$car.dis<15))

# Map the dot population density by race in DFW ---------------------------
capi <- "fef9c7528cd8685bed392e03c1e6ab78dcab41a5"  ## My Census Data API key
census_api_key(capi,install = T,overwrite = T)
readRenviron("~/.Renviron")

dec.variables <- load_variables(2020,"pl")
acs.variabels <- load_variables(2020,"acs5")

dfw.race <- get_decennial(
  geography = "tract",
  state = "TX",
  county = c("Tarrant","Dallas"),
  variables = c(
    Hispanic = "P2_002N",
    White = "P2_005N",
    Black = "P2_006N",
    Native = "P2_007N",
    Asian = "P2_008N"
  ),
  summary_var = "P2_001N",
  year = 2020,
  geometry = TRUE
) %>%
  mutate(percent = 100 * (value / summary_value))

background.tracts <- filter(dfw.race, variable == "White")

dfw.dots <- dfw.race %>%
  as_dot_density(
    value = "value",
    values_per_dot = 100,
    group = "variable"
  )

tm_shape(background.tracts) + 
  tm_polygons(col = "white", 
              border.col = "grey") + 
  tm_layout(
    #title = "Race/Ethnicity population, 2020 US Census",
    frame = F)+
  tm_shape(dfw.dots) +
  tm_dots(col = "variable", 
          palette = "Set1",
          size = 0.01, 
          title = "1 dot = 100 people") + 
  tm_layout(legend.outside = T,title = "Population by race/ethnicity,\n2020 US Census",title.size = 1.1)


# write out the final dfw zcta file ---------------------------------------
dfw.zcta.output <- dfw.zcta.full[,1:21]
dfw.zcta.output$risk <- dfw.zcta.full$X5.Year.Cumulative.ZIP.Code.Risk
dfw.zcta.output$cases <- dfw.zcta.full$X5.Year.Cumulative.ZIP.Code.Cases

dfw.zcta.output <- left_join(dfw.zcta.output,SDOH.sel2,"GEOID")

which(is.na(dfw.zcta.output$ACS_PCT_ASIAN_ZC))
which(is.na(dfw.zcta.output$ACS_PCT_BLACK_ZC))
which(is.na(dfw.zcta.output$ACS_PCT_HISPANIC_ZC))
which(is.na(dfw.zcta.output$ACS_PCT_UNINSURED_ZC))
which(is.na(dfw.zcta.output$CEN_POPDENSITY_ZC))
which(is.na(dfw.zcta.output$ACS_PCT_WHITE_ZC))

dfw.zcta.output$ACS_PCT_ASIAN_ZC[8] <- mean(dfw.zcta.output$ACS_PCT_ASIAN_ZC,na.rm=T)
dfw.zcta.output$ACS_PCT_BLACK_ZC[8] <- mean(dfw.zcta.output$ACS_PCT_BLACK_ZC,na.rm=T)
dfw.zcta.output$ACS_PCT_HISPANIC_ZC[8] <- mean(dfw.zcta.output$ACS_PCT_HISPANIC_ZC,na.rm=T)
dfw.zcta.output$ACS_PCT_UNINSURED_ZC[8] <- mean(dfw.zcta.output$ACS_PCT_UNINSURED_ZC,na.rm=T)
# dfw.zcta.output$ACS_PCT_WHITE_ZC[8] <- mean(dfw.zcta.output$ACS_PCT_WHITE_ZC,na.rm=T)

## Remove polygons with suppressed HIV data
dfw.zcta.output2 <- dfw.zcta.output[-which(dfw.zcta.output$cases<0),]

# st_write(dfw.zcta.output, "dfw-zcta21-output-Aug 14.shp")

# Plot vehicle owership ---------------------------------------------------
dfw.output <- st_read(paste0(here(),"/dfw-zcta21-output.shp"))
car.owner20 <- read.csv(paste0(here(),"/Data/Vehicle owership-2020-zcta-Texas.csv"))
dfw.car <- left_join(dfw.output,car.owner20,by = c("GEOID"="Geography.Name"))

quantile(car.owner20$Percent.Housing.Units.with.0.Vehicles.Available,probs = seq(0,1,by=1/5),na.rm = T)
mean(dfw.car$Percent.Housing.Units.with.0.Vehicles.Available,na.rm=T) ## 6.03%
sd(dfw.car$Percent.Housing.Units.with.0.Vehicles.Available,na.rm = T) ## 6.27%

class.vehicle <- 5
breaks.vehicle <- c(-1,0.79,2.74,4.55,7.48,100)
pal.vechile <- brewer.pal(class.vehicle,"YlGnBu")
labels.vehicle <- c("0-0.79","0.80-2.74","2.75-4.55","4.56-7.48","7.49-100")

## % of household without vehicles
plot.vehicle <- plot_map(dfw.car,"Percent.Housing.Units.with.0.Vehicles.Available","% no vehicle ownership",class.vehicle,breaks.vehicle,pal.vechile,labels.vehicle)+
  tm_layout(frame = F,legend.title.size = 0.5,legend.text.size = 0.5,legend.position = c(0,0))

tmap_save(plot.vehicle,filename = "Vehicle ownership-quntile.png",width = 5,height= 4, dpi = 600)

dfw.car$Percent.Housing.Units.with.0.Vehicles.Available[which(dfw.car$GEOID %in% c(76112,75236,75241,75217,75180,75240))]

dfw.car$Percent.Housing.Units.with.0.Vehicles.Available[which(dfw.car$GEOID ==76112)]
dfw.car$Percent.Housing.Units.with.0.Vehicles.Available[which(dfw.car$GEOID ==75236)]
dfw.car$Percent.Housing.Units.with.0.Vehicles.Available[which(dfw.car$GEOID ==75241)]
dfw.car$Percent.Housing.Units.with.0.Vehicles.Available[which(dfw.car$GEOID ==75217)]
dfw.car$Percent.Housing.Units.with.0.Vehicles.Available[which(dfw.car$GEOID ==75180)]
dfw.car$Percent.Housing.Units.with.0.Vehicles.Available[which(dfw.car$GEOID ==75240)]
