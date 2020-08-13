library(dtw)
library(purrr)
library(dtwclust)

weather = read.table(
  "/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/replaced_daily_weather.txt",
  header = T,
  sep = "\t"
)
weather$GS = weather$Date.Harvested - weather$Date.Planted
weather = weather[, c("Year_Exp",
                      "PRCP",
                      "TMIN",
                      "TMAX",
                      "HMEAN",
                      "solar_radiation_NASA")]

ds <-
  split(weather[,-which(colnames(weather) == "Year_Exp")], weather$Year_Exp)
for (i in 1:length(ds)) {
  name=names(ds[[i]])
  ind <- which(is.na(ds[[i]])[, "solar_radiation_NASA"])
  if (length(ind) != 0) {
    for (s in ind) {
      ds[[i]][s, "solar_radiation_NASA"] <-
        ds[[i]][s + 1, "solar_radiation_NASA"]
    }
  }
  if (length(which(is.na(ds[[i]]))) != 0) {
    print(i)
    ds[[i]] <- ds[[i]][complete.cases(ds[[i]]),]
  }
  ds[[i]]<-as.matrix(ds[[i]])
  names(ds[[i]])<-name
}




hc <- dtwclust::tsclust(
  ds,
  type = "h",
  k = 10,
  preproc = zscore,
  seed = 899,
  distance = "dtw_basic"
)
hcd <- as.dendrogram(hc)
plot(hcd, type = "rectangle", ylab = "Height")

cfgs <-
  compare_clusterings_configs(
    types = c("p", "h"),
    k = c(5:12),
    distances = pdc_configs("distance", dtw_basic = list()),
    centroids = pdc_configs("centroid", default = list()),
    controls = list(
      partitional = partitional_control(iter.max = 30L,
                                        nrep = 1L),
    hierarchical = hierarchical_control(method = "all"),
    preprocs = pdc_configs("zscore")
      
    )
  )

res <-
  compare_clusterings(ds,types = c("h","p"),
    configs = cfgs,
    #score.clus = score_fun,
    #pick.clus = pick_fun,
    return.objects = TRUE
  )


























