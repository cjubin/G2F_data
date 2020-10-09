## Top 30 variables importance for LVY CV scheme, and for each machine learning model ##
library(reshape2)
library(pacman)
library(ggbump)
library(dplyr)
library(viridis)
library(ggplot2)
library('cowplot')

# First method lightgbm
##Compare the robustness across years
features_imp_lightgbm<-list()
n=1
for (i in c(2014:2017)) {
  setwd(paste0("/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/ML_PREDICTIONS/LYOUT/lightgbm/",i))
  features_imp_lightgbm[[n]]<-as.data.frame(readRDS('useful_results_PCspredictors_G+WC+SC+Y+Lon+Lat_yld_bu_ac_tuning_method_random-CV_bootstrap_strategy_total_random.RDS')['varimps_average_bootstraps']) %>% top_n(40,varimps_average_bootstraps.average_gain)
  
  features_imp_lightgbm[[n]][,'year']<-i
  features_imp_lightgbm[[n]][,'rank']<-c(1:40)
  n<-n+1
}

df<-do.call(rbind,features_imp_lightgbm)
colnames(df)[1]<-'Feature'
## Plot variable ranking ##


p<-ggplot(df, aes(year, rank, color = Feature)) +
  geom_point(size = 7) +
  geom_text(data = df %>% filter(year == min(year)),
            aes(x = year - .1, label = Feature), size = 5, hjust = 1) +
  geom_text(data = df %>% filter(year == max(year)),
            aes(x = year + .1, label = Feature), size = 5, hjust = 0) +
  geom_bump(size = 2, smooth = 8) +
  scale_x_continuous(limits = c(2010.6, 2013.4),
                     breaks = seq(2011, 2013, 1)) +
  theme_minimal_grid(font_size = 14, line_size = 0) +
  theme(legend.position = "none",
        panel.grid.major = element_blank()) +
  labs(y = "RANK",
       x = NULL) +
  scale_y_reverse() +
  scale_color_manual(values = viridis(40))
plot(p)
ggsave(p,'variable_importances_LVY.pdf',width = 10,height = 10,device = 'pdf')

## Second option by including
