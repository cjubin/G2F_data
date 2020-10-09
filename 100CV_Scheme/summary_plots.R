## LIGHTGBM results##
library(stringr, lib.loc = "/home/uni08/jubin1/R/x86_64-redhat-linux-gnu-library/3.6")
library(ggplot2)
library(cowplot)

## Summary plots Random splits with new environments, training set 80:20 ##
setwd('/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/ML_PREDICTIONS/Heslot_scheme/lightgbm/results_index_20')
j=20
university_factor=read.table('/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/WEATHER_PROCESSING/Env_data_processing/factoruniversity.txt',header = T)


## Comparison for model G+WC+SC+Lon+Lat vs G+Y+Lon+Lat




# Read output results for G+Y+Lon+Lat

results_model_1 <- list()
for (i in 1:j) {
  results_model_1[i]<-readRDS(paste0('final_list_results_index',i,'geno_info_PCs_G+Y+Lon+Lat_yld_bu_ac_tuning_method_random-CV.RDS'))['te_corr_by_yearexp']
  
}


results_model_1final_df<-do.call(rbind, results_model_1)








# Read output results for G+WC+SC+Lon+Lat
# Need to select the best model obtained from feature selection for each split train/test

results_model_2final<- list()
winner_model <- vector()
winner_model_nb <-vector()
results_model_2a <- list()
results_model_2b <- list()
results_model_2c <- list()
results_model_2d <- list()

for (i in 1:j) {
  
  average_corr_noFS <- readRDS(paste0('final_list_results_index',i,'geno_info_PCs_G+WC+SC+Lon+Lat_yld_bu_ac_tuning_method_random-CV.RDS'))['te_mean_corr']
  average_corr_FS1 <- readRDS(paste0('FS1_final_list_results_index',i,'geno_info_PCs_G+WC+SC+Lon+Lat_yld_bu_ac_tuning_method_random-CV.RDS'))['te_mean_corr']
  average_corr_FS2 <- readRDS(paste0('FS2_final_list_results_index',i,'geno_info_PCs_G+WC+SC+Lon+Lat_yld_bu_ac_tuning_method_random-CV.RDS'))['te_mean_corr']
  average_corr_FS3 <- readRDS(paste0('FS3_final_list_results_index',i,'_geno_info_PCs_G+WC+SC+Lon+Lat_yld_bu_ac_tuning_method_random-CV.RDS'))['te_mean_corr']
  
  winner_model[i] <- c('average_corr_noFS','average_corr_FS1','average_corr_FS2','average_corr_FS3')[which.max(c(average_corr_noFS,average_corr_FS1,average_corr_FS2,average_corr_FS3))]
  
  print(winner_model[i])
  winner_model_nb[i] <- which.max(c(average_corr_noFS,average_corr_FS1,average_corr_FS2,average_corr_FS3))
  print(c(average_corr_noFS,average_corr_FS1,average_corr_FS2,average_corr_FS3)[which.max(c(average_corr_noFS,average_corr_FS1,average_corr_FS2,average_corr_FS3))])
  
  
}

table(winner_model)






## OPTION 1 : for each train/test split, the best model after feature selection is selected in each case

## We see that the model with the least included number of covariates (top 10 after recursive feature elimination) is the best most of the time.
## Therefore, we consider for further plots only the results with the model G+WC+SC+Lon+Lat including only top 10 WC+SC variables

for (i in 1:j) {
  
  results_model_2final[[i]] <- as.data.frame(readRDS(paste0('FS3_final_list_results_index',i,'_geno_info_PCs_G+WC+SC+Lon+Lat_yld_bu_ac_tuning_method_random-CV.RDS'))['te_corr_by_yearexp'])
  results_model_2final[[i]][,2] <- as.numeric(results_model_2final[[i]][,2])
}

results_model_2final_df<-do.call(rbind, results_model_2final)
colnames(results_model_2final_df)<-c('Year_Exp','COR')


## Join the two datasets 
all(results_model_1final_df$Year_Exp==results_model_2final_df$Year_Exp)
df<-cbind(results_model_1final_df,results_model_2final_df[,2])
colnames(df)[c(2,3)]<-c('x','y')
df$year<-str_sub(df$Year_Exp,1,4)
df$university<-university_factor[match(df$Year_Exp,university_factor$Year_Exp),'UnivFactor']


#Create a custom color scale
library(RColorBrewer)
library(ggplot2)
myColors <- colorRampPalette(brewer.pal(8, "Set2"))(length(unique(df$university)))
names(myColors) <- unique(as.character(df$university))
colScale <- scale_colour_manual(name = "University locations",values = myColors)

#pdf(file="option1_modif_G+WC+SC+Lon+Lat_VS_G+Y+Lon+Lat.pdf",width=17,height=9)
p <- ggplot(df,aes(x,y,colour = university,shape=year)) + geom_point()
p1 <- p + colScale + geom_abline(intercept = 0,slope = 1) + labs(x='Environment prediction accuracy, Model G+Y+Lon+Lat',y='Environment prediction accuracy, Model G+WC+SC+Lon+Lat') 
print(p1)
#dev.off()
save_plot(p1,filename='option1_plot_G+WC+SC+Lon+Lat_VS_G+Y+Lon+Lat.pdf',base_height = 9,base_width = 17)
















## OPTION 2


for (i in 1:j) {
  corr_noFS <- readRDS(paste0('final_list_results_index',i,'geno_info_PCs_G+WC+SC+Lon+Lat_yld_bu_ac_tuning_method_random-CV.RDS'))['te_corr_by_yearexp']
  corr_FS1 <- readRDS(paste0('FS1_final_list_results_index',i,'geno_info_PCs_G+WC+SC+Lon+Lat_yld_bu_ac_tuning_method_random-CV.RDS'))['te_corr_by_yearexp']
  corr_FS2 <- readRDS(paste0('FS2_final_list_results_index',i,'geno_info_PCs_G+WC+SC+Lon+Lat_yld_bu_ac_tuning_method_random-CV.RDS'))['te_corr_by_yearexp']
  corr_FS3 <- readRDS(paste0('FS3_final_list_results_index',i,'_geno_info_PCs_G+WC+SC+Lon+Lat_yld_bu_ac_tuning_method_random-CV.RDS'))['te_corr_by_yearexp']
  
  
  results_model_2final[[i]] <- as.data.frame(list(corr_noFS,corr_FS1,corr_FS2,corr_FS3)[winner_model_nb[i]])
  results_model_2final[[i]][,2] <- as.numeric(results_model_2final[[i]][,2])
}

results_model_2final_df<-do.call(rbind, results_model_2final)
colnames(results_model_2final_df)<-c('Year_Exp','COR')


## Join the two datasets 
all(results_model_1final_df$Year_Exp==results_model_2final_df$Year_Exp)
df<-cbind(results_model_1final_df,results_model_2final_df[,2])
colnames(df)[c(2,3)]<-c('x','y')
df$year<-str_sub(df$Year_Exp,1,4)
df$university<-university_factor[match(df$Year_Exp,university_factor$Year_Exp),'UnivFactor']


#Create a custom color scale
library(RColorBrewer)
library(ggplot2)
myColors <- colorRampPalette(brewer.pal(8, "Set2"))(length(unique(df$university)))
names(myColors) <- unique(as.character(df$university))
colScale <- scale_colour_manual(name = "University locations",values = myColors)

pdf(file="option2_modif_G+WC+SC+Lon+Lat_VS_G+Y+Lon+Lat.pdf",width=17,height=9)
p <- ggplot(df,aes(x,y,colour = university,shape=year)) + geom_point()
p1 <- p + colScale + geom_abline(intercept = 0,slope = 1) + labs(x='Environment prediction accuracy, Model G+Y+Lon+Lat',y='Environment prediction accuracy, Model G+WC+SC+Lon+Lat') 
print(p1)
dev.off()
ggsave(p1,filename='option2_plot_G+WC+SC+Lon+Lat_VS_G+Y+Lon+Lat.pdf',height = 9,width = 17,device = 'pdf')

