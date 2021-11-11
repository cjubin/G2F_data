args = commandArgs(trailingOnly = TRUE)

x1 = args[1]
x2 = as.numeric(args[2]) + 2013

params <-
  list(
    seed = 105L,
    geno_info = "PCs",
    phenos_file = "phenos_with_EC_covariates.txt",
    geno_file = "geno_hybrids.txt",
    nbootstrap = 50L,
    year_to_predict = as.numeric(x2),
    trait = "yld_bu_ac",
    sets_predictors = as.character(x1),
    method = "xgboost",
    tuning_hyper_parameters = "random-CV",
    bootstrap_sampling = "total_random",
    WC_features_removed = TRUE,
    `recipes::step_feature_selection` = TRUE
  )


## ----echo=FALSE---------------------------------------------------------------
start.time <- Sys.time()
print(params)


## ---- include=FALSE,echo=FALSE------------------------------------------------


source(
  'setSeeds.R'
)
`%notin%` <- Negate(`%in%`)




## ----include=FALSE,echo=FALSE-------------------------------------------------
#library(furrr)
#library(rmarkdown)
#library(checkmate)
library(data.table)
library(tidymodels)
library(vip)
#library(tidyverse)
#library(SHAPforxgboost)
library(doFuture)
library(gridExtra)
#library(future)
library(xgboost)
library(doMC)
library(doParallel)


## ---- include=FALSE,echo=FALSE------------------------------------------------
phenos = read.table(params$phenos_file, header = T, sep = '\t')
phenos$year = as.factor(phenos$year)
phenos$Year_Exp = as.factor(phenos$Year_Exp)


training = phenos[which(phenos$year %notin% params$year_to_predict), ]
test = phenos[which(phenos$year %in% params$year_to_predict), ]
training <- select(training, -'UniversityFactor')
test <- select(test, -'UniversityFactor')



print('The number of training observations is:')
print(nrow(training))
print('The number of test observations to predict is:')
print(nrow(test))
print('The number of Year_Exp to predict is:')
print(length(unique(test$Year_Exp)))



#test$P.Total=test$P.V+test$P.F+test$P.G
#training$P.Total=training$P.V+training$P.F+training$P.G

training <-
  training[, -which(
    colnames(training) %in% c(
      "parent1",
      "parent2",
      "parent1.GBS.sample",
      "parent2.GBS.sample"
    )
  )]
test <-
  test[, -which(
    colnames(test) %in% c(
      "parent1",
      "parent2",
      "parent1.GBS.sample",
      "parent2.GBS.sample"
    )
  )]




cores <- as.integer(Sys.getenv('SLURM_CPUS_PER_TASK'))

#registerDoParallel(cores)

registerDoFuture()
cl <- makeCluster(cores)
plan(cluster, workers = cl)



if (params$geno_info == 'PCs') {
  print('Genomic information reduced to PCs extracted from SNPs genotype matrix.')
  
  geno_hybrids = fread(params$geno_file)
  geno_hybrids = as.data.frame(geno_hybrids)
  
  
  geno_hybrids <-
    geno_hybrids[,-which(
      colnames(geno_hybrids) %in% c(
        "parent1",
        "parent2",
        "parent1.GBS.sample",
        "parent2.GBS.sample"
      )
    )]
  
  print('Data read')
  
  colnames(geno_hybrids)[2:ncol(geno_hybrids)] <-
    paste0('SNP', colnames(geno_hybrids)[2:ncol(geno_hybrids)])
  
  markers_tr = geno_hybrids[geno_hybrids$pedigree %in% training$pedigree, ]
  markers_te = geno_hybrids[geno_hybrids$pedigree %in% test$pedigree &
                              geno_hybrids$pedigree %notin% training$pedigree, ]
  rownames(markers_te) <- NULL
  
  rec1 <- recipes::recipe(pedigree ~ . ,
                          data = markers_tr) %>%
    recipes::step_pca(starts_with('SNP'),
                      num_comp = 275,
                      options = list(center = T, scale. = T))
  
  norm_obj <- prep(rec1, training = markers_tr)
  
  te <- bake(norm_obj, markers_te)
  te$pedigree = markers_te$pedigree
  tr <- juice(norm_obj)
  pc_values <- rbind(tr, te)
  
  
  training = merge(training, pc_values, by = 'pedigree', all.x = T)
  test = merge(test, pc_values, by = 'pedigree', all.x = T)
  
} else if (params$geno_info == 'G') {
  
    geno_hybrids = fread(params$geno_file)
    geno_hybrids = as.data.frame(geno_hybrids)
    
    geno_hybrids <-
      geno_hybrids[,-which(
        colnames(geno_hybrids) %in% c(
          "parent1",
          "parent2",
          "parent1.GBS.sample",
          "parent2.GBS.sample"
        )
      )]
    
    print('Data read')
    
    colnames(geno_hybrids)[2:ncol(geno_hybrids)] <-
      paste0('SNP', colnames(geno_hybrids)[2:ncol(geno_hybrids)])
    
    
    G = A.mat(geno_hybrids[, 2:ncol(geno_hybrids)])
    row.names(G) = as.character(geno_hybrids[, 1])
    colnames(G) = as.character(geno_hybrids[, 1])
    
      
} else{
  geno_hybrids = fread(params$geno_file)
  geno_hybrids = as.data.frame(geno_hybrids)
  
  geno_hybrids <-
    geno_hybrids[,-which(
      colnames(geno_hybrids) %in% c(
        "parent1",
        "parent2",
        "parent1.GBS.sample",
        "parent2.GBS.sample"
      )
    )]
  
  print('Data read')
  
  colnames(geno_hybrids)[2:ncol(geno_hybrids)] <-
    paste0('SNP', colnames(geno_hybrids)[2:ncol(geno_hybrids)])
  
  
  training = merge(training, geno_hybrids, by = 'pedigree', all.x = T)
  test = merge(test, geno_hybrids, by = 'pedigree', all.x = T)
}

print('PCA achieved')


## Conversion to factors
#training$year = as.factor(as.vector(training$year))
#training$counties = as.factor(as.vector(training$counties))
#test$year = as.factor(as.vector(test$year))
#test$counties = as.factor(as.vector(test$counties))

## Conversion to numeric
if (!is.numeric(training$Latitude) |
    !is.numeric(training$Longitude) |
    !is.numeric(test$Latitude) | !is.numeric(test$Latitude)) {
  training$Latitude = as.numeric(as.vector(training$Latitude))
  training$Latitude = as.numeric(as.vector(training$Longitude))
  test$Latitude = as.numeric(as.vector(test$Latitude))
  test$Longitude = as.numeric(as.vector(test$Longitude))
}

## Predictors variables included for prediction according to the sets of predictors defined by the user.
toMatch = c('.Total', '.V', '.F', '.G', 'length.growing.season')
if (params$geno_info == 'PCs')
{
  toMatch2 = 'PC'
} else if (params$geno_info == 'PLS') {
  toMatch2 = 'PLS'
} else if (params$geno_info == 'G') {
  toMatch2 = 'pedigree'
} else{
  toMatch2 = 'SNP'
}

toMatch3 = c('year')
toMatch4 = c('counties')
toMatch5 = c('.SC')
toMatch6 = c('Year_Exp')
toMatch7 = c('Longitude')
toMatch8 = c('Latitude')
toMatch9 = c('P.V', 'P.F', 'P.G')
toMatch10 = c(
  'MeanT.V',
  'MeanT.F',
  'MeanT.G',
  'MinT.V',
  'MinT.F',
  'MinT.G',
  'MaxT.V',
  'MaxT.F',
  'MaxT.G',
  'GDD.V',
  'GDD.F',
  'GDD.G',
  'FreqMaxT30.V',
  'FreqMaxT30.F',
  'FreqMaxT30.G',
  'FreqMaxT35.V',
  'FreqMaxT35.F',
  'FreqMaxT35.G',
  'Photothermal.time.Tot.V',
  'Photothermal.time.Tot.F',
  'Photothermal.time.Tot.G'
)

matches1 <-
  grep(paste(toMatch, collapse = "|"), colnames(training), value = TRUE)
matches2 <-
  grep(paste(toMatch2, collapse = "|"), colnames(training), value = TRUE)
matches3 <-
  grep(paste(toMatch3, collapse = "|"), colnames(training), value = TRUE)
matches4 <-
  grep(paste(toMatch4, collapse = "|"), colnames(training), value = TRUE)
matches5 <-
  grep(paste(toMatch5, collapse = "|"), colnames(training), value = TRUE)
matches6 <-
  grep(paste(toMatch6, collapse = "|"), colnames(training), value = TRUE)
matches7 <-
  grep(paste(toMatch7, collapse = "|"), colnames(training), value = TRUE)
matches8 <-
  grep(paste(toMatch8, collapse = "|"), colnames(training), value = TRUE)
matches9 <-
  grep(paste(toMatch9, collapse = "|"), colnames(training), value = TRUE)
matches10 <-
  grep(paste(toMatch10, collapse = "|"), colnames(training), value = TRUE)


if (params$sets_predictors == 'WC+SC') {
  predictors = c(matches1, matches5, matches3, matches6)
}
if (params$sets_predictors == 'G') {
  predictors = c(matches2, matches3, matches6)
}
if (params$sets_predictors == 'Y+L') {
  predictors = c(matches3, matches4, matches6)
}
if (params$sets_predictors == 'G+WC') {
  predictors = c(matches1, matches2, matches3, matches6)
}
if (params$sets_predictors == 'G+SC') {
  predictors = c(matches2, matches5, matches3, matches6)
}
if (params$sets_predictors == 'WC+SC+Y+L') {
  predictors = c(matches1, matches3, matches4, matches5, matches6)
}
if (params$sets_predictors == 'G+WC+SC') {
  predictors = c(matches1, matches2, matches5, matches3, matches6)
}
if (params$sets_predictors == 'G+Y+L') {
  predictors = c(matches2, matches3, matches4, matches6)
}
if (params$sets_predictors == 'G+SC+Y+Lon+Lat') {
  predictors = c(matches5, matches2, matches3, matches6, matches7, matches8)
}
if (params$sets_predictors == 'G+WC+SC+Y+L') {
  predictors = c(matches1, matches2, matches3, matches4, matches5, matches6)
}
if (params$sets_predictors == 'G+Y') {
  predictors = c(matches2, matches3, matches6)
}
if (params$sets_predictors == 'G+L') {
  predictors = c(matches2, matches4, matches3, matches6)
}
if (params$sets_predictors == 'G+WC+SC+Lon+Lat') {
  predictors = c(matches1,
                 matches2,
                 matches3,
                 matches5,
                 matches6,
                 matches7,
                 matches8)
}
if (params$sets_predictors == 'G+WC+Lon+Lat') {
  predictors = c(matches1, matches2, matches3, matches6, matches7, matches8)
}
if (params$sets_predictors == 'G+SC+Lon+Lat') {
  predictors = c(matches5, matches2, matches3, matches6, matches7, matches8)
}
if (params$sets_predictors == 'G+Y+Lon+Lat+P') {
  predictors = c(matches2, matches3, matches6, matches7, matches8, matches9)
}
if (params$sets_predictors == 'G+Y+Lon+Lat+T') {
  predictors = c(matches2, matches3, matches6, matches7, matches8, matches10)
}
if (params$sets_predictors == 'G+Y+Lon+Lat+T+P') {
  predictors = c(matches2,
                 matches3,
                 matches6,
                 matches7,
                 matches8,
                 matches9,
                 matches10)
}
if (params$sets_predictors == 'G+Lon+Lat+P') {
  predictors = c(matches2, matches3, matches6, matches7, matches8, matches9)
}
if (params$sets_predictors == 'G+Lon+Lat+T') {
  predictors = c(matches2, matches3, matches6, matches7, matches8, matches10)
}
if (params$sets_predictors == 'G+Lon+Lat+T+P') {
  predictors = c(matches2,
                 matches3,
                 matches6,
                 matches7,
                 matches8,
                 matches9,
                 matches10)
}
if (params$sets_predictors == 'G+Y+Lon+Lat') {
  predictors = c(matches2, matches3, matches6, matches7, matches8)
}
if (params$sets_predictors == 'G+WC+SC+Y+Lon+Lat') {
  predictors = c(matches1,
                 matches2,
                 matches3,
                 matches5,
                 matches6,
                 matches7,
                 matches8)
}
if (params$sets_predictors == 'G+WC+Y+Lon+Lat') {
  predictors = c(matches1, matches2, matches3, matches6, matches7, matches8)
}



#Selecting the variables based on the list of predictors defined previously
#Create dummy variables if factor variables are present
training = training[, colnames(training) %in% c(predictors, params$trait)]
test = test[, colnames(test) %in% c(predictors, params$trait)]



#Excluding some weather covariates if WC_feature_selection TRUE
if (params$WC_features_removed == TRUE &
    params$sets_predictors %in% c(
      'WC+SC',
      'WC+SC+Y+L',
      'G+WC+SC',
      'G+WC+SC+Y+L',
      'G+WC',
      'G+WC+SC+Lon+Lat',
      'G+WC+Lon+Lat',
      'G+WC+SC+Y+Lon+Lat',
      'G+WC+Y+Lon+Lat'
    )) {
  to_exclude = c(
    "Mean.Wind.V",
    "Mean.Wind.F",
    "Mean.Wind.G",
    "MaxWind.V",
    "MaxWind.F",
    "MaxWind.G",
    "FreqDaysWindSpeedSup12.V",
    "FreqDaysWindSpeedSup12.F",
    "FreqDaysWindSpeedSup12.G",
    "SumDiffP_ETP.V",
    "SumDiffP_ETP.F",
    "SumDiffP_ETP.G"
  )
  training = training[, colnames(training) %notin% c(to_exclude)]
  test = test[, colnames(test) %notin% c(to_exclude)]
  
}

if (params$sets_predictors %in% c(
  'G',
  'WC+SC',
  'G+WC+SC',
  'G+WC+SC+Lon+Lat',
  'G+WC+Lon+Lat',
  'G+SC+Lon+Lat',
  'G+Lon+Lat+T',
  'G+Lon+Lat+P',
  'G+Lon+Lat+T+P'
)) {
  rec <- recipe(yld_bu_ac ~ . ,
                data = training) %>%
    update_role(yld_bu_ac, new_role = 'outcome') %>%
    update_role(year, new_role = "id variable") %>%
    update_role(Year_Exp, new_role = "id variable") %>%
    update_role(-yld_bu_ac,-year, -Year_Exp, new_role = 'predictor') %>%
    step_nzv(all_predictors(), -starts_with('PC')) %>%
    step_normalize(all_numeric(),-all_outcomes(), -starts_with('PC'))
} else if (params$sets_predictors %in% c('G+L')) {
  rec <- recipes::recipe(yld_bu_ac ~ . ,
                         data = training) %>%
    recipes::update_role(yld_bu_ac, new_role = 'outcome') %>%
    recipes::update_role(year, new_role = "id variable") %>%
    recipes::update_role(Year_Exp, new_role = "id variable") %>%
    recipes::update_role(-yld_bu_ac,-year, new_role = 'predictor') %>%
    step_nzv(all_predictors(), -starts_with('PC')) %>%
    step_normalize(all_numeric(),-all_outcomes(), -starts_with('PC'))
  recipes::step_dummy(counties, one_hot = TRUE)
  
} else if (params$sets_predictors %in% c('G+Y')) {
  rec <- recipes::recipe(yld_bu_ac ~ . ,
                         data = training) %>%
    recipes::update_role(Year_Exp, new_role = "id variable") %>%
    recipes::update_role(yld_bu_ac, new_role = 'outcome') %>%
    recipes::update_role(-yld_bu_ac, new_role = 'predictor') %>%
    step_nzv(all_predictors(), -starts_with('PC')) %>%
    step_normalize(all_numeric(),-all_outcomes(), -starts_with('PC')) %>%
    recipes::step_dummy(year, preserve = F, one_hot = TRUE)
} else if (params$sets_predictors %in% c('G+Y+Lon+Lat')) {
  rec <- recipes::recipe(yld_bu_ac ~ . ,
                         data = training) %>%
    recipes::update_role(Year_Exp, new_role = "id variable") %>%
    recipes::update_role(yld_bu_ac, new_role = 'outcome') %>%
    recipes::update_role(-yld_bu_ac, new_role = 'predictor') %>%
    recipes::step_normalize(all_numeric(),-all_outcomes(), -starts_with('PC')) %>%
    recipes::step_dummy(year, preserve = F, one_hot = TRUE)
} else if (params$sets_predictors %in% c('G+WC')) {
  rec <- recipes::recipe(yld_bu_ac ~ . ,
                         data = training) %>%
    recipes::update_role(yld_bu_ac, new_role = 'outcome') %>%
    recipes::update_role(-yld_bu_ac, new_role = 'predictor') %>%
    recipes::update_role(year, new_role = "id variable") %>%
    recipes::update_role(Year_Exp, new_role = "id variable") %>%
    step_nzv(all_predictors(), -starts_with('PC')) %>%
    step_normalize(all_numeric(),-all_outcomes(), -starts_with('PC'))
} else if (params$sets_predictors %in% c('G+SC')) {
  rec <- recipes::recipe(yld_bu_ac ~ . ,
                         data = training) %>%
    recipes::update_role(yld_bu_ac, new_role = 'outcome') %>%
    recipes::update_role(-yld_bu_ac, new_role = 'predictor') %>%
    recipes::update_role(year, new_role = "id variable") %>%
    recipes::update_role(Year_Exp, new_role = "id variable") %>%
    step_nzv(all_predictors(), -starts_with('PC')) %>%
    step_normalize(all_numeric(),-all_outcomes(), -starts_with('PC'))
} else if (params$sets_predictors %in% c('G+Y+L', 'G+WC+Y+L', 'G+WC+SC+Y+L')) {
  rec <- recipes::recipe(yld_bu_ac ~ . ,
                         data = training) %>%
    recipes::update_role(yld_bu_ac, new_role = 'outcome') %>%
    recipes::update_role(-yld_bu_ac, new_role = 'predictor') %>%
    recipes::update_role(Year_Exp, new_role = "id variable") %>%
    step_nzv(all_predictors(), -starts_with('PC')) %>%
    step_normalize(all_numeric(),-all_outcomes(), -starts_with('PC')) %>%
    recipes::step_dummy(year, preserve = F, one_hot = TRUE) %>%
    recipes::step_dummy(counties, preserve = F, one_hot = TRUE)
} else if (params$sets_predictors %in% c(
  'G+WC+SC+Y+Lon+Lat',
  'G+Y+Lon+Lat+T',
  'G+Y+Lon+Lat+P',
  'G+Y+Lon+Lat+T+P',
  'G+SC+Y+Lon+Lat',
  'G+WC+Y+Lon+Lat'
)) {
  rec <- recipes::recipe(yld_bu_ac ~ . ,
                         data = training) %>%
    recipes::update_role(yld_bu_ac, new_role = 'outcome') %>%
    recipes::update_role(-yld_bu_ac, new_role = 'predictor') %>%
    recipes::update_role(Year_Exp, new_role = "id variable") %>%
    step_nzv(all_predictors(), -starts_with('PC')) %>%
    step_normalize(all_numeric(),-all_outcomes(), -starts_with('PC')) %>%
    recipes::step_dummy(year, preserve = F, one_hot = TRUE)
} else {
  rec <- recipes::recipe(yld_bu_ac ~ . ,
                         data = training) %>%
    recipes::update_role(yld_bu_ac, new_role = 'outcome') %>%
    recipes::update_role(-yld_bu_ac, new_role = 'predictor') %>%
    recipes::update_role(Year_Exp, new_role = "id variable") %>%
    step_nzv(all_predictors(), -starts_with('PC')) %>%
    step_normalize(all_numeric(),-all_outcomes(), -starts_with('PC')) %>%
    recipes::step_dummy(year, preserve = F, one_hot = TRUE) %>%
    recipes::step_dummy(counties, one_hot = TRUE)
}






## ---- echo=FALSE--------------------------------------------------------------
print(rec)

norm_obj <- prep(rec, training = training)

transformed_te <-
  bake(norm_obj, test)[, -which(colnames(bake(norm_obj, test)) %in% c('year', 'Year_Exp'))]
transformed_te_original <- bake(norm_obj, test)
transformed_tr <-
  juice(norm_obj)[, -which(colnames(juice(norm_obj)) %in% c('year', 'Year_Exp'))]
transformed_tr_original <- juice(norm_obj)

print('Number of columns after data pre-processing:')
print(ncol(transformed_tr))
print('Name columns after data pre-processing:')
print(colnames(transformed_tr))


print('Data prepared for hyperparameter optimization - First step')


## ---- echo=FALSE,include=FALSE------------------------------------------------


xgboost_model <-
  parsnip::boost_tree(
    mode = "regression",
    trees = tune(),
    min_n = tune(),
    tree_depth = tune(),
    learn_rate = tune(),
    mtry = tune(),
    sample_size = 1
  ) %>%
  set_engine("xgboost", objective = "reg:linear") %>%
  translate()





## ---- echo=FALSE--------------------------------------------------------------
xgb_grid <- parameters(trees(),
                       learn_rate(),
                       tree_depth(),
                       min_n(),
                       mtry())

xgb_grid <-
  xgb_grid %>% update(
    mtry = mtry(c(round(316 * 0.4), round(316 * 0.8))),
    trees = trees(c(4000, 7000)),
    min_n = min_n(c(5, 18)),
    learn_rate = learn_rate(range(c(3e-04, 0.01)), trans = NULL),
    tree_depth = tree_depth(c(2, 12))
  )
#levels = c(trees=4,min_n=4,sample_size=1,tree_depth=4,learn_rate=4,mtry=1)



## ---- echo=FALSE,message=FALSE,warning=FALSE----------------------------------
print('Starting hyperparameter optimization')
if (params$sets_predictors %in% c(
  'G',
  'WC+SC',
  'G+WC+SC',
  'G+WC+SC+Lon+Lat',
  'G+WC+Lon+Lat',
  'G+SC+Lon+Lat',
  'G+Lon+Lat+T',
  'G+Lon+Lat+P',
  'G+Lon+Lat+T+P'
)) {
  wf <- workflow() %>%
    add_model(xgboost_model) %>%
    add_formula(yld_bu_ac ~ .)
  
} else if (params$sets_predictors %in% c('G+L')) {
  wf <- workflow() %>%
    add_model(xgboost_model) %>%
    add_formula(yld_bu_ac ~ .)
  
  
} else if (params$sets_predictors %in% c('G+WC')) {
  wf <- workflow() %>%
    add_model(xgboost_model) %>%
    add_formula(yld_bu_ac ~ .)
  
  
} else if (params$sets_predictors %in% c('G+SC')) {
  wf <- workflow() %>%
    add_model(xgboost_model) %>%
    add_formula(yld_bu_ac ~ .)
  
  
} else if (params$sets_predictors %in% c('G+Y+L')) {
  wf <- workflow() %>%
    add_model(xgboost_model) %>%
    add_formula(yld_bu_ac ~ .)
  
  
} else if (params$sets_predictors %in% c(
  'G+Y',
  'G+Y+Lon+Lat',
  'G+WC+Y+L',
  'G+WC+SC+Y+L',
  'G+WC+SC+Y+Lon+Lat',
  'G+Y+Lon+Lat+T',
  'G+Y+Lon+Lat+P',
  'G+Y+Lon+Lat+T+P',
  'G+SC+Y+Lon+Lat',
  'G+WC+Y+Lon+Lat'
)) {
  wf <- workflow() %>%
    add_model(xgboost_model) %>%
    add_formula(yld_bu_ac ~ .)
  
  
} else {
  wf <- workflow() %>%
    add_model(xgboost_model) %>%
    add_formula(yld_bu_ac ~ .)
  
}



## Write the datasets used for training and test datasets

tr_dataset <-
  paste(
    '~/LYOUT_revisions/',
    params$method,
    '/',
    params$year_to_predict,
    '/training_dataset_',
    'geno_info_',
    params$geno_info,
    '_',
    params$sets_predictors,
    '.txt',
    sep = ''
  )

write.table(transformed_tr_original, file = tr_dataset)

te_dataset <-
  paste(
    '~/LYOUT_revisions/',
    params$method,
    '/',
    params$year_to_predict,
    '/test_dataset_',
    'geno_info_',
    params$geno_info,
    '_',
    params$sets_predictors,
    '.txt',
    sep = ''
  )

write.table(transformed_te_original, file = te_dataset)


registerDoFuture()
cl <- makeCluster(cores)
plan(cluster, workers = cl)

##

if (params$tuning_hyper_parameters == '4-fold-CV-year') {
  set.seed(params$seed)
  folds <- vfold_cv(transformed_tr,
                    strata = 'year',
                    repeats = 2,
                    v = 5)
  
  set.seed(params$seed)
  opt_res <- wf %>%
    tune_bayes(
      resamples = folds,
      param_info = xgb_grid,
      iter = 30,
      initial = 10,
      metrics = yardstick::metric_set(rmse),
      control = tune::control_bayes(verbose = TRUE, no_improve = 20)
    )
  
  print(opt_res)
  
}

if (params$tuning_hyper_parameters == 'random-CV') {
  set.seed(params$seed)
  folds <- vfold_cv(transformed_tr, repeats = 2, v = 5)
  set.seed(params$seed)
  opt_res <- wf %>%
    tune_bayes(
      resamples = folds,
      param_info = xgb_grid,
      iter = 30,
      initial = 10,
      metrics = yardstick::metric_set(rmse),
      control = tune::control_bayes(verbose = TRUE, no_improve = 20)
    )
  print(opt_res)
  
}







## ----echo=FALSE---------------------------------------------------------------
print('Hyperparameter optimization done.')


## ----echo=FALSE---------------------------------------------------------------
print(opt_res %>%
        collect_metrics())



## ----echo=FALSE---------------------------------------------------------------
autoplot(opt_res)


## ----echo=FALSE---------------------------------------------------------------
xgboost_best_params <- opt_res %>%
  tune::select_best("rmse", maximize = FALSE) -> highest_score

xgboost_best_params

xgboost_model_final <- finalize_workflow(wf,
                                         xgboost_best_params)

output_file <-
  paste(
    '~/LYOUT_revisions/',
    params$method,
    '/',
    params$year_to_predict,
    '/best_hyperparameters_params_',
    'geno_info_',
    params$geno_info,
    '_',
    params$sets_predictors,
    '_',
    params$trait,
    '_tuning_method_',
    params$tuning_hyper_parameters,
    '.RDS',
    sep = ''
  )

saveRDS(opt_res, file = output_file)





## ----echo=FALSE,include=FALSE-------------------------------------------------

final_xgboost <-
  xgboost_model_final %>%
  fit(data = transformed_tr)

#transformed_tr2<-list()
#transformed_tr2$data<-as.matrix(transformed_tr %>% select(-yld_bu_ac))
#transformed_tr2$label<-as.matrix(transformed_tr[, 'yld_bu_ac'])

#mod <-
#  xgboost::xgboost(
#    data = transformed_tr2$data,
#    label = transformed_tr2$label,
#    xgb_param = list(
#      objective = "reg:linear",
#      nrounds = xgboost_best_params$trees,
#      eta = xgboost_best_params$learn_rate,
#      max_depth = xgboost_best_params$tree_depth,
#      min_child_weight = xgboost_best_params$min_n,
#    ),
#    nrounds = xgboost_best_params$trees,
#    verbose = TRUE,
#    early_stopping_rounds = 50
#  )





## Across Environments

m <- final_xgboost %>% predict(new_data = transformed_tr)  %>%
  bind_cols(transformed_tr_original) %>%
  yardstick::metrics(yld_bu_ac, .pred) %>%
  mutate(.estimate = format(round(.estimate, 2), big.mark = ","))


m1 <- final_xgboost %>% predict(new_data = transformed_tr) %>%
  bind_cols(transformed_tr_original) %>%
  select(yld_bu_ac, .pred)  %>%
  cor(method = 'pearson')

m2 <- final_xgboost %>% predict(new_data = transformed_te)  %>%
  bind_cols(transformed_te_original) %>%
  yardstick::metrics(yld_bu_ac, .pred) %>%
  mutate(.estimate = format(round(.estimate, 2), big.mark = ","))


m3 <- final_xgboost %>% predict(new_data = transformed_te) %>%
  bind_cols(transformed_te_original) %>%
  select(yld_bu_ac, .pred)  %>%
  cor(method = 'pearson')

## By Environment

s <- final_xgboost %>% predict(new_data = transformed_tr)  %>%
  bind_cols(transformed_tr_original) %>%
  group_by(Year_Exp) %>%
  yardstick::metrics(yld_bu_ac, .pred) %>%
  mutate(.estimate = format(round(.estimate, 2), big.mark = ","))

s$.estimate = as.numeric(s$.estimate)
s_mean <-
  s %>% group_by(.metric) %>% mutate(mean = mean(.estimate)) %>% mutate(sd =
                                                                          sd(.estimate))
s_mean <- unique(s_mean[, c(2, 5, 6)])

s1 <- final_xgboost %>% predict(new_data = transformed_tr) %>%
  bind_cols(transformed_tr_original) %>%
  group_by(Year_Exp) %>%
  summarize(COR = cor(yld_bu_ac, .pred, method = 'pearson'))

s1_mean <- mean(s1$COR)
s1_sd <- sd(s1$COR)

predicted_values_test_set <-
  final_xgboost %>% predict(new_data = transformed_te)  %>%
  bind_cols(transformed_te_original)

s2 <- final_xgboost %>% predict(new_data = transformed_te)  %>%
  bind_cols(transformed_te_original) %>%
  group_by(Year_Exp) %>%
  yardstick::metrics(yld_bu_ac, .pred) %>%
  mutate(.estimate = format(round(.estimate, 2), big.mark = ","))

s2$.estimate = as.numeric(s2$.estimate)
s2_mean <-
  s2 %>% group_by(.metric) %>% mutate(mean = mean(.estimate)) %>% mutate(sd =
                                                                           sd(.estimate))
s2_mean <- unique(s2_mean[, c(2, 5, 6)])

s3 <- final_xgboost %>% predict(new_data = transformed_te) %>%
  bind_cols(transformed_te_original) %>%
  group_by(Year_Exp) %>%
  summarize(COR = cor(yld_bu_ac, .pred, method = 'pearson'))

s3_mean <- mean(s3$COR)
s3_sd <- sd(s3$COR)



varimps <-
  xgb.importance(model = pull_workflow_fit(final_xgboost)$fit)





final_object <-
  list(
    predicted_values_test_set,
    m,
    m1,
    m2,
    m3,
    varimps,
    s,
    s1,
    s2,
    s3,
    s_mean,
    s1_mean,
    s1_sd,
    s2_mean,
    s3_mean,
    s3_sd
  )


names(final_object) <-
  c(
    'predicted_values_test_set',
    'tr_metrics_across_all',
    'tr_corr_across_all',
    'te_metrics_across_all',
    'te_corr_across_all',
    'varimps',
    'tr_metrics_by_yearexp',
    'tr_corr_by_yearexp',
    'te_metrics_by_yearexp',
    'te_corr_by_yearexp',
    'tr_mean_sd_metrics',
    'tr_mean_corr',
    'tr_sd_corr',
    'te_mean_sd_metrics',
    'te_mean_corr',
    'te_sd_corr'
  )


output_file <-
  paste(
    '~/LYOUT_revisions/',
    params$method,
    '/',
    params$year_to_predict,
    '/final_list_results_year',
    params$year_to_predict,
    'geno_info_',
    params$geno_info,
    '_',
    params$sets_predictors,
    '_',
    params$trait,
    '_tuning_method_',
    params$tuning_hyper_parameters,
    '.RDS',
    sep = ''
  )

saveRDS(final_object, file = output_file)

output_file <-
  paste(
    '~/LYOUT_revisions/',
    params$method,
    '/',
    params$year_to_predict,
    '/fitted_model_index_',
    params$year_to_predict,
    'geno_info_',
    params$geno_info,
    '_',
    params$sets_predictors,
    '_',
    params$trait,
    '_tuning_method_',
    params$tuning_hyper_parameters,
    '.RDS',
    sep = ''
  )

saveRDS(final_xgboost, file = output_file)


## ----echo=FALSE---------------------------------------------------------------
varimps %>%
  top_n(40, Gain) %>%
  ggplot(aes(reorder(Feature, Gain), Gain)) +
  geom_col(aes(fill = Gain))  +
  scale_fill_gradient2(low = "white",
                       high = "blue",
                       midpoint = median(varimps$Gain))  +
  coord_flip()  +
  labs(x = "Feature")

varimps[order(-varimps$Gain), ]



## ----echo=FALSE---------------------------------------------------------------
print(m)
print(m1)



## ----echo=FALSE---------------------------------------------------------------
print(s)


## ----echo=FALSE---------------------------------------------------------------
print(s_mean)


## ----echo=FALSE---------------------------------------------------------------
print(s1)


## ----echo=FALSE---------------------------------------------------------------
print(s1_mean)


## ----echo=FALSE---------------------------------------------------------------
print(s1_sd)


## ----echo=FALSE---------------------------------------------------------------
print(m2)
print(m3)



## ----echo=FALSE---------------------------------------------------------------
print(s2)


## ----echo=FALSE---------------------------------------------------------------
print(s2_mean)


## ----echo=FALSE---------------------------------------------------------------
print(s3)


## ----echo=FALSE---------------------------------------------------------------
print(s3_mean)


## ----echo=FALSE---------------------------------------------------------------
print(s3_sd)


## ----echo=FALSE---------------------------------------------------------------
if (params$sets_predictors %in% c(
  'WC+SC',
  'WC+SC+Y+L',
  'G+WC+SC',
  'G+WC+SC+Y+L',
  'G+WC',
  'G+WC+Y+L',
  'G+WC+SC+Lon+Lat',
  'G+WC+Lon+Lat',
  'G+WC+SC+Y+Lon+Lat'
)) {
  draft <- TRUE
} else{
  draft <- FALSE
}
print(draft)


end.time <- Sys.time()
time.taken <- difftime(end.time, start.time, units = 'hours')
print(time.taken)



## ----echo=FALSE---------------------------------------------------------------
print(sessionInfo())
