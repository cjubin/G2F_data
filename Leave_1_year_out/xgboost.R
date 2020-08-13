#######################
#######################
args = commandArgs(trailingOnly = TRUE)
x1 = args[1]
x2 = args[2]
x3 = args[3]
x4 = args[4]
x5 = args[5]
x6 = args[6]
x7 = args[7]
x8 = args[8]
x9 = args[9]

##




xgboost_year_prediction = function(geno_info = c('snps', 'PCs'),
                                   phenos_file,
                                   geno_file,
                                   seed = 105,
                                   nbootstrap = 100,
                                   year_to_predict,
                                   trait,
                                   sets_predictors = c('G',
                                                       'WC+SC',
                                                       'WC+SC+Y+L',
                                                       'Y+L',
                                                       'G+WC+SC',
                                                       'G+Y+L',
                                                       'G+WC+SC+Y+L',
                                                       'G+Y',
                                                       'G+L'),
                                   method,
                                   tuning_hyper_parameters = c('default', '4-fold-CV-year', 'random-CV'),
                                   bootstrap_sampling = c('total_random', 'yearly_based')) {
  library(caret)
  library(doParallel)
  library(data.table)
  library(xgboost)
  source('/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/ML_PREDICTIONS/setSeeds.R')
  `%notin%` <- Negate(`%in%`)
  
  
  ### First step: pre-processing the data and split the data according to the year to predict
  ### Load dataset with all observations with environmental and genomic predictors
  
  phenos = read.table(phenos_file,header = T,sep = '\t')
  phenos<-phenos[,-which(colnames(phenos)%in%c("parent1","parent2","parent1.GBS.sample","parent2.GBS.sample"))]
  
  geno_hybrids=fread(geno_file)
  geno_hybrids=as.data.frame(geno_hybrids)
  geno_hybrids<-geno_hybrids[,-which(colnames(geno_hybrids)%in%c("parent1","parent2","parent1.GBS.sample","parent2.GBS.sample"))]
  print('Data read')
  
  colnames(geno_hybrids)[2:ncol(geno_hybrids)]<-paste0('SNP',colnames(geno_hybrids)[2:ncol(geno_hybrids)])
  phenos=merge(phenos,geno_hybrids,by='pedigree',all.x = T)
  
  
  
  phenos$Year = as.factor(phenos$year)
  
  ##Predictors variables included for prediction according to the sets of predictors defined by the user.
  toMatch = c('.V', '.F', '.G', 'length.growing.season')
  toMatch2 = c('SNP', 'UD')
  toMatch3 = c('year')
  toMatch4 = c('counties')
  toMatch5 = c('.SC')
  matches1 <-
    grep(paste(toMatch, collapse = "|"), colnames(phenos), value = TRUE)
  matches2 <-
    grep(paste(toMatch2, collapse = "|"), colnames(phenos), value = TRUE)
  matches3 <-
    grep(paste(toMatch3, collapse = "|"), colnames(phenos), value = TRUE)
  matches4 <-
    grep(paste(toMatch4, collapse = "|"), colnames(phenos), value = TRUE)
  matches5 <-
    grep(paste(toMatch5, collapse = "|"), colnames(phenos), value = TRUE)
  
  if (sets_predictors == 'WC+SC') {
    predictors = c(matches1, matches5)
  }
  if (sets_predictors == 'G') {
    predictors = c(matches2)
  }
  if (sets_predictors == 'Y+L') {
    predictors = c(matches3, matches4)
  }
  if (sets_predictors == 'WC+SC+Y+L') {
    predictors = c(matches1, matches3, matches4, matches5)
  }
  if (sets_predictors == 'G+WC+SC') {
    predictors = c(matches1, matches2, matches5)
  }
  if (sets_predictors == 'G+Y+L') {
    predictors = c(matches2, matches3, matches4)
  }
  if (sets_predictors == 'G+WC+SC+Y+L') {
    predictors = c(matches1, matches2, matches3, matches4, matches5)
  }
  if (sets_predictors == 'G+Y') {
    predictors = c(matches2, matches3)
  }
  if (sets_predictors == 'G+L') {
    predictors = c(matches2, matches4)
  }
  
  #Selecting the variables based on the list of predictors defined previously
  #Create dummy variables if factor variables are present
  phenos2 = phenos[, colnames(phenos) %in% c(predictors, trait)]
  
  ##Retain the factor variables which should not be pre-processed by standardization once they have been converted to dummies
  # year variable to split train/test afterwards
  
  year = phenos[, 'Year']
  unique_years=as.character(unique(year))
  loc=phenos[,matches4]
  unique_loc=as.character(unique(loc))
  
  unique_years=paste('Year.',unique_years,sep = '')
  unique_loc=paste(matches4,'.',unique_loc,sep = '')
  unique_loc=gsub(' ','.',unique_loc)
  unique_loc=gsub('-','.',unique_loc)  
  
  #Conversion of factor or character variables to dummy variables
  print(colnames(phenos2))
  
  if (sets_predictors%in%c('WC+SC+Y+L',
      'Y+L',
      'G+Y+L',
      'G+WC+SC+Y+L',
      'G+Y',
      'G+L')){
  converted_dummies <-
    data.frame(predict(caret::dummyVars(~ ., data = phenos2), newdata =
                         phenos2))
  
  #Splitting the data into a training and testing set based on the year to predict
  training =converted_dummies [-which(year==year_to_predict),]
  training_set_year = year[year %notin% year_to_predict]
  row.names(training_set_year) = NULL
  
  test = converted_dummies [which(year==year_to_predict),]}
  
  #Pre-processing
  training =phenos2[-which(year==year_to_predict),]
  training_set_year = year[year %notin% year_to_predict]
  row.names(training_set_year) = NULL
  
  test = phenos2[which(year==year_to_predict),]
  
  df1<-data.matrix(training[, -which(colnames(training) %in% c(trait,unique_years,unique_loc))])
  df2<-data.matrix(test[, -which(colnames(test) %in% c(trait,unique_years,unique_loc))])
  
  preProcValues <-
    caret::preProcess(df1, method = c("center", "scale"))
  
  trainTransformed <-
    predict(preProcValues, df1)
  testTransformed <-
    predict(preProcValues, df2)
  
  if (sets_predictors %in% c('WC+SC+Y+L', 'Y+L', 'G+Y+L', 'G+WC+SC+Y+L')) {
    trainTransformed <-
      as.data.frame(cbind(training[, c(trait, unique_years, unique_loc)], trainTransformed))
    colnames(trainTransformed)[1] = trait
    testTransformed <-
      as.data.frame(cbind(test[, c(trait, unique_years, unique_loc)], testTransformed))
  } else if (sets_predictors %in% c('G+Y')) {
    trainTransformed <-
      as.data.frame(cbind(training[, c(trait, unique_years)], trainTransformed))
    colnames(trainTransformed)[1] = trait
    testTransformed <-
      as.data.frame(cbind(test[, c(trait, unique_years)], testTransformed))
  } else if (sets_predictors %in% c('G+L')) {
    trainTransformed <-
      as.data.frame(cbind(training[, c(trait, unique_loc)], trainTransformed))
    colnames(trainTransformed)[1] = trait
    testTransformed <-
      as.data.frame(cbind(test[, c(trait, unique_loc)], testTransformed))
  } else {
    trainTransformed <-
      as.data.frame(cbind(training[, c(trait)], trainTransformed))
    colnames(trainTransformed)[1] = trait
    print(trainTransformed[,trait])
    print(class(trainTransformed))
    print(which(is.na(trainTransformed[,trait])))
    print(which(is.na(trainTransformed), arr.ind=TRUE))
    testTransformed <-
      as.data.frame(cbind(test[, c(trait)], testTransformed))
  }
  
  
  colnames(testTransformed)[1] = trait
  row.names(trainTransformed) = NULL
  row.names(testTransformed) = NULL
  
  print('Data prepared for hyperparameter optimization - First step')
  
  ###Second step: hyperparameter optimization according to the option set in tuning_hyper_parameters
  #Grid search: grid table
  
  tune_grid <- expand.grid(
    nrounds = c(800,1500),
    eta = c(0.01),
    max_depth = c(4),
    gamma = 0,
    colsample_bytree = 1,
    min_child_weight = c(10),
    subsample = 1
  )
  
  ##List of observations according to the year, to specify the indexes within each fold
  #index: a list with elements for each resampling iteration. Each list element is a vector of integers corresponding to the rows used for training at that iteration.
  tmp = unique(training_set_year)
  list1 = list()
  list2 = list()
  for (j in 1:length(tmp)) {
    list1[[j]] = which(training_set_year %notin% tmp[j])
    list2[[j]] = which(training_set_year == tmp[j])
  }
  
  ##Define the method for hyperparameter tuning: default parameters 'default'
  #or cv on 4 folds, each fold defined by the year of phenotypic observations: '4-fold-CV-year'
  #or random cv (5-fold cv repeated two times), independent of the year: 'random-CV'
  if (tuning_hyper_parameters %in% c('4-fold-CV-year', 'random-CV')) {
    if (tuning_hyper_parameters == '4-fold-CV-year') {
      fitControl <- trainControl(
        method = 'cv',
        # k-fold cross validation
        number = 4,
        index = list1,
        indexOut = list2,
        savePredictions = 'final'
      )
    }
    if (tuning_hyper_parameters == 'random-CV') {
      cvseeds <-
        setSeeds(
          method = 'repeatedcv',
          numbers = 4,
          repeats = 2,
          seed = 1598
        )
      fitControl <- trainControl(
        method = 'repeatedcv',
        number = 4,
        # k-fold cross validation
        repeats = 2,
        savePredictions = 'final',
        seeds = cvseeds
      )
    }
    
    #Detecting number of cores available
    cores <- as.integer(Sys.getenv('SLURM_NTASKS'))
    print(cores)
    
    
    print('Starting hyperparameter optimization')
    xgb_fit_hyperparameters <-
      caret::train(
        yld_bu_ac ~ .,
        data = trainTransformed,
        method = "xgbTree",
        metric = 'RMSE',
        tuneGrid = tune_grid,
        trControl = fitControl,
        verbose = TRUE,
        num.threads = cores 
      )
    
    print('Hyperparameter optimization done.')
    
    
    
    ##Extract the best hyperparameters in a grid table
    best_hyperparameters = expand.grid(
      nrounds = xgb_fit_hyperparameters$bestTune$nrounds,
      eta = xgb_fit_hyperparameters$bestTune$eta,
      max_depth = xgb_fit_hyperparameters$bestTune$max_depth,
      gamma = xgb_fit_hyperparameters$bestTune$gamma,
      colsample_bytree = xgb_fit_hyperparameters$bestTune$colsample_bytree,
      min_child_weight = xgb_fit_hyperparameters$bestTune$min_child_weight,
      subsample = xgb_fit_hyperparameters$bestTune$subsample
    )
    output_file <-
      paste(
        '/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/ML_PREDICTIONS/LYOUT/',
        method,
        '/best_hyperparameters',
        year_to_predict,
        '_',
        sets_predictors,
        '_',
        trait,
        '_tuning_method_',
        tuning_hyper_parameters,
        '.RDS',
        sep = ''
      )
    
    saveRDS(xgb_fit_hyperparameters, file = output_file)
    
    #Plot tuning hyperparameters
    #trellis.par.set(caretTheme())
    pdf(
      paste(
        '/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/ML_PREDICTIONS/LYOUT/',
        method,
        '/resampling_plot_hyperparameters',
        year_to_predict,
        '_',
        sets_predictors,
        '_',
        trait,
        '_tuning_method_',
        tuning_hyper_parameters,
        '.pdf',
        sep = ''
      ),
      width = 8,
      height = 8
    )
    
    print(plot(xgb_fit_hyperparameters))
    dev.off()
  }
  
  ####100 Bootstrap sample datasets generated using the TrainingTransformed to predict the test set
  
  ##Option bootstrap_sampling
  
  #Option 'yearly_based'
  #By setting a seed, we ensure that the same bootstrap samples can be generated to test fairly each method
  #Generating nbootstrap datasets: sampling same number of observations per year
  
  if (bootstrap_sampling == 'yearly_based') {
    
    set.seed(seed)
    
    list_bootstrap = list()
    for (j in 1:4) {
      length_training_set = nrow(trainTransformed)
      list_bootstrap[[j]] = lapply(1:nbootstrap, function(x) {
        sample(which(training_set_year == tmp[j]),
               length(which(training_set_year == tmp[j])),
               replace = TRUE)
      })
    }
    bootstrap_samples = list()
    for (i in 1:nbootstrap) {
      bootstrap_samples[[i]] = c(list_bootstrap[[1]][[i]],
                                 list_bootstrap[[2]][[i]],
                                 list_bootstrap[[3]][[i]],
                                 list_bootstrap[[4]][[i]])
    }
    
  }
  #Option total_random
  #More observations from a year compared to the others in the training set can be sampled - no control about the sampling procedure
  if (bootstrap_sampling == 'total_random') {
    set.seed(seed)
    for (i in 1:nbootstrap) {
      bootstrap_samples[[i]] = sample(rownames(trainTransformed),
                                      nrow(trainTransformed),
                                      replace = TRUE)
    }
    
  }
  
  output_file <-
    paste(
      '/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/ML_PREDICTIONS/LYOUT/',
      method,
      '/bootstrap_samples',
      year_to_predict,
      '_',
      sets_predictors,
      '_',
      trait,
      '_bootstrap_strategy_',
      bootstrap_sampling,
      '.RDS',
      sep = ''
    )
  
  saveRDS(bootstrap_samples, file = output_file)
  
  ###############################
  ###############################
  #Prediction of the year_to_predict using nbootstrap samples and best hyperparameters/default (option tuning_hyper_parameters)
  
  cores <- as.integer(Sys.getenv('SLURM_NTASKS'))
  
  print(cores)
  cl <- makeForkCluster(cores)
  registerDoParallel(cl)
  
  if (tuning_hyper_parameters %in% c('4-fold-CV-year', 'random-CV')) {
    
    #Training the nbootstrap samples with the best hyperparameters found
    
    #Note on varImp:
    #The function automatically scales the importance scores to be between 0 and 100. Using scale = FALSE avoids this normalization step.
    
    xgb_fit_bootstrap <- foreach(i = 1:nbootstrap) %dopar% {
      model <- caret::train(
        yld_bu_ac ~ .,
        data = trainTransformed[bootstrap_samples[[i]],],
        method = "xgbTree",
        metric = 'RMSE',
        trControl = trainControl(method = "none"),
        tuneGrid = xgb_fit_hyperparameters$bestTune,
        verbose = TRUE,
        num.threads = 7
      )
      pred <-
        predict(model, newdata = testTransformed[,-which(colnames(testTransformed) %in%
                                                           trait)])
      R2squared = caret::R2(
        pred = pred,
        obs = testTransformed[, trait] ,
        formula = 'traditional',
        na.rm = FALSE
      )#Record model metric: R2squared: R^2 = 1-\frac{\sum (y_i - hat{y}_i)^2}{\sum (y_i - \bar{y}_i)^2}
      RMSE = RMSE(pred = pred,
                  obs = testTransformed[, trait],
                  na.rm = FALSE)#RMSE
      pearson = cor(pred, testTransformed[, trait], method = 'pearson')#Pearson coefficient
      varImp = varImp(model, scale = FALSE)#varImp scores
      list(
        'model' = model,
        'pred' = pred,
        'R2squared' = R2squared,
        'RMSE' = RMSE,
        'pearson' = pearson,
        'varImp' = varImp
      )
    }
  }
  
  
  
  
  if (tuning_hyper_parameters %in% c('default')) {
    
    #Training bootstrap samples with default parameters
    
    #Note on varImp:
    #The function automatically scales the importance scores to be between 0 and 100. Using scale = FALSE avoids this normalization step.
    
    xgb_fit_bootstrap <- foreach(i = 1:nbootstrap) %dopar% {
      model <- caret::train(
        yld_bu_ac ~ .,
        data = trainTransformed[bootstrap_samples[[i]],],
        method = "xgbTree",
        metric = 'RMSE',
        trControl = trainControl(method = "none"),
        verbose = TRUE,
        num.threads = 7
      )
      pred <-
        predict(model, newdata = testTransformed[,-which(colnames(testTransformed) %in%
                                                           trait)])
      R2squared = caret::R2(
        pred = pred,
        obs = testTransformed[, trait] ,
        formula = 'traditional',
        na.rm = FALSE
      )#Record model metric: R2squared: R^2 = 1-\frac{\sum (y_i - hat{y}_i)^2}{\sum (y_i - \bar{y}_i)^2}
      RMSE = RMSE(pred = pred,
                  obs = testTransformed[, trait],
                  na.rm = FALSE)#RMSE
      pearson = cor(pred, testTransformed[, trait], method = 'pearson')#Pearson coefficient
      varImp = varImp(model, scale = FALSE)#varImp scores
      list(
        'model' = model,
        'pred' = pred,
        'R2squared' = R2squared,
        'RMSE' = RMSE,
        'pearson' = pearson,
        'varImp' = varImp
      )
    }
  }
  
  
  
  stopCluster(cl)
  
  output_file <-
    paste(
      '/home/uni08/jubin1/Data/GenomesToFields/G2F20142018/ML_PREDICTIONS/LYOUT/',
      method,
      '/',
      year_to_predict,
      '_',
      sets_predictors,
      '_',
      trait,
      '_tuning_method_',
      tuning_hyper_parameters,
      '_bootstrap_strategy_',
      bootstrap_sampling,
      '.RDS',
      sep = ''
    )
  
  saveRDS(xgb_fit_bootstrap, file = output_file)
  
  
  
  return(xgb_fit_bootstrap)
  
}

xgboost_year_prediction(
  phenos_file=as.character(x8),
  geno_file=as.character(x9),
  nbootstrap = as.numeric(x1),
  year_to_predict = x2,
  trait = as.character(x3),
  sets_predictors = as.character(x4),
  method = as.character(x5),
  tuning_hyper_parameters = as.character(x6),
  bootstrap_sampling = as.character(x7)
)
