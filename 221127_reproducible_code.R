# > sessionInfo()
# R version 4.2.0 (2022-04-22 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19043)
# 
# Matrix products: default
# 
# locale:
# [1] LC_COLLATE=German_Germany.utf8  LC_CTYPE=German_Germany.utf8    LC_MONETARY=German_Germany.utf8
# [4] LC_NUMERIC=C                    LC_TIME=German_Germany.utf8    
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] pROC_1.18.0       wesanderson_0.3.6 gbm_2.1.8         timbR_0.1.0       rpart.plot_3.1.1 
# [6] rpart_4.1.16      ranger_0.14.1    
# 
# loaded via a namespace (and not attached):
# [1] reticulate_1.25   tidyselect_1.1.2  xfun_0.30         purrr_0.3.4       splines_4.2.0    
# [6] haven_2.5.0       lattice_0.20-45   colorspace_2.0-3  vctrs_0.4.1       generics_0.1.3   
# [11] htmltools_0.5.2   yaml_2.3.5        utf8_1.2.2        survival_3.3-1    rlang_1.0.2      
# [16] pillar_1.8.0      glue_1.6.2        DBI_1.1.2         plyr_1.8.7        lifecycle_1.0.1  
# [21] munsell_0.5.0     gtable_0.3.0      evaluate_0.15     DALEX_2.4.1       knitr_1.39       
# [26] forcats_0.5.1     fastmap_1.1.0     fansi_1.0.3       Rcpp_1.0.9        scales_1.2.0     
# [31] backports_1.4.1   DALEXtra_2.2.1    checkmate_2.1.0   mlr3viz_0.5.9     jsonlite_1.8.0   
# [36] ggplot2_3.3.6     hms_1.1.1         png_0.1-7         digest_0.6.29     dplyr_1.0.9      
# [41] grid_4.2.0        cli_3.3.0         tools_4.2.0       magrittr_2.0.3    tibble_3.1.8     
# [46] mlr3misc_0.10.0   pkgconfig_2.0.3   ellipsis_0.3.2    Matrix_1.4-1      data.table_1.14.2
# [51] assertthat_0.2.1  rmarkdown_2.14    R6_2.5.1          compiler_4.2.0   

#devtools::install_github("imbs-hl/timbR", "develop")
library(ranger)
library(rpart)
library(rpart.plot)
library(timbR)
library(gbm)
library(wesanderson)
load("tvdata.Robj")


# modify function from timbR package for posterior probabilities (method == "prediction")
measure_distances <- function(rf, metric = "splitting variables", test_data = NULL, positive = NULL){
  
  ## Check inputs ----
  if (!checkmate::testClass(rf, "ranger")){
    stop("rf must be of class ranger")
  }
  if (!checkmate::testList(rf$forest)){
    stop("rf must be trained using write.forest = TRUE.")
  }
  if (!checkmate::testChoice(metric,
                             choices = c("splitting variables", "weighted splitting variables", "terminal nodes", "prediction"))){
    stop(paste("metric has to be from c('splitting variables', 'weighted splitting variables', 'terminal nodes', 'prediction')."))
  }
  if (metric %in% c("terminal nodes", "prediction")){
    if (checkmate::testNull(test_data)){
      stop("You have to provide a test data set for distance measure by terminal nodes or prediction.")
    }
    if ("try-error" %in% class(try(predict(rf, data = test_data)))){
      stop("The provided test data set does not fit to the provided ranger object")
    }
    if (nrow(test_data) < 2){
      stop("You have to provide at least two samples as a test data set")
    }
  } else {
    if (!checkmate::testNull(test_data)){
      message("You provided a test data set for a distance measure by splitting variables. This is not necessary and will be ignored.")
    }
  }
  
  ## Prepare matrix for output ----
  distances <- matrix(data = NA, nrow = rf$num.trees, ncol = rf$num.trees)
  
  ## Extract outcome id
  outcome_id <- rf$forest$dependent.varID
  
  ## Extract number of features
  num_features <- rf$num.independent.variables
  
  ## Calculation for d0 of Banerjee et al. (2012) ----
  if (metric == "splitting variables"){
    ## Simplify for each tree which features were used
    feature_usage <- lapply(X      = 1:rf$num.trees,
                            FUN    = function(x){
                              splitting_variables <- sort(unique(treeInfo(rf, x)$splitvarID))
                              fu <- rep(0, num_features)
                              fu[(splitting_variables+1)] <- 1
                              fu
                            })
    
    ## Calculate standardized pair-wise distances
    for (i in 1:rf$num.trees){
      for (j in 1:rf$num.trees){
        distances[i,j] <- sum((feature_usage[[i]] - feature_usage[[j]])^2)/num_features
      }
    }
  }
  
  ## Calculation weighted version of d0 ----
  if (metric == "weighted splitting variables"){
    ## Calculate usage score for each variable
    US <- lapply(1:rf$num.trees, function(i){
      ## initialize levels for ith tree
      split_level <- rep(1, length(rf$forest$split.varIDs[[i]]))
      
      ## Extract child node IDs for ith tree
      child.nodeIDs <- rf$forest$child.nodeIDs[[i]]
      
      ## Extract split var IDs for ith tree
      split.varIDs <- rf$forest$split.varIDs[[i]]
      
      ## Child nodes get level of parent node + 1 ----
      for (j in 1:length(rf$forest$split.varIDs[[i]])){
        if (child.nodeIDs[[1]][j] != 0){
          split_level[child.nodeIDs[[1]][j] + 1] <- split_level[j] + 1
        }
        
        if (child.nodeIDs[[2]][j] != 0){
          split_level[child.nodeIDs[[2]][j] + 1] <- split_level[j] + 1
        }
      }
      
      ## Usage score for each variable
      US <- rep(0, rf$num.independent.variables)
      
      US <- lapply(1:rf$num.independent.variables, function(j){
        sum(1/(2^(split_level[split.varIDs == j] - 1))) / (max(split_level) - 1)
      })
      
      as.numeric(do.call("cbind", US))
    })
    
    US <- do.call("rbind", US)
    
    
    distance <- lapply(1:rf$num.trees, function(x){
      distance <- lapply(1:rf$num.trees, function(y){
        #1/rf$num.independent.variables * sum((US[x,] - US[y,])^2)
        sum((US[x,] - US[y,])^2)
      })
      
      as.numeric(do.call("rbind", distance))
    })
    
    distances <- as.matrix(do.call("rbind", distance))
    
  }
  
  ## Calculation for d1 of Banerjee et al. (2012) ----
  if (metric == "terminal nodes"){
    # distances <- matrix(data = 0, nrow = rf$num.trees, ncol = rf$num.trees)
    
    ## Initialize matrix for terminal nodes
    term_node <- predict(rf, data = test_data, type = "terminalNodes")$predictions
    
    ## Calculate if observations end in same terminal node for each tree
    I <- list()
    
    for (x in 1:rf$num.trees){
      I[[x]] <- matrix(data = NA, nrow = nrow(test_data), ncol = nrow(test_data))
      
      for (i in 1:nrow(test_data)){
        for (j in 1:nrow(test_data)){
          if (term_node[i,x] == term_node[j,x]){
            I[[x]][i,j] <- 1
          } else {
            I[[x]][i,j] <- 0
          }
        }
      }
    }
    
    ## Calculate distances
    for (x in 1:rf$num.trees){
      for (y in 1:rf$num.trees){
        for (i in 1:(nrow(test_data)-1)){
          for (j in (i+1):nrow(test_data)){
            distances[x,y] <- distances[x,y] + abs(I[[x]][i,j] - I[[y]][i,j])
          }
        }
      }
    }
    ## Normailze distances
    distances <- distances / choose(nrow(test_data), 2)
    
  }
  
  ## Calculation for d2 of Banerjee et al. (2012) ----
  if (metric == "prediction"){
    
    ## Predict outcome for all test data
    pred <- predict(rf, data = test_data, predict.all = TRUE)
    
    ## Controll if outcome is factor
    if (is.factor(rf$predictions)){
      ## Calculate standardized pair-wise distances
      for (i in 1:rf$num.trees){
        for (j in 1:rf$num.trees){
          distances[i,j] <- sum(pred$predictions[,i] != pred$predictions[,j])/nrow(test_data)
        }
      }
    } else {
      ## Calculate standardized pair-wise distances
      if(!any(positive %in% attributes(pred$predictions)$dimnames[[2]])) {
        warning("Argument positive does not match any class level!")
        positive <- NULL
      }
      if(is.null(positive)){
        warning("Last level has been selected!!")
        positive <- attributes(pred$predictions)$dimnames[[2]][length(attributes(pred$predictions)$dimnames[[2]])]
      }
      for (i in 1:rf$num.trees){
        for (j in 1:rf$num.trees){
          distances[i,j] <- sum((pred$predictions[,,i] - pred$predictions[,,j])^2)/nrow(test_data)
        }
      }
    }
  }
  
  ## Return distance matrix ----
  return(distances)
}


#########################################
### create ranger & rpart baseline models  
set.seed(1896)
rf <- ranger(default ~ ., data = train, write.forest=TRUE, probability = TRUE, importance = "impurity") #, num.trees = 100)

rp <- rpart(default ~ ., data = train)

# visualize tree
rpart.plot(rp, cex = 0.55, faclen = 4, main = "Decision Tree") # clip.facs = , snip = TRUE

# create predictions
predictions <- data.frame(default = valid$default,
                          ranger  =  predict(rf, valid)$predictions[,2], # score = vorhersage good
                          rpart   = predict(rp, valid)[,2]
                          )

### complexity (# rules) 
# ranger:
nrules_rf <- 0
for(i in 1:rf$num.trees){
  tinf <- treeInfo(rf, i)
  nrules_rf <- nrules_rf + nrow(tinf) - sum(tinf$terminal) # number of rules
}
complexity <- nrules_rf # 40272

# rpart:
complexity <- c(complexity, nrow(rp$frame[rp$frame$var != "<leaf>",])) # 16
names(complexity) <- c("ranger", "rpart")


### Performance (AUC): 
library(pROC)
scores <- predictions$ranger 
curve <- roc(predictions$default, scores, levels = c("good","bad"), direction = ">")
                                        # levels = c("controls", "cases"),  direction = controls > cases
aucs <- auc(curve) # 0.7762
set.seed(42)
ci(auc(curve), method = "bootstrap") # [0.7159, 0.8271]

for(j in 3:ncol(predictions)){
  scores <- predictions[,j]
  curve <- roc(predictions$default, scores, levels = c("good","bad"), direction = ">")
  aucs <- c(aucs, auc(curve))
}
names(aucs) <- c("ranger", "rpart")

# compare auc of bnoth models on bootstrap samples
set.seed(42)
nboot  <- 1000
nvalid <- nrow(valid)
diffs  <- numeric(nboot) 
for(i in 1:nboot){
  ids <- sample(nvalid, nboot, replace = TRUE)
  curve1 <- roc(predictions$default[ids], predictions$ranger[ids], levels = c("good","bad"), direction = ">")
  curve2 <- roc(predictions$default[ids], predictions$rpart[ids], levels = c("good","bad"), direction = ">")
  diffs[i] <- auc(curve1) - auc(curve2)
}
plot(density(diffs), main = "AUC(RF) - AUC(Tree)")

# ### Fig1 paper ###
# pdf("Figs/BootTree.pdf", 12, 6)
# par(mfrow=c(1,2))
# plot(density(diffs), main = "AUC(RF) - AUC(Tree)")
# rpart.plot(rp, tweak = 2.0, varlen = 12, faclen = 2, main = "Decision Tree")#, clip.facs = T, snip = TRUE
# dev.off()
# 
# setEPS()
# postscript("Figs/BootTree.eps", width = 12, height = 6)
# par(mfrow=c(1,2))
# plot(density(diffs), main = "AUC(RF) - AUC(Tree)")
# rpart.plot(rp, tweak = 2.0, varlen = 12, faclen = 2, main = "Decision Tree")#, clip.facs = T, snip = TRUE
# dev.off()
# ##########################


### most representative tree(s) 
m <- "splitting variables"
d0 <- timbR::measure_distances(rf = rf, metric = m) 

# m <- "terminal nodes"
# d1 <- measure_distances(rf = rf, metric = m, test_data = train) 
# rf_reptree1 <- select_trees(rf = rf, num.trees = 1, distance.matrix = d1)
# # skipped for reasons of computation time

m <- "prediction"
d2 <- measure_distances(rf = rf, metric = m, test_data = train, positive = "bad") 

# ...according to https://www.biorxiv.org/content/10.1101/2022.05.15.492004v1
m <- "weighted splitting variables"
d3 <- timbR::measure_distances(rf = rf, metric = m) 


### number of rules and explainability of MRT_d0
### (REM: In the paper d0, d1, d2 and d3 are called d1, d2 d3 and d4 respectively.)
# rules
dist_score <- rowSums(d0)
table(dist_score) # 151

# ...for d0: compute average upsilon over all representative trees (no unique tree => no unique prediction)  
tids <- which(dist_score == min(dist_score)) # set.seed(1896); tid <- sample(tid, 1)
nrules   <- NULL
upsilons <- NULL
pred_all <- predict(rf, valid, predict.all = TRUE)

for(tid in tids){
  # complexity
  tinf <- treeInfo(rf, tid)
  nrules <- c(nrules, nrow(tinf) - sum(tinf$terminal))
  # explainability
  ASE <- mean((predictions$ranger - pred_all$predictions[,2,tid])^2)
  ASE0 <- mean((predictions$ranger - mean(predictions$ranger))^2)
  upsilons <- c(upsilons, 1 - ASE / ASE0)
}
mean(nrules)   # 80.54
mean(upsilons) # -2.60

complexity <- c(complexity, mean(nrules)) # 80.5 > forest tree much deeper than simple rpart but worse performance on test data (overfitting the training subset)
names(complexity)[length(complexity)] <- "MRT_d0"
upsilon <- mean(upsilons) 
names(upsilon) <- "MRT_d0"

# ...d2
# rules
dist_score <- rowSums(d2)
tid <- which(dist_score == min(dist_score))

tinf <- treeInfo(rf, tid)
nrules <- nrow(tinf) - sum(tinf$terminal) # 80

# explainability
ASE <- mean((predictions$ranger - pred_all$predictions[,2,tid])^2)
ASE0 <- mean((predictions$ranger - mean(predictions$ranger))^2)
upsilons <- 1 - ASE / ASE0

complexity <- c(complexity, nrules)
names(complexity)[length(complexity)] <- "MRT_d2"
upsilon <- c(upsilon, upsilons)
names(upsilon)[length(upsilon)] <- "MRT_d2"

predictions$MRT_d2 <- pred_all$predictions[,2,tid]

plot(predictions$ranger, predictions$MRT_d2, xlim = c(0,1), ylim = c(0,1),
     xlab = "Predictions Ranger", ylab = "Predictions MRT", pch = 16, cex = 0.6)
abline(h = mean(predictions$ranger), lty = "dotted")
lines(c(0,1), c(0,1), col = "grey")


#...d3
# rules
dist_score <- rowSums(d3)
tid <- which(dist_score == min(dist_score))

tinf <- treeInfo(rf, tid)
nrules <- nrow(tinf) - sum(tinf$terminal) # 78
# explainability
ASE <- mean((predictions$ranger - pred_all$predictions[,2,tid])^2)
ASE0 <- mean((predictions$ranger - mean(predictions$ranger))^2)
upsilons <- 1 - ASE / ASE0

complexity <- c(complexity, nrules)
names(complexity)[length(complexity)] <- "MRT_d3"
upsilon <- c(upsilon, upsilons)
names(upsilon)[length(upsilon)] <- "MRT_d3"

predictions$MRT_d3 <- pred_all$predictions[,2,tid]


### ...most representative GROVEs (cluster of trees)
# create distance object for hierarchical clustering
distance.matrix <- d2
distances <- NULL
for (i in 1:(ncol(distance.matrix) - 1)) distances <- c(distances, distance.matrix[(i + 1):nrow(distance.matrix), i])
attr(distances, "Labels") <- colnames(distance.matrix)
attr(distances, "Size")   <- ncol(distance.matrix)
attr(distances, "Metric") <- "tree dissimilarity"
class(distances) <- "dissimilarity"

dendro <- hclust(distances, method = "ward.D")

plot(dendro, main = "Clusters of Trees", xlab = "Trees", ylab = "Dissimilarity")
plot(length(dendro$height):1, y=dendro$height, main = "Screeplot", 
     xlab = "Number of Representative Trees", ylab = "Dissimilarity", 
     xlim = c(0,40), type = "both", pch = 16)
# 3, 10, 24

# ### Fig3 paper ###
# pdf("Figs/clust.pdf", 12, 6)
# par(mfrow=c(1,2))
# plot(dendro, main = "Clusters of Trees", xlab = "Trees", ylab = "Dissimilarity")
# plot(length(dendro$height):1, y=dendro$height, main = "Screeplot", 
#      xlab = "Number of Representative Trees", ylab = "Dissimilarity", 
#      xlim = c(0,40), type = "both", pch = 16)
# for(j in c(3,10,24)) lines(c(j,j), c(0,rev(dendro$height)[j]), col = "grey", lty = "dotted")
# dev.off()
# 
# setEPS()
# postscript("Figs/clust.eps", width = 12, height = 6)
# par(mfrow=c(1,2))
# plot(dendro, main = "Clusters of Trees", xlab = "Trees", ylab = "Dissimilarity")
# plot(length(dendro$height):1, y=dendro$height, main = "Screeplot", 
#      xlab = "Number of Representative Trees", ylab = "Dissimilarity", 
#      xlim = c(0,40), type = "both", pch = 16)
# for(j in c(3,10,24)) lines(c(j,j), c(0,rev(dendro$height)[j]), col = "grey", lty = "dotted")
# dev.off()
# ##########################

for(num.trees in c(3,10,24)){
  # assign trees to clusters
  clusts <- cutree(dendro, num.trees)
  
  dmats      <- list()
  tree.ids   <- list()
  for(i in seq(along.with = sort(unique(clusts)))){
    tree.ids[[i]] <- which(clusts == i)
    dmats[[i]]    <- distance.matrix[clusts == i, clusts == i]
  }
  
  # select one tree per cluster
  within.cl.dists <- lapply(dmats, function(z){if(length(dim(z))>1){z <- rowSums(z)}; return(z)})
  within.ids <- sapply(within.cl.dists, function(z) which.min(z)[1]) # REM: if several candidates, the first one is picked
  
  ids <- sapply(seq(along.with = within.ids), function(j) tree.ids[[j]][within.ids[j]])
  
  grove <- rf
  grove$num.trees             <- num.trees
  grove$forest$num.trees      <- num.trees
  grove$forest$child.nodeIDs  <- grove$forest$child.nodeIDs[ids]
  grove$forest$split.varIDs   <- grove$forest$split.varIDs[ids]
  grove$forest$split.values   <- grove$forest$split.values[ids]
  # 
  if (length(grove$inbag.counts) > 0) grove$inbag.counts <- grove$inbag.counts[ids]
  grove$predictions      <- NULL
  grove$prediction.error <- NULL
  # #grove
  
  # # the subsequent code does make Rstudio abort due to a fatal error:
  # predictions[[paste0("MRgrove",num.trees)]] =  predict(object = grove, data = valid)$predictions[,2]
  
  # ...instead: compute predictions by hand (by averaging over all trees' predictions):
  predictions_all_grove =  predict(object = rf, data = valid, predict.all = TRUE)$predictions   
  # str(predictions_all_grove)
  predictions_grove <- predictions_all_grove[,,ids]
  
  probs <- predictions_grove[,2,1]
  for(k in 2:dim(predictions_grove)[3]) probs <- probs + predictions_grove[,2,k]
  probs <- probs / dim(predictions_grove)[3]
  predictions[[paste0("MRgrove",num.trees)]] <- probs
  
  nrules_grove <- 0
  for(i in 1:grove$num.trees){
    tinf <- treeInfo(rf, i)
    nrules_grove <- nrules_grove + nrow(tinf) - sum(tinf$terminal) # number of rules
  }
  complexity <- c(complexity, nrules_grove) # 340
  names(complexity)[length(complexity)] <- paste0("MRgrove",num.trees)
  
  ASE <- mean((predictions$ranger - predictions[[paste0("MRgrove",num.trees)]])^2)
  ASE0 <- mean((predictions$ranger - mean(predictions$ranger))^2)
  upsilons <- 1 - ASE / ASE0
  upsilon <- c(upsilon, upsilons)
  names(upsilon)[length(upsilon)] <- paste0("MRgrove",num.trees)
  
}

### visualize upsilon
plot(predictions$ranger, predictions$MRT_d2, xlim = c(0,1), ylim = c(0,1),
     xlab = "Predictions Ranger", ylab = "Predictions MRT", pch = 16, cex = 0.6)
abline(h = mean(predictions$ranger), lty = "dotted")
lines(c(0,1), c(0,1), col = "grey")


# ### Fig2 paper ###
# pdf("Figs/ppplots.pdf", 12, 6)
# par(mfrow=c(1,2))
# plot(predictions$ranger, predictions$MRT_d2, xlim = c(0,1), ylim = c(0,1),
#      xlab = "Predictions Ranger", ylab = "Predictions MRT", pch = 16, cex = 0.6, 
#      main = "Most Representative Tree")
# abline(h = mean(predictions$ranger), lty = "dotted")
# lines(c(0,1), c(0,1), col = "grey")
# 
# plot(predictions$ranger, predictions$MRgrove24, xlim = c(0,1), ylim = c(0,1),
#      xlab = "Predictions Ranger", ylab = "Predictions MR Grove 24", pch = 16, cex = 0.6,
#      main = "Grove of 24 Most Representative Trees")
# abline(h = mean(predictions$ranger), lty = "dotted")
# lines(c(0,1), c(0,1), col = "grey")      
# dev.off()
# 
# setEPS()
# postscript("Figs/ppplots.eps", width = 12, height = 6)
# par(mfrow=c(1,2))
# plot(predictions$ranger, predictions$MRT_d2, xlim = c(0,1), ylim = c(0,1),
#      xlab = "Predictions Ranger", ylab = "Predictions MRT", pch = 16, cex = 0.6, 
#      main = "Most Representative Tree")
# abline(h = mean(predictions$ranger), lty = "dotted")
# lines(c(0,1), c(0,1), col = "grey")
# 
# plot(predictions$ranger, predictions$MRgrove24, xlim = c(0,1), ylim = c(0,1),
#      xlab = "Predictions Ranger", ylab = "Predictions MR Grove 24", pch = 16, cex = 0.6,
#      main = "Grove of 24 Most Representative Trees")
# abline(h = mean(predictions$ranger), lty = "dotted")
# lines(c(0,1), c(0,1), col = "grey")      
# dev.off()


### surrogate tree
train_surrogate         <-  train
train_surrogate$default <- predict(rf, train)$predictions[,2]
rf_surrogate            <- rpart(default~., data = train_surrogate)

rpart.plot(rf_surrogate, cex = 0.55, faclen = 4, main = "Decision Tree") # clip.facs = , snip = TRUE

predictions$surrogate <- predict(rf_surrogate, valid) # auc: 0.7195 ...sogar leicht besser als rp (0.6785)

nrules_rf_surrogate <- nrow(rf_surrogate$frame[rf_surrogate$frame$var != "<leaf>",])
complexity <- c(complexity, nrules_rf_surrogate) # 15 ...entsprechend wieder einfacher
names(complexity)[length(complexity)] <- "surrogate tree"

ASE <- mean((predictions$ranger - predictions$surrogate)^2)
ASE0 <- mean((predictions$ranger - mean(predictions$ranger))^2)
upsilons <- 1 - ASE / ASE0
upsilon <- c(upsilon, upsilons)
names(upsilon)[length(upsilon)] <- "surrogate tree"

plot(predictions$ranger, predictions$surrogate, xlim = c(0,1), ylim = c(0,1),
     xlab = "Predictions Ranger", ylab = "Predictions Surrogate Tree", pch = 16, cex = 0.6)
abline(h = mean(predictions$ranger), lty = "dotted")
lines(c(0,1), c(0,1), col = "grey")      


### surrogate grove (via gbm)
train_surrogate         <-  train 
train_surrogate$default <- predict(rf, train)$predictions[,2]

for(num.trees in c(3,5,7,9,20,30,50,100)){
  set.seed(42)
  surrogate_grove <- gbm(default ~., data = train_surrogate, n.trees = num.trees)
  predictions[[paste0("surrogate_grove", num.trees)]] <- predict(surrogate_grove, valid) 
 
  rules <- NULL 
  for(tid in 1:num.trees){
    tinf  <- pretty.gbm.tree(surrogate_grove, i.tree = tid)
    rules <- rbind(rules, tinf[tinf$SplitVar != -1,])
  }
  if(num.trees == 9) rules9 <- rules # store it for analysis below
  
  complexity <- c(complexity, nrow(rules)) #
  names(complexity)[length(complexity)] <- paste0("surrogate_grove", num.trees)
  
  ASE <- mean((predictions$ranger - predictions[[paste0("surrogate_grove", num.trees)]])^2)
  ASE0 <- mean((predictions$ranger - mean(predictions$ranger))^2)
  upsilons <- 1 - ASE / ASE0
  upsilon <- c(upsilon, upsilons)
  names(upsilon)[length(upsilon)] <- paste0("surrogate_grove", num.trees)
}


### summarize results
ups <- upsilon[-(1:3)] 
complex <- complexity[-(1:5)]
rbind(ups, complex)

cp <- wes_palette("Darjeeling1", 3, type = c("discrete"))
plot(complex, ups, log = "x", pch = 15, ylim = c(0,1),
     xlab = "Complexity", ylab = "Explainability",
     col = cp[c(rep(1,3), 3, rep(2,8))]
     )
legend(x = 2.5, y = 1, legend = c("Surrogate Tree", "MRT Grove 3/10/24", "Surrogate Grove"), fill = cp[c(3,1,2)], bty = "n", cex = 0.75)
abline(v = 10, col = "grey", lty = "dotted")
abline(h = 0.8, col = "grey", lty = "dotted")


# ### Fig4 paper ###
# pdf("Figs/surrogates.pdf", 12, 6)
# par(mfrow=c(1,2))
# rpart.plot(rf_surrogate, cex = 0.55, faclen = 4, main = "Decision Tree") # clip.facs = , snip = TRUE
# 
# cp <- wes_palette("Darjeeling1", 3, type = c("discrete"))
# plot(complex, ups, log = "x", pch = 15, ylim = c(0,1),
#      xlab = "Complexity", ylab = "Explainability",
#      col = cp[c(rep(1,3), 3, rep(2,8))]
# )
# legend(x = 2.5, y = 1, legend = c("Surrogate Tree", "MRT Grove 3/10/24", "Surrogate Grove"), fill = cp[c(3,1,2)], bty = "n", cex = 0.75)
# abline(v = 10, col = "grey", lty = "dotted")
# abline(h = 0.8, col = "grey", lty = "dotted")
# dev.off()
# 
# setEPS()
# postscript("Figs/surrogates.eps", width = 12, height = 6)
# par(mfrow=c(1,2))
# rpart.plot(rf_surrogate, cex = 0.55, faclen = 4, main = "Decision Tree") # clip.facs = , snip = TRUE
# 
# cp <- wes_palette("Darjeeling1", 3, type = c("discrete"))
# plot(complex, ups, log = "x", pch = 15, ylim = c(0,1),
#      xlab = "Complexity", ylab = "Explainability",
#      col = cp[c(rep(1,3), 3, rep(2,8))]
# )
# legend(x = 2.5, y = 1, legend = c("Surrogate Tree", "MRT Grove 3/10/24", "Surrogate Grove"), fill = cp[c(3,1,2)], bty = "n", cex = 0.75)
# abline(v = 10, col = "grey", lty = "dotted")
# abline(h = 0.8, col = "grey", lty = "dotted")
# dev.off()
# ####################


# results <- rbind(ups, complex)
# save(results, rf, rp, rf_surrogate, d0, d2, d3, predictions, train, valid, file = "r_objects.Robj")

### interpret surrogate grove

set.seed(42)
num.trees <- 9
surrogate_grove <- gbm(default ~., data = train_surrogate, n.trees = num.trees)

rules <- NULL 
for(tid in 1:num.trees){
  tinf  <- pretty.gbm.tree(surrogate_grove, i.tree = tid)
  newrule <- tinf[tinf$SplitVar != -1,]
  newrule <- data.frame(newrule, pleft = tinf$Prediction[rownames(tinf) == newrule$LeftNode], pright = tinf$Prediction[rownames(tinf) == newrule$RightNode])
  rules <- rbind(rules, newrule)
}
# results in the same tree as above: rules; rules9

vars   <- NULL
splits <- NULL
csplits_left <- NULL
pleft  <- NULL
pright <- NULL

for(i in 1:nrow(rules)){
  vars   <- c(vars,   names(train_surrogate)[rules$SplitVar[i]+1])
  if(is.numeric(train_surrogate[,rules$SplitVar[i]+1])){
    splits       <- c(splits, rules$SplitCodePred[i])
    csplits_left <- c(csplits_left, NA)
  }
  if(is.factor(train_surrogate[,rules$SplitVar[i]+1])){
    levs <- levels(train_surrogate[,(rules$SplitVar[i]+1)])
    lids <- surrogate_grove$c.splits[[(rules$SplitCodePred[i] +1)]] == -1
    if(sum(lids) == 1) levs <- levs[lids]
    if(sum(lids) > 1)  levs <- paste(levs[lids], sep = "|")
    csl <- levs[1]
    for(j in 2:length(levs)) csl <- paste(csl, levs[j], sep = " | ")
    splits       <- c(splits, "")
    csplits_left <- c(csplits_left, csl)
  }
  
  pleft  <- c(pleft,  rules$pleft[i])
  pright <- c(pright, rules$pright[i])
}

basepred <- surrogate_grove$initF

df <- data.frame(vars, splits, left = csplits_left, pleft = round(pleft, 4), pright = round(pright,4))
# #write.table(df, "221120_export_rules9_gbm.csv"  ,sep = ";")
# xtable::xtable(df, digits = 3)


######################################################################################################

# SplitCodePred	  if the split variable is continuous then this component is the split point. 
#                 If the split variable is categorical then this component contains the index of object$c.split that describes the categorical split. 
#                 If the node is a terminal node then this is the prediction.

# c.splits	A list of all the categorical splits in the collection of trees. 
#           If the trees[[i]] component of a gbm object describes a categorical split then the splitting value will refer to a component of c.splits. 
#           That component of c.splits will be a vector of length equal to the number of levels in the categorical split variable. 
#           -1 indicates left, +1 indicates right, and 0 indicates that the level was not present in the training data.


# https://stackoverflow.com/questions/31296541/understanding-tree-structure-in-r-gbm-package

# SplitVar:	index of which variable is used to split. -1 indicates a terminal node.
# SplitCodePred:	if the split variable is continuous then this component is the split point. If the split variable is categorical then this component contains the index of object$c.split that describes the categorical split. If the node is a terminal node then this is the prediction.
# LeftNode: the index of the row corresponding to the left node.
# RightNode: the index of the row corresponding to the right node.

# In your example:
# Id SplitVar SplitCodePred LeftNode RightNode MissingNode ErrorReduction Weight   Prediction
# 0         9  6.250000e+01        1         2          21      0.6634681   5981  0.005000061

# - means that the root node (indicated by the row number 0) is split by the 9-th split variable 
#   (the numbering of the split variable here starts from 0, 
#   so the split variable is the 10th column in the training set x).

# - SplitCodePred of 6.25 denotes that all points less than 6.25 went to the LeftNode 1 
#   and all points greater than 6.25 went to RightNode 2. 
#   All points that had a missing value in this column were assigned to the MissingNode 21.

# - The ErrorReduction was 0.6634 due to this split and there were 5981 (Weight) in the root node.

# - Prediction of 0.005 denotes the value assigned to all values at this node 
#   before the point was split. 
#   In the case of terminal nodes (or leaves) denoted by -1 in SplitVar, LeftNode, RightNode, and MissingNode, 
#   the Prediction denotes the value predicted for all the points belonging to this leaf node 
#   ...adjusted (times) times the shrinkage. 
