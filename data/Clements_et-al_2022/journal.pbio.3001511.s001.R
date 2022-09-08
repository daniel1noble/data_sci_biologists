
  ---------------------------------------------------------------------------------------
  ---------------------------------------------------------------------------------------
    
  ----------------------------- Supplementary information to ----------------------------
 
  ------------------ META-ANALYSIS REVEALS AN EXTREME "DECLINE EFFECT" ------------------
  ------------------ IN OCEAN ACIDIFICATION IMPACTS ON FISH BEHAVIOUR -------------------
  
  -------- Jeff C. Clements, Josefin Sundin, Timothy D. Clark, Fredrik Jutfelt ----------
    
  ---------------------------------------------------------------------------------------
  ---------------------------------------------------------------------------------------

    
##### Get packages #####
library(pacman)
pacman::p_load(metafor, MCMCglmm, tidyverse, rotl, magrittr, kableExtra, rmarkdown,gridExtra, psych, bindrcpp, pander)
library(BiocManager)
library(ggplot2)
library(viridis)
library(patchwork)

##### Run base code for metaAidR #####
#' @title I^2 Function 
#' @description Function for calculating commonly reported  I^2 measures from MCMCglmm model objects.  See Nakagawa and Santos (2012) for detailed explanation on the various types of I^2
#' @param model The MCMCglmm or metafor model object. Note that if using a metafor model object an observation level random effect must be used to calculate the residual variance. This should be input as "~1|obs" in the random effect list. 
#' @param v The vector of sampling variance for each effect size. 
#' @param sims The number of simulations used for calculating confidence intervals on I^2 estimates for metafor objects.
#' @param phylo A character string with the name of the phylogenetic random effect. Defaults to FALSE meaning that no phylogenetic heritability is calculated. 
#' @param obs A character string with the name of the observation-level random effect in metafor rma.mv models (e.g. "obs", "effectid", "rowid" etc.). The I^2 value returned for this effect refers to the residual among-effect size heterogeneity.
#' @param ME A character string with the name of the sampling error random effect. This is important if one wishes to enter the sampling variance matrix in as a sparse matrix (i.e. 'ginverse' argument) for MCMCglmm. Otherwise, assumed that the 'mev' argument is used. 
#' @return A data.frame containing the relevant I^2 measures along with the 95 percent confidence / credible intervals.
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @references Nakagawa, S. and Santos, E.S.A. (2012) Methodological issues and advances in biological meta-analysis. Evolutionary Ecology, 26:1253-1274.
#' @export

I2 <- function(model, v, ME = FALSE, sims = 1500, phylo = FALSE, obs = FALSE){
  
  if(class(model) != "MCMCglmm" && class(model) != "rma.mv" && class(model) != "rma"){
    stop("The model object is not of class 'MCMCglmm' or 'metafor'")
  }
  
  wi <- 1/v  #weight
  Vw <- sum((wi) * (length(wi) - 1))  / (((sum(wi)^2) - sum((wi)^2)))
  
  if("MCMCglmm" %in% class(model)){
    # Get posterior distribution
    # TO DO : NEED TO MAKE THIS WORK WITH ginverse in MCMCglmm. Added a bit, but needs work.
    if(ME == FALSE){
      post <- model$VCV[,-match(c("sqrt(mev):sqrt(mev).meta"), colnames(model$VCV))]
    } else{
      post <- model$VCV[,-match(ME, colnames(model$VCV))]
    }
    #Calculate total variance
    VT <- rowSums(cbind(post, Vw))
    Vt <- rowSums(post)  # remove Vw
    
    # For each variance component divide by the total variance. Note this needs to be fixed for phylo, but does deal with variable random effects.
    I2_re <- post / VT
    I2_total  <- Vt / VT
    
    if(phylo == FALSE){
      tmpMatrix <- cbind(I2_re, total = I2_total)
    }else{
      I2_phylo <- post[,match(phylo, colnames(post))] / Vt
      tmpMatrix <- cbind(I2_re, I2_phylo, total = I2_total)
    }
    
    mode <- MCMCglmm::posterior.mode(coda::as.mcmc(tmpMatrix))
    CI <- coda::HPDinterval(coda::as.mcmc(tmpMatrix))
    colnames(CI) <- c("2.5% CI", "97.5% CI")
    
    I2_Table <- as.data.frame(cbind(I2_Est. = mode, CI))
    class(I2_Table) <- c("metaAidR", "data.frame")
    
    return(round_df(I2_Table, digits = 4))
  }
  
  if("rma.mv" %in% class(model) | "rma" %in% class(model)){
    # Monte Carlo Simulations
    # From metafor extract the important statistics
    sigma2 <- matrix(model$sigma2, nrow = 1, ncol = length(model$sigma2))
    colnames(sigma2) <- model$s.names
    sigmaN <- model$s.nlevels
    
    if(obs == FALSE){
      stop("Please add the name of the observation-level random effect, obs. If models do not include this, re-run models including (~1|obs) in the random effect list")
    }
    
    #For each variance estimate use Monte Carlo simulation of data
    Sims <- data.frame(mapply(function(x,y) simMonteCarlo(x, y, sims = sims), x = sigma2, y = sigmaN))
    colnames(Sims) <- colnames(sigma2) 
    
    #Calculate total variance
    VT <- rowSums(cbind(Sims, Vw))
    Vt <- rowSums(Sims)  # remove Vw
    
    # For each variance component divide by the total variance. Note this needs to be fixed for phylo, but does deal with variable random effects.
    I2_re <- Sims / VT
    I2_total <- data.frame(Vt / VT)
    
    if(phylo == FALSE){
      tmpMatrix <- data.frame(I2_re[, -match("obs", colnames(I2_re))], total = I2_total)
      names(tmpMatrix) = c(colnames(I2_re)[!colnames(I2_re) %in% 'obs'], 'total')
      
    }else{
      I2_phylo <- Sims[, match(phylo, colnames(sigma2))] / Vt
      tmpMatrix <- cbind(I2_re, phylo = I2_phylo, total = I2_total)
    }
    
    CI <- lapply(tmpMatrix, function(x) stats::quantile(x, c(0.025, 0.975), na.rm = TRUE))
    I_CI <- as.data.frame(do.call(rbind, CI))
    colnames(I_CI) <- c("2.5% CI", "97.5% CI")
    I2_table <- cbind(I2_Est. = colMeans(tmpMatrix), I_CI )
    
    class(I2_table) <- c("metaAidR", "data.frame")
    
    return(round_df(I2_table, digits = 4))
  }
  
}


#' @title Parametric simulation 
#' @description Function for calculating I2 estimates using parametric simulations of model estimates taken from metafor. Note that the effectiveness of these simulations depends on the accuracy of model variance estimates.
#' @param estimate The estimate (i.e. variance) from a metafor model
#' @param sims The number of simulations 
#' @param n The sample size used in estimating the variance 
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @export
simMonteCarlo <- function(estimate, n, sims){
  tmp <- data.frame(num = base::rep(1:sims, each = n), 
                    y = stats::rnorm(n*sims, 0, base::sqrt(estimate)))
  Var <- tmp %>% dplyr::group_by(num) %>% dplyr::summarise(Mean_var = stats::var(y))
  return(as.numeric(Var$Mean_var))
}

## NOTE about PIPE: Run usethis::use_pipe() in the console. The package usethis will add what you need to import the pipe to your NAMESPACE and it will also drop warnings in checks

# Function for  rounding a data frame
round_df <- function(x, digits) {
  numeric_columns <- sapply(x, class) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}


#' @title Folded normal distribution
#' @description Applying the folded normal distribution using the posterior parameter estimates of effect sizes that are normally distributed. This function will allow one to understand the overall magnitude of effect regardless of effect size direction. 
#' @param mu The posterior distribution of mean estimates from an MCMCglmm object. Alternatively, you can give it the mean and sd of a normal distribution and it will provide the mean (only used if type = "raw").
#' @param sd The standard deviation of the posterior distribution of variance-covariance matrix from an MCMCglmm object or just the sd of a normal distribution. 
#' @param type Indicates whether the posterior mean (i.e. "mean") or mode (i.e. "mode") should be returned. Alternatively, if one knows the mean and the sd of you can use type = "raw" to get the mean and variance of the folded normal distribution.
#' @return Posterior mean or mode along with the 95 percent credible intervals (i.e. the highest posterior density interval - HPDinterval). Alternatively if type = "raw" it will return the mean and variance for the folded normal distribution.
#' @references Morrisey, M.B. 2016. Meta-analysis of magnitudes, differences and variation in evolutionary parameters. Journal of Evolutionary Biology, 29:1882-1904
#' @references Morrisey, M.B. 2016. Rejoinder: Further considerations for meta-analysis of transformed quantities such as absolute values. Journal of Evolutionary Biology, 29:1922-1931
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @examples
#' set.seed(60)
#' x <- rnorm(1000, 1, 1)
#' mean(x)
#' sd(x)
#' # Calculate the mean of the folded norm from raw x
#' foldnorm(mu = mean(x), sd = sd(x), type = "raw")
#' #Should be close to: 
#' mean(abs(x))
#' var(abs(x))
#' @export

foldnorm <-  function(mu, sd, type = c("mean", "mode", "raw")){
  
  postfnorm <- stats::dnorm(mu, 0, sd)*2*sd^2 + mu*(2*stats::pnorm(mu, 0, sd) -1)
  
  if(type == "raw"){
    est <- postfnorm
    var.fnorm <- mu^2 + sd^2 - (sd*sqrt(2/pi)*exp((-1*mu^2)/(2*sd^2)) + mu*(1-2*stats::pnorm(-1*mu/sd, 0, 1)))^2
    est <- data.frame(Mean=est, Variance = var.fnorm)
  }
  
  if(type == "mean"){
    mean <- mean(postfnorm)
    CI <- coda::HPDinterval(postfnorm)
    est <- Matrix::cBind(mean = mean, CI)
  }
  
  if(type == "mode"){
    mode <- MCMCglmm::posterior.mode(postfnorm)
    CI <- coda::HPDinterval(postfnorm)
    est <- Matrix::cBind(mode = mode, CI)
  }
  return(est)
}

#' @title Covariance and correlation matrix function basing on shared level ID
#' @description Function for generating simple covariance and correlation matrices 
#' @param data Dataframe object containing effect sizes, their variance, unique IDs and clustering variable
#' @param V Name of the variable (as a string – e.g, "V1") containing effect size variances variances
#' @param cluster Name of the variable (as a string – e.g, "V1") indicating which effects belong to the same cluster. Same value of 'cluster' are assumed to be nonindependent (correlated).
#' @param obs Name of the variable (as a string – e.g, "V1") containing individual IDs for each value in the V (Vector of variances). If this parameter is missing, label will be labelled with consecutive integers starting from 1.
#' @param rho Known or assumed correlation value among effect sizes sharing same 'cluster' value. Default value is 0.5.
#' @param type Optional logical parameter indicating whether a full variance-covariance matrix (default or "vcv") is needed or a correlation matrix ("cor") for the non-independent blocks of variance values.
#' @export

make_VCV_matrix <- function(data, V, cluster, obs, type=c("vcv", "cor"), rho=0.5){
  type <- match.arg(type)
  if (missing(data)) {
    stop("Must specify dataframe via 'data' argument.")
  }
  if (missing(V)) {
    stop("Must specify name of the variance variable via 'V' argument.")
  }
  if (missing(cluster)) {
    stop("Must specify name of the clustering variable via 'cluster' argument.")
  }
  if (missing(obs)) {
    obs <- 1:length(V)   
  }
  if (missing(type)) {
    type <- "vcv" 
  }
  
  new_matrix <- matrix(0,nrow = dim(data)[1],ncol = dim(data)[1]) #make empty matrix of the same size as data length
  rownames(new_matrix) <- data[ ,obs]
  colnames(new_matrix) <- data[ ,obs]
  # find start and end coordinates for the subsets
  shared_coord <- which(data[ ,cluster] %in% data[duplicated(data[ ,cluster]), cluster]==TRUE)
  # matrix of combinations of coordinates for each experiment with shared control
  combinations <- do.call("rbind", tapply(shared_coord, data[shared_coord,cluster], function(x) t(utils::combn(x,2))))
  
  if(type == "vcv"){
    # calculate covariance values between  values at the positions in shared_list and place them on the matrix
    for (i in 1:dim(combinations)[1]){
      p1 <- combinations[i,1]
      p2 <- combinations[i,2]
      p1_p2_cov <- rho * sqrt(data[p1,V]) * sqrt(data[p2,V])
      new_matrix[p1,p2] <- p1_p2_cov
      new_matrix[p2,p1] <- p1_p2_cov
    }
    diag(new_matrix) <- data[ ,V]   #add the diagonal
  }
  
  if(type == "cor"){
    # calculate covariance values between  values at the positions in shared_list and place them on the matrix
    for (i in 1:dim(combinations)[1]){
      p1 <- combinations[i,1]
      p2 <- combinations[i,2]
      p1_p2_cov <- rho
      new_matrix[p1,p2] <- p1_p2_cov
      new_matrix[p2,p1] <- p1_p2_cov
    }
    diag(new_matrix) <- 1   #add the diagonal of 1
  }
  
  return(new_matrix)
}


#Value: Labelled full variance-covariance or correlation matrice of the size and labels matching initial dataframe will be returned 

## TO DO: Implement a shared control equations for SMD and lnRR, maybe similar to eho = 0.5, but we have exact calculations so probably worth implementing

#' @title es_stat: Function for calculating effect size statistics
#' @description Function for calculating common effect size statistics.  Effect size statistics include both the effect size itself and its sampling error. This function estimates both commonly used effect sizes (Hedges' d & g, Zr, lnOR, lnRR) along with an ability to convert between the different effect size statistics. For two-group comparisons, directionally is dependent on the group used for m1 or p1. In all cases, m1 - m2 or p1 - p2 is contrasted.
#' @param m1 Mean of group 1
#' @param m2 Mean of group 2
#' @param sd1 Standard deviation of group 1
#' @param sd2 Standard deviation of group 2
#' @param n1 Sample size of group 1
#' @param n2 Sample size of group 2
#' @param p1 Proportion in group 1. Used with type = "lnOR" only.
#' @param p2 Proportion in group 2. Used with type = "lnOR" only.
#' @param r Correlation coefficient. Used with type = "r" only .
#' @param nr Sample size used for estimating the correlation coefficient. Used with type "r" only.
#' @param type The type specifies the specific effect statistics one wishes to calculate. Types include: "d" = Hedges' d, "g" = bias corrected Hedges' d, "Zr" = Fisher's z-transformed correlation coefficient,"lnOR" = log odds ratio.
#' @return Function returns the effect size and its sampling variance in a matrix (two column, n rows). The arguments can be vectors. 
#' @author Daniel Noble - daniel.noble@unsw.edu.au
#' @export

es_stat <- function(m1, m2, sd1, sd2, n1, n2, p1, p2, r, nr, type = c("d", "g", "Zr","lnOR")){
  
  if(type == "g"){
    g     <- hedge(m1, m2, sd1, sd2, n1, n2, type = "g")
    v_g <- var_hedge(g, n1, n2, type = "g")
    return(cbind(g, v_g))
  }
  
  if(type == "d"){
    d     <- hedge(m1, m2, sd1, sd2, n1, n2, type = "d")
    v_d <- var_hedge(d, n1, n2, type = "d")	
    return(cbind(d, v_d))
  }
  
  if(type == "Zr"){
    return(Zr_es(r, nr))
  }
  
  if(type == "lnOR"){
    return(lOR_es(p1, p2, n1, n2))
  }
  
}

#' @title es_ratio: Function for calculating ratio-based effect size statistics
#' @description  This function estimates less commonly used (lnHR) or newly developed effect sizes (i.e. for variance: lnCVR, lnVR, lnSD). Effect size statistics include both the effect size itself and its sampling error. For two-group comparisons, directionally is dependent on the group used for m1 or d1. In all cases, m1 - m2 or d1-d2 is contrasted.
#' @param m1 Mean of group 1
#' @param m2 Mean of group 2
#' @param sd1 Standard deviation of group 1
#' @param sd2 Standard deviation of group 2
#' @param n1 Sample size of group 1
#' @param n2 Sample size of group 2
#' @param d1 Number of events (e.g. deaths) in time interval t-1 to t for group 1. Used for type "lnHR" only.
#' @param d2 Number of events (e.g. deaths) in time interval t-1 to t for group 2. Used for type "lnHR" only.
#' @param n1HR Number at risk for group 1 during time interval t-1 to t. Used for type "lnHR" only.
#' @param n2HR Number at risk for group 2 during time interval t-1 to t. Used for type "lnHR" only.
#' @param cor Assumed correlation between mean and variance. Used when a vector of statistics of length 1 is provided. Used only with type = "lnCVR".
#' @param equal.corr Logical indicating whether the correlation between log mean and log sd are assumed to be the same (TRUE) or different (FALSE). Used only with type = "lnCVR".
#' @param type The type specifies the specific effect statistics one wishes to calculate. Types include: "lnHR" = log hazards ratio, "lnRR" = log response ratio, "lnCVR" = log coefficient of variation ratio, "lnVR" = log variance ratio. 
#' @return Function returns the effect size and its sampling variance in a matrix (two column, n rows). The arguments can be vectors. 
#' @author Daniel Noble - daniel.noble@unsw.edu.au
#' @export
es_ratio <- function(m1, m2, sd1, sd2, n1, n2, d1, d2, n1HR, n2HR, cor = cor, equal.corr=TRUE, type = c("lnRR", "lnCVR", "lnVR", "lnHR")){
  if(type == "lnRR"){
    return(lnRR_es(m1, m2, sd1, sd2, n1, n2))
  }
  
  if(type == "lnHR"){
    return(lnHR_es(d1, d2, n1, n2))
  }
  
  if(type == "lnVR"){
    return(lnVR_es(sd1, sd2, n1, n2))
  }
  
  if(type == "lnCVR"){
    lnCVR <- lnCVR_es(m1, m2, sd1, sd2, n1, n2)
    v_lnCVR <-var_lnCVR(m1, m2, sd1, sd2, n1, n2, cor = cor, Equal.E.C.Corr = equal.corr)
    return(cbind(lnCVR, v_lnCVR))
  }
}


#' @title hedge
#' @description Function for calculating Hedges' d or biased corrected Hedges' g. m1 is subtracted from m2.
#' @param m1 Mean of group 1
#' @param m2 Mean of group 2
#' @param sd1 Standard deviation of group 1
#' @param sd2 Standard deviation of group 2
#' @param n1 Sample size of group 1
#' @param n2 Sample size of group 2
#' @param type  Whether Hedges' d (type = "d") and biased corrected Hedges' d or Hedges' g (type = "g") is the be used. 
#' @references Borenstein et al. 2009. Effect sizes based on means. In: Introduction to Meta-Analysis (eds. Borenstein, M., Hedges, L.V., Higgins, J.P.T. and Rothstein, H.R.). pg 21-32.
#' @author Daniel Noble - daniel.noble@unsw.edu.au
#' @export

hedge <- function(m1, m2,  sd1, sd2, n1, n2, type = c("d", "g")){
  if(type == "g"){
    J <- 1 - (3 / ((4 * (n1 + n2 - 2)) -1))
    sd_pool<- (((n1 - 1) * (sd1^2)) + ((n2 - 1) * (sd2^2)))/(n1 + n2 - 2)
    h_g <- ((m1 - m2) / sqrt(sd_pool))*J
  }
  
  if(type == "d"){
    sd_pool<- (((n1 - 1) * (sd1^2)) + ((n2 - 1) * (sd2^2)))/(n1 + n2 - 2)
    h_g <- ((m1 - m2) / sqrt(sd_pool))
  }
  return(h_g)
}

#' @title var_hedge
#' @description Function for calculating the sampling variance for Hedges' d or biased corrected Hedges' g. 
#' @param hedge Calculation of Hedges g or d
#' @param n1 Sample size of group 1
#' @param n2 Sample size of group 2
#' @param type  Whether Hedges' d (type = "d") and biased corrected Hedges' d or Hedges' g (type = "g") is the be used. 
#' @author Daniel Noble - daniel.noble@unsw.edu.au
#' @export
var_hedge <- function(hedge, n1, n2, type = c("d", "g")){
  if(type == "g"){
    J <- 1 - (3 / ((4 * (n1 + n2 - 2)) -1))
    vg <- ((n1 + n2)/(n1 * n2)) + ((hedge^2)/(2 * (n1 + n2 - 2)))
    v <- vg * (J^2)
  }
  
  if(type == "d"){
    v <- ((n1 + n2)/(n1 * n2)) + ((hedge^2)/(2 * (n1 + n2 - 2)))
  }
  return(v)
}	

#' @title Zr_es
#' @description Function for calculating Fisher's z-transformed correlation. 
#' @param r Correlation coefficient
#' @param n Sample size
#' @author Daniel Noble - daniel.noble@unsw.edu.au
#' @export

Zr_es <- function(r, n){
  z    <- 0.5*(log((1+r)/(1-r)))
  Vz <- 1 / (n -3)
  return(cbind(z, Vz))
}

#' @title lnOR
#' @description Function for calculating log odds ratio effect size statistics
#' @param p1 Proportion in group 1
#' @param p2 Proportion in group 2
#' @param n1 Sample size in group 1
#' @param n2 Sample size in group 2
#' @author Daniel Noble - daniel.noble@unsw.edu.au
#' @export

lOR_es <- function(p1, p2, n1, n2){
  lnOR    <- log(p1 / (1-p1)) - log(p2 / (1-p2))
  v_lnOR   <- (1 / (p1*n1)) + (1/ ((1-p1)*n1)) + (1 / (p2*n2)) + (1/ ((1-p2)*n2))
  return(cbind(lnOR, v_lnOR))
}

#' @title lnRR
#' @description Function for calculating log response ratio
#' @param m1 Mean of group 1
#' @param m2 Mean of group 2
#' @param sd1 Standard deviation of group 1
#' @param sd2 Standard deviation of group 2
#' @param n1 Sample size of group 1
#' @param n2 Sample size of group 2
#' @author Daniel Noble - daniel.noble@unsw.edu.au
#' @export

lnRR_es <- function(m1, m2,  sd1, sd2, n1, n2){
  sd_pool<- (((n1 - 1) * (sd1^2)) + ((n2 - 1) * (sd2^2)))/(n1 + n2 - 2)
  lnRR <- log(m1) - log(m2)
  VlnRR <- (sd_pool^2)*((1/(n1*(m1^2)) + 1/(n2*(m2^2))))
  return(cbind(lnRR,VlnRR))
}

#' @title lnCVR
#' @description Function for calculating log coefficient of variation ratio
#' @param m1 Mean of group 1
#' @param m2 Mean of group 2
#' @param sd1 Standard deviation of group 1
#' @param sd2 Standard deviation of group 2
#' @param n1 Sample size of group 1
#' @param n2 Sample size of group 2
#' @references Nakagawa et al. (2015) Meta-analysis of variation: ecological and evolutionary applications and beyond. Methods in Ecology and Evolution, 6:143-152.
#' @author Alistair Senior - Alistair.senior@sydney.edu.au
#' @export

lnCVR_es<-function(m1, m2, sd1, sd2, n1, n2){
  lnCVR <-(log(sd1) - log(m1) + 1 / (2*(n1 - 1))) - (log(sd2) - log(m2) + 1 / (2*(n2 - 1)))
  return(lnCVR)	
}


#' @title Sampling variance for lnCVR
#' @description Function for calculating sampling variance of log coefficient of variation ratio under different assumptions about mean-variance relationships. 
#' @param m1 Mean of group 1
#' @param m2 Mean of group 2
#' @param sd1 Standard deviation of group 1
#' @param sd2 Standard deviation of group 2
#' @param n1 Sample size of group 1
#' @param n2 Sample size of group 2
#' @param cor Assumed correlation between mean and variance. Used when a vector of statistics of length 1 is provided. 
#' @param Equal.E.C.Corr Logical indicating whether the the correlation between mean and variance is equal in both group groups.
#' @author Alistair Senior - Alistair.senior@sydney.edu.au
#' @references Nakagawa et al. (2015) Meta-analysis of variation: ecological and evolutionary applications and beyond. Methods in Ecology and Evolution, 6:143-152.
#' @export
var_lnCVR<-function(m1, m2, sd1, sd2, n1, n2, cor, Equal.E.C.Corr=TRUE){
  #To Do: Function depends on a vector being provided. But situations were only a single valued vector used. Need to solve correlation between log(sd) and log(mean) in these situations.
  if(length(m1) <= 1 & length(sd1) <= 1){
    #Assume correlation is 0.5
    S2<- sd1^2 / (n1 * (m1^2)) + 1 / (2 * (n1 - 1)) - 2 * cor * sqrt((sd1^2 / (n1 * (m1^2))) * (1 / (2 * (n1 - 1)))) + sd2^2 / (n2 * (m2^2)) + 1 / (2 * (n2 - 1)) - 2 * cor * sqrt((sd2^2 / (n2 * (m2^2))) * (1 / (2 * (n2 - 1))))
  }
  
  if(Equal.E.C.Corr==TRUE & length(m1) > 1 & length(sd1) > 1){
    mvcorr<-stats::cor.test(log(c(m1, m2)), log(c(sd1, sd2)))$estimate
    
    S2<- sd1^2 / (n1 * (m1^2)) + 1 / (2 * (n1 - 1)) - 2 * mvcorr * sqrt((sd1^2 / (n1 * (m1^2))) * (1 / (2 * (n1 - 1)))) + sd2^2 / (n2 * (m2^2)) + 1 / (2 * (n2 - 1)) - 2 * mvcorr * sqrt((sd2^2 / (n2 * (m2^2))) * (1 / (2 * (n2 - 1))))
  }
  
  if(Equal.E.C.Corr==FALSE & length(m1) > 1 & length(sd1) > 1){
    Cmvcorr<-stats::cor.test(log(m1), log(sd1))$estimate
    Emvcorr<-stats::cor.test(log(m2), (sd2))$estimate
    
    S2<- sd1^2 / (n1 * (m1^2)) + 1 / (2 * (n1 - 1)) - 2 * Cmvcorr * sqrt((sd1^2 / (n1 * (m1^2))) * (1 / (2 * (n1 - 1)))) + sd2^2 / (n2 * (m2^2)) + 1 / (2 * (n2 - 1)) - 2 * Emvcorr * sqrt((sd2^2 / (n2 * (m2^2))) * (1 / (2 * (n2 - 1))))
  }	
  return(S2)	
}

#' @title lnSD_es
#' @description Function for calculating log standard deviation effect size statistics
#' @param sd Standard deviation of group
#' @param n Sample size of group
#' @references Nakagawa et al. (2015) Meta-analysis of variation: ecological and evolutionary applications and beyond. Methods in Ecology and Evolution, 6:143-152. 
#' @author Alistair Senior - Alistair.senior@sydney.edu.au
#' @export
lnSD_es<-function(sd, n){
  lnSD <- log(sd) + (1 / (2 * (n - 1)))
  v_lnSD <- (1 / (2 * (n - 1)))	
  return(cbind(lnSD, v_lnSD))
}

#' @title lnVR_es
#' @description Function for calculating log variation ratio effect size statistics.
#' @param sd1 Standard deviation of group 1
#' @param sd2 Standard deviation of group 2
#' @param n1 Sample size of group 1
#' @param n2 Sample size of group 2
#' @author Alistair Senior - Alistair.senior@sydney.edu.au
#' @references Nakagawa et al. (2015) Meta-analysis of variation: ecological and evolutionary applications and beyond. Methods in Ecology and Evolution, 6:143-152.
#' @export
lnVR_es<-function(sd1, sd2, n1, n2){
  lnVR <- log(sd1 / sd2) + (1 / (2 * (n1 - 1))) - (1 / (2 * (n2 - 1)))
  v_lnVR <- (1 / (2 * (n1 - 1))) + (1 / (2 * (n2 - 1)))
  return(cbind(lnVR, v_lnVR))
}


#' @title lnHR_es
#' @description Function for calculating log hazards ratio effect size statistics, for a specified period in time-to-event or survival curve data.
#' @param d1 Number of events (e.g. deaths) in time interval t-1 to t for group 1
#' @param d2 Number of events (e.g. deaths) in time interval t-1 to t for group 2
#' @param n1HR Number at risk for group 1 during time interval t-1 to t
#' @param n2HR Number at risk for group 2 during time interval t-1 to t
#' @author Daniel Noble - daniel.noble@unsw.edu.au
#' @references Williamson PR, Smith CT, Hutton JL,  Marson AG (2002). Aggregate data meta-analysis with time-to-event outcomes. Stat. Med. 21, 3337-3351.
#' @references Parmar MKB, Torri V,  Stewart L (1998). Extracting summary statistics to perform meta-analyses of the published literature for survival endpoints. Stat. Med. 17, 2815-2834.
#' @export
lnHR_es <- function(d1, d2, n1HR, n2HR){
  et <- (d1 + d2) * ((n1HR*n2HR) / (n1HR + n2HR))
  wt <- (d1 + d2) * ((n1HR*n2HR) / (n1HR + n2HR)^2)
  
  lnHRt <- (d1 - et) / wt
  var_lnHRt <- 1 / wt
  return(cbind(lnHRt, var_lnHRt))
}

##### META-ANALYSIS - YEAR ONLINE - FULL DATASET #####
##attach dataset
decline<-read.csv(file.choose()) ##use dataset "S5 Data"
attach(decline)

##set factors
decline$year.online<-as.factor(decline$year.online)
decline$year.print<-as.factor(decline$year.print)
decline$obs<-as.factor(decline$obs)
decline$study<-as.factor(decline$study)

##view summary
summary(decline)

##subset by year
y2009 <- filter(decline, year.online == "2009")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline))]
y2009$obs <- 1:nrow(y2009)
y2010 <- filter(decline, year.online == "2010")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline))]
y2010$obs <- 1:nrow(y2010)
y2011 <- filter(decline, year.online == "2011")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline))]
y2011$obs <- 1:nrow(y2011)
y2012 <- filter(decline, year.online == "2012")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline))]
y2012$obs <- 1:nrow(y2012)
y2013 <- filter(decline, year.online == "2013")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline))]
y2013$obs <- 1:nrow(y2013)
y2014 <- filter(decline, year.online == "2014")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline))]
y2014$obs <- 1:nrow(y2014)
y2015 <- filter(decline, year.online == "2015")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline))]
y2015$obs <- 1:nrow(y2015)
y2016 <- filter(decline, year.online == "2016")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline))]
y2016$obs <- 1:nrow(y2016)
y2017 <- filter(decline, year.online == "2017")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline))]
y2017$obs <- 1:nrow(y2017)
y2018 <- filter(decline, year.online == "2018")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline))]
y2018$obs <- 1:nrow(y2018)
y2019 <- filter(decline, year.online == "2019")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline))]
y2019$obs <- 1:nrow(y2019)


##comupte effect sizes for each year
lnRR2009 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2009,append=TRUE)
lnRR2010 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2010,append=TRUE)
lnRR2011 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2011,append=TRUE)
lnRR2012 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2012,append=TRUE)
lnRR2013 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2013,append=TRUE)
lnRR2014 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2014,append=TRUE)
lnRR2015 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2015,append=TRUE)
lnRR2016 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2016,append=TRUE)
lnRR2017 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2017,append=TRUE)
lnRR2018 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2018,append=TRUE)
lnRR2019 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2019,append=TRUE)

##remove NAs
lnRR2009clean<-na.omit(lnRR2009)
lnRR2010clean<-na.omit(lnRR2010)
lnRR2011clean<-na.omit(lnRR2011)
lnRR2012clean<-na.omit(lnRR2012)
lnRR2013clean<-na.omit(lnRR2013)
lnRR2014clean<-na.omit(lnRR2014)
lnRR2015clean<-na.omit(lnRR2015)
lnRR2016clean<-na.omit(lnRR2016)
lnRR2017clean<-na.omit(lnRR2017)
lnRR2018clean<-na.omit(lnRR2018)
lnRR2019clean<-na.omit(lnRR2019)

##view mean-variance relationship
pp1<-ggplot(decline,aes(x=log(ctrl.mean),y=log(ctrl.sd),col=year.online))+ geom_point(size=2,na.rm = TRUE)
pp2<-ggplot(decline,aes(x=log(oa.mean),y=log(oa.sd),col=year.online))+ geom_point(size=2,na.rm = TRUE)
grid.arrange(pp1,pp2, nrow =1)

##look at lnRR by year
MLMA_2009_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2009)
summary(MLMA_2009_lnRR)
MLMA_2010_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2010)
summary(MLMA_2010_lnRR)
MLMA_2011_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2011)
summary(MLMA_2011_lnRR)
MLMA_2012_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2012)
summary(MLMA_2012_lnRR)
MLMA_2013_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2013)
summary(MLMA_2013_lnRR)
MLMA_2014_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2014)
summary(MLMA_2014_lnRR)
MLMA_2015_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2015)
summary(MLMA_2015_lnRR)
MLMA_2016_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2016)
summary(MLMA_2016_lnRR)
MLMA_2017_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2017)
summary(MLMA_2017_lnRR)
MLMA_2018_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2018)
summary(MLMA_2018_lnRR)
MLMA_2019_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2019)
summary(MLMA_2019_lnRR)

##set prior
prior <- list(R=list(V = 1, nu =0.002), G = list(G = list(V=1, nu = 0.002)))

##run bayesian MLMA models
model_magnitude_bayes_2009 <- MCMCglmm(yi ~ 1, mev = lnRR2009clean$vi, random = ~obs, data = lnRR2009clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2010 <- MCMCglmm(yi ~ 1, mev = lnRR2010clean$vi, random = ~obs, data = lnRR2010clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2011 <- MCMCglmm(yi ~ 1, mev = lnRR2011clean$vi, random = ~obs, data = lnRR2011clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2012 <- MCMCglmm(yi ~ 1, mev = lnRR2012clean$vi, random = ~obs, data = lnRR2012clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2013 <- MCMCglmm(yi ~ 1, mev = lnRR2013clean$vi, random = ~obs, data = lnRR2013clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2014 <- MCMCglmm(yi ~ 1, mev = lnRR2014clean$vi, random = ~obs, data = lnRR2014clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2015 <- MCMCglmm(yi ~ 1, mev = lnRR2015clean$vi, random = ~obs, data = lnRR2015clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2016 <- MCMCglmm(yi ~ 1, mev = lnRR2016clean$vi, random = ~obs, data = lnRR2016clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2017 <- MCMCglmm(yi ~ 1, mev = lnRR2017clean$vi, random = ~obs, data = lnRR2017clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2018 <- MCMCglmm(yi ~ 1, mev = lnRR2018clean$vi, random = ~obs, data = lnRR2018clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2019 <- MCMCglmm(yi ~ 1, mev = lnRR2019clean$vi, random = ~obs, data = lnRR2019clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)

##get model summaries
summary(model_magnitude_bayes_2009)
summary(model_magnitude_bayes_2010)
summary(model_magnitude_bayes_2011)
summary(model_magnitude_bayes_2012)
summary(model_magnitude_bayes_2013)
summary(model_magnitude_bayes_2014)
summary(model_magnitude_bayes_2015)
summary(model_magnitude_bayes_2016)
summary(model_magnitude_bayes_2017)
summary(model_magnitude_bayes_2018)
summary(model_magnitude_bayes_2019)

##extract posteriors
sol2009 <- model_magnitude_bayes_2009$Sol
VCV2009 <- model_magnitude_bayes_2009$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2009$VCV))]
sol2010 <- model_magnitude_bayes_2010$Sol 
VCV2010 <- model_magnitude_bayes_2010$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2010$VCV))]
sol2011 <- model_magnitude_bayes_2011$Sol 
VCV2011 <- model_magnitude_bayes_2011$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2011$VCV))]
sol2012 <- model_magnitude_bayes_2012$Sol 
VCV2012 <- model_magnitude_bayes_2012$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2012$VCV))]
sol2013 <- model_magnitude_bayes_2013$Sol 
VCV2013 <- model_magnitude_bayes_2013$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2013$VCV))]
sol2014 <- model_magnitude_bayes_2014$Sol 
VCV2014 <- model_magnitude_bayes_2014$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2014$VCV))]
sol2015 <- model_magnitude_bayes_2015$Sol 
VCV2015 <- model_magnitude_bayes_2015$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2015$VCV))]
sol2016 <- model_magnitude_bayes_2016$Sol 
VCV2016 <- model_magnitude_bayes_2016$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2016$VCV))]
sol2017 <- model_magnitude_bayes_2017$Sol 
VCV2017 <- model_magnitude_bayes_2017$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2017$VCV))]
sol2018 <- model_magnitude_bayes_2018$Sol 
VCV2018 <- model_magnitude_bayes_2018$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2018$VCV))]
sol2019 <- model_magnitude_bayes_2019$Sol 
VCV2019 <- model_magnitude_bayes_2019$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2019$VCV))]

##get folded normal function
mu.fnorm <- function(mu, sigma){dnorm(mu, 0, sigma)*2*sigma^2 + mu*(2*pnorm(mu, 0, sigma) -1)}

##get magnitude means + variance
magnitude_mean_2009 <- mu.fnorm(sol2009[,1], sqrt(rowSums(VCV2009)))
magnitude_2009_mean <- data.frame(mean_mag = mean(magnitude_mean_2009), L_CI = HPDinterval(magnitude_mean_2009)[1], U_CI = HPDinterval(magnitude_mean_2009)[2])
magnitude_mean_2010 <- mu.fnorm(sol2010[,1], sqrt(rowSums(VCV2010)))
magnitude_2010_mean <- data.frame(mean_mag = mean(magnitude_mean_2010), L_CI = HPDinterval(magnitude_mean_2010)[1], U_CI = HPDinterval(magnitude_mean_2010)[2])
magnitude_mean_2011 <- mu.fnorm(sol2011[,1], sqrt(rowSums(VCV2011)))
magnitude_2011_mean <- data.frame(mean_mag = mean(magnitude_mean_2011), L_CI = HPDinterval(magnitude_mean_2011)[1], U_CI = HPDinterval(magnitude_mean_2011)[2])
magnitude_mean_2012 <- mu.fnorm(sol2012[,1], sqrt(rowSums(VCV2012)))
magnitude_2012_mean <- data.frame(mean_mag = mean(magnitude_mean_2012), L_CI = HPDinterval(magnitude_mean_2012)[1], U_CI = HPDinterval(magnitude_mean_2012)[2])
magnitude_mean_2013 <- mu.fnorm(sol2013[,1], sqrt(rowSums(VCV2013)))
magnitude_2013_mean <- data.frame(mean_mag = mean(magnitude_mean_2013), L_CI = HPDinterval(magnitude_mean_2013)[1], U_CI = HPDinterval(magnitude_mean_2013)[2])
magnitude_mean_2014 <- mu.fnorm(sol2014[,1], sqrt(rowSums(VCV2014)))
magnitude_2014_mean <- data.frame(mean_mag = mean(magnitude_mean_2014), L_CI = HPDinterval(magnitude_mean_2014)[1], U_CI = HPDinterval(magnitude_mean_2014)[2])
magnitude_mean_2015 <- mu.fnorm(sol2015[,1], sqrt(rowSums(VCV2015)))
magnitude_2015_mean <- data.frame(mean_mag = mean(magnitude_mean_2015), L_CI = HPDinterval(magnitude_mean_2015)[1], U_CI = HPDinterval(magnitude_mean_2015)[2])
magnitude_mean_2016 <- mu.fnorm(sol2016[,1], sqrt(rowSums(VCV2016)))
magnitude_2016_mean <- data.frame(mean_mag = mean(magnitude_mean_2016), L_CI = HPDinterval(magnitude_mean_2016)[1], U_CI = HPDinterval(magnitude_mean_2016)[2])
magnitude_mean_2017 <- mu.fnorm(sol2017[,1], sqrt(rowSums(VCV2017)))
magnitude_2017_mean <- data.frame(mean_mag = mean(magnitude_mean_2017), L_CI = HPDinterval(magnitude_mean_2017)[1], U_CI = HPDinterval(magnitude_mean_2017)[2])
magnitude_mean_2018 <- mu.fnorm(sol2018[,1], sqrt(rowSums(VCV2018)))
magnitude_2018_mean <- data.frame(mean_mag = mean(magnitude_mean_2018), L_CI = HPDinterval(magnitude_mean_2018)[1], U_CI = HPDinterval(magnitude_mean_2018)[2])
magnitude_mean_2019 <- mu.fnorm(sol2019[,1], sqrt(rowSums(VCV2019)))
magnitude_2019_mean <- data.frame(mean_mag = mean(magnitude_mean_2019), L_CI = HPDinterval(magnitude_mean_2019)[1], U_CI = HPDinterval(magnitude_mean_2019)[2])

##view ES magnitudes and uncertainty 
magnitude_2009_mean
magnitude_2010_mean
magnitude_2011_mean
magnitude_2012_mean
magnitude_2013_mean
magnitude_2014_mean
magnitude_2015_mean
magnitude_2016_mean
magnitude_2017_mean
magnitude_2018_mean
magnitude_2019_mean


##### META-ANALYSIS - YEAR ONLINE - WARM WATER ONLY #####
##attach dataset
decline_climate<-read.csv(file.choose()) ##use dataset "S6 Data"
attach(decline_climate)

##set factors
decline_climate$year.online<-as.factor(decline_climate$year.online)
decline_climate$year.print<-as.factor(decline_climate$year.print)
decline_climate$obs<-as.factor(decline_climate$obs)
decline_climate$study<-as.factor(decline_climate$study)

##view summary
summary(decline_climate)

##subset by year
y2009 <- filter(decline_climate, year.online == "2009")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline_climate))]
y2009$obs <- 1:nrow(y2009)
y2010 <- filter(decline_climate, year.online == "2010")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline_climate))]
y2010$obs <- 1:nrow(y2010)
y2011 <- filter(decline_climate, year.online == "2011")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline_climate))]
y2011$obs <- 1:nrow(y2011)
y2012 <- filter(decline_climate, year.online == "2012")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline_climate))]
y2012$obs <- 1:nrow(y2012)
y2013 <- filter(decline_climate, year.online == "2013")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline_climate))]
y2013$obs <- 1:nrow(y2013)
y2014 <- filter(decline_climate, year.online == "2014")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline_climate))]
y2014$obs <- 1:nrow(y2014)
y2015 <- filter(decline_climate, year.online == "2015")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline_climate))]
y2015$obs <- 1:nrow(y2015)
y2016 <- filter(decline_climate, year.online == "2016")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline_climate))]
y2016$obs <- 1:nrow(y2016)
y2017 <- filter(decline_climate, year.online == "2017")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline_climate))]
y2017$obs <- 1:nrow(y2017)
y2018 <- filter(decline_climate, year.online == "2018")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline_climate))]
y2018$obs <- 1:nrow(y2018)
y2019 <- filter(decline_climate, year.online == "2019")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline_climate))]
y2019$obs <- 1:nrow(y2019)


##comupte effect sizes for each year
lnRR2009 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2009,append=TRUE)
lnRR2010 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2010,append=TRUE)
lnRR2011 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2011,append=TRUE)
lnRR2012 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2012,append=TRUE)
lnRR2013 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2013,append=TRUE)
lnRR2014 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2014,append=TRUE)
lnRR2015 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2015,append=TRUE)
lnRR2016 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2016,append=TRUE)
lnRR2017 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2017,append=TRUE)
lnRR2018 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2018,append=TRUE)
lnRR2019 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2019,append=TRUE)

##remove NAs
lnRR2009clean<-na.omit(lnRR2009)
lnRR2010clean<-na.omit(lnRR2010)
lnRR2011clean<-na.omit(lnRR2011)
lnRR2012clean<-na.omit(lnRR2012)
lnRR2013clean<-na.omit(lnRR2013)
lnRR2014clean<-na.omit(lnRR2014)
lnRR2015clean<-na.omit(lnRR2015)
lnRR2016clean<-na.omit(lnRR2016)
lnRR2017clean<-na.omit(lnRR2017)
lnRR2018clean<-na.omit(lnRR2018)
lnRR2019clean<-na.omit(lnRR2019)

##view mean-variance relationship
pp1<-ggplot(decline,aes(x=log(ctrl.mean),y=log(ctrl.sd),col=year.online))+ geom_point(size=2,na.rm = TRUE)
pp2<-ggplot(decline,aes(x=log(oa.mean),y=log(oa.sd),col=year.online))+ geom_point(size=2,na.rm = TRUE)
grid.arrange(pp1,pp2, nrow =1)

##look at lnRR by year
MLMA_2009_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2009)
summary(MLMA_2009_lnRR)
MLMA_2010_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2010)
summary(MLMA_2010_lnRR)
MLMA_2011_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2011)
summary(MLMA_2011_lnRR)
MLMA_2012_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2012)
summary(MLMA_2012_lnRR)
MLMA_2013_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2013)
summary(MLMA_2013_lnRR)
MLMA_2014_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2014)
summary(MLMA_2014_lnRR)
MLMA_2015_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2015)
summary(MLMA_2015_lnRR)
MLMA_2016_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2016)
summary(MLMA_2016_lnRR)
MLMA_2017_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2017)
summary(MLMA_2017_lnRR)
MLMA_2018_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2018)
summary(MLMA_2018_lnRR)
MLMA_2019_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2019)
summary(MLMA_2019_lnRR)

##set prior
prior <- list(R=list(V = 1, nu =0.002), G=list(G = list(V=1, nu = 0.002)))

##run bayesian MLMA models
model_magnitude_bayes_2009 <- MCMCglmm(yi ~ 1, mev = lnRR2009clean$vi, random = ~obs, data = lnRR2009clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2010 <- MCMCglmm(yi ~ 1, mev = lnRR2010clean$vi, random = ~obs, data = lnRR2010clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2011 <- MCMCglmm(yi ~ 1, mev = lnRR2011clean$vi, random = ~obs, data = lnRR2011clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2012 <- MCMCglmm(yi ~ 1, mev = lnRR2012clean$vi, random = ~obs, data = lnRR2012clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2013 <- MCMCglmm(yi ~ 1, mev = lnRR2013clean$vi, random = ~obs, data = lnRR2013clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2014 <- MCMCglmm(yi ~ 1, mev = lnRR2014clean$vi, random = ~obs, data = lnRR2014clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2015 <- MCMCglmm(yi ~ 1, mev = lnRR2015clean$vi, random = ~obs, data = lnRR2015clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2016 <- MCMCglmm(yi ~ 1, mev = lnRR2016clean$vi, random = ~obs, data = lnRR2016clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2017 <- MCMCglmm(yi ~ 1, mev = lnRR2017clean$vi, random = ~obs, data = lnRR2017clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2018 <- MCMCglmm(yi ~ 1, mev = lnRR2018clean$vi, random = ~obs, data = lnRR2018clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2019 <- MCMCglmm(yi ~ 1, mev = lnRR2019clean$vi, random = ~obs, data = lnRR2019clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)

##get model summaries
summary(model_magnitude_bayes_2009)
summary(model_magnitude_bayes_2010)
summary(model_magnitude_bayes_2011)
summary(model_magnitude_bayes_2012)
summary(model_magnitude_bayes_2013)
summary(model_magnitude_bayes_2014)
summary(model_magnitude_bayes_2015)
summary(model_magnitude_bayes_2016)
summary(model_magnitude_bayes_2017)
summary(model_magnitude_bayes_2018)
summary(model_magnitude_bayes_2019)

##extract posteriors
sol2009 <- model_magnitude_bayes_2009$Sol
VCV2009 <- model_magnitude_bayes_2009$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2009$VCV))]
sol2010 <- model_magnitude_bayes_2010$Sol 
VCV2010 <- model_magnitude_bayes_2010$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2010$VCV))]
sol2011 <- model_magnitude_bayes_2011$Sol 
VCV2011 <- model_magnitude_bayes_2011$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2011$VCV))]
sol2012 <- model_magnitude_bayes_2012$Sol 
VCV2012 <- model_magnitude_bayes_2012$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2012$VCV))]
sol2013 <- model_magnitude_bayes_2013$Sol 
VCV2013 <- model_magnitude_bayes_2013$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2013$VCV))]
sol2014 <- model_magnitude_bayes_2014$Sol 
VCV2014 <- model_magnitude_bayes_2014$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2014$VCV))]
sol2015 <- model_magnitude_bayes_2015$Sol 
VCV2015 <- model_magnitude_bayes_2015$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2015$VCV))]
sol2016 <- model_magnitude_bayes_2016$Sol 
VCV2016 <- model_magnitude_bayes_2016$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2016$VCV))]
sol2017 <- model_magnitude_bayes_2017$Sol 
VCV2017 <- model_magnitude_bayes_2017$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2017$VCV))]
sol2018 <- model_magnitude_bayes_2018$Sol 
VCV2018 <- model_magnitude_bayes_2018$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2018$VCV))]
sol2019 <- model_magnitude_bayes_2019$Sol 
VCV2019 <- model_magnitude_bayes_2019$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2019$VCV))]

##get folded normal function
mu.fnorm <- function(mu, sigma){dnorm(mu, 0, sigma)*2*sigma^2 + mu*(2*pnorm(mu, 0, sigma) -1)}

##get magnitude means + variance
magnitude_mean_2009 <- mu.fnorm(sol2009[,1], sqrt(rowSums(VCV2009)))
magnitude_2009_mean <- data.frame(mean_mag = mean(magnitude_mean_2009), L_CI = HPDinterval(magnitude_mean_2009)[1], U_CI = HPDinterval(magnitude_mean_2009)[2])
magnitude_mean_2010 <- mu.fnorm(sol2010[,1], sqrt(rowSums(VCV2010)))
magnitude_2010_mean <- data.frame(mean_mag = mean(magnitude_mean_2010), L_CI = HPDinterval(magnitude_mean_2010)[1], U_CI = HPDinterval(magnitude_mean_2010)[2])
magnitude_mean_2011 <- mu.fnorm(sol2011[,1], sqrt(rowSums(VCV2011)))
magnitude_2011_mean <- data.frame(mean_mag = mean(magnitude_mean_2011), L_CI = HPDinterval(magnitude_mean_2011)[1], U_CI = HPDinterval(magnitude_mean_2011)[2])
magnitude_mean_2012 <- mu.fnorm(sol2012[,1], sqrt(rowSums(VCV2012)))
magnitude_2012_mean <- data.frame(mean_mag = mean(magnitude_mean_2012), L_CI = HPDinterval(magnitude_mean_2012)[1], U_CI = HPDinterval(magnitude_mean_2012)[2])
magnitude_mean_2013 <- mu.fnorm(sol2013[,1], sqrt(rowSums(VCV2013)))
magnitude_2013_mean <- data.frame(mean_mag = mean(magnitude_mean_2013), L_CI = HPDinterval(magnitude_mean_2013)[1], U_CI = HPDinterval(magnitude_mean_2013)[2])
magnitude_mean_2014 <- mu.fnorm(sol2014[,1], sqrt(rowSums(VCV2014)))
magnitude_2014_mean <- data.frame(mean_mag = mean(magnitude_mean_2014), L_CI = HPDinterval(magnitude_mean_2014)[1], U_CI = HPDinterval(magnitude_mean_2014)[2])
magnitude_mean_2015 <- mu.fnorm(sol2015[,1], sqrt(rowSums(VCV2015)))
magnitude_2015_mean <- data.frame(mean_mag = mean(magnitude_mean_2015), L_CI = HPDinterval(magnitude_mean_2015)[1], U_CI = HPDinterval(magnitude_mean_2015)[2])
magnitude_mean_2016 <- mu.fnorm(sol2016[,1], sqrt(rowSums(VCV2016)))
magnitude_2016_mean <- data.frame(mean_mag = mean(magnitude_mean_2016), L_CI = HPDinterval(magnitude_mean_2016)[1], U_CI = HPDinterval(magnitude_mean_2016)[2])
magnitude_mean_2017 <- mu.fnorm(sol2017[,1], sqrt(rowSums(VCV2017)))
magnitude_2017_mean <- data.frame(mean_mag = mean(magnitude_mean_2017), L_CI = HPDinterval(magnitude_mean_2017)[1], U_CI = HPDinterval(magnitude_mean_2017)[2])
magnitude_mean_2018 <- mu.fnorm(sol2018[,1], sqrt(rowSums(VCV2018)))
magnitude_2018_mean <- data.frame(mean_mag = mean(magnitude_mean_2018), L_CI = HPDinterval(magnitude_mean_2018)[1], U_CI = HPDinterval(magnitude_mean_2018)[2])
magnitude_mean_2019 <- mu.fnorm(sol2019[,1], sqrt(rowSums(VCV2019)))
magnitude_2019_mean <- data.frame(mean_mag = mean(magnitude_mean_2019), L_CI = HPDinterval(magnitude_mean_2019)[1], U_CI = HPDinterval(magnitude_mean_2019)[2])

##view ES magnitudes and uncertainty 
magnitude_2009_mean
magnitude_2010_mean
magnitude_2011_mean
magnitude_2012_mean
magnitude_2013_mean
magnitude_2014_mean
magnitude_2015_mean
magnitude_2016_mean
magnitude_2017_mean
magnitude_2018_mean
magnitude_2019_mean


##### META-ANALYSIS - YEAR ONLINE - OLFACTORY CUES ONLY #####
##attach dataset
decline_cues<-read.csv(file.choose()) ##use dataset "S7 Data"
attach(decline_cues)

##set factors
decline_cues$year.online<-as.factor(decline_cues$year.online)
decline_cues$year.print<-as.factor(decline_cues$year.print)
decline_cues$obs<-as.factor(decline_cues$obs)
decline_cues$study<-as.factor(decline_cues$study)

##view summary
summary(decline_cues)

##subset by year
y2009 <- filter(decline_cues, year.online == "2009")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline_cues))]
y2009$obs <- 1:nrow(y2009)
y2010 <- filter(decline_cues, year.online == "2010")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline_cues))]
y2010$obs <- 1:nrow(y2010)
y2011 <- filter(decline_cues, year.online == "2011")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline_cues))]
y2011$obs <- 1:nrow(y2011)
y2012 <- filter(decline_cues, year.online == "2012")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline_cues))]
y2012$obs <- 1:nrow(y2012)
y2013 <- filter(decline_cues, year.online == "2013")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline_cues))]
y2013$obs <- 1:nrow(y2013)
y2014 <- filter(decline_cues, year.online == "2014")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline_cues))]
y2014$obs <- 1:nrow(y2014)
y2015 <- filter(decline_cues, year.online == "2015")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline_cues))]
y2015$obs <- 1:nrow(y2015)
y2016 <- filter(decline_cues, year.online == "2016")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline_cues))]
y2016$obs <- 1:nrow(y2016)
y2017 <- filter(decline_cues, year.online == "2017")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline_cues))]
y2017$obs <- 1:nrow(y2017)
y2018 <- filter(decline_cues, year.online == "2018")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline_cues))]
y2018$obs <- 1:nrow(y2018)
y2019 <- filter(decline_cues, year.online == "2019")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline_cues))]
y2019$obs <- 1:nrow(y2019)


##comupte effect sizes for each year
lnRR2009 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2009,append=TRUE)
lnRR2010 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2010,append=TRUE)
lnRR2011 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2011,append=TRUE)
lnRR2012 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2012,append=TRUE)
lnRR2013 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2013,append=TRUE)
lnRR2014 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2014,append=TRUE)
lnRR2015 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2015,append=TRUE)
lnRR2016 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2016,append=TRUE)
lnRR2017 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2017,append=TRUE)
lnRR2018 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2018,append=TRUE)
lnRR2019 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2019,append=TRUE)

##remove NAs
lnRR2009clean<-na.omit(lnRR2009)
lnRR2010clean<-na.omit(lnRR2010)
lnRR2011clean<-na.omit(lnRR2011)
lnRR2012clean<-na.omit(lnRR2012)
lnRR2013clean<-na.omit(lnRR2013)
lnRR2014clean<-na.omit(lnRR2014)
lnRR2015clean<-na.omit(lnRR2015)
lnRR2016clean<-na.omit(lnRR2016)
lnRR2017clean<-na.omit(lnRR2017)
lnRR2018clean<-na.omit(lnRR2018)
lnRR2019clean<-na.omit(lnRR2019)

##view mean-variance relationship
pp1<-ggplot(decline,aes(x=log(ctrl.mean),y=log(ctrl.sd),col=year.online))+ geom_point(size=2,na.rm = TRUE)
pp2<-ggplot(decline,aes(x=log(oa.mean),y=log(oa.sd),col=year.online))+ geom_point(size=2,na.rm = TRUE)
grid.arrange(pp1,pp2, nrow =1)

##look at lnRR by year
MLMA_2009_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2009)
summary(MLMA_2009_lnRR)
MLMA_2010_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2010)
summary(MLMA_2010_lnRR)
MLMA_2011_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2011)
summary(MLMA_2011_lnRR)
MLMA_2012_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2012)
summary(MLMA_2012_lnRR)
MLMA_2013_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2013)
summary(MLMA_2013_lnRR)
MLMA_2014_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2014)
summary(MLMA_2014_lnRR)
MLMA_2015_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2015)
summary(MLMA_2015_lnRR)
MLMA_2016_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2016)
summary(MLMA_2016_lnRR)
MLMA_2017_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2017)
summary(MLMA_2017_lnRR)
MLMA_2018_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2018)
summary(MLMA_2018_lnRR)
MLMA_2019_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2019)
summary(MLMA_2019_lnRR)

##set prior
prior <- list(R=list(V = 1, nu =0.002), G=list(G = list(V=1, nu = 0.002)))

##run bayesian MLMA models
model_magnitude_bayes_2009 <- MCMCglmm(yi ~ 1, mev = lnRR2009clean$vi, random = ~obs, data = lnRR2009clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2010 <- MCMCglmm(yi ~ 1, mev = lnRR2010clean$vi, random = ~obs, data = lnRR2010clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2011 <- MCMCglmm(yi ~ 1, mev = lnRR2011clean$vi, random = ~obs, data = lnRR2011clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2012 <- MCMCglmm(yi ~ 1, mev = lnRR2012clean$vi, random = ~obs, data = lnRR2012clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2013 <- MCMCglmm(yi ~ 1, mev = lnRR2013clean$vi, random = ~obs, data = lnRR2013clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2014 <- MCMCglmm(yi ~ 1, mev = lnRR2014clean$vi, random = ~obs, data = lnRR2014clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2015 <- MCMCglmm(yi ~ 1, mev = lnRR2015clean$vi, random = ~obs, data = lnRR2015clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2016 <- MCMCglmm(yi ~ 1, mev = lnRR2016clean$vi, random = ~obs, data = lnRR2016clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2017 <- MCMCglmm(yi ~ 1, mev = lnRR2017clean$vi, random = ~obs, data = lnRR2017clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2018 <- MCMCglmm(yi ~ 1, mev = lnRR2018clean$vi, random = ~obs, data = lnRR2018clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2019 <- MCMCglmm(yi ~ 1, mev = lnRR2019clean$vi, random = ~obs, data = lnRR2019clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)


##get model summaries
summary(model_magnitude_bayes_2009)
summary(model_magnitude_bayes_2010)
summary(model_magnitude_bayes_2011)
summary(model_magnitude_bayes_2012)
summary(model_magnitude_bayes_2013)
summary(model_magnitude_bayes_2014)
summary(model_magnitude_bayes_2015)
summary(model_magnitude_bayes_2016)
summary(model_magnitude_bayes_2017)
summary(model_magnitude_bayes_2018)
summary(model_magnitude_bayes_2019)

##extract posteriors
sol2009 <- model_magnitude_bayes_2009$Sol
VCV2009 <- model_magnitude_bayes_2009$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2009$VCV))]
sol2010 <- model_magnitude_bayes_2010$Sol 
VCV2010 <- model_magnitude_bayes_2010$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2010$VCV))]
sol2011 <- model_magnitude_bayes_2011$Sol 
VCV2011 <- model_magnitude_bayes_2011$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2011$VCV))]
sol2012 <- model_magnitude_bayes_2012$Sol 
VCV2012 <- model_magnitude_bayes_2012$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2012$VCV))]
sol2013 <- model_magnitude_bayes_2013$Sol 
VCV2013 <- model_magnitude_bayes_2013$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2013$VCV))]
sol2014 <- model_magnitude_bayes_2014$Sol 
VCV2014 <- model_magnitude_bayes_2014$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2014$VCV))]
sol2015 <- model_magnitude_bayes_2015$Sol 
VCV2015 <- model_magnitude_bayes_2015$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2015$VCV))]
sol2016 <- model_magnitude_bayes_2016$Sol 
VCV2016 <- model_magnitude_bayes_2016$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2016$VCV))]
sol2017 <- model_magnitude_bayes_2017$Sol 
VCV2017 <- model_magnitude_bayes_2017$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2017$VCV))]
sol2018 <- model_magnitude_bayes_2018$Sol 
VCV2018 <- model_magnitude_bayes_2018$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2018$VCV))]
sol2019 <- model_magnitude_bayes_2019$Sol 
VCV2019 <- model_magnitude_bayes_2019$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2019$VCV))]

##get folded normal function
mu.fnorm <- function(mu, sigma){dnorm(mu, 0, sigma)*2*sigma^2 + mu*(2*pnorm(mu, 0, sigma) -1)}

##get magnitude means + variance
magnitude_mean_2009 <- mu.fnorm(sol2009[,1], sqrt(rowSums(VCV2009)))
magnitude_2009_mean <- data.frame(mean_mag = mean(magnitude_mean_2009), L_CI = HPDinterval(magnitude_mean_2009)[1], U_CI = HPDinterval(magnitude_mean_2009)[2])
magnitude_mean_2010 <- mu.fnorm(sol2010[,1], sqrt(rowSums(VCV2010)))
magnitude_2010_mean <- data.frame(mean_mag = mean(magnitude_mean_2010), L_CI = HPDinterval(magnitude_mean_2010)[1], U_CI = HPDinterval(magnitude_mean_2010)[2])
magnitude_mean_2011 <- mu.fnorm(sol2011[,1], sqrt(rowSums(VCV2011)))
magnitude_2011_mean <- data.frame(mean_mag = mean(magnitude_mean_2011), L_CI = HPDinterval(magnitude_mean_2011)[1], U_CI = HPDinterval(magnitude_mean_2011)[2])
magnitude_mean_2012 <- mu.fnorm(sol2012[,1], sqrt(rowSums(VCV2012)))
magnitude_2012_mean <- data.frame(mean_mag = mean(magnitude_mean_2012), L_CI = HPDinterval(magnitude_mean_2012)[1], U_CI = HPDinterval(magnitude_mean_2012)[2])
magnitude_mean_2013 <- mu.fnorm(sol2013[,1], sqrt(rowSums(VCV2013)))
magnitude_2013_mean <- data.frame(mean_mag = mean(magnitude_mean_2013), L_CI = HPDinterval(magnitude_mean_2013)[1], U_CI = HPDinterval(magnitude_mean_2013)[2])
magnitude_mean_2014 <- mu.fnorm(sol2014[,1], sqrt(rowSums(VCV2014)))
magnitude_2014_mean <- data.frame(mean_mag = mean(magnitude_mean_2014), L_CI = HPDinterval(magnitude_mean_2014)[1], U_CI = HPDinterval(magnitude_mean_2014)[2])
magnitude_mean_2015 <- mu.fnorm(sol2015[,1], sqrt(rowSums(VCV2015)))
magnitude_2015_mean <- data.frame(mean_mag = mean(magnitude_mean_2015), L_CI = HPDinterval(magnitude_mean_2015)[1], U_CI = HPDinterval(magnitude_mean_2015)[2])
magnitude_mean_2016 <- mu.fnorm(sol2016[,1], sqrt(rowSums(VCV2016)))
magnitude_2016_mean <- data.frame(mean_mag = mean(magnitude_mean_2016), L_CI = HPDinterval(magnitude_mean_2016)[1], U_CI = HPDinterval(magnitude_mean_2016)[2])
magnitude_mean_2017 <- mu.fnorm(sol2017[,1], sqrt(rowSums(VCV2017)))
magnitude_2017_mean <- data.frame(mean_mag = mean(magnitude_mean_2017), L_CI = HPDinterval(magnitude_mean_2017)[1], U_CI = HPDinterval(magnitude_mean_2017)[2])
magnitude_mean_2018 <- mu.fnorm(sol2018[,1], sqrt(rowSums(VCV2018)))
magnitude_2018_mean <- data.frame(mean_mag = mean(magnitude_mean_2018), L_CI = HPDinterval(magnitude_mean_2018)[1], U_CI = HPDinterval(magnitude_mean_2018)[2])
magnitude_mean_2019 <- mu.fnorm(sol2019[,1], sqrt(rowSums(VCV2019)))
magnitude_2019_mean <- data.frame(mean_mag = mean(magnitude_mean_2019), L_CI = HPDinterval(magnitude_mean_2019)[1], U_CI = HPDinterval(magnitude_mean_2019)[2])

##view ES magnitudes and uncertainty 
magnitude_2009_mean
magnitude_2010_mean
magnitude_2011_mean
magnitude_2012_mean
magnitude_2013_mean
magnitude_2014_mean
magnitude_2015_mean
magnitude_2016_mean
magnitude_2017_mean
magnitude_2018_mean
magnitude_2019_mean


##### META-ANALYSIS - YEAR ONLINE - LARVAE ONLY #####
##attach dataset
decline_larvae<-read.csv(file.choose()) ##use dataset "S8 Data"
attach(decline_larvae)

##set factors
decline_larvae$year.online<-as.factor(decline_larvae$year.online)
decline_larvae$year.print<-as.factor(decline_larvae$year.print)
decline_larvae$obs<-as.factor(decline_larvae$obs)
decline_larvae$study<-as.factor(decline_larvae$study)

##view summary
summary(decline_larvae)

##subset by year
y2009 <- filter(decline_larvae, year.online == "2009")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline_larvae))]
y2009$obs <- 1:nrow(y2009)
y2010 <- filter(decline_larvae, year.online == "2010")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline_larvae))]
y2010$obs <- 1:nrow(y2010)
y2011 <- filter(decline_larvae, year.online == "2011")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline_larvae))]
y2011$obs <- 1:nrow(y2011)
y2012 <- filter(decline_larvae, year.online == "2012")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline_larvae))]
y2012$obs <- 1:nrow(y2012)
y2013 <- filter(decline_larvae, year.online == "2013")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline_larvae))]
y2013$obs <- 1:nrow(y2013)
y2014 <- filter(decline_larvae, year.online == "2014")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline_larvae))]
y2014$obs <- 1:nrow(y2014)
y2015 <- filter(decline_larvae, year.online == "2015")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline_larvae))]
y2015$obs <- 1:nrow(y2015)
y2016 <- filter(decline_larvae, year.online == "2016")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline_larvae))]
y2016$obs <- 1:nrow(y2016)
y2017 <- filter(decline_larvae, year.online == "2017")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline_larvae))]
y2017$obs <- 1:nrow(y2017)
y2018 <- filter(decline_larvae, year.online == "2018")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline_larvae))]
y2018$obs <- 1:nrow(y2018)
y2019 <- filter(decline_larvae, year.online == "2019")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline_larvae))]
y2019$obs <- 1:nrow(y2019)


##comupte effect sizes for each year
lnRR2009 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2009,append=TRUE)
lnRR2010 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2010,append=TRUE)
lnRR2011 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2011,append=TRUE)
lnRR2012 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2012,append=TRUE)
lnRR2013 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2013,append=TRUE)
lnRR2014 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2014,append=TRUE)
lnRR2015 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2015,append=TRUE)
lnRR2016 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2016,append=TRUE)
lnRR2017 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2017,append=TRUE)
lnRR2018 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2018,append=TRUE)
lnRR2019 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2019,append=TRUE)

##remove NAs
lnRR2009clean<-na.omit(lnRR2009)
lnRR2010clean<-na.omit(lnRR2010)
lnRR2011clean<-na.omit(lnRR2011)
lnRR2012clean<-na.omit(lnRR2012)
lnRR2013clean<-na.omit(lnRR2013)
lnRR2014clean<-na.omit(lnRR2014)
lnRR2015clean<-na.omit(lnRR2015)
lnRR2016clean<-na.omit(lnRR2016)
lnRR2017clean<-na.omit(lnRR2017)
lnRR2018clean<-na.omit(lnRR2018)
lnRR2019clean<-na.omit(lnRR2019)

##look at lnRR by year
MLMA_2009_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2009)
summary(MLMA_2009_lnRR)
MLMA_2010_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2010)
summary(MLMA_2010_lnRR)
MLMA_2011_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2011)
summary(MLMA_2011_lnRR)
MLMA_2012_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2012)
summary(MLMA_2012_lnRR)
MLMA_2013_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2013)
summary(MLMA_2013_lnRR)
MLMA_2014_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2014)
summary(MLMA_2014_lnRR)
MLMA_2015_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2015)
summary(MLMA_2015_lnRR)
MLMA_2016_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2016)
summary(MLMA_2016_lnRR)
MLMA_2017_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2017)
summary(MLMA_2017_lnRR)
MLMA_2018_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2018)
summary(MLMA_2018_lnRR)
MLMA_2019_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2019)
summary(MLMA_2019_lnRR)

##set prior
prior <- list(R=list(V = 1, nu =0.002), G=list(G = list(V=1, nu = 0.002)))

##run bayesian MLMA models
model_magnitude_bayes_2009 <- MCMCglmm(yi ~ 1, mev = lnRR2009clean$vi, random = ~obs, data = lnRR2009clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2010 <- MCMCglmm(yi ~ 1, mev = lnRR2010clean$vi, random = ~obs, data = lnRR2010clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2011 <- MCMCglmm(yi ~ 1, mev = lnRR2011clean$vi, random = ~obs, data = lnRR2011clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2012 <- MCMCglmm(yi ~ 1, mev = lnRR2012clean$vi, random = ~obs, data = lnRR2012clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2013 <- MCMCglmm(yi ~ 1, mev = lnRR2013clean$vi, random = ~obs, data = lnRR2013clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2014 <- MCMCglmm(yi ~ 1, mev = lnRR2014clean$vi, random = ~obs, data = lnRR2014clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2015 <- MCMCglmm(yi ~ 1, mev = lnRR2015clean$vi, random = ~obs, data = lnRR2015clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2016 <- MCMCglmm(yi ~ 1, mev = lnRR2016clean$vi, random = ~obs, data = lnRR2016clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2017 <- MCMCglmm(yi ~ 1, mev = lnRR2017clean$vi, random = ~obs, data = lnRR2017clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2018 <- MCMCglmm(yi ~ 1, mev = lnRR2018clean$vi, random = ~obs, data = lnRR2018clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2019 <- MCMCglmm(yi ~ 1, mev = lnRR2019clean$vi, random = ~obs, data = lnRR2019clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)


##get model summaries
summary(model_magnitude_bayes_2009)
summary(model_magnitude_bayes_2010)
summary(model_magnitude_bayes_2011)
summary(model_magnitude_bayes_2012)
summary(model_magnitude_bayes_2013)
summary(model_magnitude_bayes_2014)
summary(model_magnitude_bayes_2015)
summary(model_magnitude_bayes_2016)
summary(model_magnitude_bayes_2017)
summary(model_magnitude_bayes_2018)
summary(model_magnitude_bayes_2019)

##extract posteriors
sol2009 <- model_magnitude_bayes_2009$Sol
VCV2009 <- model_magnitude_bayes_2009$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2009$VCV))]
sol2010 <- model_magnitude_bayes_2010$Sol 
VCV2010 <- model_magnitude_bayes_2010$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2010$VCV))]
sol2011 <- model_magnitude_bayes_2011$Sol 
VCV2011 <- model_magnitude_bayes_2011$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2011$VCV))]
sol2012 <- model_magnitude_bayes_2012$Sol 
VCV2012 <- model_magnitude_bayes_2012$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2012$VCV))]
sol2013 <- model_magnitude_bayes_2013$Sol 
VCV2013 <- model_magnitude_bayes_2013$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2013$VCV))]
sol2014 <- model_magnitude_bayes_2014$Sol 
VCV2014 <- model_magnitude_bayes_2014$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2014$VCV))]
sol2015 <- model_magnitude_bayes_2015$Sol 
VCV2015 <- model_magnitude_bayes_2015$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2015$VCV))]
sol2016 <- model_magnitude_bayes_2016$Sol 
VCV2016 <- model_magnitude_bayes_2016$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2016$VCV))]
sol2017 <- model_magnitude_bayes_2017$Sol 
VCV2017 <- model_magnitude_bayes_2017$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2017$VCV))]
sol2018 <- model_magnitude_bayes_2018$Sol 
VCV2018 <- model_magnitude_bayes_2018$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2018$VCV))]
sol2019 <- model_magnitude_bayes_2019$Sol 
VCV2019 <- model_magnitude_bayes_2019$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2019$VCV))]

##get folded normal function
mu.fnorm <- function(mu, sigma){dnorm(mu, 0, sigma)*2*sigma^2 + mu*(2*pnorm(mu, 0, sigma) -1)}

##get magnitude means + variance
magnitude_mean_2009 <- mu.fnorm(sol2009[,1], sqrt(rowSums(VCV2009)))
magnitude_2009_mean <- data.frame(mean_mag = mean(magnitude_mean_2009), L_CI = HPDinterval(magnitude_mean_2009)[1], U_CI = HPDinterval(magnitude_mean_2009)[2])
magnitude_mean_2010 <- mu.fnorm(sol2010[,1], sqrt(rowSums(VCV2010)))
magnitude_2010_mean <- data.frame(mean_mag = mean(magnitude_mean_2010), L_CI = HPDinterval(magnitude_mean_2010)[1], U_CI = HPDinterval(magnitude_mean_2010)[2])
magnitude_mean_2011 <- mu.fnorm(sol2011[,1], sqrt(rowSums(VCV2011)))
magnitude_2011_mean <- data.frame(mean_mag = mean(magnitude_mean_2011), L_CI = HPDinterval(magnitude_mean_2011)[1], U_CI = HPDinterval(magnitude_mean_2011)[2])
magnitude_mean_2012 <- mu.fnorm(sol2012[,1], sqrt(rowSums(VCV2012)))
magnitude_2012_mean <- data.frame(mean_mag = mean(magnitude_mean_2012), L_CI = HPDinterval(magnitude_mean_2012)[1], U_CI = HPDinterval(magnitude_mean_2012)[2])
magnitude_mean_2013 <- mu.fnorm(sol2013[,1], sqrt(rowSums(VCV2013)))
magnitude_2013_mean <- data.frame(mean_mag = mean(magnitude_mean_2013), L_CI = HPDinterval(magnitude_mean_2013)[1], U_CI = HPDinterval(magnitude_mean_2013)[2])
magnitude_mean_2014 <- mu.fnorm(sol2014[,1], sqrt(rowSums(VCV2014)))
magnitude_2014_mean <- data.frame(mean_mag = mean(magnitude_mean_2014), L_CI = HPDinterval(magnitude_mean_2014)[1], U_CI = HPDinterval(magnitude_mean_2014)[2])
magnitude_mean_2015 <- mu.fnorm(sol2015[,1], sqrt(rowSums(VCV2015)))
magnitude_2015_mean <- data.frame(mean_mag = mean(magnitude_mean_2015), L_CI = HPDinterval(magnitude_mean_2015)[1], U_CI = HPDinterval(magnitude_mean_2015)[2])
magnitude_mean_2016 <- mu.fnorm(sol2016[,1], sqrt(rowSums(VCV2016)))
magnitude_2016_mean <- data.frame(mean_mag = mean(magnitude_mean_2016), L_CI = HPDinterval(magnitude_mean_2016)[1], U_CI = HPDinterval(magnitude_mean_2016)[2])
magnitude_mean_2017 <- mu.fnorm(sol2017[,1], sqrt(rowSums(VCV2017)))
magnitude_2017_mean <- data.frame(mean_mag = mean(magnitude_mean_2017), L_CI = HPDinterval(magnitude_mean_2017)[1], U_CI = HPDinterval(magnitude_mean_2017)[2])
magnitude_mean_2018 <- mu.fnorm(sol2018[,1], sqrt(rowSums(VCV2018)))
magnitude_2018_mean <- data.frame(mean_mag = mean(magnitude_mean_2018), L_CI = HPDinterval(magnitude_mean_2018)[1], U_CI = HPDinterval(magnitude_mean_2018)[2])
magnitude_mean_2019 <- mu.fnorm(sol2019[,1], sqrt(rowSums(VCV2019)))
magnitude_2019_mean <- data.frame(mean_mag = mean(magnitude_mean_2019), L_CI = HPDinterval(magnitude_mean_2019)[1], U_CI = HPDinterval(magnitude_mean_2019)[2])

##view ES magnitudes and uncertainty 
magnitude_2009_mean
magnitude_2010_mean
magnitude_2011_mean
magnitude_2012_mean
magnitude_2013_mean
magnitude_2014_mean
magnitude_2015_mean
magnitude_2016_mean
magnitude_2017_mean
magnitude_2018_mean
magnitude_2019_mean


##### META-ANALYSIS - INVESTIGATOR EFFECTS #####
##attach dataset
decline_authors<-read.csv(file.choose()) ##use dataset "S9 Data"
attach(decline_authors)

##set factors
decline_authors$year.online<-as.factor(decline_authors$year.online)
decline_authors$year.print<-as.factor(decline_authors$year.print)
decline_authors$obs<-as.factor(decline_authors$obs)
decline_authors$study<-as.factor(decline_authors$study)

##view summary
summary(decline_authors)

##subset by year
y2012 <- filter(decline_authors, year.online == "2012")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline_authors))]
y2012$obs <- 1:nrow(y2012)
y2013 <- filter(decline_authors, year.online == "2013")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline_authors))]
y2013$obs <- 1:nrow(y2013)
y2014 <- filter(decline_authors, year.online == "2014")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline_authors))]
y2014$obs <- 1:nrow(y2014)
y2015 <- filter(decline_authors, year.online == "2015")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline_authors))]
y2015$obs <- 1:nrow(y2015)
y2016 <- filter(decline_authors, year.online == "2016")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline_authors))]
y2016$obs <- 1:nrow(y2016)
y2017 <- filter(decline_authors, year.online == "2017")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline_authors))]
y2017$obs <- 1:nrow(y2017)
y2018 <- filter(decline_authors, year.online == "2018")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline_authors))]
y2018$obs <- 1:nrow(y2018)
y2019 <- filter(decline_authors, year.online == "2019")[,-match(c("avg.n","cue.type","if.at.pub","X2017.if","if.group"), colnames (decline_authors))]
y2019$obs <- 1:nrow(y2019)


##comupte effect sizes for each year
lnRR2012 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2012,append=TRUE)
lnRR2013 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2013,append=TRUE)
lnRR2014 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2014,append=TRUE)
lnRR2015 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2015,append=TRUE)
lnRR2016 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2016,append=TRUE)
lnRR2017 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2017,append=TRUE)
lnRR2018 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2018,append=TRUE)
lnRR2019 <- escalc(measure= "ROM", m1i=oa.mean,sd1i=oa.sd,n1i=oa.n,m2i=ctrl.mean,sd2i=ctrl.sd,n2i=ctrl.n,data=y2019,append=TRUE)

##remove NAs
lnRR2012clean<-na.omit(lnRR2012)
lnRR2013clean<-na.omit(lnRR2013)
lnRR2014clean<-na.omit(lnRR2014)
lnRR2015clean<-na.omit(lnRR2015)
lnRR2016clean<-na.omit(lnRR2016)
lnRR2017clean<-na.omit(lnRR2017)
lnRR2018clean<-na.omit(lnRR2018)
lnRR2019clean<-na.omit(lnRR2019)

##view mean-variance relationship
pp1<-ggplot(decline_authors,aes(x=log(ctrl.mean),y=log(ctrl.sd),col=year.online))+ geom_point(size=2,na.rm = TRUE)
pp2<-ggplot(decline_authors,aes(x=log(oa.mean),y=log(oa.sd),col=year.online))+ geom_point(size=2,na.rm = TRUE)
grid.arrange(pp1,pp2, nrow =1)

##look at lnRR by year
MLMA_2012_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2012)
summary(MLMA_2012_lnRR)
MLMA_2013_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2013)
summary(MLMA_2013_lnRR)
MLMA_2014_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2014)
summary(MLMA_2014_lnRR)
MLMA_2015_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2015)
summary(MLMA_2015_lnRR)
MLMA_2016_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2016)
summary(MLMA_2016_lnRR)
MLMA_2017_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2017)
summary(MLMA_2017_lnRR)
MLMA_2018_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2018)
summary(MLMA_2018_lnRR)
MLMA_2019_lnRR <- rma.mv(yi ~ 1, V = vi, random=~1|obs, data=lnRR2019)
summary(MLMA_2019_lnRR)

##set prior
prior <- list(R=list(V = 1, nu =0.002), G = list(G = list(V=1, nu = 0.002)))

##run bayesian MLMA models
model_magnitude_bayes_2012 <- MCMCglmm(yi ~ 1, mev = lnRR2012clean$vi, random = ~obs, data = lnRR2012clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2013 <- MCMCglmm(yi ~ 1, mev = lnRR2013clean$vi, random = ~obs, data = lnRR2013clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2014 <- MCMCglmm(yi ~ 1, mev = lnRR2014clean$vi, random = ~obs, data = lnRR2014clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2015 <- MCMCglmm(yi ~ 1, mev = lnRR2015clean$vi, random = ~obs, data = lnRR2015clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2016 <- MCMCglmm(yi ~ 1, mev = lnRR2016clean$vi, random = ~obs, data = lnRR2016clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2017 <- MCMCglmm(yi ~ 1, mev = lnRR2017clean$vi, random = ~obs, data = lnRR2017clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2018 <- MCMCglmm(yi ~ 1, mev = lnRR2018clean$vi, random = ~obs, data = lnRR2018clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)
model_magnitude_bayes_2019 <- MCMCglmm(yi ~ 1, mev = lnRR2019clean$vi, random = ~obs, data = lnRR2019clean, prior = prior, burnin = 10000, nitt = 1000000, thin = 100, verbose = FALSE)

##get model summaries
summary(model_magnitude_bayes_2012)
summary(model_magnitude_bayes_2013)
summary(model_magnitude_bayes_2014)
summary(model_magnitude_bayes_2015)
summary(model_magnitude_bayes_2016)
summary(model_magnitude_bayes_2017)
summary(model_magnitude_bayes_2018)
summary(model_magnitude_bayes_2019)

##extract posteriors
sol2012 <- model_magnitude_bayes_2012$Sol 
VCV2012 <- model_magnitude_bayes_2012$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2012$VCV))]
sol2013 <- model_magnitude_bayes_2013$Sol 
VCV2013 <- model_magnitude_bayes_2013$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2013$VCV))]
sol2014 <- model_magnitude_bayes_2014$Sol 
VCV2014 <- model_magnitude_bayes_2014$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2014$VCV))]
sol2015 <- model_magnitude_bayes_2015$Sol 
VCV2015 <- model_magnitude_bayes_2015$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2015$VCV))]
sol2016 <- model_magnitude_bayes_2016$Sol 
VCV2016 <- model_magnitude_bayes_2016$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2016$VCV))]
sol2017 <- model_magnitude_bayes_2017$Sol 
VCV2017 <- model_magnitude_bayes_2017$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2017$VCV))]
sol2018 <- model_magnitude_bayes_2018$Sol 
VCV2018 <- model_magnitude_bayes_2018$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2018$VCV))]
sol2019 <- model_magnitude_bayes_2019$Sol 
VCV2019 <- model_magnitude_bayes_2019$VCV[,-match("sqrt(mev):sqrt(mev).meta", colnames(model_magnitude_bayes_2019$VCV))]

##get folded normal function
mu.fnorm <- function(mu, sigma){dnorm(mu, 0, sigma)*2*sigma^2 + mu*(2*pnorm(mu, 0, sigma) -1)}

##get magnitude means + variance
magnitude_mean_2012 <- mu.fnorm(sol2012[,1], sqrt(rowSums(VCV2012)))
magnitude_2012_mean <- data.frame(mean_mag = mean(magnitude_mean_2012), L_CI = HPDinterval(magnitude_mean_2012)[1], U_CI = HPDinterval(magnitude_mean_2012)[2])
magnitude_mean_2013 <- mu.fnorm(sol2013[,1], sqrt(rowSums(VCV2013)))
magnitude_2013_mean <- data.frame(mean_mag = mean(magnitude_mean_2013), L_CI = HPDinterval(magnitude_mean_2013)[1], U_CI = HPDinterval(magnitude_mean_2013)[2])
magnitude_mean_2014 <- mu.fnorm(sol2014[,1], sqrt(rowSums(VCV2014)))
magnitude_2014_mean <- data.frame(mean_mag = mean(magnitude_mean_2014), L_CI = HPDinterval(magnitude_mean_2014)[1], U_CI = HPDinterval(magnitude_mean_2014)[2])
magnitude_mean_2015 <- mu.fnorm(sol2015[,1], sqrt(rowSums(VCV2015)))
magnitude_2015_mean <- data.frame(mean_mag = mean(magnitude_mean_2015), L_CI = HPDinterval(magnitude_mean_2015)[1], U_CI = HPDinterval(magnitude_mean_2015)[2])
magnitude_mean_2016 <- mu.fnorm(sol2016[,1], sqrt(rowSums(VCV2016)))
magnitude_2016_mean <- data.frame(mean_mag = mean(magnitude_mean_2016), L_CI = HPDinterval(magnitude_mean_2016)[1], U_CI = HPDinterval(magnitude_mean_2016)[2])
magnitude_mean_2017 <- mu.fnorm(sol2017[,1], sqrt(rowSums(VCV2017)))
magnitude_2017_mean <- data.frame(mean_mag = mean(magnitude_mean_2017), L_CI = HPDinterval(magnitude_mean_2017)[1], U_CI = HPDinterval(magnitude_mean_2017)[2])
magnitude_mean_2018 <- mu.fnorm(sol2018[,1], sqrt(rowSums(VCV2018)))
magnitude_2018_mean <- data.frame(mean_mag = mean(magnitude_mean_2018), L_CI = HPDinterval(magnitude_mean_2018)[1], U_CI = HPDinterval(magnitude_mean_2018)[2])
magnitude_mean_2019 <- mu.fnorm(sol2019[,1], sqrt(rowSums(VCV2019)))
magnitude_2019_mean <- data.frame(mean_mag = mean(magnitude_mean_2019), L_CI = HPDinterval(magnitude_mean_2019)[1], U_CI = HPDinterval(magnitude_mean_2019)[2])

##view ES magnitudes and uncertainty 
magnitude_2012_mean
magnitude_2013_mean
magnitude_2014_mean
magnitude_2015_mean
magnitude_2016_mean
magnitude_2017_mean
magnitude_2018_mean
magnitude_2019_mean




##### CREATE SCATTERPLOT FIGURE TO VISUALIZE MEAN EFFECT SIZE MAGNITUDE FOR EACH OBSERVATION OVER TIME (FIG 1A) #####
##attach dataset
decline_allobs<-read.csv(file.choose()) ##use dataset "S10 Data"
attach(decline_allobs)
summary(decline_allobs)

#Create plot
Decline_studies_loess<-ggplot(decline_allobs,aes(x=year.online, y=lnrr.mag,color=study)) + geom_smooth(method="loess", se=TRUE, fullrange=TRUE, level=0.95,color="black")+geom_point(size=ctrl.n*0.03,alpha=0.6) + geom_smooth(method="loess", se=TRUE, fullrange=TRUE, level=0.95,color="black") + scale_size(range = c(1, 2), name="Sample size")+ scale_color_viridis(discrete=TRUE)+ xlab("Year")+ylab("Effect size magnitude (lnRR)")+ scale_x_continuous(breaks = round(seq(min(study$Year), max(study$Year), by = 1),1))+scale_y_continuous(breaks = round(seq(min(study$lnrr), 15, by = 1),1)) + theme_minimal(12)+theme(legend.position = "none")+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
Decline_studies_loess

tiff("Decline all obs loess", width = 4.2, height = 4.6, units = "in", res = 800,compression="lzw")
Decline_studies_loess
dev.off()

##### CREATE SCATTERPLOT FIGURE TO VISUALIZE SAMPLE SIZE BIAS (FIG 3) #####
##attach dataset
decline_n<-read.csv(file.choose()) ##use dataset "S11 Data"
attach(decline_n)
summary(decline_n)

lnrr_samplesize<-ggplot(decline_n,aes(x=avg.n, y=mean.lnrr.mag, color=study))+geom_point(size=2,alpha=0.6)+ geom_smooth(method="lm", se=TRUE, fullrange=TRUE, level=0.95,color="black")+ scale_color_viridis(discrete=TRUE)+ xlab("Mean sample size")+ylab("Effect size magitude (lnRR)")+ scale_x_continuous(breaks = round(seq(0, 600, by = 50),1))+scale_y_continuous(breaks = round(seq(min(mean.lnrr.mag), max(mean.lnrr.mag), by = 0.5),1))+geom_hline(yintercept = 1,linetype="dashed",color="red")+geom_vline(xintercept = 30,linetype="dashed",color="red")+theme_bw(12)+theme(legend.position = "NONE")+ theme(panel.grid.minor = element_blank())
lnrr_samplesize


tiff("Sample size decline", width = 6, height = 4, units = "in", res = 800,compression="lzw")
lnrr_samplesize
dev.off()


##### CREATE SCATTERPLOT FIGURE TO VISUALIZE PUBLICATION BIAS (FIG 4C-E) #####
##attach dataset
decline_impact<-read.csv(file.choose()) ##use dataset "S12 Data"
attach(decline_impact)
summary(decline_impact)


#Create plot 
lnrr_impact_lm<-ggplot(decline_impact,aes(x=if.at.pub, y=mean.lnrr.mag, color=study))+geom_point(size=avg.n*0.07,alpha=0.6) + geom_smooth(method="glm", se=TRUE, fullrange=TRUE, level=0.95,color="black") + scale_size(range = c(1, 3), name="Sample size")+ scale_color_viridis(discrete=TRUE)+ xlab("Impact factor at time of publication online")+ylab("Effect size magnitude (lnRR)")+ scale_x_continuous(breaks = round(seq(0, 20, by = 2),1))+scale_y_continuous(breaks = round(seq(min(mean.lnrr.mag), max(mean.lnrr.mag), by = 0.5),1))+theme_bw(12)+theme(legend.position = "NONE")+ theme(panel.grid.minor = element_blank())

citations_impact_lm<-ggplot(decline_impact,aes(x=if.at.pub, y=citations.per.year, color=study))+geom_point(size=avg.n*0.07,alpha=0.6) + geom_smooth(method="glm", se=TRUE, fullrange=TRUE, level=0.95,color="black") + scale_size(range = c(1, 3), name="Sample size")+ scale_color_viridis(discrete=TRUE)+ xlab("Impact factor at time of publication online")+ylab("Citations per year")+ scale_x_continuous(breaks = round(seq(0, 20, by = 2),1))+scale_y_continuous(breaks = round(seq(min(citations.per.year), max(citations.per.year), by = 5),1))+theme_bw(12)+theme(legend.position = "NONE")+ theme(panel.grid.minor = element_blank())

citations_lnrr_lm<-ggplot(decline_impact,aes(x=mean.lnrr.mag, y=citations.per.year, color=study))+geom_point(size=avg.n*0.07,alpha=0.6) + geom_smooth(method="glm", se=TRUE, fullrange=TRUE, level=0.95,color="black") + scale_size(range = c(1, 3), name="Sample size")+ scale_color_viridis(discrete=TRUE)+ xlab("Effect size magnitude (lnRR)")+ylab("Citations per year")+ scale_x_continuous(breaks = round(seq(min(mean.lnrr.mag), max(mean.lnrr.mag), by = 0.5),1))+scale_y_continuous(breaks = round(seq(min(citations.per.year), max(citations.per.year), by = 5),1))+theme_bw(12)+theme(legend.position = "NONE")+ theme(panel.grid.minor = element_blank())

impactfactor<-lnrr_impact_lm + citations_impact_lm + citations_lnrr_lm

impactfactor

tiff("Decline impact", width = 12, height = 5, units = "in", res = 800,compression="lzw")
impactfactor
dev.off()



##### CREATE SCATTERPLOT FIGURE TO VISUALIZE INVESTIGATOR EFFECTS (FIG 5) #####
##attach dataset
nodixsonmunday<-read.csv(file.choose()) ##use dataset "S13 Data"
attach(nodixsonmunday)
summary(nodixsonmunday)

#Create plot
Decline_nodixsonmunday_loess<-ggplot(nodixsonmunday,aes(x=year.online, y=lnrr.mag,color=study)) + geom_smooth(method="loess", se=TRUE, fullrange=TRUE, level=0.95,color="black")+geom_point(size=ctrl.n*0.03,alpha=0.6) + geom_smooth(method="loess", se=TRUE, fullrange=TRUE, level=0.95,color="black") + scale_size(range = c(1, 2), name="Sample size")+ scale_color_viridis(discrete=TRUE)+ xlab("Year")+ylab("Effect size magnitude (lnRR)")+scale_x_continuous(breaks = round(seq(min(nodixsonmunday$year.online), max(nodixsonmunday$year.online), by=1),1)) + scale_y_continuous(limits = c(-0.3, 14), breaks = seq(0, 14, by = 1))+geom_hline(yintercept = 0,color="black")+ theme_bw(12)+theme(legend.position = "none")+ theme(panel.grid.minor = element_blank())
Decline_nodixsonmunday_loess

tiff("Decline no dixson munday loess", width = 4.2, height = 4.6, units = "in", res = 800,compression="lzw")
Decline_nodixsonmunday_loess
dev.off()
