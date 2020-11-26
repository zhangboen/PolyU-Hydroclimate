rm(list=ls())
library(pracma)
library(foreach)
library(doParallel)
library(VineCopula)
library(cubature)
library(Emcdf)

### ------------------univariate cdf and pdf---------------------------
# Compute empirical univariate probability distribution
Calc_Emp_Prob <- function(D) {
  # Length of data
  n = length(D)
  # Pre-assign probability array
  P = zeros(n,1)
  
  # Loop through the data
  for (i in 1:n) {
    P[i,1] = sum( D <= D[i] )
  }
  
  # Gringorten plotting position
  Y = (P - 0.44) / (n + 0.12)
  
}

# fit univariate cdf, select from 13 cdf candidate (no error)
fitalldist <- function(data) {
  
  require(fitdistrplus)
  require(extRemes)
  require(ismev)
  require(SCI) 
  require(goft)
  require(gPdtest)
  require(actuar)
  require(evd)
  
  # change the same data value slightly
  dup_data <- duplicated(data)
  for (i in 1:length(dup_data)) {
    if (dup_data[i] == 0) {
      data[data=dup_data[i]] = data[data=dup_data[i]] + runif(data[data=dup_data[i]],0,0.001)
    } else {
      data[data=dup_data[i]] = data[data=dup_data[i]] + runif(data[data=dup_data[i]],-0.001,0.001)
    }
  }
  
  # calculate empirical CDF
  emp <- Calc_Emp_Prob(data)
  
  # total 13 distributions
  distnames <- c("gamma","exponential","weibull","normal","logistic",
                 "log-normal","log-logistic","cauchy","gumbel",
                 "generalized extreme value","generalized pareto","pearson iii",
                 "inverse gaussian")
  
  results <- data.frame('aic' = numeric(), 'aicc' = numeric(), 
                        'bic' = numeric(), 'p' = numeric())
  
  ep_mat <- matrix( NA, length(data), 13 )
  pdf_mat <- matrix(NA, length(data), 13 )
  para_list <- list()
  
  loop_id <- 1
  # distribution not for negative values
  for (i in c("gamma","exp","weibull","norm","logis","lnorm","llogis","cauchy")) {
    if (any(data<0)) {
      if (i %in% c("gamma","exp","weibull","lnorm","llogis")) {
        results[loop_id,] <- c(NA,NA,NA,NA)
        loop_id <- loop_id + 1
        next
      } 
    }
    possibleError <- tryCatch({fit <- fitdist(data,i)},error = function(e) e)
    if (!inherits(possibleError, "error")) {
      fit <- fitdist(data,i)
      para <- c(as.list(fit$estimate), as.list(fit$fix.arg))
      para_list[[loop_id]] = para
      name <- fit$distname
      pdistname <- paste("p", name, sep = "")
      ddistname <- paste("d", name, sep = "")
      ep <- do.call(pdistname, c(list(data), as.list(para)))
      pdf <- do.call(ddistname, c(list(data), as.list(para)))
      ep_mat[,loop_id] <- ep
      pdf_mat[,loop_id] <- pdf
      
      aic <- fit$aic
      aicc <- aic + 2*length(para)*(length(para)+1)/(fit$n-length(para)-1)
      bic <- fit$bic
      p <- ks.test(emp,ep)$p
      
      results[loop_id,] <- c(aic,aicc,bic,p)
      loop_id <- loop_id + 1
    } else {
      results[loop_id,] <- c(NA,NA,NA,NA)
      loop_id <- loop_id + 1
    }
  }
  
  # function to calcuate good-of-fitness 
  calc_gof <- function(ep, emp, loop_id, results, parlen, name) {
    res <- ep - emp
    # sample size
    k <- length(ep)
    # number of parameter
    m <- parlen
    aic <- k * log(sum(res^2)/k) + 2*m
    aicc <- aic + 2*m*(m+1)/(k-m-1)
    bic <- k * log(sum(res^2)/k) + m * log(k)
    p <- suppressWarnings(ks.test(emp,ep)$p)
    results[loop_id,] <<- c(aic, aicc, bic, p)
    
    loop_id <<- loop_id + 1
  }
  
  # gumbel
  possibleError <- tryCatch({fitgumb <- dist.start(data,"gumbel")},error = function(e) e)
  if (!inherits(possibleError, "error")) {
    fitgumb <- dist.start(data,"gumbel")
    ep <- actuar::pgumbel(data, alpha = fitgumb$loc, scale = fitgumb$scale)
    pdf <- actuar::dgumbel(data, alpha = fitgumb$loc, scale = fitgumb$scale)
    para_list[[loop_id]] = list( alpha = fitgumb$loc, scale = fitgumb$scale )
    ep_mat[,loop_id] <- ep
    pdf_mat[,loop_id] <- pdf
    name <- "gumbel"
    calc_gof(ep, emp, loop_id, results, 2, name)
  } else {
    results[loop_id,] <- c(NA,NA,NA,NA)
    loop_id <- loop_id + 1
  }
  
  # generalized extreme value
  possibleError <- tryCatch({fitgev <- extRemes::fevd(data, threshold = min(data)-1, type = 'GEV',method = 'Lmoments')},
                            error = function(e) e)
  if (!inherits(possibleError, "error")) {
    fitgev <- extRemes::fevd(data, threshold = min(data)-1, type = 'GEV',method = 'Lmoments')
    gevpar <- fitgev$results
    ep <- evd::pgev(data, shape = as.numeric(gevpar[3]), loc =as.numeric(gevpar[1]), scale = as.numeric(gevpar[2]))
    pdf <- evd::dgev(data, shape = as.numeric(gevpar[3]), loc =as.numeric(gevpar[1]), scale = as.numeric(gevpar[2]))
    para_list[[loop_id]] = list(shape = as.numeric(gevpar[3]), loc =as.numeric(gevpar[1]), scale = as.numeric(gevpar[2]))
    ep_mat[,loop_id] <- ep
    pdf_mat[,loop_id] <- pdf
    name <- "gev"
    calc_gof(ep, emp, loop_id, results, 3, name)
  } else {
    results[loop_id,] <- c(NA,NA,NA,NA)
    loop_id <- loop_id + 1
  }
  
  ## generalized pareto, use min(data) as the threshold of gpd
  possibleError <- tryCatch({fitgpd <- extRemes::fevd(data, threshold = min(data)-1, type = 'GP',method = 'Lmoments')},
                            error = function(e) e)
  if (!inherits(possibleError, "error")) {
    fitgpd <- extRemes::fevd(data, threshold = min(data)-1, type = 'GP',method = 'Lmoments')
    gpdpar <- fitgpd$results
    ep <- evd::pgpd(data, loc = min(data)-1, scale = as.numeric(gpdpar[1]), shape = as.numeric(gpdpar[2]))
    pdf <- evd::dgpd(data, loc = min(data)-1, scale = as.numeric(gpdpar[1]), shape = as.numeric(gpdpar[2]))
    para_list[[loop_id]] = list(loc = min(data)-1, scale = as.numeric(gpdpar[1]), shape = as.numeric(gpdpar[2]))
    ep_mat[,loop_id] <- ep
    pdf_mat[,loop_id] <- pdf
    name <- "gpd"
    calc_gof(ep, emp, loop_id, results, 3, name)
  } else {
    results[loop_id,] <- c(NA,NA,NA,NA)
    loop_id <- loop_id + 1
  }
  
  # Three-parameter Gamma (Pearson Type III)
  possibleError <- tryCatch({fitpe3 <- dist.start(data,"pe3")},error = function(e) e)
  if (!inherits(possibleError, "error")) {
    fitpe3 <- dist.start(data,"pe3")
    ep <- SCI::ppe3(data, shape=fitpe3$shape, scale=fitpe3$scale, location=fitpe3$location)
    pdf <- SCI::dpe3(data, shape=fitpe3$shape, scale=fitpe3$scale, location=fitpe3$location)
    para_list[[loop_id]] = list(shape=fitpe3$shape, scale=fitpe3$scale, location=fitpe3$location)
    ep_mat[,loop_id] <- ep
    pdf_mat[,loop_id] <- pdf
    name <- "pe3"
    calc_gof(ep, emp, loop_id, results, 3, name)
  } else {
    results[loop_id,] <- c(NA,NA,NA,NA)
    loop_id <- loop_id + 1
  }
  
  # inverse gaussian
  possibleError <- tryCatch({fitig <- goft::ig_fit(data)},error = function(e) e)
  if (!inherits(possibleError, "error")) {
    fitig <- goft::ig_fit(data)
    ep <- actuar::pinvgauss(data, mean = fitig[1,1],shape = fitig[2,1])
    pdf <- actuar::dinvgauss(data, mean = fitig[1,1],shape = fitig[2,1])
    para_list[[loop_id]] = list( mean = fitig[1,1],shape = fitig[2,1] )
    ep_mat[,loop_id] <- ep
    pdf_mat[,loop_id] <- pdf
    name <- "invgauss"
    calc_gof(ep, emp, loop_id, results, 2, name)
  } else {
    results[loop_id,] <- c(NA,NA,NA,NA)
    loop_id <- loop_id + 1
  }
  
  results['family'] = c("gamma","exp","weibull","norm","logis","lnorm","llogis",
                        "cauchy","gumbel","gev","gpd","pe3","invgauss")
  
  # reorder results according to p-value and AIC
  results2 <- results[!is.na(results$p),]
  results2 <- results2[order(rank(results2$p), -rank(results2$aic), decreasing = T), ]
  rownames(results2) <- NULL
  
  # output
  PD_name <- results2[results2$p>0.1,"family"][1]
  if (is.na(PD_name)) {
    warning('The fitting of univariate distribution is not significant !')
    PD_name <- results2[1,'family']
  }
  PD_id <- which( results$family == PD_name )
  ep <- ep_mat[,PD_id]
  pdf <- pdf_mat[,PD_id]
  para <- para_list[[PD_id]]
  
  # check whether any NaN in ep and pdf
  loop_id = 1
  while( any(is.na(ep)) | any(is.na(pdf)) ) {
    PD_name <- results2[results2$p>0.1,"family"][loop_id + 1]
    PD_id <- which( results$family == PD_name )
    ep <- ep_mat[,PD_id]
    pdf <- pdf_mat[,PD_id]
  }
  para <- para_list[[PD_id]]
  
  # parameteric bootstrap for uncertainty
  if (PD_name %in% c("gamma","exp","weibull","norm","logis","lnorm","llogis","cauchy")) {
    f1 <- fitdist(data,PD_name)
    bnor <- bootdist(f1)
    resboot <- bnor$estim
  } else if (PD_name=='gumbel') {
    rdata <- actuar::rgumbel(1001*length(data), fitgumb$loc, fitgumb$scale)
    dim(rdata) <- c(length(data), 1001)
    func <- function(iter) {
      res <- dist.start(rdata[, iter],"gumbel")
      return(c(res$loc, res$scale))
    }
    resboot <- t(sapply(1:1001, func))
    resboot <- as.data.frame(resboot); colnames(resboot) = c("loc","scale")
  } else if (PD_name=='gev') {
    gevpar <- fitgev$results
    rdata <- evd::rgev(1001*length(data),shape = as.numeric(gevpar[3]),loc =as.numeric(gevpar[1]),scale = as.numeric(gevpar[2]))
    dim(rdata) <- c(length(data), 1001)
    func <- function(iter) {
      res <- extRemes::fevd(rdata[, iter], threshold = min(rdata)-1, type = 'GEV',method = 'Lmoments')
      par <- res$results
      return(c(as.numeric(par[3]),as.numeric(par[1]),as.numeric(par[2])))
    }
    resboot <- t(sapply(1:1001, func))
    resboot <- as.data.frame(resboot); colnames(resboot) = c("shape","loc","scale")
  } else if (PD_name=='gpd') {
    gpdpar <- fitgpd$results
    rdata <- evd::rgpd(1001*length(data),loc = min(data)-1, scale = as.numeric(gpdpar[1]), shape = as.numeric(gpdpar[2]))
    dim(rdata) <- c(length(data), 1001)
    func <- function(iter) {
      res <- extRemes::fevd(rdata[, iter], threshold = min(rdata)-1, type = 'GP',method = 'Lmoments')
      par <- res$results
      return(c(min(rdata[,iter])-1, as.numeric(par[1]), as.numeric(par[2])))
    }
    resboot <- t(sapply(1:1001, func))
    resboot <- as.data.frame(resboot); colnames(resboot) = c("loc","scale","shape")
  } else if (PD_name=='pe3') {
    rdata <- SCI::rpe3(1001*length(data), shape=fitpe3$shape, scale=fitpe3$scale, location=fitpe3$location)
    dim(rdata) <- c(length(data), 1001)
    func <- function(iter) {
      res <- dist.start(rdata[,iter],"pe3")
      return(c(res$shape, res$scale, res$location))
    }
    resboot <- t(sapply(1:1001, func))
    resboot <- as.data.frame(resboot); colnames(resboot) = c("shape","scale","location")
  } else if (PD_name=='invgauss') {
    rdata <- actuar::rinvgauss(1001*length(data), mean = fitig[1,1],shape = fitig[2,1])
    dim(rdata) <- c(length(data), 1001)
    func <- function(iter) {
      res <- goft::ig_fit(rdata[,iter])
      return(c(res[1,1],res[2,1]))
    }
    resboot <- t(sapply(1:1001, func))
    resboot <- as.data.frame(resboot); colnames(resboot) = c("mean","shape")
  }
  
  return(list("PD_name"=PD_name,"ep"=ep, "pdf"=pdf, "emp"=emp, "STAT"=results2, "para"=para, parboot=resboot))
  
}

# univariate pdf
unipdf <- function(x, dname, par) {
  # total 13 distributions
  
  distnames <- c("norm","logis","lnorm","llogis","weibull","gumbel","exp","gpd","cauchy","gamma","invgauss","gev","pe3")
  
  if (dname=="norm") {
    f <- function(x) {1/(par[2]*sqrt(2*pi))*exp(-(x-par[1])^2/(2*par[2]^2))}
  } else if (dname == "logis") {
    f <- function(x) {exp(-(x-par[1])/par[2])/(par[2]*(1+exp(-(x-par[1])/par[2]))^2)}
  } else if (dname == "lnorm") {
    f <- function(x) {exp(-(log(x)-par[1])^2/(2*par[2]^2))/(x*par[2]*sqrt(2*pi))} 
  } else if (dname == "llogis") {
    f <- function(x) {(par[2]/par[1])*(x/par[1])^(par[2]-1)/(1+(x/par[1])^par[2])^2}
  } else if (dname == "weibull") {
    f <- function(x) {par[2]/par[1]*(x/par[1])^(par[1]-1)*exp(-(x/par[1])^par[2])}
  } else if (dname == "gumbel") {
    f <- function(x) {exp(-(exp(-(x-par[1])/par[2])+(x-par[1])/par[2]))/par[2]}
  } else if (dname == "exp") {
    f <- function(x) {par[1]*exp(-par[1]*x)}
  } else if (dname == "gpd") {
    f <- function(x) {1/par[2]*(1+par[3]*(x-par[1])/par[2])^(-(1+1/par[3]))}
  } else if (dname == "cauchy") {
    f <- function(x) {1/(pi*par[2]*(1+((x-par[1])/par[2])^2))}
  } else if (dname == "gamma") {
    f <- function(x) {dgamma(x,par[1],par[2])}
  } else if (dname == "invgauss") {
    f <- function(x) {(par[1]/(2*pi*x^3))^(1/2)*exp(-par[1]*(x-par[2])^2/(2*(par[2])^2*x))}
  } else if (dname == "gev") {
    if (par[3] == 0) {
      f <- function(x) {1/par[2]*(exp(-(x-par[1])/par[2]))^(par[3]+1)*exp(-exp(-(x-par[1])/par[2]))}
    } else {
      f <- function(x) {1/par[2]*((1+par[3]*((x-par[1])/par[2]))^(-1/par[3]))^(par[3]+1)*exp(-(1+par[3]*((x-par[1])/par[2]))^(-1/par[3]))}
    }
  } else if (dname == "pe3") {
    if (par[3] == 0) {
      f <- function(x) {exp(-((x-par[1])/par[2])^2/2)/sqrt(2*pi)}
    } else {
      if (par[3] > 0) {
        (dgamma((x-(par[1]-2*par[2]/par[3]))/(0.5*par[2]*abs(par[3])), 4/(par[3])^2))/(0.5*par[2]*abs(par[3])) 
      } else {
        (dgamma(((par[1]-2*par[2]/par[3])-x)/(0.5*par[2]*abs(par[3])), 4/(par[3])^2))/(0.5*par[2]*abs(par[3])) 
      } 
    }
  }
  
  return (f)
}

### ----------------MCMC-based bivariate copula estimation---------------
# MCMC-based Bivariate Copula estimation and return tpar, bestfamily, and crdf
ProBiCopSelect <- function(u, v, familyset, show = T) {
  # number of iteration
  IT = 2000
  EP = cbind(u,v)
  
  # calculate empirical bivariate cdf
  EBVP = Calc_Emp_BiVar_Prob(EP)
  
  # save select copula criterion 
  AIC <- as.numeric()
  BIC <- as.numeric()
  AICc <- as.numeric()
  NSE <- as.numeric()
  RMSE <- as.numeric()
  
  # save all tpar, best par, and id_best_par for each family
  tpar_matrix <- array(NA, dim=c(2500, 5, length(familyset)))
  best_par <- matrix(NA,length(familyset),2)
  id_best_par <- numeric()
  
  # LOOP estimate copula parameter using MCMC
  loop_id = 1
  for (i in familyset) {
    PAR_RANGE = list()
    if (i==0) {
      PAR_RANGE$min = NA; PAR_RANGE$max = NA  # Independence Copula
      BiCDF = u * v
      D = 0; n = length(u)
      RES <- BiCDF - EBVP
      AIC[loop_id] <- n * log(sum(RES*RES)/n) + 2 * D
      AICc[loop_id] <- AIC[loop_id] + (2*D^2 + 2*D)/(n-D-1)
      BIC[loop_id] <- n * log(sum(RES*RES)/n) + D * log(n)
      RMSE[loop_id] <- sqrt(mean(RES*RES))
      NSE[loop_id] <- 1 - sum( RES*RES )/sum( (EBVP - mean(EBVP))^2 )
      loop_id <- loop_id + 1
      next
    } else if (i==1) {
      PAR_RANGE$min = -0.9999; PAR_RANGE$max = 0.9999 # Gaussian Copula
    } else if (i==2) {
      PAR_RANGE$min = c(-0.9999,0.0001); PAR_RANGE$max = c(0.9999,35) # t Copula
    } else if (i==3) {
      PAR_RANGE$min = 0.0001; PAR_RANGE$max = 35 # Clayton Copula
    } else if (i==4) {
      PAR_RANGE$min = 1.0001; PAR_RANGE$max = 35 # Gumbel Copula
    } else if (i==5) {
      PAR_RANGE$min = -35; PAR_RANGE$max = 35    # Frank Copula
    } else if (i==6) {
      PAR_RANGE$min = 1+1e-4; PAR_RANGE$max = 35      # Joe Copula
    } else if (i==7) {
      PAR_RANGE$min = c(1e-4,1); PAR_RANGE$max = c(35,35)# BB1 Copula
    } else if (i==8) {
      PAR_RANGE$min = c(1,1); PAR_RANGE$max = c(35,35)# BB6 Copula
    } else if (i==9) {
      PAR_RANGE$min = c(1,1e-4); PAR_RANGE$max = c(35,35)# BB7 Copula
    } else if (i==10) {
      PAR_RANGE$min = c(1,1e-4); PAR_RANGE$max = c(35,1)# BB8 Copula
    } else if (i==13) {
      PAR_RANGE$min = 1e-4; PAR_RANGE$max = 35     # 180_Clayton
    } else if (i==14) {
      PAR_RANGE$min = 1+1e-4; PAR_RANGE$max = 35   # 180_Gumbel
    } else if (i==16) {
      PAR_RANGE$min = 1+1e-4; PAR_RANGE$max = 35   # 180_Joe
    } else if (i==17) {
      PAR_RANGE$min = c(1e-4,1); PAR_RANGE$max = c(35,35)# 180_BB1
    } else if (i==18) {
      PAR_RANGE$min = c(1,1); PAR_RANGE$max = c(35,35) # 180_BB6
    } else if (i==19) {
      PAR_RANGE$min = c(1,1e-4); PAR_RANGE$max = c(35,35) # 180_BB7
    } else if (i==20) {
      PAR_RANGE$min = c(1,1e-4); PAR_RANGE$max = c(35,1) # 180_BB8
    } else if (i==23) {
      PAR_RANGE$min = -35; PAR_RANGE$max = -1e-4    # 90_Clayton
    } else if (i==24) {
      PAR_RANGE$min = -35; PAR_RANGE$max = -1       # 90_Gumbel
    } else if (i==26) {
      PAR_RANGE$min = -35; PAR_RANGE$max = -1-1e-4  # 90_Joe
    } else if (i==27) {
      PAR_RANGE$min = c(-35,-35); PAR_RANGE$max = c(-1e-4,-1) # 90_BB1
    } else if (i==28) {
      PAR_RANGE$min = c(-35,-35); PAR_RANGE$max = c(1,1)      # 90_BB6
    } else if (i==29) {
      PAR_RANGE$min = c(-35,-35); PAR_RANGE$max = c(-1,-1e-4) # 90_BB7
    } else if (i==30) {
      PAR_RANGE$min = c(-35,-1); PAR_RANGE$max = c(-1,-1e-4)  # 90_BB8
    } else if (i==33) {
      PAR_RANGE$min = -35; PAR_RANGE$max = -1e-4    # 270_Clayton
    } else if (i==34) {
      PAR_RANGE$min = -35; PAR_RANGE$max = -1       # 270_Gumbel
    } else if (i==36) {
      PAR_RANGE$min = -35; PAR_RANGE$max = -1-1e-4  # 270_Joe
    } else if (i==37) {
      PAR_RANGE$min = c(-35,-35); PAR_RANGE$max = c(-1e-4,-1) # 270_BB1
    } else if (i==38) {
      PAR_RANGE$min = c(-35,-35); PAR_RANGE$max = c(1,1)      # 270_BB6
    } else if (i==39) {
      PAR_RANGE$min = c(-35,-35); PAR_RANGE$max = c(-1,-1e-4) # 270_BB7
    } else if (i==40) {
      PAR_RANGE$min = c(-35,-1); PAR_RANGE$max = c(-1,-1e-4)  # 270_BB8
    } else if (i==104) {
      PAR_RANGE$min = c(1,1e-4); PAR_RANGE$max = c(35,1)  # Tawn1
    } else if (i==114) {
      PAR_RANGE$min = c(1,1e-4); PAR_RANGE$max = c(35,1)  # 180_Tawn1
    } else if (i==124) {
      PAR_RANGE$min = c(1,1e-4); PAR_RANGE$max = c(35,1)  # 90_Tawn1
    } else if (i==134) {
      PAR_RANGE$min = c(1,1e-4); PAR_RANGE$max = c(35,1)  # 270_Tawn1
    } else if (i==204) {
      PAR_RANGE$min = c(1,1e-4); PAR_RANGE$max = c(35,1)  # Tawn2
    } else if (i==214) {
      PAR_RANGE$min = c(1,1e-4); PAR_RANGE$max = c(35,1)  # 180_Tawn2
    } else if (i==224) {
      PAR_RANGE$min = c(1,1e-4); PAR_RANGE$max = c(35,1)  # 90_Tawn2
    } else if (i==234) {
      PAR_RANGE$min = c(1,1e-4); PAR_RANGE$max = c(35,1)  # 270_Tawn2
    }
    
    D = length(PAR_RANGE$min)
    n <- nrow(EP)
    
    # estimate parameter using MCMC 
    Family = i
    MCMC_out = MCMC(PAR_RANGE, IT, Family, EP, EBVP, 0, show)
    
    # only save last 2500 samples
    tpar <- MCMC_out$TPAR
    tpar <- tpar[(nrow(tpar)-2499):nrow(tpar),1:(D+2), drop=F]
    tpar_matrix[,1:(D+2),loop_id] <- tpar[,1:(D+2)]
    
    # find the best sample with the maximum likelihood
    id_best_par[loop_id] = which.max(tpar[,D+2])
    best_par[loop_id,1:D] = tpar[id_best_par[loop_id],1:D, drop=F]
    
    # calculate fitted bivariate cdf using the best parameter
    if (D==1) {
      BiCDF <- BiCopCDF(EP[,1],EP[,2], i, par = as.numeric(best_par[loop_id,1]), par2 = 0)
    } else {
      BiCDF <- BiCopCDF(EP[,1],EP[,2], i, par = as.numeric(best_par[loop_id,1]), par2 = as.numeric(best_par[loop_id,2]))
    }
    
    # calculate RES, AIC, BIC, RMSE, NSE
    RES <- BiCDF - EBVP
    
    AIC[loop_id] <- n * log(sum(RES*RES)/n) + 2 * D
    AICc[loop_id] <- AIC[loop_id] + (2*D^2 + 2*D)/(n-D-1)
    BIC[loop_id] <- n * log(sum(RES*RES)/n) + D * log(n)
    RMSE[loop_id] <- sqrt(mean(RES*RES))
    NSE[loop_id] <- 1 - sum( RES*RES )/sum( (EBVP - mean(EBVP))^2 )
    loop_id <- loop_id + 1
  }
  
  # if the criterion is inconsistent, select based on AIC
  BestFamily <- familyset[which(AIC==min(AIC))]
  BestFamilyId <- which(familyset == BestFamily)
  
  crdf <- data.frame(familyset,AIC,AICc,BIC,RMSE,NSE)
  
  OUT <- list()
  if (BestFamily != 0) {
    OUT$tpar <- tpar_matrix[, 1:(D+2), BestFamilyId]
    OUT$id_best_par <- id_best_par[BestFamilyId]
  } else {
    OUT$tpar <- NULL
    OUT$id_best_par <- NULL
  }
  OUT$bf <- BestFamily
  OUT$crdf <- crdf
  return (OUT)
  
}

# Hybrid MCMC algorithm
MCMC <- function(PAR_RANGE,IT,Family,EP,EBVP,Parallel,show){
  # return  [CH, TPAR, OLR, Z, Rhat]
  # Notations
  # CH: Chain
  # PRP: Propsal sample
  # PAR.min, PAR.max: Min & Max range of parameters
  # N: Number of chains (min 5)
  # D: Dimension of the problem
  # IT: Number of iterations
  # LHS: Latin Hypercube Samples
  # LN: Number of LHS
  # PDV: Probability density value
  # L_PDV: Log of PDV
  # MR: Mtropolis ratio
  # Z: Past Samples
  
  ## Note
  # This MCMC code is transformed from the Matlab code of Mojtaba Sadegh et al. (2017)
  # Detailed reference is below
  # Sadegh, Mojtaba, Elisa Ragno, and Amir AghaKouchak. "Multivariate C opula A nalysis T oolbox (MvCAT): Describing dependence and underlying uncertainty using a B ayesian framework." Water Resources Research 53.6 (2017): 5166-5183.
  
  require(DiceDesign)
  require(MASS)
  require(pracma)
  require(geometry)
  # Initialize Rhat (convergence statistics)
  Rhat = as.numeric()
  
  # Number of parameters (problem dimension)
  D = length(PAR_RANGE$min)
  # Number of chains
  N = max(2*D, 5)
  # Initialize Chains with NAN
  CH = array(NA, dim = c(N, D+2, IT)) # rows: different chains; columns: parameter set; 3rd dim: iterations
  # Assure LN is multiplicative of N and at least 30*D
  LN = max(N, 30*D)
  LN = ceiling(LN/N) 
  LN = N * LN
  
  ## Start chains with Latin Hypercube Sampling: LN samples
  LHS = matrix( rep(PAR_RANGE$min, LN), LN,length(PAR_RANGE$min),byrow = T) + 
    lhsDesign(LN, D)$design * matrix(rep(PAR_RANGE$max - PAR_RANGE$min, LN),LN,length(PAR_RANGE$min),byrow = T)
  
  # Cohort of Past Samples
  Z = LHS
  # Compute probabilty density value
  LHS = cbind(LHS,Calc_PDV(LHS,Family,EP,EBVP,Parallel))
  
  # Total Parameter Sets
  TPAR = LHS
  
  # Covariance of LHS
  COV = cov(as.matrix(TPAR[,1:D]))
  
  ## Divide LHS into 5 random comlexes: Based on SCE of Duan et al, 1993
  # Create randomly distributed numbers between 1 & LN
  c_idx = sample(1:LN,LN)
  # Divide indices into N complexes
  c_idx = matrix(c_idx, nrow = LN/N)
  
  # create matrix to save outliers for each iteration
  OLR <- zeros(IT,N)
  
  # Select best candidate in each complex to be the starting point of each chain
  for (i in 1:N) {
    # Find indices of max likelihood in each complex
    id = which.max(LHS[c_idx[,i], ncol(LHS)])
    # Assign max likelihood samples to chains
    CH[i, 1:(D+2), 1] = LHS[c_idx[id,i],]
  }
  
  ## Go through time: Based on ter Braak 2006, 2008
  
  for (t in 2:IT) {
    
    # show process
    if (show) {
      if (t %% 100 ==0) {
        print(paste0('Process iteration of ',t))
      }
    }
    
    # Initialize jumps for each chain
    dPRP = zeros(N, D)
    PRP = matrix(NA, N, D+2)
    id_CH = matrix(NA, N, 3)
    
    # Use snooker update with a probability of 0.1: remark 5, page 439 of ter Braak 2008
    Snooker = runif(1,0,1) < 0.1
    
    # Update covariance each 20 iterations: remark 3, page 226 of Haario 2001
    if (t %% 20 == 0) {
      # Covariance of LHS based on the last 50% of chains
      S_T = nrow(TPAR)
      COV = cov(TPAR[(S_T/2):S_T, 1:D,drop=F])
    }
    
    ## Go through the loop for different chains: "Algorithm for updating population X" page 438 of ter Braak 2008
    for (i in 1:N) {
      # Determine indices of the chains used for DE: Eq 1 of ter Braak 2008, page 437
      id_CH[i,1:3] = randsample(1:nrow(Z) , 3, replacement = F)
      
      ## Parallel update
      if (!Snooker) {
        ## Perform Subspace Sampling
        # Randomly select the subspace for sampling
        SS = which( runif(D,0,1) <= runif(1, 0, 1) )
        
        # Shouldn't be empty
        if (length(SS)==0) { SS = sample(1:D,1) }
        # Size of the subspace
        d = length(SS)
        
        #First chain will follow adaptive metropolis
        if (i <= 1) {
          #Create a proposal candidate: Eq.3, page 3 of Roberts 2009
          # Small positive constant
          Beta = runif(1, 0, 0.1)
          dPRP[i, SS] = (1 - Beta) * mvrnorm(mu = zeros(1, d), Sigma = (2.38^2/d)*COV[SS,SS] ) + 
            Beta * mvrnorm( mu = zeros(1, d), Sigma = (0.1^2/d)*eye(d) )
        }
        
        ## Next chains based on DE-MC 
        if (i >= 2) {
          # Select gamma with 10% probability of gamma = 1 or 0.98 for direct jumps: first line on the right in page 242 & first
          # paragraph on the right in page 248 of ter Braak 2006
          gamma = randsample(c(2.38/sqrt(2*d), 0.98), 1, c(0.9, 0.1), replacement = T)
          
          # The difference between the two chains: Eq 2, page 241 of ter Braak 2006
          dCH = Z[id_CH[i, 1], 1:D,drop=F] - Z[id_CH[i, 2], 1:D, drop=F]
          
          # Select jump: adopted from Eq 2, page 241 of ter Braak 2006 and E4, page 274 of Vrugt 2009
          dPRP[i, SS] =  runif(1, gamma-0.1, gamma+0.1) * dCH[SS] + rnorm(d, 0, 1e-12)
        }
        
        ## Snooker Update
      } 
      else {
        # Find the direction of the update: "Algorithm of DE Snooker update (Fig. 3)", page 438 of ter Braak 2008
        DIR = CH[i, 1:D, t-1] - Z[id_CH[i, 1], 1:D]
        
        # Project vector a onto b (https://en.wikipedia.org/wiki/Vector_projection):
        # a.b/b.b * b
        
        # Difference between z1 and z2 and its length on DIR
        DIF = ( (Z[id_CH[i, 2], 1:D] - Z[id_CH[i, 3], 1:D]) %*% DIR ) / ( DIR %*% DIR )
        
        DIF = pmax(DIF, 0, na.rm=T)
        
        # Resize DIR
        dCH = DIR %*% DIF
        
        # Select jump: page 439 of ter Braak, 2008
        dPRP[i, 1:D] = runif(1,1.2,2.2) * dCH
        
      }
      
    }
    
    ## Create a proposal candidate: Current + Jump 
    PRP[1:N, 1:D] = CH[1:N, 1:D, t-1] + dPRP[1:N, 1:D,drop=F]
    
    ## Boundary handling: fold back in if fall outside: General Knowledge 
    PRP[1:N , 1:D] = Bound_handle(PRP[1:N, 1:D], PAR_RANGE)
    
    ## Snooker correction for the metropolis ratio: Eq 4, page 439 of ter Braak, 2008
    C_sn = ones(N,1) # No need to correct if parallel updating
    if (Snooker) {
      DSS = PRP[1:N , 1:D] - Z[id_CH[1:N,1], 1:D] # nominator difference of Eq 4
      DS  = CH[1:N, 1:D, t-1] - Z[id_CH[1:N,1], 1:D] # denominator difference of Eq 4
      # Dot function yields norm of each chain in matrix of all chains
      C_sn = ( geometry::dot(DSS,DSS,2) / geometry::dot(DS,DS,2) )^((D-1)/2)
      C_sn = matrix(C_sn,N,1)
    }
    
    ## Calculate Likelihood 
    
    PRP = cbind(PRP[1:N,1:D],Calc_PDV(PRP[1:N, 1:D], Family, EP, EBVP, Parallel))
    
    # Compute Metropolis ratio: a/b = exp( log(a) - log(b) )
    MR = pmin(1, as.matrix(C_sn * exp(PRP[1:N, ncol(PRP)] - CH[1:N, dim(CH)[2], t-1])) )
    
    # Accept/reject the proposal point
    id_accp = which( MR >= runif(N,0,1) )
    
    id_rjct = setdiff(1:N, id_accp)
    
    ## Accept proposal and update chain
    if (length(id_accp)>0) {
      CH[id_accp, 1:(D+2), t] = PRP[id_accp, 1:(D+2)]
    }
    
    # Reject proposal and remain at the previous state
    if (length(id_rjct)>0) {
      CH[id_rjct, 1:(D+2), t] = CH[id_rjct, 1:(D+2), t-1]
      dPRP[id_rjct, 1:D] = 0
    }
    
    # Check for outliers: remark 3, page 275 of Vrugt 2009
    olr = OUTLIER(CH, t)
    if (length(olr) > 0) {
      OLR[t,olr] = 1
    }
    if (D > 3 & length(olr) > 0) {
      id_olr = order(CH[,dim(CH)[2],t-1],decreasing = T)
      CH[olr,,t] = CH[id_olr[1:length(olr)],,t-1]
    }
    
    ## Total Parameter Sets
    TPAR = rbind(TPAR, CH[1:N, 1:(D+2), t])
    
    ## Check convergence, every 10 steps
    if (t %% 10 == 0) {
      dummy = Convergence_GR( CH[1:N, 1:D, 1:t,drop=F] )
      Rhat = rbind(Rhat, cbind(t, dummy))
    }
    
    ## Update Z each 10 iterations
    if (t %% 10 == 0) {
      Z = rbind(Z, as.matrix(CH[1:N, 1:D, t]))
    }
    
  }
  
  # create output list
  OUT = list()
  OUT$CH <- CH
  OUT$TPAR <- TPAR
  OUT$OLR <- OLR
  OUT$Z <- Z
  OUT$Rhat <- Rhat
  
  return (OUT)
  
}

## Calculate probabilty density value
# Log likelihood of residuals, e, if they are uncorrelated, homoskedastic, gaussian distributed with zero mean
# log(L) = -n/2 * log(2pi) - n/2 * log(sigma^2) - 1/(2*sigma^2) * sum(e^2)
# sigma^2 = sum(e^2)/n ==> 1/(2*sigma^2) * sum(e^2) = n/2
# ==> log(L) = -n/2 - n/2 * log(2pi) - n/2 * log(sigma^2)
# ==> log(L) = Constant              - n/2 * log(sum(e^2)/n) 
Calc_PDV <- function(X,Family,EP,EBVP,Parallel) {
  require(CDVine)
  
  X = as.matrix(X)
  PDV <- as.numeric()
  L_PDV <- as.numeric()
  
  if (dim(X)[2]==1) {
    par1 = as.numeric(X[,1])
    par2 = rep(0,length(par1))
  } else {
    par1 = as.numeric(X[,1])
    par2 = as.numeric(X[,2])
  }
  
  for ( i in 1:(dim(X)[1]) ) {
    # Simulate copula probabilities
    BiCDF <- BiCopCDF(EP[,1],EP[,2],Family, par = par1[i], par2 = par2[i])
    RES <- BiCDF - EBVP
    # Assign RES vector to the PDV values (trick to decrease simuations later)
    PDV[i] = NA
    # Compute log-likelihood
    L_PDV[i] = -length(RES)/2 * log( sum(RES^2)/length(RES) )
  }
  
  OUT <- as.matrix(cbind(PDV,L_PDV))
  
  return (OUT)
  
}

# Compute empirical bivariate probability distribution
Calc_Emp_BiVar_Prob <- function(D) {
  # Length of data
  D = as.matrix(D)
  n = dim(D)[1]
  # Pre-assign bivariate probability array
  BVP = zeros(n,1)
  
  # Loop through the data
  for (i in 1:n) {
    # Pre-assign CD in each time step
    CD = zeros(n,3)
    # Find count of data points for which p(X < x)
    CD[D[,1] <= D[i,1], 1]  = 1
    # Find count of data points for which p(Y < y)
    CD[ D[,2] <= D[i,2], 2 ] = 1
    # Find count of data points for which p(X < x & Y < y)
    CD[,3] = CD[,1] * CD[,2]
    BVP[i,1] = sum( CD[,3] )
  }
  
  # Gringorten plotting position (bivariate)
  Y = (BVP - 0.44) / (n + 0.12)
  
  return (Y)
}

## Boundary handling: General Knowledge
Bound_handle <- function(PRP, PAR_RANGE) {
  PRP <- as.matrix(PRP)
  for ( i in 1:(dim(PRP)[1]) ) {
    # Find parameter values less than minimum range
    id = which(PRP[i,] < PAR_RANGE$min)
    # Fold them back in the range       make sure not to fall outside of the other limit!
    if (!is.null(id) & length(id)>0) {
      PRP[i, id] = PAR_RANGE$min[id] + pmin( abs(PAR_RANGE$min[id] - PRP[i, id]), PAR_RANGE$max[id] - PAR_RANGE$min[id] )
    }
    # Find parameter values larger than maximum range
    id = which(PRP[i ,] > PAR_RANGE$max)
    # Fold them back in the range       make sure not to fall outside of the other limit!
    if (!is.null(id) & length(id)>0) {
      PRP[i, id] = PAR_RANGE$max[id] - pmin( abs(PAR_RANGE$max[id] - PRP[i, id]), PAR_RANGE$max[id] - PAR_RANGE$min[id] )
    }
  }
  return (PRP)
}

## Check for outliers: remark 3, page 275 of Vrugt 2009
OUTLIER <- function(CH, t) {
  # Compute mean of the loglikelihood for each chain
  L = apply( as.array(CH[, dim(CH)[2], (floor(t/2)):(t-1)]), 1, mean)
  # Compute interquartile
  iqr = IQR(L,type=5)
  # First Quartile
  Q1 = quantile(L, 0.25, type=5)
  # Find out if any chain is outlier
  OLR = which( CH[, dim(CH)[2], t] < (Q1 - 2*iqr) )
  return (OLR)
  
}

Convergence_GR <- function(CH) {
  # References:
  # Gelman, Andrew, and Donald B. Rubin. "Inference from iterative simulation using multiple sequences." Statistical science (1992): 457-472.
  # Brooks, Stephen P., and Andrew Gelman. "General methods for monitoring convergence of iterative simulations." Journal of computational and graphical statistics 7.4 (1998): 434-455.
  
  # Determine number of iterations and take the second half of them
  n_iter = dim(CH)[3]
  
  CH = CH[ , , (round(n_iter/2)):n_iter,drop=F]
  
  
  # Update n_iter
  n_chain = dim(CH)[1]
  n_par = dim(CH)[2]
  n_iter = dim(CH)[3]
  
  # Compute mean and variances of each chain, as well variance of means
  MEAN = apply(CH, c(1,2), mean)
  MEAN_VAR = apply(MEAN,2,var)
  VAR = apply(CH,c(1,2),var)
  
  # Compute B (between) and W (within chain variances); following page 436 of ref 2
  B = n_iter * MEAN_VAR
  W = mean(VAR)
  
  # Compute sigma^2: first eq of page 437 of ref 2
  S2 = (n_iter-1)/n_iter * W + B/n_iter
  
  # Compute rhat for univariate analysis
  Rhat = sqrt( (n_chain+1)/n_chain * S2/W - (n_iter-1)/(n_chain*n_iter) )
  
  return (Rhat)
}

### --------------deterministic 3-D copula function-----------------------
# calculate joint cdf values for samples x1, x2, and x3
DVine3CDF <- function(x1, x2, x3, VarOrder=VarOrder, familyset = familyset) {
  
  # determine univariate pdf
  P1 <- fitalldist(x1)
  P2 <- fitalldist(x2)
  P3 <- fitalldist(x3)
  
  # extract cdf parameter
  para1 = P1$para
  para2 = P2$para
  para3 = P3$para
  
  # calculate cdf value
  cdf1name <- paste("p",P1$PD_name,sep = "")
  cdf2name <- paste("p",P2$PD_name,sep = "")
  cdf3name <- paste("p",P3$PD_name,sep = "")
  
  cdf1 <- P1$ep
  cdf2 <- P2$ep
  cdf3 <- P3$ep
  
  # pdf name
  pdf1name <- paste("d",P1$PD_name,sep = "")
  pdf2name <- paste("d",P2$PD_name,sep = "")
  pdf3name <- paste("d",P3$PD_name,sep = "")
  
  invcdf1name <- paste("q",P1$PD_name,sep = "")
  invcdf2name <- paste("q",P2$PD_name,sep = "")
  invcdf3name <- paste("q",P3$PD_name,sep = "")
  
  # determine bivariate pdf
  data <- cbind(cdf1,cdf2,cdf3)
  data <- data[,VarOrder]
  
  # construct RVine structure matrix
  Matrix <- c(1,3,2,
              0,3,2,
              0,0,2)
  Matrix <- matrix(Matrix,3,3)
  PPP <- RVineCopSelect(data,familyset,Matrix)
  
  tpar12 = PPP$par[3,1]
  tpar23 = PPP$par[3,2]
  tpar132 = PPP$par[2,1]
  # sample ID of tpar12 and tpar23
  
  family12 = PPP$family[3,1]
  family23 = PPP$family[3,2]
  family132 = PPP$family[2,1]
  
  ## construct joint PDF
  pdf1 <- function(w) {do.call(pdf1name,c(list(w),as.list(para1)))}
  pdf2 <- function(w) {do.call(pdf2name,c(list(w),as.list(para2)))}
  pdf3 <- function(w) {do.call(pdf3name,c(list(w),as.list(para3)))}
  
  cdf1 <- function(w) {do.call(cdf1name,c(list(w),as.list(para1)))}
  cdf2 <- function(w) {do.call(cdf2name,c(list(w),as.list(para2)))}
  cdf3 <- function(w) {do.call(cdf3name,c(list(w),as.list(para3)))}
  
  x123 <- cbind(x1,x2,x3)
  
  # calculate cdf i is the row of x123
  calc_cdf <- function(i) {
    
    par12 = tpar12
    par23 = tpar23
    par132 = tpar132
    upper <- as.vector(x123[i,])
    
    # construct decomposition pdf function
    coppdf12 <- function(a,b) {
      cop <- BiCop(family = family12,par = par12)
      return(VineCopula::BiCopPDF(cdf1(a),cdf2(b),cop))}
    coppdf23 <- function(a,b) {
      cop <- BiCop(family = family23,par = par23)
      return(VineCopula::BiCopPDF(cdf2(a),cdf3(b),cop))}
    
    h12 <- function(a,b) {VineCopula::BiCopHfunc2(cdf1(a), cdf2(b), family = family12, par = par12)}
    h23 <- function(a,b) {VineCopula::BiCopHfunc1(cdf2(a), cdf3(b), family = family23, par = par23)}
    coppdf132 <- function(a,b,c) {
      cop <- BiCop(family = family132,par = par132)
      return(VineCopula::BiCopPDF(h12(a,b),h23(b,c),cop))}
    
    # pair-copula decomposition of joint pdf
    JointPdf <- function(arg) { 
      x <- arg[1]
      y <- arg[2]
      z <- arg[3]
      ff <- pdf1(x) * pdf2(y) * pdf3(z) * coppdf12(x,y) * coppdf23(y,z) * coppdf132(x,y,z) 
      return (ff)
    }
    
    cdf <- hcubature(JointPdf, lowerLimit = rep(0, 3), upperLimit = upper, tol = 1e-3, maxEval = 10^4)
    
    cdf = cdf$integral
    
    return (cdf)
    
  }
  
  # parallel calculate cdf
  # Calculate the number of cores
  no_cores <- detectCores() - 1
  cl <- makeCluster( no_cores )
  registerDoParallel(cl)
  
  t1 = Sys.time()
  allcdf = foreach(i=1:length(x1), .combine='c',.packages = c("VineCopula","cubature","actuar",
                                                              "SCI","evd")) %dopar% calc_cdf(i)
  t2 = Sys.time()
  print(t2-t1)
  
  stopCluster(cl)
  
  allcdf = data.frame(id=1:length(x1),cdf=allcdf)
  
  return (allcdf)
  
}

# calcualte all cdf and pdf for x,y,z in 0.1 step between [0,1]
DVine3CdfPdf <- function(x1, x2, x3, VarOrder=VarOrder, familyset = familyset) {
  
  # determine univariate pdf
  P1 <- fitalldist(x1)
  P2 <- fitalldist(x2)
  P3 <- fitalldist(x3)
  
  # calculate cdf value
  cdf1name <- paste("p",P1$PD_name,sep = "")
  cdf2name <- paste("p",P2$PD_name,sep = "")
  cdf3name <- paste("p",P3$PD_name,sep = "")
  PAR1 <- P1$para; PAR2 <- P2$para; PAR3 <- P3$para
  cdf1 <- do.call(cdf1name,c(list(x1),as.list(PAR1)))
  cdf2 <- do.call(cdf2name,c(list(x2),as.list(PAR2)))
  cdf3 <- do.call(cdf3name,c(list(x3),as.list(PAR3)))
  
  # pdf name
  pdf1name <- paste("d",P1$PD_name,sep = "")
  pdf2name <- paste("d",P2$PD_name,sep = "")
  pdf3name <- paste("d",P3$PD_name,sep = "")
  
  invcdf1name <- paste("q",P1$PD_name,sep = "")
  invcdf2name <- paste("q",P2$PD_name,sep = "")
  invcdf3name <- paste("q",P3$PD_name,sep = "")
  
  # determine bivariate pdf
  data <- cbind(cdf1,cdf2,cdf3)
  data <- data[,VarOrder]
  
  # construct RVine structure
  Matrix <- c(1,3,2,
              0,3,2,
              0,0,2)
  Matrix <- matrix(Matrix,3,3)
  PPP <- RVineCopSelect(data,familyset,Matrix)
  
  tpar12 = PPP$par[3,1]
  tpar23 = PPP$par[3,2]
  tpar132 = PPP$par[2,1]
  # sample ID of tpar12 and tpar23
  
  family12 = PPP$family[3,1]
  family23 = PPP$family[3,2]
  family132 = PPP$family[2,1]
  
  ## construct joint PDF
  pdf1 <- function(w) {do.call(pdf1name,c(list(w),as.list(PAR1)))}
  pdf2 <- function(w) {do.call(pdf2name,c(list(w),as.list(PAR2)))}
  pdf3 <- function(w) {do.call(pdf3name,c(list(w),as.list(PAR3)))}
  
  cdf1 <- function(w) {do.call(cdf1name,c(list(w),as.list(PAR1)))}
  cdf2 <- function(w) {do.call(cdf2name,c(list(w),as.list(PAR2)))}
  cdf3 <- function(w) {do.call(cdf3name,c(list(w),as.list(PAR3)))}
  
  coppdf12 <- function(a,b) {
    cop <- BiCop(family = family12,par = tpar12)
    return(BiCopPDF(cdf1(a),cdf2(b),cop))}
  coppdf23 <- function(a,b) {
    cop <- BiCop(family = family23,par = tpar23)
    return(BiCopPDF(cdf2(a),cdf3(b),cop))}
  h12 <- function(a,b) {BiCopHfunc2(cdf1(a),cdf2(b),family = family12, par = tpar12)}
  h23 <- function(a,b) {BiCopHfunc1(cdf2(a),cdf3(b),family = family23, par = tpar23)}
  coppdf132 <- function(a,b,c) {
    cop <- BiCop(family = family132,par = tpar132)
    return(BiCopPDF(h12(a,b),h23(b,c),cop))}
  
  JointPdf <- function(arg) { 
    x <- arg[1]
    y <- arg[2]
    z <- arg[3]
    ff <- pdf1(x) * pdf2(y) * pdf3(z) * coppdf12(x,y) * coppdf23(y,z) * coppdf132(x,y,z) 
    return (ff)
  }
  
  # calculate cdf i is the row of tpar and x123
  calc_cdf <- function(i) {
    upper <- as.vector(x123[i,])
    cdf <- hcubature(JointPdf, lowerLimit = rep(0, 3), upperLimit = upper, tol = 1e-3, maxEval = 10^4)
    cdf = cdf$integral
    return (cdf)
  }
  
  # integrate to cdf
  a <-  rep( rep( seq(0.05,0.95,0.1), each=10) , times=10 )
  b <- rep( seq(0.05,0.95,0.1), times=100 )
  c <- rep( seq(0.05,0.95,0.1), each=100 )
  
  raw1 <- do.call(invcdf1name, c(list(a),as.list(PAR1)))
  raw2 <- do.call(invcdf2name, c(list(b),as.list(PAR2)))
  raw3 <- do.call(invcdf3name, c(list(c),as.list(PAR3)))
  
  x123 <- cbind(raw1,raw2,raw3)
  
  # parallel computing cdf
  no_cores <- detectCores() - 1
  registerDoSEQ()
  cl <- makeCluster( no_cores )
  registerDoParallel(cl)
  allcdf = foreach(i=1:nrow(x123), .combine='c',.packages = c("VineCopula","cubature","actuar",
                                                              "SCI","evd")) %dopar% calc_cdf(i)
  stopCluster(cl)
  
  # calculate pdf
  Matrix <- c(1,3,2,
              0,3,2,
              0,0,2)
  Matrix <- matrix(Matrix,3,3)
  family <- c(0, family132, family12,
              0,        0,  family23,
              0,        0,         0)
  family <- matrix(family, 3,3)
  par <- c(0, tpar132, tpar12,
           0,      0, tpar23,
           0,      0,     0)
  par <- matrix(par, 3,3)
  par2 <- matrix(0, 3,3)
  RVM <- RVineMatrix(Matrix = Matrix, family = family, par = par, par2 = par2,
                     names = c("V1", "V2", "V3"))
  allpdf = RVinePDF(cbind(a,b,c), RVM)
  
  out = cbind(a,b,c,raw1,raw2,raw3,allcdf,allpdf)
  
  return (out)
  
}

### --------------MCMC-based 3-D Copula function-----------------------
# Compute empirical bivariate probability distribution
Calc_Emp_TriVar_Prob <- function(D) {
  # Length of data
  D = as.matrix(D)
  n = dim(D)[1]
  
  # register parallel backend
  registerDoSEQ()
  no_cores <- detectCores() - 1
  cl <- makeCluster( no_cores )
  registerDoParallel(cl)
  
  # Loop through the data
  emp <- function(i) {
    # Pre-assign CD in each time step
    CD = zeros(n,4)
    # Find count of data points for which p(X < x)
    CD[D[,1] <= D[i,1], 1]  = 1
    # Find count of data points for which p(Y < y)
    CD[D[,2] <= D[i,2], 2] = 1
    # Find count of data points for which p(Z < z)
    CD[D[,3] <= D[i,3], 3]  = 1
    # Find count of data points for which p(X < x & Y < y & Z < z)
    CD[,4] = CD[,1] * CD[,2] * CD[,3]
    return (sum( CD[,4] ))
  }
  TVP <- foreach(i=1:n, .combine = 'c', .packages = c("pracma")) %dopar% emp(i)
  
  # stop parallel backend
  stopCluster(cl)
  
  # Gringorten plotting position (bivariate)
  Y = (TVP - 0.44) / (n + 0.12)
  
  return (Y)
}

# Probabilistic 3-variable pair-copula selection and estimation for R-Vine Models using sequential MCMC
# SampleNum is number of parameter estimates
ProDVine3Est <- function(x1, x2, x3, familyset,SampleNum) {
  
  # determine univariate pdf
  P1 <- fitalldist(x1)
  P2 <- fitalldist(x2)
  P3 <- fitalldist(x3)
  
  # calculate cdf value
  cdf1name <- paste("p",P1$PD_name,sep = "")
  cdf2name <- paste("p",P2$PD_name,sep = "")
  cdf3name <- paste("p",P3$PD_name,sep = "")
  
  cdf1 <- P1$ep
  cdf2 <- P2$ep
  cdf3 <- P3$ep
  
  # pdf name
  pdf1name <- paste("d",P1$PD_name,sep = "")
  pdf2name <- paste("d",P2$PD_name,sep = "")
  pdf3name <- paste("d",P3$PD_name,sep = "")
  
  invcdf1name <- paste("q",P1$PD_name,sep = "")
  invcdf2name <- paste("q",P2$PD_name,sep = "")
  invcdf3name <- paste("q",P3$PD_name,sep = "")
  
  # reorder data
  data <- cbind(cdf1,cdf2,cdf3)
  D <- ncol(data)
  n <- nrow(data)
  
  # estimate tree 1 sequentially
  UncProBi12 <- ProBiCopSelect(data[,1], data[,2], familyset)
  UncProBi23 <- ProBiCopSelect(data[,2], data[,3], familyset)
  
  ## calculate Hfun for conditional copula fitting of var1 and var2
  tpar12 <- UncProBi12$tpar
  family12 <- UncProBi12$bf
  
  # single or double paramter
  if ( is.na(tpar12[1,2])) {
    par1 = tpar12[,1]
    par2 = rep(0,nrow(tpar12))
  } else {
    par1 = tpar12[,1]
    par2 = tpar12[,2]
  }
  
  # register parallel backend
  registerDoSEQ()
  no_cores <- detectCores() - 1
  cl <- makeCluster( no_cores )
  registerDoParallel(cl)
  
  # loop calcualte Hfun2 i.e., H2(data1|data2,theta)
  Hfunc2 <- function(i) {BiCopHfunc2(data[,1],data[,2],family12,par=par1[i],par2 = par2[i])}
  ProHfunc12 <- foreach(i=1:2500, .combine='rbind',.packages = c("VineCopula")) %dopar% Hfunc2(i)
  
  ## calculate Hfun for conditional copula fitting of var2 and var3
  tpar23 <- UncProBi23$tpar
  family23 <- UncProBi23$bf
  
  # single or double paramter
  if ( is.na(tpar23[1,2]) ) {
    par1 = tpar23[,1]
    par2 = rep(0,nrow(tpar23))
  } else {
    par1 = tpar23[,1]
    par2 = tpar23[,2]
  }
  
  # loop calcualte Hfun1 i.e., H1(data3|data2,theta)
  Hfunc1 <- function(i) {BiCopHfunc1(data[,2],data[,3],family23,par=par1[i],par2 = par2[i])}
  ProHfunc32 <- foreach(i=1:2500, .combine='rbind',.packages = c("VineCopula")) %dopar% Hfunc1(i)
  
  # stop parallel backend
  stopCluster(cl)
  
  ## MCMC-based estimate unconditional copula fitting of C_13|2(v1|v2,v3|v2)
  # determine the best family using the best parameter
  BestHfunc12 <- ProHfunc12[ UncProBi12$id_best_par , ]
  BestHfunc32 <- ProHfunc32[ UncProBi23$id_best_par , ]
  # calculate the empirical bivariate cdf of x1|x2 and x3|x2
  EBVP132 = Calc_Emp_BiVar_Prob(cbind(BestHfunc12,BestHfunc32))
  # determine best family of copula 132
  ConProBi132 <- ProBiCopSelect( BestHfunc12, BestHfunc32, familyset)
  bestFamily132 <- ConProBi132$bf
  
  # select SampleNum samples from ProHfunc12 and ProHfunc23 during CI95
  q <- seq(0, 1, length.out = SampleNum)
  sampleID12 <- unique(sapply(quantile(ProHfunc12[,1],q), function(y)which.min(abs(ProHfunc12[,1] - y))))
  sampleID23 <- unique(sapply(quantile(ProHfunc32[,1],q), function(y)which.min(abs(ProHfunc32[,1] - y))))
  
  if (nrow(ProHfunc12) > SampleNum) {
    ProHfunc12 <- ProHfunc12[sampleID12,]
    Par12 <- tpar12[sampleID12,]
  }
  if (nrow(ProHfunc32) > SampleNum) {
    ProHfunc32 <- ProHfunc32[sampleID23,]
    Par23 <- tpar23[sampleID23,]
  }
  
  ProHfunc12 <- t(ProHfunc12)
  ProHfunc32 <- t(ProHfunc32)
  # resample and add 4 cases: Par12.min, Par12.max, Par23.min, Par23.max
  sample12 <- sample.int(SampleNum)
  ProHfunc12 <- ProHfunc12[,c(sample12, rep(which.min(Par12[,1]),2), rep(which.max(Par12[,1]),2))]
  Par12 <- Par12[c(sample12, rep(which.min(Par12[,1]),2), rep(which.max(Par12[,1]),2)),]
  
  sample23 <- sample.int(SampleNum)
  ProHfunc32 <- ProHfunc32[,c(sample23,rep(c(which.min(Par23[,1]), which.max(Par23[,1])),2))]
  Par23 <- Par23[c(sample23,rep(c(which.min(Par23[,1]), which.max(Par23[,1])),2)),]
  
  ## loop estimate copula C_13|2 using MCMC and return 100x100x2500 parameter samples
  # create array to save parameter samples
  familyset <- bestFamily132
  
  # MCMC-based Bivariate Copula estimation and return tpar, bestfamily, and crdf
  tpar132 <- matrix(NA,(SampleNum+4)*2500, 4)
  tpar12_rep = matrix(NA,(SampleNum+4)*2500, 4)
  tpar23_rep = matrix(NA,(SampleNum+4)*2500, 4)
  for (i in 1:(SampleNum+4)) {
    ProBi132 <- ProBiCopSelect(ProHfunc12[,i], ProHfunc32[,i], familyset, show = F)
    tpar <- ProBi132$tpar
    tpar132[(2500*(i-1)+1):(i*2500), 1:(dim(tpar)[2])] = tpar
    tpar12_rep[(2500*(i-1)+1):(i*2500),1:(dim(Par12)[2])] = matrix(rep(Par12[i,],2500),nrow=2500,byrow = T)
    tpar23_rep[(2500*(i-1)+1):(i*2500),1:(dim(Par23)[2])] = matrix(rep(Par23[i,],2500),nrow=2500,byrow = T)
    
    print(paste0("The ",i,"th MCMC estimation is finished"))
  }
  
  print('Finish estimating tree 2')
  
  lst <- list(tpar12 = tpar12,
              tpar23 = tpar23,
              tpar12_rep = tpar12_rep,
              tpar23_rep = tpar23_rep,
              tpar132 = tpar132,
              bestpar12 = tpar12[UncProBi12$id_best_par],
              bestpar23 = tpar23[UncProBi23$id_best_par],
              bestpar132 = tpar132[which.max(tpar132[,3]),1],
              sampleID12 = sampleID12,
              sampleID23 = sampleID23,
              family12 = family12,
              family23 = family23,
              family132 = bestFamily132,
              crdf12 = UncProBi12$crdf,
              crdf23 = UncProBi23$crdf,
              crdf132 = ConProBi132$crdf )
  
  return (lst)
  
}

# calculate probabilistic cdf for all samples 
# Pro3Est is results of estimation, EstNum is number of sample in uncertainty
ProDVine3CDF <- function(x1, x2, x3, Pro3Est, EstNum) {
  
  # determine univariate pdf
  P1 <- fitalldist(x1)
  P2 <- fitalldist(x2)
  P3 <- fitalldist(x3)
  
  # calculate cdf value
  cdf1value <- P1$ep
  cdf2value <- P2$ep
  cdf3value <- P3$ep
  
  # extract cdf parameter
  para1 = P1$para
  para2 = P2$para
  para3 = P3$para
  
  # pdf, cdf, invcdf name
  pdf1name <- paste("d",P1$PD_name,sep = "")
  pdf2name <- paste("d",P2$PD_name,sep = "")
  pdf3name <- paste("d",P3$PD_name,sep = "")
  
  cdf1name <- paste("p",P1$PD_name,sep = "")
  cdf2name <- paste("p",P2$PD_name,sep = "")
  cdf3name <- paste("p",P3$PD_name,sep = "")
  
  invcdf1name <- paste("q",P1$PD_name,sep = "")
  invcdf2name <- paste("q",P2$PD_name,sep = "")
  invcdf3name <- paste("q",P3$PD_name,sep = "")
  
  ## construct joint PDF
  pdf1 <- function(w) {do.call(pdf1name,c(list(w),as.list(para1)))}
  pdf2 <- function(w) {do.call(pdf2name,c(list(w),as.list(para2)))}
  pdf3 <- function(w) {do.call(pdf3name,c(list(w),as.list(para3)))}
  
  cdf1 <- function(w) {do.call(cdf1name,c(list(w),as.list(para1)))}
  cdf2 <- function(w) {do.call(cdf2name,c(list(w),as.list(para2)))}
  cdf3 <- function(w) {do.call(cdf3name,c(list(w),as.list(para3)))}
  
  # extract estimation results: family, parameter set, sampleID
  family12 = Pro3Est$family12
  family23 = Pro3Est$family23
  family132 = Pro3Est$family132
  
  tpar12 = Pro3Est$tpar12[, 1]
  tpar23 = Pro3Est$tpar23[, 1]
  tpar132 = Pro3Est$tpar132[,1]
  
  # generate sub-parameter set with length of EstNum
  tpar12 = tpar12[round(seq(1,length(tpar12),length.out = EstNum))]
  tpar23 = tpar23[round(seq(1,length(tpar12),length.out = EstNum))]
  tpar132 = tpar132[round(seq(1,length(tpar12),length.out = EstNum))]
  
  # # meshgrid tpar12,tpar23,tpar132
  # n = length(tpar12)
  # tpar12 = rep( rep( tpar12, each = n),times = n )
  # tpar23 = rep( tpar23, times = n*n )
  # tpar132 = rep( tpar132, each = n*n )
  
  # subset x1 x2 x3
  x1 <- x1[(length(x1)-179):length(x1)]
  x2 <- x2[(length(x2)-179):length(x2)]
  x3 <- x3[(length(x3)-179):length(x3)]
  
  x123 <- cbind(rep(x1, each=length(tpar12)),
                rep(x2, each=length(tpar12)),
                rep(x3, each=length(tpar12)))
  
  tpar12 = rep(tpar12, times = length(x1) )
  tpar23 = rep(tpar23, times = length(x1) )
  tpar132 = rep(tpar132, times = length(x1) )
  
  # calculate cdf i is the row of tpar and x123
  calc_cdf <- function(i) {
    
    par12 = tpar12[i]
    par23 = tpar23[i]
    par132 = tpar132[i]
    upper <- as.vector(x123[i,])
    
    # construct decomposition pdf function
    coppdf12 <- function(a,b) {
      cop <- BiCop(family = family12,par = par12)
      return(VineCopula::BiCopPDF(cdf1(a),cdf2(b),cop))}
    coppdf23 <- function(a,b) {
      cop <- BiCop(family = family23,par = par23)
      return(VineCopula::BiCopPDF(cdf2(a),cdf3(b),cop))}
    
    h12 <- function(a,b) {VineCopula::BiCopHfunc2(cdf1(a), cdf2(b), family = family12, par = par12)}
    h23 <- function(a,b) {VineCopula::BiCopHfunc1(cdf2(a), cdf3(b), family = family23, par = par23)}
    coppdf132 <- function(a,b,c) {
      cop <- BiCop(family = family132,par = par132)
      return(VineCopula::BiCopPDF(h12(a,b),h23(b,c),cop))}
    
    # pair-copula decomposition of joint pdf
    JointPdf <- function(arg) { 
      x <- arg[1]
      y <- arg[2]
      z <- arg[3]
      ff <- pdf1(x) * pdf2(y) * pdf3(z) * coppdf12(x,y) * coppdf23(y,z) * coppdf132(x,y,z) 
      return (ff)
    }
    
    cdf <- hcubature(JointPdf, lowerLimit = rep(0, 3), upperLimit = upper, tol = 1e-3, maxEval = 10^4)
    
    cdf = cdf$integral
    
    return (cdf)
    
  }
  
  # parallel calculate cdf
  # Calculate the number of cores
  no_cores <- detectCores() - 1
  registerDoSEQ()
  cl <- makeCluster( no_cores )
  registerDoParallel(cl)
  
  t1 = Sys.time()
  allcdf = foreach(i=1:length(tpar12), .combine='c',.packages = c("VineCopula","cubature","actuar",
                                                                  "SCI","evd")) %dopar% calc_cdf(i)
  t2 = Sys.time()
  print(t2-t1)
  
  stopCluster(cl)
  
  allcdf = matrix(allcdf, round(length(allcdf)/length(x1)), length(x1))
  allcdf = t(allcdf)
  # transform to normal space
  allcdf = qnorm(allcdf)
  # add total uncertainty by adding rnorm(0,0.15)
  allcdf2 = allcdf + matrix( rnorm( nrow(allcdf)*ncol(allcdf), 0, 0.15), nrow(allcdf), ncol(allcdf) )
  
  par_unc = apply( allcdf, 1, FUN=quantile, probs=c(.025,.975) )
  tot_unc = apply( allcdf2, 1, FUN=quantile, probs=c(.025,.975) )
  
  par_unc = t(par_unc)
  tot_unc = t(tot_unc)
  
  allcdf = data.frame(id=1:length(x1),
                      par_unc1 = as.numeric(par_unc[,1] ),
                      par_unc2 = as.numeric(par_unc[,2] ),
                      tot_unc1 = as.numeric(tot_unc[,1] ),
                      tot_unc2 = as.numeric(tot_unc[,2] )
  )
  
  return (allcdf)
  
}

# deterministic bivariate cdf
BiCopDI <- function(x1,x2,familyset) {
  # fit univariate distribution
  P1 <- fitalldist(x1)
  P2 <- fitalldist(x2)
  
  u1 <- P1$ep
  u2 <- P2$ep
  
  cop = ProBiCopSelect(u1, u2, familyset, show = F)
  
  family <- cop$bf
  best_par <- cop$tpar[cop$id_best_par,1]
  
  cdf = BiCopCDF(u1, u2, par=best_par, family = family)
  
  TDI = qnorm(cdf)
  TDI = c(rep(NA,360-length(TDI)),TDI)
  return (TDI)
  
}

# estimate flood return periods for past and future climates
MCMC_RP <- function(MCMC_out, P1, P2, P3, PVD.hist, PVD.futu, mu, SampleNum=20) {
  
  # calculate cdf value
  cdf1value <- P1$ep
  cdf2value <- P2$ep
  cdf3value <- P3$ep
  
  # extract cdf parameter
  para1 = P1$para
  para2 = P2$para
  para3 = P3$para
  
  # pdf, cdf, invcdf name
  pdf1name <- paste("d",P1$PD_name,sep = "")
  pdf2name <- paste("d",P2$PD_name,sep = "")
  pdf3name <- paste("d",P3$PD_name,sep = "")
  
  cdf1name <- paste("p",P1$PD_name,sep = "")
  cdf2name <- paste("p",P2$PD_name,sep = "")
  cdf3name <- paste("p",P3$PD_name,sep = "")
  
  invcdf1name <- paste("q",P1$PD_name,sep = "")
  invcdf2name <- paste("q",P2$PD_name,sep = "")
  invcdf3name <- paste("q",P3$PD_name,sep = "")
  
  pdf1 <- function(w) {do.call(pdf1name,c(list(w),as.list(para1)))}
  pdf2 <- function(w) {do.call(pdf2name,c(list(w),as.list(para2)))}
  pdf3 <- function(w) {do.call(pdf3name,c(list(w),as.list(para3)))}
  
  cdf1 <- function(w) {do.call(cdf1name,c(list(w),as.list(para1)))}
  cdf2 <- function(w) {do.call(cdf2name,c(list(w),as.list(para2)))}
  cdf3 <- function(w) {do.call(cdf3name,c(list(w),as.list(para3)))}
  
  ## ----------------- univariate cdf without uncertainty --------------------
  cdf1value2 <- c( cdf1(PVD.hist[,1]), cdf1(PVD.futu[,1]) )
  cdf2value2 <- c( cdf2(PVD.hist[,2]), cdf2(PVD.futu[,2]) )
  cdf3value2 <- c( cdf3(PVD.hist[,3]), cdf3(PVD.futu[,3]) )
  data <- cbind(cdf1value2, cdf2value2, cdf3value2)
  ## ----------------- univariate cdf with uncertainty ----------------------
  cdf1_unc <- foreach(i=1:nrow(P1$parboot), .combine = 'rbind') %do% do.call(cdf1name,c(list(c(PVD.hist[,1], PVD.futu[,1] )),as.list(P1$parboot[i,])))
  cdf2_unc <- foreach(i=1:nrow(P2$parboot), .combine = 'rbind') %do% do.call(cdf2name,c(list(c(PVD.hist[,2], PVD.futu[,2] )),as.list(P2$parboot[i,])))
  cdf3_unc <- foreach(i=1:nrow(P3$parboot), .combine = 'rbind') %do% do.call(cdf3name,c(list(c(PVD.hist[,3], PVD.futu[,3] )),as.list(P3$parboot[i,])))
  # retain q025 and q975 cdf
  cdf1_unc <- as.numeric(t(apply(cdf1_unc, 2, function(x) quantile(x,probs = c(.025,.975)))))
  cdf2_unc <- as.numeric(t(apply(cdf2_unc, 2, function(x) quantile(x,probs = c(.025,.975)))))
  cdf3_unc <- as.numeric(t(apply(cdf3_unc, 2, function(x) quantile(x,probs = c(.025,.975)))))
  data <- cbind(cdf1_unc, cdf2_unc, cdf3_unc)
  data[data<=0] = 1e-5
  
  familyset = MCMC_out$familyset
  family12 = familyset[1]
  family23 = familyset[2]
  family132 = familyset[3]
  
  tpar = MCMC_out$tpar[,1:(ncol(MCMC_out$tpar)-2)]
  loglik = MCMC_out$tpar[,ncol(MCMC_out$tpar)]
  q <- unlist(sapply(seq(min(loglik), max(loglik), length.out = SampleNum), function(x) which.min(abs(x-loglik))))
  tpar0 = tpar[q,]
  tpar = cbind(
    rep(tpar0[,1] , nrow(data) ),
    rep(tpar0[,2] , nrow(data) ),
    rep(tpar0[,3] , nrow(data) ),
    rep(tpar0[,4] , nrow(data) )
  )
  
  raw1 <- rep( c(PVD.hist[,1], PVD.futu[,1]), each = SampleNum )
  raw2 <- rep( c(PVD.hist[,2], PVD.futu[,2]), each = SampleNum )
  raw3 <- rep( c(PVD.hist[,3], PVD.futu[,3]), each = SampleNum )
  x123 <- cbind(raw1,raw2,raw3)
  
  # calculate cdf
  calc_cdf <- function(i) {
    
    par12 = tpar[i,1]
    par23_1 = tpar[i,2]
    par23_2 = tpar[i,3]
    par132 = tpar[i,4]
    upper <- as.vector(x123[i,])
    
    # construct decomposition pdf function
    coppdf12 <- function(a,b) {
      cop <- BiCop(family = family12, par = par12)
      return(VineCopula::BiCopPDF(cdf1(a),cdf2(b),cop))}
    coppdf23 <- function(a,b) {
      cop <- BiCop(family = family23, par = par23_1, par2 = par23_2)
      return(VineCopula::BiCopPDF(cdf2(a),cdf3(b),cop))}
    
    h12 <- function(a,b) {VineCopula::BiCopHfunc2(cdf1(a), cdf2(b), family = family12, par = par12)}
    h23 <- function(a,b) {VineCopula::BiCopHfunc1(cdf2(a), cdf3(b), family = family23, par = par23_1, par2 = par23_2)}
    coppdf132 <- function(a,b,c) {
      cop <- BiCop(family = family132,par = par132)
      return(VineCopula::BiCopPDF(h12(a,b),h23(b,c),cop))}
    
    # pair-copula decomposition of joint pdf
    JointPdf <- function(arg) { 
      x <- arg[1]
      y <- arg[2]
      z <- arg[3]
      ff <- pdf1(x) * pdf2(y) * pdf3(z) * coppdf12(x,y) * coppdf23(y,z) * coppdf132(x,y,z) 
      return (ff)
    }
    
    cdf <- hcubature(JointPdf, lowerLimit = rep(0, 3), upperLimit = upper, tol = 1e-4, maxEval = 10^4)
    
    cdf = cdf$integral
    
    return (cdf)
    
  }
  
  # Calculate the number of cores
  no_cores <- detectCores() - 1
  cl <- makeCluster( no_cores )
  registerDoParallel(cl)
  # parallel estimate joint cdf
  t1 = Sys.time()
  allcdf = foreach(i=1:nrow(x123), .combine='c', .packages = c("VineCopula","cubature",
                                                               "actuar","SCI","evd")) %dopar% calc_cdf(i)
  t2 = Sys.time()
  print(t2-t1)
  stopCluster(cl)
  
  tricdf = t ( matrix( allcdf, SampleNum, nrow(data)) )
  
  # estimate tpar of var1 and var3
  local13 = VineCopula::BiCopSelect(cdf1value, cdf3value, familyset = c(0,1,2,3,4,5,6,7))
  mcmc13 = ProBiCopEst(cbind(cdf1value, cdf3value), family = local13$family, show = T)
  tpar13 = rep(mcmc13$tpar[q,1], nrow(data) )
  
  data1 = apply(data, 2, function(x) repmat(x, SampleNum,1))
  
  # calculate bicdf
  bicdf13 = sapply(1:nrow(data1), function(x) VineCopula::BiCopCDF(data1[x,1], data1[x,3], family = local13$family, par = tpar13[x]))
  bicdf12 = sapply(1:nrow(data1), function(x) VineCopula::BiCopCDF(data1[x,1], data1[x,2], family = family12, par = tpar[x,1]))
  bicdf23 = sapply(1:nrow(data1), function(x) VineCopula::BiCopCDF(data1[x,2], data1[x,3], family = family23, par = tpar[x,2], par2 = tpar[x,3]))
  
  bicdf13 = t ( matrix( bicdf13, SampleNum, nrow(data)) )
  bicdf12 = t ( matrix( bicdf12, SampleNum, nrow(data)) )
  bicdf23 = t ( matrix( bicdf23, SampleNum, nrow(data)) )
  
  cdf1 = t ( matrix( data1[,1], SampleNum, nrow(data)) )
  cdf2 = t ( matrix( data1[,2], SampleNum, nrow(data)) )
  cdf3 = t ( matrix( data1[,3], SampleNum, nrow(data)) )
  
  out = list(tricdf = tricdf,
             bicdf13 = bicdf13,
             bicdf12 = bicdf12,
             bicdf23 = bicdf23,
             cdf1 = cdf1,
             cdf2 = cdf2,
             cdf3 = cdf3)
  
  return(out)
}

# --------------- application ------------------
# load guadalupe flood dataset
df = read.csv('Guadalupe_70year_PVD.csv')
P = df$P; V = df$V; D = df$D
# randomization procedure
D = D + rand(length(D),1)
V = V + rand(length(V),1)
P = P + rand(length(P),1)

x1 = P; x2 = V; x3 = D
# Independence, Gaussian, Clayton, Gumbel, Frank, Joe, BB1, BB6, BB7, Tawn1 
familyset <- 0:7
# Bayesian 3-D vine copula simulation
SampleNum = 100
DVine = ProDVine3Est(x1, x2, x3, familyset, SampleNum)
# estimate flood return period for past and future climates
P1 = readRDS('guadalupe_cdf_fit_P.rds')
P2 = readRDS('guadalupe_cdf_fit_V.rds')
P3 = readRDS('guadalupe_cdf_fit_D.rds')

hist = read.csv('guadalupe_flood_event_1981_1995.csv')
futu = read.csv('guadalupe_flood_event_2085_2099.csv')
PVD.hist = hist[c('P','V','D')]
PVD.futu = futu[c('P','V','D')]
mu = mean(df[2:nrow(df),'start']-df[1:(nrow(df)-1),'start'])/365
RP = MCMC_RP( DVine, P1, P2, P3, PVD.hist, PVD.futu, mu)
RP.or = mu / (1 - RP$tricdf)
RP.and = mu / (1 - RP$cdf1 - RP$cdf2 - RP$cdf3 + RP$bicdf12 + RP$bicdf23 + RP$bicdf13 - RP$tricdf )
