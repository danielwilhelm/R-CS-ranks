test_that("summary does not raise errors", {
  mod <- lmranks(r(mpg) ~ r(cyl) + disp, data=mtcars)
  expect_warning(summary(mod), "degrees of freedom")
})

test_that("vcov passes shallow checks", {
  model <- lmranks(r(mpg) ~ r(cyl) + disp, data=mtcars)
  V <- vcov(model)
  
  expect_true(isSymmetric(V))
  vals <- eigen(V, only.values = TRUE)
  expect_true(all(vals$values> 0))
  expect_equal(colnames(V),
               c("(Intercept)", "r(cyl)", "disp"))
  expect_equal(rownames(V),
               c("(Intercept)", "r(cyl)", "disp"))
})

test_that("vcov passes shallow checks in no ranked regressors case", {
  model <- lmranks(r(mpg) ~ cyl + disp, data=mtcars)
  V <- vcov(model)
  
  expect_true(isSymmetric(V))
  vals <- eigen(V, only.values = TRUE)
  expect_true(all(vals$values> 0))
})

test_that("vcov works for singular model matrix", {
  data(mtcars)
  W <- cbind(mtcars$disp, mtcars$disp)
  mod <- lmranks(r(mpg) ~ r(cyl) + W, data=mtcars)
  cov_mat <- vcov(mod)
  
  expect_true(all(!is.na(cov_mat[1:3, 1:3])))
  expect_true(all(is.na(cov_mat[4,])))
  expect_true(all(is.na(cov_mat[,4])))
})

test_that("vcov works for singular model matrix, complete=FALSE", {
  data(mtcars)
  W <- cbind(mtcars$disp, mtcars$disp)
  mod <- lmranks(r(mpg) ~ r(cyl) + W, data=mtcars)
  cov_mat <- vcov(mod, complete=FALSE)
  
  expect_true(all(!is.na(cov_mat)))
  expect_equal(c(nrow(cov_mat), ncol(cov_mat)),
               c(3,3))
})

test_that("get_and_separate_regressors works",{
  data(mtcars)
  model_1 <- lm(mpg ~ disp + cyl + hp, data=mtcars)
  model_1$rank_terms_indices <- 1
  expected_W <- cbind(rep(1, nrow(mtcars)),
                      as.matrix(mtcars[,c("cyl", "hp")]))
  colnames(expected_W) <- c("(Intercept)", "cyl", "hp")
  expected_RX <- mtcars$disp
  names(expected_RX) <- rownames(mtcars)
  expected_out_1 <- list(RX = expected_RX,
                         W = expected_W,
                         rank_column_index = 2) # Intercept
  expect_equal(get_and_separate_regressors(model_1),
               expected_out_1)
  
  W <- as.matrix(mtcars[,c("cyl", "hp")])
  RX <- mtcars$disp
  Y <- mtcars$mpg
  model_2 <- lm(Y ~ W + RX)
  model_2$rank_terms_indices <- 2
  expected_out_2 <- expected_out_1
  names(expected_out_2$RX) <- 1:length(RX)
  rownames(expected_out_2$W) <- 1:length(RX)
  colnames(expected_out_2$W) <- c("(Intercept)", "Wcyl", "Whp")
  expected_out_2$rank_column_index <- 4
  expect_equal(get_and_separate_regressors(model_2),
               expected_out_2)
})
test_that("get_and_separate_regressors works for no ranked regressors",{
  data(mtcars)
  model <- lm(mpg ~ disp + cyl + hp, data=mtcars)
  model$rank_terms_indices <- numeric(0)
  
  expected_W <- cbind(rep(1, nrow(mtcars)),
                      as.matrix(mtcars[,c("disp", "cyl", "hp")]))
  colnames(expected_W) <- c("(Intercept)", "disp", "cyl", "hp")
  expected_out <- list(RX = integer(0),
                       W = expected_W,
                       rank_column_index = integer(0))
  
  expect_equal(get_and_separate_regressors(model),
               expected_out)
})

test_that("get_projection_model works for ranked target", {
  data(mtcars)
  original_model <- lmranks(r(mpg) ~ r(disp) + cyl + hp - 1, data=mtcars)
  
  RX <- original_model$model$`r(disp)`
  names(RX) <- rownames(mtcars)
  W <- as.matrix(mtcars[,c("cyl", "hp")])
  expected_model <- lm(RX ~ W - 1)
  expected_model$rank_terms_indices <- integer(0)
  expected_model$ranked_response <- TRUE
  
  expect_equal(get_projection_model(original_model, 1),
               expected_model)
})

test_that("get_projection_model works for usual target with ranked regressor present", {
  data(mtcars)
  original_model <- lmranks(r(mpg) ~ r(disp) + cyl + hp - 1, data=mtcars)
  
  RX <- original_model$model$`r(disp)`
  names(RX) <- rownames(mtcars)
  W_l <- mtcars$cyl
  names(W_l) <- rownames(mtcars)
  W_minus_l <- as.matrix(mtcars[,"hp"])
  colnames(W_minus_l) <- "hp"
  rownames(W_minus_l) <- rownames(mtcars)
  expected_model <- lm(W_l ~ RX + W_minus_l - 1)
  expected_model$rank_terms_indices <- 1
  expected_model$ranked_response <- FALSE
  
  expect_equal(get_projection_model(original_model, 2),
               expected_model)
})

test_that("get_projection_model works for usual target with ranked regressor present", {
  data(mtcars)
  original_model <- lmranks(r(mpg) ~ r(disp) + cyl + hp - 1, data=mtcars)
  
  RX <- original_model$model$`r(disp)`
  names(RX) <- rownames(mtcars)
  W_l <- mtcars$cyl
  names(W_l) <- rownames(mtcars)
  W_minus_l <- as.matrix(mtcars[,"hp"])
  colnames(W_minus_l) <- "hp"
  rownames(W_minus_l) <- rownames(mtcars)
  expected_model <- lm(W_l ~ RX + W_minus_l - 1)
  expected_model$rank_terms_indices <- 1
  expected_model$ranked_response <- FALSE
  
  expect_equal(get_projection_model(original_model, 2),
               expected_model)
})

test_that("get_projection_model works with no ranked regressors", {
  data(mtcars)
  original_model <- lmranks(r(mpg) ~ disp + cyl + hp - 1, data=mtcars)
  
  W_l <- mtcars$cyl
  names(W_l) <- rownames(mtcars)
  W_minus_l <- as.matrix(mtcars[,c("disp", "hp")])
  colnames(W_minus_l) <- c("disp", "hp")
  rownames(W_minus_l) <- rownames(mtcars)
  expected_model <- lm(W_l ~ W_minus_l - 1)
  expected_model$rank_terms_indices <- integer(0)
  expected_model$ranked_response <- FALSE
  
  expect_equal(get_projection_model(original_model, 2),
               expected_model)
})

test_that("get_projection_model works for intercept", {
  data(mtcars)
  original_model <- lmranks(r(mpg) ~ r(disp) + cyl + hp, data=mtcars)
  
  RX <- original_model$model$`r(disp)`
  names(RX) <- rownames(mtcars)
  W_l <- rep(1, nrow(mtcars))
  names(W_l) <- rownames(mtcars)
  W_minus_l <- as.matrix(mtcars[,c("cyl", "hp")])
  expected_model <- lm(W_l ~ RX + W_minus_l - 1)
  expected_model$rank_terms_indices <- 1
  expected_model$ranked_response <- FALSE
  
  expect_equal(get_projection_model(original_model, 1),
               expected_model)
})

test_that("get_projection_model works for model with 1 usual regressor", {
  data(mtcars)
  original_model <- lmranks(r(mpg) ~ r(disp), data=mtcars)
  
  RX <- original_model$model$`r(disp)`
  names(RX) <- rownames(mtcars)
  W_l <- rep(1, nrow(mtcars))
  names(W_l) <- rownames(mtcars)
  expected_model <- lm(W_l ~ RX - 1)
  expected_model$rank_terms_indices <- 1
  expected_model$ranked_response <- FALSE
  
  expect_equal(get_projection_model(original_model, 1),
               expected_model)
})

test_that("calculate_g_l_3 works for no ranked regressors case", {
  data(mtcars)
  I_X <- NULL
  W <- as.matrix(mtcars[, c("cyl", "hp")])
  W <- cbind(rep(1, nrow(W)),
             W)
  
  original_model <- lmranks(r(mpg) ~ cyl + hp, data=mtcars)
  RY <- original_model$model$`r(mpg)`
  W_l <- W[,2]
  W_minus_l <- W[,-2]
  proj_model <- lm(W_l ~ W_minus_l - 1)
  proj_model$rank_terms_indices <- integer(0)
  proj_model$ranked_response <- FALSE
  
  expected_out <- rep(0, nrow(mtcars))
  
  expect_equal(calculate_g_l_3(original_model, proj_model, I_X=I_X),
               expected_out)
})

test_that("calculate_g_l_3 works for proj_model with ranked response", {
  data(mtcars)
  original_model <- lmranks(r(mpg) ~ r(cyl) + hp, data=mtcars)
  RX <- original_model$model$`r(cyl)`
  I_X <- compare(mtcars$cyl)
  RY <- original_model$model$`r(mpg)`
  W <- as.matrix(mtcars[, c("hp")])
  W <- cbind(rep(1, nrow(W)),
             W)
  
  proj_model <- lm(RX ~ W - 1)
  proj_model$rank_terms_indices <- integer(0)
  proj_model$ranked_response <- TRUE
  
  original_resid <- resid(original_model)
  predictor <- proj_model$fitted.values
  
  expected_out <- sapply(1:length(RY), function(i){
    mean(
      sapply(1:length(RY), function(j){
        original_resid[j] * (I_X[i,j] - predictor[j])
      })
    )
  })
  
  expect_equal(calculate_g_l_3(original_model, proj_model, I_X=I_X),
               expected_out)
})

test_that("calculate_g_l_3 works for proj_model with usual response", {
  data(mtcars)
  original_model <- lmranks(r(mpg) ~ r(cyl) + hp, data=mtcars)
  RX <- original_model$model$`r(cyl)`
  I_X <- compare(mtcars$cyl)
  RY <- frank(mtcars$mpg)
  I_Y <- original_model$model$`r(disp)`
  W_l <- mtcars$hp
  W_minus_l <- rep(1, length(W_l))
  
  
  proj_model <- lm(W_l ~ RX + W_minus_l - 1)
  proj_model$rank_terms_indices <- 1
  proj_model$ranked_response <- FALSE
  
  original_resid <- resid(original_model)
  coefs <- coef(proj_model)
  proj_usual_predictor <- W_minus_l * coefs[2]
  
  expected_out <- sapply(1:length(RY), function(i){
    mean(
      sapply(1:length(RY), function(j){
        original_resid[j] * (W_l[j] - coefs[1]*I_X[i,j] - proj_usual_predictor[j])
      })
    )
  })
  
  expect_equal(calculate_g_l_3(original_model, proj_model, I_X=I_X),
               expected_out)
})

test_that("replace_ranks_with_ineq_indicator_and_calculate_residuals works 
          for ranked response, usual regressors", {
            data(mtcars)
            RY <- frank(mtcars$mpg)
            I_Y <- compare(mtcars$mpg)
            W <- as.matrix(mtcars[, c("disp", "hp")])
            W <- cbind(rep(1, nrow(W)),
                       W)
            original_model <- lm(RY ~ W - 1)
            original_model$rank_terms_indices <- integer(0)
            original_model$ranked_response <- TRUE
            
            coefs <- coef(original_model)
            predictor <- W %*% coefs
            expected_resid <- sapply(1:length(RY), function(j){
              I_Y[,j] - predictor[j]
            })
            expected_resid <- t(expected_resid)
            
            expect_equal(replace_ranks_with_ineq_indicator_and_calculate_residuals(
              original_model, I_Y = I_Y
            ), expected_resid)
            
          })

test_that("replace_ranks_with_ineq_indicator_and_calculate_residuals works 
          for usual response, ranked regressors", {
            data(mtcars)
            RX <- frank(mtcars$cyl)
            I_X <- compare(mtcars$cyl)
            Y <- frank(mtcars$mpg)
            W <- as.matrix(mtcars[, c("hp")])
            W <- cbind(rep(1, nrow(W)),
                       W)
            original_model <- lm(Y ~ RX + W - 1)
            original_model$rank_terms_indices <- 1
            original_model$ranked_response <- FALSE
            
            coefs <- coef(original_model)
            predictor <- W %*% coefs[-1]
            expected_resid <- sapply(1:length(Y), function(j){
              Y[j] - coefs[1] * I_X[,j] - predictor[j]
            })
            expected_resid <- t(expected_resid)
            
            expect_equal(replace_ranks_with_ineq_indicator_and_calculate_residuals(
              original_model, I_X = I_X
            ), expected_resid)
          })

test_that("replace_ranks_with_ineq_indicator_and_calculate_residuals works 
          for ranked response, ranked regressors", {
            data(mtcars)
            RX <- frank(mtcars$cyl)
            I_X <- compare(mtcars$cyl)
            RY <- frank(mtcars$mpg)
            I_Y <- compare(mtcars$mpg)
            W <- as.matrix(mtcars[, c("hp")])
            W <- cbind(rep(1, nrow(W)),
                       W)
            original_model <- lm(RY ~ RX + W - 1)
            original_model$rank_terms_indices <- 1
            original_model$ranked_response <- TRUE
            
            coefs <- coef(original_model)
            predictor <- W %*% coefs[-1]
            expected_resid <- sapply(1:length(RY), function(j){
              I_Y[,j] - coefs[1] * I_X[,j] - predictor[j]
            })
            expected_resid <- t(expected_resid)
            
            expect_equal(replace_ranks_with_ineq_indicator_and_calculate_residuals(
              original_model, I_X = I_X, I_Y = I_Y
            ), expected_resid)
          })

test_that("replace_ranks_with_ineq_indicator_and_calculate_residuals works 
          for singular fit", {
            data(mtcars)
            RX <- frank(mtcars$cyl)
            I_X <- compare(mtcars$cyl)
            RY <- frank(mtcars$mpg)
            I_Y <- compare(mtcars$mpg)
            W <- cbind(mtcars[, c("hp")],mtcars[, c("hp")])
            W <- cbind(rep(1, nrow(W)),
                       W)
            original_model <- lm(RY ~ RX + W - 1)
            original_model$rank_terms_indices <- 1
            original_model$ranked_response <- TRUE
            
            coefs <- coef(original_model)
            predictor <- W[,-3] %*% coefs[c(2:3)]
            expected_resid <- sapply(1:length(RY), function(j){
              I_Y[,j] - coefs[1] * I_X[,j] - predictor[j]
            })
            expected_resid <- t(expected_resid)
            
            expect_equal(replace_ranks_with_ineq_indicator_and_calculate_residuals(
              original_model, I_X = I_X, I_Y = I_Y
            ), expected_resid)
          })



test_that("vcov produces correct asymptotic variance estimate of rank-rank slope", {
  set.seed(100)
  
  for (covariates in c(TRUE,FALSE)) {
    
    # draw data
    n <- 10000
    if (covariates) {
      X <- rnorm(n)
      W <- matrix(rnorm(n*2), n, 2)
      Y <- X + rowSums(W) + rnorm(n,0,0.5)  
      W <- cbind(1,W)  
    } else {
      X <- rnorm(n)
      Y <- X + rnorm(n,0,0.5) 
      W <- matrix(1,n,1)
    }
    
    # compute ranks
    RY <- frank(Y, increasing=TRUE)
    RX <- frank(X, increasing=TRUE)
    
    
    # ------- compute asymptotic variance "by hand"
    
    Ifn <- function(u, v) return( u<=v )
    
    # first stage
    res <- lm(RX ~ W-1)
    Wgammahat <- predict(res)
    nuhat <- resid(res)
    gammahat <- coef(res)
    
    # outcome equation
    res <- lm(RY~RX+W-1)
    rhohat <- coef(res)[1]
    betahat <- coef(res)[-1]
    epsilonhat <- resid(res)
    
    # construct h1
    h1 <- epsilonhat * nuhat
    
    # construct h2
    h2fn <- function(xy) mean((Ifn(xy[2],Y)-rhohat*Ifn(xy[1],X)-c(W%*%betahat)) * nuhat)
    h2 <- apply(cbind(X,Y), 1, h2fn)
    
    # construct h3
    h3fn <- function(x) mean(epsilonhat * (Ifn(x,X)-Wgammahat))
    h3 <- sapply(X, h3fn)
    
    # compute asymptotic variance
    sigma2hat <- mean((h1+h2+h3)^2) / var(nuhat)^2
    
    
    # ------- compute asymptotic variance using lmranks
    
    if (covariates) { 
      res <- lmranks(r(Y) ~ r(X) + W - 1)
      sigma2hat.lmranks <- vcov(res)[1,1]*n
    } else {
      res <- lmranks(r(Y) ~ r(X))
      sigma2hat.lmranks <- vcov(res)[2,2]*n
    }
    
    
    # ------- test equality
    
    expect_equal(sigma2hat, sigma2hat.lmranks, tolerance=1e-5)
  }
})



test_that("vcov produces asymptotic variance estimate of rank-rank slope equal to that of Hoeffding (1948)", {
  set.seed(100)
  
  # draw data
  n <- 1000
  X <- rnorm(n)
  Y <- X + rnorm(n,0,0.5) 
  
  # ------- compute Hoeffding's variance "by hand"
  
  n <- length(Y)
  rhohat <- cor(frank(Y, omega=1/2, increasing=TRUE),frank(X, omega=1/2, increasing=TRUE))
  oX <- outer(X,X,'<=')
  oY <- outer(Y,Y,'<=')
  FhatX <- colMeans(oX)
  FhatY <- colMeans(oY)
  FhatAvgX <- function(x) mean(colMeans(outer(X,rep(x,n),'<=')*oY))
  FhatAvgY <- function(y) mean(colMeans(outer(Y,rep(y,n),'<=')*oX))
  psihatX <- sapply(X, FhatAvgX) - FhatX*mean(FhatY)
  psihatY <- sapply(Y, FhatAvgY) - mean(FhatX)*FhatY
  
  W <- 3*( (2*FhatX-1)*(2*FhatY-1) + 4*psihatX + 4*psihatY )
  
  # compute asymptotic variance
  sigma2hat <- var(W)
  
  
  # ------- compute asymptotic variance using lmranks
  
  res <- lmranks(r(Y) ~ r(X))
  sigma2hat.lmranks <- vcov(res)[2,2]*n
  
  
  # ------- test equality
  
  expect_equal(sigma2hat, sigma2hat.lmranks, tolerance=1e-3)
})



