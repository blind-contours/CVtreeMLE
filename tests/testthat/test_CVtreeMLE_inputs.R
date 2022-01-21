library(CVtreeMLE)
library(testthat)
library(sl3)

# Set multicore-compatible seed.
set.seed(1, "L'Ecuyer-CMRG")

# Simulation code from PH 295, Fall 2016 --
# https://github.com/wilsoncai1992/PH295-lab/blob/master/lab3/Lab3_Lecture.Rmd

# structural equation for W_1
# takes as input a vector U_W1 and returns a vector evaluating
# f_{W,1}(U_W1)
f_W1 <- function(U_W1){
  return(U_W1)
}

# structural equation for W_2
# takes as input a vector U_W2 and returns a vector evaluating
# f_{W,2}(U_W2)
f_W2 <- function(U_W2){
  return(U_W2)
}

# structural equation for A
f_A <- function(W_1, W_2, U_A){
  return(as.numeric(plogis(W_1 - W_2 + U_A) > 0.5))
}

# structural equation for Y
f_Y <- function(W_1, W_2, A, U_Y){
  return(-W_1 + W_2 + A - U_Y)
}

# function to draw n observations from an scm
# n = the number of observations to draw
# returns a data.frame with named columns
simObsSCM <- function(n){
  ## first we draw the errors
  # draw U_{W,1}
  U_W1 <- rbinom(n,1,0.5)
  # draw U_{W,2}
  U_W2 <- rbinom(n,1,0.5)
  # draw U_A
  U_A <- rnorm(n,0,1)
  # draw U_Y
  U_Y <- rnorm(n,0,1)

  ## now we can evaluate the observations sequentially
  # evaluate W_1
  W_1 <- f_W1(U_W1)
  #evaluate W_2
  W_2 <- f_W2(U_W2)
  # evaluate A
  A <- f_A(W_1 = W_1, W_2 = W_2, U_A = U_A)
  # evaluate Y
  Y <- f_Y(W_1 = W_1, W_2 = W_2, A = A, U_Y = U_Y)

  ## return a data.frame object
  out <- data.frame(W_1 = W_1, W_2 = W_2, A = A, Y = Y)
  return(out)
}

data = simObsSCM(100)
summary(data)

#################################
# set up learners

Q1_stack <- "Just filler here for test"


# we will arbitrarily investigate cell qualities as a mixture and the other variables as covariates

W <- c("W_1", "W_2")
A <- c("A")

expect_error(CVtreeMLE_results <- CVtreeMLE(data = data,
                                  W = W,
                                  Y = "Y",
                                  A = A,
                                  back_iter_SL = Q1_stack,
                                  n_folds = 2,
                                  family = "gaussian",
                                  H.AW_trunc_lvl = 10,
                                  max_iter = 10,
                                  minsize = 20,
                                  parallel = TRUE,
                                  verbose = FALSE))

W <- c("W_1", "W_2")
A <- c("A_1", "A_2")

NA_indices <- sample(1:length(data[,"Y"]), 10)

data[NA_indices, "Y"] <- NA

colnames(data)[3] <- "A_1"
data$A_2 <- rbinom(dim(data)[1], 1, prob = 0.4)

expect_error(CVtreeMLE_results <- CVtreeMLE(data = data,
                                            W = W,
                                            Y = "Y",
                                            A = A,
                                            back_iter_SL = Q1_stack,
                                            n_folds = 2,
                                            family = "gaussian",
                                            H.AW_trunc_lvl = 10,
                                            max_iter = 10,
                                            minsize = 20,
                                            parallel = TRUE,
                                            verbose = FALSE))


