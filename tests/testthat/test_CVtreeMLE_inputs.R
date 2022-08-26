library(CVtreeMLE)
library(testthat)
library(sl3)


# simulation code from ph 295, fall 2016 --
# https://github.com/wilsoncai1992/ph295-lab/blob/master/lab3/lab3_lecture.rmd

# structural equation for w_1
# takes as input a vector u_w1 and returns a vector evaluating
# f_{w,1}(u_w1)
f_w1 <- function(u_w1) {
  return(u_w1)
}

# structural equation for w_2
# takes as input a vector u_w2 and returns a vector evaluating
# f_{w,2}(u_w2)
f_w2 <- function(u_w2) {
  return(u_w2)
}

# structural equation for a
f_a <- function(w_1, w_2, u_a) {
  return(as.numeric(plogis(w_1 - w_2 + u_a) > 0.5))
}

# structural equation for y
f_y <- function(w_1, w_2, a, u_y) {
  return(-w_1 + w_2 + a - u_y)
}

# function to draw n observations from an scm
# n = the number of observations to draw
# returns a data.frame with named columns
simobsscm <- function(n) {
  ## first we draw the errors
  # draw u_{w,1}
  u_w1 <- rbinom(n, 1, 0.5)
  # draw u_{w,2}
  u_w2 <- rbinom(n, 1, 0.5)
  # draw u_a
  u_a <- rnorm(n, 0, 1)
  # draw u_y
  u_y <- rnorm(n, 0, 1)

  ## now we can evaluate the observations sequentially
  # evaluate w_1
  w_1 <- f_w1(u_w1)
  # evaluate w_2
  w_2 <- f_w2(u_w2)
  # evaluate a
  a <- f_a(w_1 = w_1, w_2 = w_2, u_a = u_a)
  # evaluate y
  y <- f_y(w_1 = w_1, w_2 = w_2, a = a, u_y = u_y)

  ## return a data.frame object
  out <- data.frame(w_1 = w_1, w_2 = w_2, a = a, y = y)
  return(out)
}

data <- simobsscm(100)
summary(data)

#################################

na_indices <- sample(seq_along(data[, "y"]), 10)

data[na_indices, "y"] <- na

colnames(data)[3] <- "a_1"
data$a_2 <- rbinom(dim(data)[1], 1, prob = 0.4)

expect_error(cvtreemle_results <- CVtreeMLE(
  data = data,
  w = w,
  y = "y",
  a = a,
  n_folds = 2,
  family = "gaussian",
  h.aw_trunc_lvl = 10,
  max_iter = 10,
  parallel = TRUE,
  verbose = FALSE
))
