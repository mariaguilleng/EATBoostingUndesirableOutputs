library(boostingDEA)
library(dplyr)
library(lpSolveAPI)
source("eff_eatboosting.R")

set.seed(4321)

data <- read.csv2("apcb.csv")
columns <- c('raw', 'machine_error'	, 'manual_error', 'other_error', 'error_free')
data <- data %>%
  subset(select = c('raw', 'machine_error'	, 'manual_error', 'other_error', 'error_free'))

# Train model using both original inputs and undesitable outputs as inputs vars.
# and desirable outputs as output vars.
EATBoost_model <- EATBoost(data, x = 1:4, y = 5, 
                           num.iterations = 12, learning.rate = 0.25, num.leaves = 12)
pred <- predict(EATBoost_model, data, x = 1:4)

# Calculate the efficiency scores for each DMU
scores <- eff_eatboost_undesriable_out(data, x = 1, y_desirable = 5, y_undesirable = 2:4, model = EATBoost_model)
