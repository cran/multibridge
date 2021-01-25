## ----load_data----------------------------------------------------------------
library(multibridge)
data(peas)
peas

## ----specify_hr---------------------------------------------------------------
x          <- peas$counts
# Test the following restricted Hypothesis:
# Hr: roundYellow > wrinkledYellow = roundGreen > wrinkledGreen
Hr   <- c('roundYellow > wrinkledYellow = roundGreen > wrinkledGreen')

# Prior specification 
# We assign a uniform Dirichlet distribution, that is, we set all concentration parameters to 1
a <- c(1, 1, 1, 1)
categories <- peas$peas

## ----compute_results----------------------------------------------------------
results <- multibridge::mult_bf_informed(x=x,Hr=Hr, a=a, factor_levels=categories, 
                                       bf_type = 'BFre', seed = 2020)

## ----show_summary-------------------------------------------------------------
m1 <- summary(results)
m1

## ----show_bf------------------------------------------------------------------
bayes_list <- bayes_factor(results)
bayes_list$bf_table
# Bayes factors in favor for informed hypothesis
bfre <- bayes_list$bf_table[bayes_list$bf_table$bf_type=='BFre', ]

