## ----load_data----------------------------------------------------------------
# load the package and data
library(multibridge)
data(journals)
journals

## ----specify_hr---------------------------------------------------------------
# since percentages are rounded to two decimal values, we round the articles
# with an error to receive integer values (step 1)
x <- round(journals$articles_with_NHST  * (journals$perc_articles_with_errors/100))
# total number of articles (step 2)
n <- journals$articles_with_NHST 

# Specifying the informed Hypothesis (step 3)
Hr <- c('JAP , PS , JCCP , PLOS , DP , FP , JEPG < JPSP')

# Prior specification (step 4 and 5)
# We assign a uniform beta distribution to each binomial propotion
a <- rep(1, 8)
b <- rep(1, 8)

# categories of interest (step 6)
journal_names <- journals$journal

## ----compute_results----------------------------------------------------------
ineq_results <- multibridge::binom_bf_informed(x=x, n=n, Hr=Hr, a=a, b=b, 
                                             factor_levels=journal_names, 
                                             bf_type = 'BFre',
                                             seed = 2020)
m1 <- summary(ineq_results)
m1

## ----percentage_error---------------------------------------------------------
ineq_results$bf_list$error_measures

## ----compute_results2---------------------------------------------------------
ineq_results <- multibridge::binom_bf_informed(x=x, n=n, Hr=Hr, a=a, b=b, 
                                             factor_levels=journal_names, 
                                             bf_type = 'BFre',
                                             seed = 2020,
                                             niter = 2e4)
m2 <- summary(ineq_results)

## ----percentage_error2--------------------------------------------------------
ineq_results$bf_list$error_measures

## ----show_bfe0----------------------------------------------------------------
eq_results      <- multibridge::binom_bf_informed(x=x, n=n, Hr=Hr, a=a, b=b, 
                                             factor_levels=journal_names, 
                                             bf_type = 'BFr0',
                                             seed = 2020,
                                             niter = 2e4)
m3 <- summary(eq_results)
m3

