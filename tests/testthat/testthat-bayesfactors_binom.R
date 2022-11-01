context("evaluate Bayes factors for informed hypotheses - Binomial")

factor_levels <- c('binom1', 'binom2', 'binom3', 'binom4')
a             <- c(1, 1, 1, 1)
b             <- c(1, 1, 1, 1)

test_that("yields equal BF estimates of 1", {
  # when no data is provided
  x <- c(0, 0, 0, 0, 0)
  n <- c(0, 0, 0, 0, 0)
  a <- c(1, 1, 1, 1, 1)
  b <- c(1, 1, 1, 1, 1)
  factor_levels <- c('theta1', 'theta2', 'theta3', 'theta4', 'theta5')
  output_total  <- binom_bf_informed(x=x, n=n, 
                                     a=a, b=b, 
                                     Hr= c('theta1', '<',  'theta2', '<', 'theta3', '&', 'theta4', ',', 'theta5'),
                                     factor_levels=factor_levels, 
                                     seed=2020, 
                                     bf_type = 'BF0r')
  expect_equal(output_total$bf_list$bf, structure(list(LogBFr0 = -0.00274743849807657, 
                                                       BF0r = 1.00275121616607, BFr0 = 0.997256332256977), class = "data.frame", row.names = c(NA, 
                                                                                                                                               -1L)))
  # when encomassing and informed hypothesis resemble each other
  x <- c(5, 10, 15, 20, 25)
  n <- c(30, 30, 30, 30, 30)
  output_total  <- binom_bf_informed(x=x, n=n, 
                                     a=a, b=b, 
                                     Hr= c('theta1', ',',  'theta2', ',', 'theta3', ',', 'theta4', ',', 'theta5'),
                                     factor_levels=factor_levels, 
                                     seed=2020, 
                                     bf_type = 'BFer') 
  expect_equal(output_total$bf_list$bf, structure(list(LogBFer = 0, 
                                                       BFer = 1, BFre = 1), class = "data.frame", row.names = c(NA, 
                                                                                                                -1L)))
  
  # when null and informed hypothesis resemble each other
  output_total  <- binom_bf_informed(x=x, n=n, 
                                     a=a, b=b, 
                                     Hr= c('theta1 = theta2 = theta3 = theta4 = theta5'),
                                     factor_levels=factor_levels, 
                                     seed=2020, 
                                     bf_type = 'BF0r')
  expect_equal(output_total$bf_list$bf, structure(list(LogBFr0 = 0, 
                                                       BF0r = 1, BFr0 = 1), class = "data.frame", row.names = c(NA, 
                                                                                                                -1L)))
})

test_that("yields equal BF estimates for Marsmans example", {
  
  # Maarten Marsmans example
  x <- c(5, 10, 15, 14)
  n <- c(17, 16, 16, 16)
  Hr            <- c('binom1', '<',  'binom2', '<', 'binom3', '<', 'binom4')
  output_total  <- binom_bf_informed(factor_levels=factor_levels, Hr=Hr, a=a, 
                                      b=b, x=x, n=n, niter = 5e3, bf_type = 'LogBFer', seed=2020)
  
  expect_equal(output_total$bf_list$bf, structure(list(LogBFer = -1.88123614919558, 
                                                       BFer = 0.152401598156163, BFre = 6.56161098110872), class = "data.frame", row.names = c(NA, 
                                                                                                                                               -1L)))
})

test_that("yields equal marginal likelihoods estimates for Gronau example", {
  
  # Gronau example
  x <- c(2, 2)
  n <- c(10, 10)
  a <- c(1, 1)
  b <- c(1, 1)
  factor_levels <- c('theta1', 'theta2')
  Hr <- c('theta1 < theta2')
  output_total  <- binom_bf_informed(x=x, n=n, Hr=Hr, a=a, b=b, factor_levels, seed=2020, niter=2e3, bf_type = 'BFre')
  summary(output_total)$logmlHe # -4.795791
  
  expect_equal(summary(output_total)$logmlHe, -4.79579054559673)
})

test_that("yields equal marginal likelihoods estimates for Sarafolgou example", {
  
  # Sarafoglou example
  x <- c(2, 2)
  n <- c(10, 10)
  a <- c(1, 1)
  b <- c(1, 1)
  K <- length(a)
  factor_levels <- c('theta1', 'theta2')
  Hr <- c('theta1 < theta2')
  output_total  <- binom_bf_informed(x=x, n=n, Hr=Hr, a=a, b=b, factor_levels, seed=2020, niter=2e3, bf_type = 'BFre')
  
  expect_equal(summary(output_total)$logmlH0, -3.91689998250764)
  
  lcoef    <- sum(lchoose(n, x))
  lbeta.xa <- lbeta(sum(x) + sum(a) - (K - 1), sum(n) - sum(x) + sum(b) - (K - 1))
  lbeta.a  <- lbeta(sum(a) - (K - 1), sum(b) - (K - 1))
  analytic_computation <- lcoef + lbeta.xa - lbeta.a # -3.9169
  expect_equal(summary(output_total)$logmlH0, analytic_computation)
  
})

test_that("yields equal BF estimates for Nuijten et al. 2016 example", {
  
  data(journals)
  a <- rep(1, 8)  
  b <- rep(1, 8)  
  counts <- round(journals$articles_with_NHST  * (journals$perc_articles_with_errors/100))
  total  <- journals$articles_with_NHST 
  factor_levels <- levels(journals$journal)
  
  Hr1 <- c('JAP , PS , JCCP , PLOS , DP , FP , JEPG < JPSP')
  Hr2 <- c('JCCP < DP < JPSP')
  Hr3 <- c('JAP < PS < JCCP < PLOS < DP < FP < JEPG < JPSP')
  output_total  <- binom_bf_informed(factor_levels=factor_levels, Hr=Hr1, a=a, 
                                      b=b, x=counts, n=total, 
                                      niter = 5e3, bf_type = 'LogBFer', seed=2020)
  expect_equal(output_total$bf_list$bf, structure(list(LogBFer = -2.00537413225981, 
                                                       BFer = 0.13460992435813, BFre = 7.42887275784734), 
                                                  class = "data.frame", row.names = c(NA, -1L)))
  output_total  <- binom_bf_informed(factor_levels=factor_levels, Hr=Hr2, a=a, 
                                      b=b, x=counts, n=total, 
                                      niter = 5e3, bf_type = 'LogBFer', seed=2020)
  expect_equal(output_total$bf_list$bf, structure(list(LogBFer = -1.05965988707527, 
                                                       BFer = 0.346573664469722, BFre = 2.88538946411309), 
                                                  class = "data.frame", row.names = c(NA, -1L)))
  output_total  <- binom_bf_informed(factor_levels=factor_levels, Hr=Hr3, a=a, 
                                      b=b, x=counts, n=total, 
                                      niter = 5e3, bf_type = 'LogBFer', seed=2020)
  expect_equal(output_total$bf_list$bf, structure(list(LogBFer = -8.94683860069302, 
                                                       BFer = 0.000130147960022851, 
                                                       BFre = 7683.56261461512), 
                                                  class = "data.frame", row.names = c(NA, -1L)))
})

test_that("BF estimate matches convergence expectations", {
  
  # yields equal BF estimates of 1 when no data is provided
  counts        <- c(0, 0, 0, 0)
  total         <- c(0, 0, 0, 0)
  Hr            <- c('binom1', '<',  'binom2', '<', 'binom3', '<', 'binom4')
  output_total  <- binom_bf_informed(factor_levels=factor_levels, Hr=Hr, a=a, 
                                      b=b, x=counts, n=total, 
                                      niter = 5e3, bf_type = 'LogBFer', seed=2020)
  
  expect_equal(output_total$bf_list$bf, structure(list(LogBFer = 0.00487829502509207, 
                                                       BFer = 1.00489021327864, BFre = 0.995133584530907), class = "data.frame", row.names = c(NA, 
                                                                                                                                               -1L)))
  # Very Large Values: BFre should converge to factorial(4)
  counts        <- c(3, 100, 800, 900)
  total         <- c(1e3, 1e3, 1e3, 1e3)
  Hr            <- c('binom1', '<',  'binom2', '<', 'binom3', '<', 'binom4')
  output_total  <- binom_bf_informed(factor_levels=factor_levels, Hr=Hr, a=a, 
                                      b=b, x=counts, n=total, 
                                      niter = 5e3, bf_type = 'LogBFer', seed=2020)
  
  expect_equal(output_total$bf_list$bf, structure(list(LogBFer = -3.23721127565167, 
                                                       BFer = 0.0392732648369703, BFre = 25.4626144312463), class = "data.frame", row.names = c(NA, 
                                                                                                                                                -1L)))
  # Very Small Bayes factor
  Hr            <- c('binom1', '>',  'binom2', '>', 'binom3', '>', 'binom4')
  output_total  <- binom_bf_informed(factor_levels=factor_levels, Hr=Hr, a=a, 
                                      b=b, x=counts, n=total, 
                                      niter = 5e3, bf_type = 'LogBFer', seed=2020)
  expect_equal(output_total$bf_list$bf, structure(list(LogBFer = 1641.78003934897, 
                                                       BFer = Inf, BFre = 0), class = "data.frame", row.names = c(NA, 
                                                                                                                  -1L)))
})

test_that("yields equal BF estimates for costraints with free parameters", {
  
  counts        <- c(3, 4, 10, 11)
  total         <- c(15, 12, 12, 12)
  Hr            <- c('binom1', '<',  'binom2', ',', 'binom3', ',', 'binom4')
  output_total  <- binom_bf_informed(factor_levels=factor_levels, Hr=Hr, a=a, 
                                      b=b, x=counts, n=total, 
                                      niter = 5e3, bf_type = 'LogBFer', seed=2020)
  
  expect_equal(output_total$bf_list$bf, structure(list(LogBFer = -1.12851829614506, 
                                                       BFer = 0.323512250819858, BFre = 3.09107305045098), class = "data.frame", row.names = c(NA, 
                                                                                                                                               -1L)))
})

test_that("BF works properly when user provides a matrix as input", {
  
  factor_levels <- c('vincent', 'maarten', 'tim', 'alex', 'felix', 'georg')
  # data
  x <- cbind(c(3, 21, 5, 20, 26, 1), c(7, 29, 5, 9, 4, 4))
  row.names(x) <- factor_levels
  # priors
  a <- c(1, 1, 1, 1, 1, 1)
  b <- c(1, 1, 1, 1, 1, 1)
  # informed hypothesis
  Hr            <- c('vincent < maarten < tim < alex < felix < georg')
  output_total  <- binom_bf_informed(x, n=NULL, Hr, a, b, niter=2e3, 
                                     factor_levels, seed=2020, bf_type = 'BFer')
  expect_equal(summary(output_total), structure(list(hyp = "vincent < maarten < tim < alex < felix < georg", 
                                                     bf = 3.7619, logmlHe = -17.3545602336964, logmlH0 = -23.8018377142112, 
                                                     logmlHr = -18.6794868357989, re2 = 0.00046584336900814, bf_type = "BFer", 
                                                     cred_level = 0.95, prior = list(a = c(1, 1, 1, 1, 1, 1), 
                                                                                     b = c(1, 1, 1, 1, 1, 1)), data = list(x = c(vincent = 3, 
                                                                                                                                 maarten = 21, tim = 5, alex = 20, felix = 26, georg = 1), 
                                                                                                                           n = c(vincent = 10, maarten = 50, tim = 10, alex = 29, 
                                                                                                                                 felix = 30, georg = 5)), nr_equal = 0L, nr_inequal = 1L, 
                                                     estimates = structure(list(factor_level = c("vincent", "maarten", 
                                                                                                 "tim", "alex", "felix", "georg"), alpha = c(4, 4, 4, 4, 4, 
                                                                                                                                             4), beta = c(8, 8, 8, 8, 8, 8), lower = c(0.109263443819098, 
                                                                                                                                                                                       0.293454873840595, 0.233793597659345, 0.506040958537446, 
                                                                                                                                                                                       0.701664170992203, 0.0432718682927417), median = c(0.323804462585187, 
                                                                                                                                                                                                                                          0.422083303632833, 0.5, 0.681278841818006, 0.850945073531973, 
                                                                                                                                                                                                                                          0.26444998329566), upper = c(0.609742559572421, 0.558307172232539, 
                                                                                                                                                                                                                                                                       0.766206402340655, 0.827125778473961, 0.945475673784916, 
                                                                                                                                                                                                                                                                       0.641234578997675)), row.names = c(NA, -6L), class = "data.frame")), class = c("summary.bmult", 
                                                                                                                                                                                                                                                                                                                                                      "list")))
  })

test_that("marginal likelihoods are factored correctly", {
  # Hr: t1 = t2 = t3 & t4 , t5
  x <- c(5, 10, 15, 20, 25)
  n <- c(30, 30, 30, 30, 30)
  a <- c(1, 1, 1, 1, 1)
  b <- c(1, 1, 1, 1, 1)
  factor_levels <- paste0('theta', 1:5)
  output_total1  <- binom_bf_informed(x=x, n=n, 
                                      a=a, b=b,
                                      Hr= c('theta1', '=',  'theta2', '=', 'theta3', '&', 'theta4', ',', 'theta5'),
                                      factor_levels=factor_levels, 
                                      seed=2020, 
                                      bf_type = 'BF0r') 
  logmlHr_combined <- summary(output_total1)$logmlHr
  
  # H0: t1 = t2 = t3
  output_total1  <- binom_bf_informed(x=x[1:3], n=n[1:3], 
                                      a=a[1:3], b=b[1:3],
                                      Hr= c('theta1', '<', 'theta2', '<', 'theta3'),
                                      factor_levels=c('theta1', 'theta2', 'theta3'), 
                                      seed=2020, 
                                      bf_type = 'BF0r') 
  logmlH0 <- summary(output_total1)$logmlH0 
  
  # He: t4 , t5
  output_total1  <- binom_bf_informed(x=x[4:5], n=n[4:5], 
                                      a=a[4:5], b=b[4:5],
                                      Hr= c('theta4', '<', 'theta5'),
                                      factor_levels=c('theta4', 'theta5'), 
                                      seed=2020, 
                                      bf_type = 'BF0r') 
  logmlHe <- summary(output_total1)$logmlHe
  
  expect_equal(logmlHr_combined, -18.3000234839776)
  expect_equal(logmlHr_combined, (logmlH0 + logmlHe))
  
})
