context("evaluate Bayes factors for informed hypotheses - Multinomial")

test_that("yields equal BF estimates of 1", {
  # when no data is provided
  a             <- c(1, 1, 1, 1, 1)
  counts        <- c(0, 0, 0, 0, 0)
  factor_levels <- paste0('g', 1:5)
  Hr            <- c('g1 < g2 = g3 = g4 < g5')
  output_total  <-
    mult_bf_informed(
      factor_levels = factor_levels,
      Hr = Hr,
      a = a,
      x = counts,
      niter = 5e3,
      bf_type = 'BFre',
      seed = 4
    )
  
  expect_equal(output_total$bf_list$bf,
               structure(
                 list(
                   LogBFer = 0,
                   BFer = 1,
                   BFre = 1
                 ),
                 class = "data.frame",
                 row.names = c(NA,-1L)
               ))
  
  # when encomassing and informed hypothesis resemble each other
  x <- c(5, 10, 15, 20, 25)
  output_total  <- mult_bf_informed(
    x = x,
    a = a,
    Hr = c('g1', ',',  'g2', ',', 'g3', ',', 'g4', ',', 'g5'),
    factor_levels = factor_levels,
    seed = 2020,
    bf_type = 'BFer'
  )
  expect_equal(output_total$bf_list$bf,
               structure(
                 list(
                   LogBFer = 0,
                   BFer = 1,
                   BFre = 1
                 ),
                 class = "data.frame",
                 row.names = c(NA,-1L)
               ))
  
  # when null and informed hypothesis resemble each other
  output_total  <- mult_bf_informed(
    x = x,
    a = a,
    Hr = c('g1', '=',  'g2', '=', 'g3', '=', 'g4', '=', 'g5'),
    factor_levels = factor_levels,
    seed = 2020,
    bf_type = 'BF0r'
  )
  expect_equal(output_total$bf_list$bf,
               structure(
                 list(
                   LogBFr0 = 0,
                   BF0r = 1,
                   BFr0 = 1
                 ),
                 class = "data.frame",
                 row.names = c(NA,-1L)
               ))
  
})

test_that("yields equal BF estimates for Haberman example", {
  data(lifestresses)
  a             <- rep(1, 18)
  counts        <- lifestresses$stress.freq
  factor_levels <- lifestresses$month
  Hr            <- paste0(1:18, collapse = ">")
  output_total  <-
    mult_bf_informed(
      factor_levels = factor_levels,
      Hr = Hr,
      a = a,
      x = counts,
      niter = 5e3,
      bf_type = 'BFre',
      seed = 4
    )
  
  expect_equal(
    output_total$bf_list$bf,
    data.frame(
      LogBFer = -5.14220,
      BFer = 0.00584,
      BFre = 171.09105
    ),
    tolerance = 0.002
  )
  
})

test_that("yields equal BF estimates for Mendelian Peas example", {
  data(peas)
  a             <- c(1, 1, 1, 1)
  counts        <- peas$counts
  factor_levels <- peas$peas
  
  ## specify hypothesis explicitely
  Hr1                   <-
    c('roundYellow > wrinkledYellow = roundGreen > wrinkledGreen')
  output_total_explicit <-
    mult_bf_informed(
      factor_levels = factor_levels,
      Hr = Hr1,
      a = a,
      x = counts,
      niter = 5e3,
      bf_type = 'BFre',
      seed = 4
    )
  # total BF
  expect_equal(
    output_total_explicit$bf_list$bf,
    data.frame(
      LogBFer = -4.071621,
      BFer = 0.01704974,
      BFre = 58.65193
    ),
    tolerance = 0.002
  )
  # inequality BF
  expect_equal(
    output_total_explicit$bf_list$logBFe_inequalities,
    data.frame(
      logBFe_inequalities = -1.73936,
      logml_prior = -1.79188,
      logml_post = -0.05252
    ),
    tolerance = 0.002
  )
  # equality BF
  expect_equal(
    output_total_explicit$bf_list$logBFe_equalities,
    data.frame(logBFe_equalities = -2.33227),
    tolerance = 0.002
  )
  
  ## specify hypothesis using indeces
  Hr2                <- c('1 > 2 = 3 > 4')
  output_total_index <-
    mult_bf_informed(
      factor_levels = factor_levels,
      Hr = Hr2,
      a = a,
      x = counts,
      niter = 5e3,
      bf_type = 'BFre',
      seed = 4
    )
  # total BF
  expect_equal(
    output_total_index$bf_list$bf,
    data.frame(
      LogBFer = -4.071621,
      BFer = 0.01704974,
      BFre = 58.65193
    ),
    tolerance = 0.002
  )
  # inequality BF
  expect_equal(
    output_total_index$bf_list$logBFe_inequalities,
    data.frame(
      logBFe_inequalities = -1.73936,
      logml_prior = -1.79188,
      logml_post = -0.05252
    ),
    tolerance = 0.002
  )
  # equality BF
  expect_equal(
    output_total_index$bf_list$logBFe_equalities,
    data.frame(logBFe_equalities = -2.33227),
    tolerance = 0.002
  )
  
  ## specify hypothesis using mix of explicit names and indeces
  Hr3              <- c('roundYellow > 2 = 3 > wrinkledGreen')
  output_total_mix <-
    mult_bf_informed(
      factor_levels = factor_levels,
      Hr = Hr3,
      a = a,
      x = counts,
      niter = 5e3,
      bf_type = 'BFre',
      seed = 4
    )
  # total BF
  expect_equal(
    output_total_mix$bf_list$bf,
    data.frame(
      LogBFer = -4.071621,
      BFer = 0.01704974,
      BFre = 58.65193
    ),
    tolerance = 0.002
  )
  # inequality BF
  expect_equal(
    output_total_mix$bf_list$logBFe_inequalities,
    data.frame(
      logBFe_inequalities = -1.73936,
      logml_prior = -1.79188,
      logml_post = -0.05252
    ),
    tolerance = 0.002
  )
  # equality BF
  expect_equal(
    output_total_mix$bf_list$logBFe_equalities,
    data.frame(logBFe_equalities = -2.33227),
    tolerance = 0.002
  )
})

test_that("yields equal BF estimates for costraints with free parameters", {
  data(peas)
  a             <- c(1, 1, 1, 1)
  counts        <- peas$counts
  factor_levels <- levels(peas$peas)
  
  Hr           <-
    c('roundYellow > wrinkledYellow , roundGreen > wrinkledGreen')
  output_total <-
    mult_bf_informed(
      factor_levels = factor_levels,
      Hr = Hr,
      a = a,
      x = counts,
      niter = 5e3,
      bf_type = 'BFre',
      seed = 4
    )
  # total BF
  expect_equal(
    output_total$bf_list$bf,
    data.frame(
      LogBFer = -2.407356,
      BFer = 0.09005307,
      BFre = 11.10456
    ),
    tolerance = 0.002
  )
  # inequality BF
  expect_equal(
    output_total$bf_list$logBFe_inequalities,
    data.frame(
      logBFe_inequalities = -2.407356,
      logml_prior = -2.481806,
      logml_post = -0.07445036
    ),
    tolerance = 0.002
  )
})

test_that("yields equal BF estimates for costraints with free and equal parameters",
          {
            data(peas)
            a             <- c(1, 1, 1, 1)
            counts        <- peas$counts
            factor_levels <- levels(peas$peas)
            
            Hr           <-
              c('roundYellow > wrinkledYellow , roundGreen = wrinkledGreen')
            output_total <-
              mult_bf_informed(
                factor_levels = factor_levels,
                Hr = Hr,
                a = a,
                x = counts,
                niter = 5e3,
                bf_type = 'BFre',
                seed = 4
              )
            # total BF
            expect_equal(
              output_total$bf_list$bf,
              data.frame(
                LogBFer = 18.19921,
                BFer = 80133886,
                BFre = 1.247912e-08
              ),
              tolerance = 0.002
            )
            # inequality BF
            expect_equal(
              output_total$bf_list$logBFe_inequalities,
              data.frame(
                logBFe_inequalities = -1.16125,
                logml_prior = -0.8742136,
                logml_post = 0.2870361
              ),
              tolerance = 0.002
            )
            # equality BF
            expect_equal(
              output_total$bf_list$logBFe_equalities,
              data.frame(logBFe_equalities = 19.36046),
              tolerance = 0.002
            )
          })

test_that("yields equal BF estimates for example with multiple equality constraints",
          {
            data(peas)
            a             <- c(1, 1, 1, 1)
            counts        <- peas$counts
            factor_levels <- levels(peas$peas)
            Hr            <-
              c('roundYellow = wrinkledYellow & roundGreen = wrinkledGreen')
            output_total  <-
              mult_bf_informed(
                factor_levels = factor_levels,
                Hr = Hr,
                a = a,
                x = counts,
                niter = 5e3,
                bf_type = 'BFre',
                seed = 4
              )
            # equality BF
            expect_equal(
              output_total$bf_list$logBFe_equalities,
              structure(
                list(logBFe_equalities = c(54.8269578439159, 19.36045897615)),
                row.names = c(NA,-2L),
                class = "data.frame"
              )
            )
          })

test_that("marginal likelihoods are computed correctly in Haberman example",
          {
            x <-
              c(15, 11, 14, 17,  5, 11, 10,  4,  8, 10,  7,  9, 11,  3,  6,  1,  1,  4)
            a <- rep(1, 18)
            factor_levels <- paste0('theta', 1:18)
            Hr <-
              paste0(paste(factor_levels[-18], collapse = ' > '), ' > theta18')
            output_total  <-
              mult_bf_informed(
                x,
                Hr,
                a,
                factor_levels,
                seed = 2020,
                niter = 2e3,
                bf_type = 'BFre'
              )
            
            expect_equal(summary(output_total)$logmlHe,-52.3340886139821)
            expect_equal(summary(output_total)$logmlH0,-55.6338529642158)
            expect_equal(summary(output_total)$bf, 166.87)
            expect_equal(exp(
              summary(output_total)$logmlHe - summary(output_total)$logmlH0
            ),
            27.1062505863656)
            
          })

test_that("marginal likelihoods are the same for some cases of He, H0 and Hr", {
  a <- c(1, 1, 1, 1, 1)
  x <- c(5, 10, 15, 20, 25)
  factor_levels <- paste0('theta', 1:5)
  output_total  <- mult_bf_informed(
    x = x,
    a = a,
    Hr = c(
      'theta1',
      ',',
      'theta2',
      ',',
      'theta3',
      ',',
      'theta4',
      ',',
      'theta5'
    ),
    factor_levels = factor_levels,
    seed = 2020,
    bf_type = 'BFer'
  )
  expect_equal(summary(output_total)$logmlHe,
               summary(output_total)$logmlHr)
  
  output_total  <- mult_bf_informed(
    x = x,
    a = a,
    Hr = c(
      'theta1',
      '=',
      'theta2',
      '=',
      'theta3',
      '=',
      'theta4',
      '=',
      'theta5'
    ),
    factor_levels = factor_levels,
    seed = 2020,
    bf_type = 'BF0r'
  )
  expect_equal(summary(output_total)$logmlH0,
               summary(output_total)$logmlHr)
  
})


test_that("marginal likelihoods are factored correctly", {
  # Hr: t1 = t2 = t3 & t4 , t5
  x <- c(5, 10, 15, 20, 25)
  a <- c(1, 1, 1, 1, 1)
  factor_levels <- paste0('theta', 1:5)
  output_total1  <- mult_bf_informed(
    x = x,
    a = a,
    Hr = c(
      'theta1',
      '=',
      'theta2',
      '=',
      'theta3',
      '&',
      'theta4',
      ',',
      'theta5'
    ),
    factor_levels = factor_levels,
    seed = 2020,
    bf_type = 'BF0r'
  )
  logmlHr_combined <- summary(output_total1)$logmlHr
  
  
  # H0: t1 = t2 = t3
  output_total1  <- mult_bf_informed(
    x = x[1:3],
    a = a[1:3],
    Hr = c('theta1', '<', 'theta2', '<', 'theta3'),
    factor_levels = c('theta1', 'theta2', 'theta3'),
    seed = 2020,
    bf_type = 'BF0r'
  )
  logmlH0 <- summary(output_total1)$logmlH0 # -6.091308
  
  # He: t4 , t5 | t1 = t2 = t3
  x <- c(30, 20, 25)
  a <- c(3, 1, 1)
  output_total1  <- mult_bf_informed(
    x = x,
    a = a,
    Hr = c('theta1', '<', 'theta2', '<', 'theta3'),
    factor_levels = c('theta1', 'theta2', 'theta3'),
    seed = 2020,
    bf_type = 'BF0r'
  )
  logmlHe_given_h0 <- summary(output_total1)$logmlHe # -8.016066
  
  expect_equal(logmlHr_combined,-14.1073736951353)
  expect_equal(logmlHr_combined, (logmlH0 + logmlHe_given_h0))
  
})
