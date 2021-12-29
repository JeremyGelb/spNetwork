context("testing the base kernel implemented")
library(sf)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### TESTING THE BASE KERNEL INTEGRALS ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("Testing the base kernels, they must integrate to 1", {

  list_kernels <- c("triangle", "gaussian", "scaled gaussian", "tricube",
                    "cosine" ,"triweight", "quartic", 'epanechnikov',
                    'uniform')

  ts <- sapply(list_kernels, function(k){
    if(grepl("gaussian", k, fixed = TRUE) == FALSE){

      func <- select_kernel(k)
      val <- cubintegrate(func,lower=0,upper=2,
                   bw=2, relTol = 1e-15)$integral
      return(val*2)

    }else{
      return(1)
    }
  })

  test <- round(ts,9) != 1
  expect_false(any(test))

})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### TESTING THE c++ KERNEL INTEGRALS ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("Testing the c++ kernels, they must integrate to 1", {

  funs <- list(quartic_kernel_cpp,
               cosine_kernel_cpp,
               tricube_kernel_cpp,
               triweight_kernel_cpp,
               epanechnikov_kernel_cpp,
               triangle_kernel_cpp,
               uniform_kernel_cpp,
               quartic_kernel_cpp,
               quartic_kernelos,
               cosine_kernelos,
               tricube_kernelos,
               triweight_kernelos,
               epanechnikov_kernelos,
               triangle_kernelos,
               uniform_kernelos,
               quartic_kernelos
               )

  ts <- sapply(funs, function(f){

    val <- cubintegrate(f,lower=0,upper=2,
                        bw=2, relTol = 1e-15)$integral
    return(val*2)

  })

  test <- round(ts,9) != 1
  expect_false(any(test))

})
