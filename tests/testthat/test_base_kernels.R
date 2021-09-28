context("testing the base kernel implemented")

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
      val <- cubintegrate(triangle_kernel,lower=0,upper=2,
                   bw=2, relTol = 1e-15)$integral
      return(val*2)

    }else{
      return(1)
    }
  })

  test <- ts != 1
  expect_false(any(test))

})
