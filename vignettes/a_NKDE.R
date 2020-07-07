## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----message=FALSE, warning=FALSE---------------------------------------------
library(ggplot2)
library(reshape2)
library(kableExtra)
library(spNetwork)

x <- seq(-15,15,by = 0.01)
df <- data.frame(
  x = x,
  gaussian = gaussian_kernel(x,15),
  epanechnikov = epanechnikov_kernel(x,15),
  quartic = quartic_kernel(x,15),
  triangle = triangle_kernel(x,15),
  tricube = tricube_kernel(x,15),
  triweight = triweight_kernel(x,15),
  cosine = cosine_kernel(x,15))

df2 <- melt(df, id.vars = "x")
names(df2) <- c("x","kernel","y")

ggplot(df2) + 
  geom_line(aes(x=x,y=y,color=kernel))


## ----message=FALSE, warning=FALSE---------------------------------------------
Funs <- c(gaussian_kernel,epanechnikov_kernel,
          quartic_kernel,triangle_kernel,
          tricube_kernel,triweight_kernel,
          cosine_kernel)
Names <- c("gaussian", "epanechnikov", "quartic",
           "triangle", "tricube", "triweight", "cosine")
integrals <- sapply(Funs,function(f){
  return(round(integrate(f,upper=15,lower=-15,bw=15)$value,3))
         })
df <- data.frame("kernel"=Names,
                 "integrals" = integrals)

kable(df)


