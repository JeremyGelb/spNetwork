context("testing function used to perform NKDE")
library(sf)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### TEST FOR BW SELECTION WITH CV LIKELIHOOD ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("Testing the clean_events function", {

  ## creating the simple situation
  # definition of some events
  event <- data.frame(x=c(2,0,3,0,5,10),
                      y=c(2,3,0,-3,5,10),
                      weight = c(1,1,1,1,1,1))
  event <- st_as_sf(event, coords = c("x","y"))

  # applying the function
  new_events <- clean_events(event,digits = 1, agg = 6)
  obs_coords <- st_coordinates(new_events)

  # checking
  exp_coords <- data.frame(
    x = c(mean(c(2,0,3,0,5)),10),
    y = c(mean(c(2,3,0,-3,5)),10)
  )

  expect_equal(sum(obs_coords - exp_coords),0)


})
