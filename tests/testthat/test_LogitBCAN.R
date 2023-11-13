test_that("OVL.LogitBCAN() returns a list", {
  output_list <- OVL.LogitBCAN(x=test_data$controls,y=test_data$cases)
  expect_type(output_list, "list")
})

test_that("OVL.LogitBCAN() returns a valid result", {
  output_list <- OVL.LogitBCAN(x=test_data$controls,y=test_data$cases)
  expect_true( OVL.LogitBCAN(x=test_data$controls,y=test_data$cases)[[1]]>= 0.3& OVL.LogitBCAN(x=test_data$controls,y=test_data$cases)[[2]] <= 1)
})

test_that("OVL.LogitBCAN() returns a valid result", {
  set.seed(1)
  output_list <- OVL.LogitBCAN(x=rnorm(50,0,1),y=rnorm(50,0.5,0.5))
  expect_true( OVL.LogitBCAN(x=test_data$controls,y=test_data$cases)[[1]]>= 0.3& OVL.LogitBCAN(x=test_data$controls,y=test_data$cases)[[2]] <= 1)
})

test_that("OVL.LogitBCAN() returns a list", {
  set.seed(1)
  output_list <- OVL.LogitBCAN(x=rnorm(50,80,10),y=rnorm(50,0,11))
  expect_type(output_list, "list")
})


test_that("OVL.LogitBCAN returns a list", {
  set.seed(1)
  output_list <- OVL.LogitBCAN(x=rnorm(50,80,10),y=rnorm(50,80.2,11))
  expect_type(output_list, "list")
})
