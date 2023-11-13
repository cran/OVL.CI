test_that("kernel.e.density() returns a valid result", {
  x=test_data$controls
  gridd = seq(-20,20,length.out=1000)
  h = (4/3)^(1/5)*sd(x)*length(x)^(-1/5)
  output_list <- kernel.e.density(x,gridd,h)
  expect_type(output_list,"double" )
  expect_true(all(output_list>=0&output_list<=1))
})

test_that("kernel.g.density() returns a valid result", {
  x=test_data$controls
  gridd = seq(-20,20,length.out=1000)
  h = (4/3)^(1/5)*sd(x)*length(x)^(-1/5)
  output_list <- kernel.g.density(x,gridd,h)
  expect_type(output_list,"double" )
  expect_true(all(output_list>=0&output_list<=1))
})
