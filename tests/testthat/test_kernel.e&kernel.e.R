test_that("kernel.e() returns a double", {
  controls=test_data$controls
  output_list <- kernel.e(controls)
  expect_type(output_list,"double" )
})

test_that("kernel.g() returns a numeric", {
  controls=test_data$controls
  output_list <- kernel.g(controls)
  expect_type(output_list,"double" )
})
