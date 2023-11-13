test_that("ssdd() returns a valid result", {
  x=test_data$controls
  output_list <- ssdd(x)
  expect_type(output_list,"double" )
  expect_true(output_list>0)
  expect_false(output_list==Inf)

})
