test_that("likbox() returns a valid result", {
  h=-1.6
  controls=test_data$controls
  cases=test_data$cases
  output_list <- likbox(h,c(controls,cases),n=length(controls))
  expect_type(output_list,"double" )
  expect_false(is.na(output_list), info = "Output should not be NA")
  expect_false(any(output_list == -Inf), info = "Output should not be negative infinity")
  expect_false(any(output_list == Inf), info = "Output should not be positive infinity")
})

test_that("likbox() returns a valid result", {
  h=0
  controls=test_data$controls
  cases=test_data$cases
  output_list <- likbox(h,c(controls,cases),n=length(controls))
  expect_type(output_list,"double" )
  expect_false(is.na(output_list), info = "Output should not be NA")
  expect_false(any(output_list == -Inf), info = "Output should not be negative infinity")
  expect_false(any(output_list == Inf), info = "Output should not be positive infinity")
})
