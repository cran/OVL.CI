test_that("U() returns a valid result", {
  controls=test_data$controls
  cases=test_data$cases
  mu1=mean(controls)
  mu2=mean(cases)
  sigma1=ssdd(controls)
  sigma2=ssdd(cases)
  output_list <- U(mu1,mu2,sigma1,sigma2)
  expect_type(output_list,"double" )
  expect_false(is.na(output_list), info = "Output should not be NA")
  expect_false(any(output_list == -Inf), info = "Output should not be negative infinity")
  expect_false(any(output_list == Inf), info = "Output should not be positive infinity")
})
