context("test_versioning")

test_that("test typical inputs", {
  expect_true(Version("1.0.0") < Version("2.0.0"))
  expect_true(Version("1.0.0") < Version("1.1.0"))
  expect_true(Version("1.0.0") < Version("1.0.1"))
  expect_false(Version("1.0.0") < Version("1.0.0"))
  expect_true(Version("1.0.0") <= Version("1.0.0"))
  
  expect_false(Version("1.0.0") >= Version("2.0.0"))
  expect_false(Version("1.0.0") >= Version("1.1.0"))
  expect_false(Version("1.0.0") >= Version("1.0.1"))
  expect_true(Version("1.0.0") >= Version("1.0.0"))
  expect_false(Version("1.0.0") > Version("1.0.0"))
})

test_that("test tuple construction", {
  expect_true(Version("0100.2.03") == Version(c(100, 2, 3)))
  expect_true(Version("1.2.3") == Version(c(1, 2, 3, NA)))
  expect_true(Version("1.2.3") == Version(c(1, 2, 3, "junk")))
  expect_true(Version("1.2.3") == Version(c(1, 2, 3, -1)))
})

test_that("test local version identifiers", {
  expect_true(Version("1.0.0") == Version("1.0.0+1"))
  expect_true(Version("1.0.0") == Version("1.0.0+xxx"))
  expect_true(Version("1.0.0") == Version("1.0.0+x.y.z"))
})
