source('common-fixtures.R')

fd.small <- fd[[seq(from=1, to=nrow(fd), by=11)]]
context('Testing filtering')
test_that('Can filter',{
ff <- filter(fd.small, groups='Number.of.Cells')
expect_that(ff, is_a('FluidigmAssay'))
ff2 <- filter(fd.small, groups='Number.of.Cells', filt_control=list(nOutlier=1, sigmaContinuous=3)) #can't set sigmaContinuous to zero because we can't construct and empty object
expect_true(nrow(ff2)<nrow(ff))
ff2 <- filter(fd.small, groups='Number.of.Cells', filt_control=list(nOutlier=1, sigmaContinuous=3, filter=FALSE), apply_filter=TRUE)
expect_that(ff2, is_a('list'))
ff3 <- filter(fd.small, apply_filter=FALSE)
expect_that(ff3, is_a('data.frame'))
})
