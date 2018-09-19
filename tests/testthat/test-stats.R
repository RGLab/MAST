context('Testing stats APIs')
test_that('gene stats',{
  res <- freq(fd)
  expect_equal(fivenum(res), c(MKI67 = 0.00657894736842105, LIF = 0.0888157894736842, PECAM1 = 0.160087719298246, PRDM1 = 0.422149122807018, GAPDH = 0.857456140350877))
  res <- condmean(fd)
  expect_equal(fivenum(res), c(CXCL13 = 13.74116337, CD109 = 17.8044699829261, BCL2 = 18.859895658, GATA3 = 20.2697726769207, IL3 = 22.7434740788889))
  res <- condSd(fd)
  expect_equal(fivenum(res), c(B3GAT1 = 0.799034361476095, RORC = 2.14860685194272, CD45 = 2.57447658007055, CCR4 = 3.03002336654144, CXCL13 = 5.49024005052408))
  res <- numexp(fd)
  expect_equal(fivenum(res), c(MKI67 = 3, LIF = 40.5, PECAM1 = 73, PRDM1 = 192.5, GAPDH = 391))
  
})

# dd <- fivenum(res)
# deparse(dd, width.cutoff = 500L)
