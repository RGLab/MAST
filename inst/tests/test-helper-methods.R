wfunction <- function(x, warn=TRUE){
    if(warn) warning("I warned you")
    if(!missing(x)) return(x) else return(999)
}

context("Warning/Error Handlers behave")
test_that('hushWarning hushes',{
    expect_that(wfunction(warn=TRUE), gives_warning('you'))
    expect_silent(hushWarning(wfunction(warn=TRUE),'warned'))
})

test_that("hushWarning doesn't hush", {
    expect_that(hushWarning(wfunction(warn=TRUE),'foo'),  gives_warning('you'))
})

test_that("hushWarning returns regardless", {
    suppressWarnings({
        nothushed <- hushWarning(wfunction(warn=TRUE), 'foo')
        expect_equal(nothushed, 999)
        nothushed <- hushWarning(wfunction('xxx', warn=TRUE), 'foo')
        expect_equal(nothushed, 'xxx')
    })

    hushed <- hushWarning(wfunction(warn=TRUE), 'warn')
    expect_equal(hushed, 999)
    hushed <- hushWarning(wfunction('xxx', warn=TRUE), 'warn')
    expect_equal(hushed, 'xxx')
    notwarn <- hushWarning(wfunction('xxx', warn=FALSE), 'warn')
    expect_equal(notwarn, 'xxx')
    
})
