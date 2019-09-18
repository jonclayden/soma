options(reportrPrefixFormat="")
reportr::setOutputLevel(Warning)

# A very simple, convex function for testing
result <- soma(function(x) sum(x^2), list(min=c(-5,-5),max=c(5,5)))

# Check the optimum has been found (but with a wide tolerance)
expect_equal(result$cost[result$leader], 0, tol=0.1)
expect_equal(result$population[,result$leader], c(0,0), tol=0.1)
