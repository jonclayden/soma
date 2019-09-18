options(reportrPrefixFormat="")
reportr::setOutputLevel(Warning)

# Rastrigin's function, which contains many local minima
rastrigin <- function (a) 20 + a[1]^2 + a[2]^2 - 10*(cos(2*pi*a[1])+cos(2*pi*a[2]))

# Find the global minimum over the range -5 to 5 in each parameter
result <- soma(rastrigin, list(min=c(-5,-5),max=c(5,5)))

expect_equal(result$cost[result$leader], 0, tol=0.01)
expect_equal(result$population[,result$leader], c(0,0), tol=0.1)
