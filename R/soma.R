#' @import reportr
#' @importFrom graphics plot
#' @importFrom stats runif
NULL

#' Options for the available SOMA variants
#' 
#' These functions generate option lists (and provide defaults) for the SOMA
#' algorithm variants available in the package, which control how the algorithm
#' will proceed and when it will terminate. Each function corresponds to a
#' different top-level strategy, described in a different reference.
#' 
#' All To One (the \code{all2one} function) is the original SOMA strategy. At
#' each ``migration'', the cost function is evaluated for all individuals in
#' the population, and the one with the lowest value is designated the
#' ``leader''. All other individuals migrate towards the leader's position in
#' some or all dimensions of the parameter space, with a fixed probability of
#' perturbation in each dimension. Each migration is evaluated against the cost
#' function at several points on the line towards the leader, and the location
#' with the lowest value becomes the individual's starting position for the
#' next migration.
#' 
#' The Team To Team Adaptive (T3A) strategy (Diep, 2019) differs in that only a
#' random subset of individuals are selected into a migrant pool and a leader
#' pool for any given migration. A subset of most optimal migrants are then
#' migrated towards the single most optimal individual from the leader pool.
#' The perturbation probability and step length along the trajectory towards
#' the leader also vary according to formulae given by the strategy author as
#' the algorithm progresses through the migrations.
#' 
#' In the Pareto strategy (Diep et al., 2019), all individuals are sorted by
#' cost function value at the start of each migration. The leader is selected
#' randomly from the top 4\% (20\% of 20\%) of most optimal individuals, and
#' a single migrant is chosen at random from between the 20th and the 36th
#' percentiles of the population (the top 20\% of the bottom 80\%). The
#' perturbation probability and the step length again vary across migrations,
#' but this time in a sinusoidal fashion, and the migrant is updated in all
#' dimensions, but some more slowly than others.
#' 
#' @param populationSize The number of individuals in the population. It is
#'   recommended that this be somewhat larger than the number of parameters
#'   being optimised over, and it should not be less than 2. The default varies
#'   by strategy.
#' @param nMigrations The maximum number of migrations to complete.
#' @param pathLength The distance towards the leader that individuals may
#'   migrate. A value of 1 corresponds to the leader's position itself, and
#'   values greater than one (recommended) allow for some overshoot.
#' @param stepLength The granularity at which potential steps are evaluated.
#'   It is recommended that the \code{pathLength} not be a whole multiple of
#'   this value.
#' @param perturbationChance The probability that individual parameters are
#'   changed on any given step.
#' @param minAbsoluteSep The smallest absolute difference between the maximum
#'   and minimum cost function values. If the difference falls below this
#'   minimum, the algorithm will terminate. The default is 0, meaning that this
#'   termination criterion will never be met.
#' @param minRelativeSep The smallest relative difference between the maximum
#'   and minimum cost function values. If the difference falls below this
#'   minimum, the algorithm will terminate.
#' @param nSteps The number of candidate steps towards the leader per migrating
#'   individual. This option is used instead of \code{pathLength} and
#'   \code{stepLength} under the T3A strategy, since the step length is
#'   variable.
#' @param migrantPoolSize,leaderPoolSize The number of randomly selected
#'   individuals to include in the migrant and leader pools, respectively,
#'   under the T3A strategy.
#' @param nMigrants The number of individuals that will migrate, at each
#'   migration, under the T3A strategy.
#' @return A list of class \code{"soma.options"}.
#' @references
#'   I. Zelinka (2004). SOMA - self-organizing migrating algorithm. In G.C.
#'   Onwubolu & B.V. Babu, eds, New optimization techniques in engineering.
#'   Volume 141 of ``Studies in Fuzziness and Soft Computing'', pp. 167-217.
#'   Springer.
#'   
#'   Q. B. Diep (2019). Self-Organizing Migrating Algorithm Team To Team
#'   Adaptive – SOMA T3A. 
#'   
#' @author Jon Clayden <code@@clayden.org>
#' @aliases soma.options
#' @rdname soma.options
#' @export
all2one <- function (populationSize = 10L, nMigrations = 20L, pathLength = 3, stepLength = 0.11, perturbationChance = 0.1, minAbsoluteSep = 0, minRelativeSep = 1e-3)
{
    assert(populationSize >= 2L, "Population size must be at least 2")
    options <- list(strategy="all2one", populationSize=populationSize, nMigrations=nMigrations, pathLength=pathLength, stepLength=stepLength, perturbationChance=perturbationChance, minAbsoluteSep=minAbsoluteSep, minRelativeSep=minRelativeSep)
    class(options) <- "soma.options"
    return (options)
}

#' @rdname soma.options
#' @export
t3a <- function (populationSize = 30L, nMigrations = 20L, nSteps = 45L, migrantPoolSize = 10L, leaderPoolSize = 10L, nMigrants = 4L, minAbsoluteSep = 0, minRelativeSep = 1e-3)
{
    assert(populationSize >= 2L, "Population size must be at least 2")
    assert(nMigrants <= migrantPoolSize, "The number of migrants can't be larger than the migrant pool size")
    options <- list(strategy="t3a", populationSize=populationSize, nMigrations=nMigrations, nSteps=nSteps, migrantPoolSize=migrantPoolSize, leaderPoolSize=leaderPoolSize, nMigrants=nMigrants, minAbsoluteSep=minAbsoluteSep, minRelativeSep=minRelativeSep)
    class(options) <- "soma.options"
    return (options)
}

#' @rdname soma.options
#' @export
pareto <- function (populationSize = 100L, nMigrations = 20L, nSteps = 10L, perturbationFrequency = 1, stepFrequency = 1, minAbsoluteSep = 0, minRelativeSep = 1e-3)
{
    assert(populationSize >= 2L, "Population size must be at least 2")
    options <- list(strategy="pareto", populationSize=populationSize, nMigrations=nMigrations, nSteps=nSteps, perturbationFrequency=perturbationFrequency, stepFrequency=stepFrequency, minAbsoluteSep=minAbsoluteSep, minRelativeSep=minRelativeSep)
    class(options) <- "soma.options"
    return (options)
}

"%||%" <- function (X, Y) { if (is.null(X) || length(X)==0) Y else X }

#' The Self-Organising Migrating Algorithm
#' 
#' The Self-Organising Migrating Algorithm is a general-purpose, stochastic
#' optimisation algorithm. The approach is similar to that of genetic
#' algorithms, although it is based on the idea of a series of ``migrations''
#' by a fixed set of individuals, rather than the development of successive
#' generations. It can be applied to any cost-minimisation problem with a
#' bounded parameter space, and is robust to local minima.
#' 
#' @param costFunction A cost function which takes a numeric vector of
#'   parameters as its first argument, and returns a numeric scalar
#'   representing the associated cost value.
#' @param bounds A list with elements \code{min} and \code{max}, each a numeric
#'   vector giving the upper and lower bounds for each parameter, respectively.
#' @param min,max Vectors of minimum and maximum bound values for each
#'   parameter to the \code{costFunction}.
#' @param options A list of options for the SOMA algorithm itself, usually
#'   generated by functions like \code{\link{all2one}}.
#' @param init An optional matrix giving the starting population's positions in
#'   parameter space, one per column. If omitted, initialisation is random (as
#'   is usual for SOMA), but specifying a starting state can be helpful when
#'   running the algorithm in stages or investigating the consistency of
#'   solutions.
#' @param \dots Additional parameters to \code{costFunction} (for \code{soma})
#'   or the default plotting method (for \code{plot.soma}).
#' @param x An object of class \code{"soma"}.
#' @param y Ignored.
#' @return A list of class \code{"soma"}, containing the following elements.
#'   \describe{
#'     \item{leader}{The index of the ``leader'', the individual in the
#'       population with the lowest cost.}
#'     \item{population}{A matrix whose columns give the parameter values for
#'       each individual in the population at convergence.}
#'     \item{cost}{A vector giving the cost function values for each individual
#'       at convergence.}
#'     \item{history}{A vector giving the cost of the leader for each migration
#'       during the optimisation. This should be nonincreasing.}
#'     \item{migrations}{The number of migrations completed.}
#'     \item{evaluations}{The number of times the \code{costFunction} was
#'       evaluated.}
#'   }
#'   A \code{plot} method is available for this class, which shows the history
#'     of leader cost values during the optimisation.
#' @examples
#' # Rastrigin's function, which contains many local minima
#' rastrigin <- function (a) 20 + a[1]^2 + a[2]^2 - 10*(cos(2*pi*a[1])+cos(2*pi*a[2]))
#' 
#' # Find the global minimum over the range -5 to 5 in each parameter
#' x <- soma(rastrigin, bounds(c(-5,-5), c(5,5)))
#' 
#' # Find the location of the leader - should be near the true minimum of c(0,0)
#' print(x$population[,x$leader])
#' 
#' # Plot the cost history of the leaders
#' plot(x)
#' @references
#'   I. Zelinka (2004). SOMA - self-organizing migrating algorithm. In G.C.
#'   Onwubolu & B.V. Babu, eds, New optimization techniques in engineering.
#'   Volume 141 of ``Studies in Fuzziness and Soft Computing'', pp. 167-217.
#'   Springer.
#' @author R implementation by Jon Clayden <code@@clayden.org>.
#' @seealso \code{\link{soma.options}} for setting options. \code{\link{optim}}
#'   implements other general-purpose optimisation methods.
#' @export
soma <- function (costFunction, bounds, options = list(), init = NULL, ...)
{
    # Check bounds and options
    if (!inherits(bounds, "soma.bounds"))
        bounds <- do.call(soma::bounds, bounds)
    if (!inherits(options, "soma.options"))
        options <- do.call(soma::all2one, options)
    
    nParams <- length(bounds$min)
    nParamsTotal <- nParams * options$populationSize
    if (!is.null(options$pathLength) && !is.null(options$stepLength))
        steps <- seq(0, options$pathLength, options$stepLength)
    nSteps <- options$nSteps %||% length(steps)
    perturbationChance <- options$perturbationChance %||% 1
    
    leaderPoolIndices <- migrantPoolIndices <- NULL
    if (options$strategy == "pareto")
    {
        leaderPoolIndices <- seq_len(ceiling(options$populationSize / 25))
        migrantPoolIndices <- seq_len(ceiling(options$populationSize / 6.25)) + round(options$populationSize / 5)
    }
    
    # Create the population
    if (!is.null(init))
        population <- init
    else
    {
        population <- matrix(runif(nParamsTotal), nrow=nParams, ncol=options$populationSize)
        population <- population * (bounds$max-bounds$min) + bounds$min
    }
    
    # Calculate initial costs
    evaluations <- NULL
    evaluationCount <- 0
    costFunctionWrapper <- function(x) { evaluationCount <<- evaluationCount + 1; costFunction(x, ...) }
    costFunctionValues <- apply(population, 2, costFunctionWrapper)
    
    migrationCount <- 0
    leaderCostHistory <- numeric(0)
    
    report(OL$Info, "Starting SOMA optimisation")
    
    repeat
    {
        progress <- migrationCount / options$nMigrations
        
        if (options$strategy == "all2one")
        {
            # Find the current leader
            leader <- which.min(costFunctionValues)
            migrants <- seq_len(options$populationSize)
        }
        else if (options$strategy == "t3a")
        {
            perturbationChance <- 0.05 + 0.9 * progress
            steps <- (seq_len(nSteps)-1) * (0.15 - 0.08 * progress)
            leaderPool <- sample(options$populationSize, options$leaderPoolSize)
            leader <- leaderPool[which.min(costFunctionValues[leaderPool])]
            migrantPool <- sample(options$populationSize, options$migrantPoolSize)
            migrants <- migrantPool[order(costFunctionValues[migrantPool])[seq_len(options$nMigrants)]]
        }
        else if (options$strategy == "pareto")
        {
            perturbationChance <- 0.5 + 0.45 * cos(options$perturbationFrequency * pi * progress + pi)
            steps <- (seq_len(nSteps)-1) * (0.35 + 0.15 * cos(options$stepFrequency * pi * progress))
            order <- order(costFunctionValues)
            leader <- sample(order[leaderPoolIndices], 1)
            migrants <- sample(order[migrantPoolIndices], 1)
        }
        
        leaderValue <- costFunctionValues[leader]
        separationOfExtremes <- max(costFunctionValues) - min(costFunctionValues)
        sumOfExtremes <- max(costFunctionValues) + min(costFunctionValues)
        
        # Check termination criteria
        if (migrationCount == options$nMigrations)
        {
            report(OL$Info, "Migration limit (#{options$nMigrations}) reached - stopping")
            break
        }
        if (separationOfExtremes < options$minAbsoluteSep)
        {
            report(OL$Info, "Absolute cost separation (#{separationOfExtremes}) is below threshold (#{options$minAbsoluteSep}) - stopping", signif=3)
            break
        }
        # isTRUE() needed here in case extremes are infinite: Inf/Inf => NaN
        if (isTRUE(abs(separationOfExtremes/sumOfExtremes) < options$minRelativeSep))
        {
            report(OL$Info, "Relative cost separation (#{separationOfExtremes/sumOfExtremes}) is below threshold (#{options$minRelativeSep}) - stopping", signif=3)
            break
        }
        
        # Establish which parameters will be changed, and which individuals will migrate
        toPerturb <- array(FALSE, dim=dim(population))
        toPerturb[,migrants] <- runif(nParams * length(migrants)) < perturbationChance
        if (options$strategy == "pareto")
            toPerturb[,migrants] <- ifelse(toPerturb[,migrants], 1, progress)
        toMigrate <- colSums(toPerturb) > 0
        toMigrate[leader] <- FALSE
        nMigrating <- sum(toMigrate)
        
        # No-op migrations aren't counted
        if (nMigrating == 0)
        {
            report(OL$Verbose, "No parameters to perturb - skipping migration")
            next
        }
        else
            report(OL$Verbose, "Migration ##{migrationCount+1}: #{nMigrating} individuals moving towards leader ##{leader} with cost #{leaderValue}", signif=3)
        
        # Find the migration direction for each individual
        # NB: This form exploits recycling, which depends on the shape of "population" (params x individuals)
        directionsFromLeader <- population[,toMigrate] - population[,leader]
        
        # Second line here has a minus because directions are away from leader
        populationSteps <- array(rep(population[,toMigrate],nSteps), dim=c(nParams,nMigrating,nSteps))
        populationSteps <- populationSteps - rep(steps, each=nParams*nMigrating) * rep(directionsFromLeader * toPerturb[,toMigrate], nSteps)
        
        # Replace out-of-bounds parameters with random valid values
        outOfBounds <- which(populationSteps < bounds$min | populationSteps > bounds$max)
        randomSteps <- array(runif(nParamsTotal*nSteps), dim=c(nParams,options$populationSize,nSteps))
        randomSteps <- randomSteps * (bounds$max-bounds$min) + bounds$min
        populationSteps[outOfBounds] <- randomSteps[outOfBounds]
        
        # Values over potential locations
        allCostFunctionValues <- apply(populationSteps, 2:3, costFunctionWrapper)
        individualBestLocs <- apply(allCostFunctionValues, 1, which.min)
        
        # Migrate each individual to its best new location, and update costs
        indexingMatrix <- cbind(seq_len(nMigrating), individualBestLocs)
        subpopulation <- t(apply(populationSteps, 1, "[", indexingMatrix))
        population[,toMigrate] <- subpopulation
        bestCostFunctionValues <- allCostFunctionValues[indexingMatrix]
        costFunctionValues[toMigrate] <- bestCostFunctionValues
        
        migrationCount <- migrationCount + 1
        evaluations <- c(evaluations, evaluationCount)
        leaderCostHistory <- c(leaderCostHistory, min(costFunctionValues))
    }
    
    leader <- which.min(costFunctionValues)
    report(OL$Info, "Final leader is ##{leader}, with cost #{costFunctionValues[leader]}", signif=3)
    
    returnValue <- list(leader=leader, population=population, cost=costFunctionValues, history=leaderCostHistory, migrations=migrationCount, evaluations=evaluations)
    class(returnValue) <- "soma"
    
    return (returnValue)
}

#' @rdname soma
#' @export
bounds <- function (min, max)
{
    min <- as.numeric(min)
    max <- as.numeric(max)
    assert(length(min) == length(max), "Bounds are not of equal length")
    assert(all(is.finite(min)) && all(is.finite(max)), "Bounds must be finite")
    assert(all(min <= max), "At least one minimum is greater than the equivalent maximum", level=OL$Warning)
    assert(all(min != max), "At least one minimum and maximum bound are equal - consider making these fixed arguments to the cost function", level=OL$Warning)
    return (structure(list(min=min, max=max), class="soma.bounds"))
}

#' @rdname soma
#' @export
plot.soma <- function (x, y = NULL, ...)
{
    plot(seq_along(x$history), x$history, xlab="Migration number", ylab="Leader cost value", type="b", ...)
}
