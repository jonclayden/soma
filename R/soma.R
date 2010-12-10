soma <- function (costFunction, bounds, options = list(), strategy = "all2one", ...)
{
    if (!(all(c("min","max") %in% names(bounds))))
        output(OL$Error, "Bounds list must contain \"min\" and \"max\" vector elements")
    if (length(bounds$min) != length(bounds$max))
        output(OL$Error, "Bounds are not of equal length")
    if (strategy != "all2one")
        output(OL$Error, "Only the \"all2one\" strategy is currently supported")
    
    defaultOptions <- list(pathLength=3, stepLength=0.11, perturbationChance=0.1, minAbsoluteSep=0, minRelativeSep=1e-3, nMigrations=20, populationSize=10)
    defaultsNeeded <- setdiff(names(defaultOptions), names(options))
    options[defaultsNeeded] <- defaultOptions[defaultsNeeded]
    
    nParams <- length(bounds$min)
    nParamsTotal <- nParams * options$populationSize
    steps <- seq(0, options$pathLength, options$stepLength)
    nSteps <- length(steps)
    steps <- rep(steps, each=nParamsTotal)
    
    population <- matrix(runif(nParamsTotal), nrow=nParams, ncol=options$populationSize)
    population <- population * (bounds$max-bounds$min) + bounds$min
    
    migrationCount <- 0
    leaderCostHistory <- numeric(0)
    
    repeat
    {
        # Values over actual locations
        costFunctionValues <- apply(population, 2, costFunction, ...)
        leaderIndex <- which.min(costFunctionValues)
        leaderValue <- costFunctionValues[leaderIndex]
        separationOfExtremes <- max(costFunctionValues) - leaderValue
        sumOfExtremes <- max(costFunctionValues) + leaderValue
        
        leaderCostHistory <- c(leaderCostHistory, leaderValue)
        
        # Check termination criteria
        if (migrationCount == options$nMigrations)
        {
            output(OL$Info, "Migration limit (", options$nMigrations, ") reached - stopping")
            break
        }
        if (separationOfExtremes < options$minAbsoluteSep)
        {
            output(OL$Info, "Absolute cost separation (", signif(separationOfExtremes,3), ") is below threshold (", signif(options$minAbsoluteSep,3), ") - stopping")
            break
        }
        if (separationOfExtremes/sumOfExtremes < options$minRelativeSep)
        {
            output(OL$Info, "Relative cost separation (", signif(separationOfExtremes/sumOfExtremes,3), ") is below threshold (", signif(options$minRelativeSep,3), ") - stopping")
            break
        }
        
        directionsFromLeader <- apply(population, 2, "-", population[,leaderIndex])
        toPerturb <- runif(nParamsTotal) < options$perturbationChance
        
        # Second line here has a minus because directions are away from leader
        populationSteps <- array(rep(population,nSteps), dim=c(nParams,options$populationSize,nSteps))
        populationSteps <- populationSteps - steps * rep(directionsFromLeader * toPerturb, nSteps)
        
        outOfBounds <- which(populationSteps < bounds$min | populationSteps > bounds$max)
        randomSteps <- array(runif(nParamsTotal*nSteps), dim=c(nParams,options$populationSize,nSteps))
        randomSteps <- randomSteps * (bounds$max-bounds$min) + bounds$min
        populationSteps[outOfBounds] <- randomSteps[outOfBounds]
        
        # Values over potential locations
        costFunctionValues <- apply(populationSteps, 2:3, costFunction, ...)
        individualBestLocs <- apply(costFunctionValues, 1, which.min)
        
        indexingMatrix <- cbind(seq_len(options$populationSize), individualBestLocs)
        population <- t(apply(populationSteps, 1, "[", indexingMatrix))
        
        migrationCount <- migrationCount + 1
        if (migrationCount %% 10 == 0)
            output(OL$Verbose, "Completed ", migrationCount, " migrations")
    }
    
    output(OL$Info, "Leader is #", leaderIndex, ", with cost ", signif(costFunctionValues[leaderIndex],3))
    
    returnValue <- list(leader=leaderIndex, population=population, cost=costFunctionValues, history=leaderCostHistory, migrations=migrationCount)
    class(returnValue) <- "soma"
    
    return (returnValue)
}
