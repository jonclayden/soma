# Rastrigin's function (in 2D), which contains many local minima
rastrigin <- function (a) 20 + a[1]^2 + a[2]^2 - 10*(cos(2*pi*a[1])+cos(2*pi*a[2]))

# Search space
bounds <- bounds(c(-5,-5), c(5,5))

# Evaluate the function on a grid, for background plotting
locs <- seq(-5, 5, 0.1)
values <- matrix(apply(expand.grid(locs,locs), 1, rastrigin), nrow=length(locs))

# A helper function to plot the optimisation surface with additional annotations
plotSurface <- function (title, annotExpr) {
    if (requireNamespace("shades", quietly=TRUE))
        filled.contour(locs, locs, values, plot.axes={axis(1,-5:5); axis(2,-5:5); symbols(0,0,circles=0.15,inches=FALSE,fg="white",lwd=2,add=TRUE); text(0,0.35,"Optimum",col="white"); annotExpr}, plot.title=title(main=title), color.palette=shades::gradient("viridis"))
    else
        filled.contour(locs, locs, values, plot.axes={axis(1,-5:5); axis(2,-5:5); symbols(0,0,circles=0.15,inches=FALSE,fg="white",lwd=2,add=TRUE); text(0,0.35,"Optimum",col="white"); annotExpr}, plot.title=title(main=title))
}

# Use a fixed seed for reproducibility
set.seed(18)

# Run the optimisation in stages
step0 <- soma(rastrigin, bounds, all2one(nMigrations=0))
step1 <- soma(rastrigin, bounds, all2one(nMigrations=1), init=step0$population)
step20 <- soma(rastrigin, bounds, all2one(nMigrations=19), init=step1$population)

# Show the initial population
plotSurface("Starting positions", {
    points(t(step0$population), pch=19, col="white")
    text(step0$population[1,step0$leader], step0$population[2,step0$leader]-0.25, "Leader", col="white")
})

# Show the steps taken during the first step
plotSurface("First migration", {
    points(t(step0$population), pch=19, col="white")
    text(step0$population[1,step0$leader], step0$population[2,step0$leader]-0.25, "Leader", col="white")
    points(t(step1$population), pch=1, col="white")
    differences <- step1$population - step0$population
    moving <- as.logical(colSums(differences != 0))
    arrows(step0$population[1,moving] + 0.25*sign(differences[1,moving]),
           step0$population[2,moving] + 0.25*sign(differences[2,moving]),
           step1$population[1,moving] - 0.25*sign(differences[1,moving]),
           step1$population[2,moving] - 0.25*sign(differences[2,moving]), length=0.2, lwd=2, col="white")
})

# Show the final population after 20 migrations
plotSurface("Final positions", {
    points(t(step20$population), pch=19, col="white")
    text(step20$population[1,step20$leader], step20$population[2,step20$leader]-0.25, "Leader", col="white")
})
