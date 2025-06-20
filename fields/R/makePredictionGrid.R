makePredictionGridList <- function(mKrigObject, nx, ny, NNSize) {
  #
  # NNSize used to add additional points so offGridWeights
  # has enough neighbors.
  #
  sDimension <- ncol(mKrigObject$x)
  xr <- range(mKrigObject$x[, 1])
  dx <- (xr[2] - xr[1]) / (nx - 2 * NNSize)
  xg <- 0:(nx - 1) * dx +  (xr[1] - dx * (NNSize - 1 / 2))
  predictionGridList <- list(x = xg)
  # add y grid if 2D
  if (sDimension == 2) {
    yr <- range(mKrigObject$x[, 2])
    dy <- (yr[2] - yr[1]) / (ny - 2 * NNSize)
    yg <- 0:(ny - 1) * dy +  (yr[1] - dy * (NNSize - 1 / 2))
    predictionGridList$y <- yg
  }
  return(predictionGridList)
}

checkPredictGrid <- function(predictionGridList) {
  testX <-
    sd(diff(predictionGridList$x)) / mean(diff(predictionGridList$x))
  if (testX > 1e-8) {
    stop("predictionGridList$x must be equally spaced")
  }
  dx <- predictionGridList$x[2] - predictionGridList$x[1]
  if (!is.null(predictionGridList$y)) {
    testY <-
      sd(diff(predictionGridList$x)) / mean(diff(predictionGridList$x))
    if (testY > 1e-8) {
      stop("predictionGridList$y must be equally spaced")
    }
    dy <- predictionGridList$y[2] - predictionGridList$y[1]
  }
}

makeSimulationGrid <- function(predictionGridList, gridRefinement) {
  dx <- predictionGridList$x[2] - predictionGridList$x[1]
  simulationGridList <- list(x = seq(
    min(predictionGridList$x),
    max(predictionGridList$x),
    dx / gridRefinement
  ))
  if (!is.null(predictionGridList$y)) {
    dy <- predictionGridList$y[2] - predictionGridList$y[1]
    simulationGridList$y <-
      seq(min(predictionGridList$y),
          max(predictionGridList$y),
          dy / gridRefinement)
  }
  return( simulationGridList )
}
