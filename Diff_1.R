# Setting up a 2D grid of diffusing particles moving one step per cycle and changing direction only after collision with a fixed barrier

# fundamental variables that will be changed depending on the experiment
particleNumber <- 200
barrierPercentage <- 0
Matrix_X <- 100
Matrix_Y <- 100
particleDistribution <- 1  # 0 for random particle placement within the grid; 1, for particles all placed in the center of the grid
useBarriers <- FALSE # if barriers are used (TRUE) each particles moves along the same direction for each step until a barrier is hit when . The next movment in in a different direction.
                    # If barriers are not (FALSE) used each particle step is in a random direction. Use low % barriers for ballistic movements.

timeSeriesLst <- c(0, 10, 30) # the experiment will step through the listed step sizes

# Misc. fixed values that are defined by the model or calculated from the experimental variables above
dataSetD <- 1:length(timeSeriesLst)
centerX <- round(Matrix_X/2, digits = 0)
centerY <- round(Matrix_Y/2, digits = 0)
directionPaths <- 8
directionLst <- 1:8
stepLst <- list(X = c(0,1,1,1,0,-1,-1,-1), Y = c( 1,1,0,-1,-1,-1,0,1)) # movment in x and y for movment direction value
gridFociLst <- as.list(1:(Matrix_X*Matrix_Y))
particle <- data.frame(X=1, Y=1, D=1, location=1) ## define a particle located at x, y coordinates,  with a movement direction D. defaults 1 for all

if (useBarriers) {barrierNumber <- round(barrierPercentage/100*Matrix_X*Matrix_Y, digits = 0)
                  barrierLst <- sample(gridFociLst, barrierNumber, replace = FALSE, prob = NULL )}
# set up obstacles in the grid space listed in barrierLst. Nonmoving, single point barriers that change the direction (must change direction) of the moving particle.

# Misc functions

locationCalc <- function(part){(part[1]-1)*Matrix_Y+part[2]}
coordCalc <- function(part){xPos <- ceiling(part[4])
                           yPos <- part[4]-(xPos-1)*Matrix_Y
                           list(xPos,yPos)}
twoDDistanceCalc <- function(x,y,x2,y2){sqrt((x-x2)*(x-x2) +(y-y2)*(y-y2))}

# create a set of diffusing particles
ParticleSet <- list(particle)
for (i in 2:particleNumber) {ParticleSet<-c(ParticleSet,list(particle))} 

cntr <- 1 
for (stepNumber in timeSeriesLst) {
# set values for the particle coordinates within the grid (0, for random; 1 for grid centered), and the randomly set particle movement direction
if (particleDistribution == 0){
for (i in 1:particleNumber) {ParticleSet[[i]][1] <- sample.int(Matrix_X, size=1, replace = TRUE, prob = NULL)
                             ParticleSet[[i]][2] <- sample.int(Matrix_Y, size=1, replace = TRUE, prob = NULL)
                             ParticleSet[[i]][3] <- sample.int(directionPaths, size = 1, replace = TRUE, prob = NULL)}}

if (particleDistribution == 1){
for (i in 1:particleNumber) {ParticleSet[[i]][1] <- centerX
                             ParticleSet[[i]][2] <- centerY
                             ParticleSet[[i]][3] <- sample.int(directionPaths, size = 1, replace = TRUE, prob = NULL)}}

for (i in 1:particleNumber) {ParticleSet[[i]][4] <- locationCalc(ParticleSet[[i]])}

# ParticleSet # print the particles' starting positions and movement direction

# The main iterations of the particle diffuision experiment
for (j in 1:stepNumber){ 
for (i in 1:particleNumber) {ParticleSet[[i]][1] <- ParticleSet[[i]][1] + stepLst[[1]][unlist(ParticleSet[[i]][3])]
                            ParticleSet[[i]][2] <- ParticleSet[[i]][2] + stepLst[[2]][unlist(ParticleSet[[i]][3])] 
                      
                      # calculate if particle i has escaped the defined grid area and bring them back through the opposite walls
                      if (ParticleSet[[i]][1] > Matrix_X) {ParticleSet[[i]][1] <- ParticleSet[[i]][1]-Matrix_X} 
                      if (ParticleSet[[i]][1] < 1) {ParticleSet[[i]][1] <- ParticleSet[[i]][1]+Matrix_X}
                      if (ParticleSet[[i]][2] > Matrix_Y) {ParticleSet[[i]][2] <- ParticleSet[[i]][2]-Matrix_Y}
                      if (ParticleSet[[i]][2] < 1) {ParticleSet[[i]][2] <- ParticleSet[[i]][2]+Matrix_Y}
                      
                      # recalculate the position number from the new x,y coordinates for particle i
                      ParticleSet[[i]][4] <- locationCalc(ParticleSet[[i]]) 
                      # change direction of particle i movement if a barrier is in the same location
                      if (useBarriers) {if (any(barrierLst == ParticleSet[[i]]$location)) {tmpLst <- directionLst[-unlist(ParticleSet[[i]]$D)]
                                                              ParticleSet[[i]]$D <- sample(tmpLst, size = 1, replace = TRUE, prob = NULL)}}
                      else {ParticleSet[[i]]$D <- sample(directionLst, size = 1, replace = TRUE, prob = NULL)}
}
}
#ParticleSet

# pulling all the distances, x & y, for all particles into distX and distY then calculating the mean change in distance from center when 
# particles are initially positioned at the center.

distX <- ParticleSet[[1]][1]
distY <- ParticleSet[[1]][2]
for (m in 2:particleNumber) { distX <-c(distX,ParticleSet[[m]][1])
                              distY <-c(distY,ParticleSet[[m]][2])}

if (particleDistribution == 1){distLst <- twoDDistanceCalc(unlist(ParticleSet[[1]][1]), unlist(ParticleSet[[1]][2]), centerX,centerY)
                               for (n in 2:particleNumber) {distLst <- c(distLst, twoDDistanceCalc(unlist(ParticleSet[[n]][1]), 
                                                                                      unlist(ParticleSet[[n]][2]), centerX,centerY))}
                               dataSetD[cntr] <- mean(unlist(distLst))}
cntr <- cntr+1
}

if (particleDistribution == 1){plot(timeSeriesLst, dataSetD, main = "Avg Distance Traveled (vs) Step Number", xlab = "Step Number", ylab = "Average Particle Distance from Center", )}
 plot(distX,distY, xlim = c(0,101), ylim = c(0,101),main = "Particle Distribution", xlab = "Grid Position X", ylab ="Grid Position Y" asp = Matrix_Y/Matrix_X)