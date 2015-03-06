# Setting up a 2D grid of diffusing particles moving one step per cycle and changing direction only after collision with a fixed barrier

# fundamental variables that will be changed depending on the experiment
particleNumber <- 1200
barrierPercentage <- 20
Matrix_X <- 200
Matrix_Y <- 200
particleDistribution <- 1  # 0 for random particle placement within the grid; 1, for particles all placed in the center of the grid
useBarriers <- FALSE # if barriers are used (TRUE) each particles moves along the same direction for each step until a barrier is hit when . The next movment in in a different direction.
                    # If barriers are not used (FALSE) each particle step is in a new random direction. Use low % barriers for ballistic movements but set the matrix size large.

timeSeriesLst <- c(5,10,15,20,30,40,60,80,100, 120,150,200) # each value is the step number for an individual diffusion experiment. The
# particle distribution plot will be for the last experiment. The mean distance for the particles for each experiment is saved in "dataSetD"
# saved as DataB with DataA holding the timeSeriesLst

# Misc. fixed values that are defined by the model or calculated from the experimental variables above
dataSetD <- 1:length(timeSeriesLst)
centerX <- round(Matrix_X/2, digits = 0)
centerY <- round(Matrix_Y/2, digits = 0)
directionPaths <- 8
directionLst <- 1:8
stepVectorX <- c(0,1,1,1,0,-1,-1,-1)
stepVectorY <- c(1,1,0,-1,-1,-1,0,1) # movment in x and y for movment direction value
gridBarrierMatrix <- matrix(FALSE, Matrix_X*Matrix_Y, 1) #actually a Boolean vector of grid locations with 
# TRUE / FALSE indicating the presence or absend of a barrier. Set up initially with no barriers.

# set up obstacles in the grid space listed in barrierLst. Non-moving, single point barriers are placed that change the direction (must change
# the direction) of the moving particles when they are in the same space.
if (useBarriers) {barrierNumber <- round(barrierPercentage/100*Matrix_X*Matrix_Y, digits = 0)
                  barrierVec <- sample(Matrix_X*Matrix_Y, barrierNumber, replace = FALSE, prob = NULL)
                  for (i in barrierVec) {gridBarrierMatrix[i] <- TRUE}}

# Misc functions
locationCalc <- function(x,y){(x-1)*Matrix_Y+y}
twoDDistanceCalc <- function(x,y,x2,y2){sqrt((x-x2)*(x-x2)+(y-y2)*(y-y2))}

# create a matrix of data for the diffusing particles with each row a different particle: 1st & 2nd columns are the X & Y positions, 3rd column is the movment direction, 
# and the 4th is the single number location.
ParticleSet <- matrix(0, particleNumber, 4)

cntr <- 1 
for (stepNumber in timeSeriesLst) {
# set values for the particle coordinates within the grid (0, for random; 1 for grid centered), and the randomly set particle movement direction
if (particleDistribution == 0){
  {ParticleSet[,1] <- sample.int(Matrix_X, size=particleNumber, replace = TRUE, prob = NULL)
   ParticleSet[,2] <- sample.int(Matrix_Y, size=particleNumber, replace = TRUE, prob = NULL)
   ParticleSet[,3] <- sample.int(directionPaths, size = particleNumber, replace = TRUE, prob = NULL)}
   ParticleSet[,4] <- locationCalc(ParticleSet[,1],ParticleSet[,2])}

if (particleDistribution == 1){
  ParticleSet[,1] <- rep(centerX,particleNumber)
  ParticleSet[,2] <- rep(centerY,particleNumber)
  ParticleSet[,3] <- sample.int(directionPaths, size = particleNumber, replace = TRUE, prob = NULL)
  tempVal <- locationCalc(centerX,centerY)
  ParticleSet[,4] <- tempVal}

# ParticleSet # print the particles' starting positions and movement direction
# The main iterations of the particle diffuision experiment
for (j in 1:stepNumber){
    ParticleSet[,1] <- ParticleSet[,1] + stepVectorX[ParticleSet[,3]]
    ParticleSet[,2] <- ParticleSet[,2] + stepVectorY[ParticleSet[,3]]
               
# calculate if particles have escaped the defined grid area and bring them back through the far opposite walls
    tempVec <- which(ParticleSet[,1] > Matrix_X)
    for(pq in tempVec) {ParticleSet[pq,1] <- ParticleSet[pq,1]-Matrix_X}
    tempVec <- which(ParticleSet[,1] < 1)
    for(pq in tempVec) {ParticleSet[pq,1] <- ParticleSet[pq,1]+Matrix_X}
    tempVec <- which(ParticleSet[,2] > Matrix_Y)
    for(pq in tempVec) {ParticleSet[pq,2] <- ParticleSet[pq,2]-Matrix_Y}
    tempVec <- which(ParticleSet[,2] < 1)
    for(pq in tempVec) {ParticleSet[pq,2] <- ParticleSet[pq,2]+Matrix_Y} 
    
# calculate the position number from the new x,y coordinates for particle i
    ParticleSet[,4] <- locationCalc(ParticleSet[,1],ParticleSet[,2])

# change direction of a particle's movement if a barrier is found in the same location as the particle in the last movement step
    if (useBarriers) {for (k in 1:particleNumber) 
        {if (gridBarrierMatrix[ParticleSet[k,4]]) {ParticleSet[k,3] <- sample(directionLst[-ParticleSet[k,3]], size = 1, replace = TRUE, prob = NULL)}}}    
    else {ParticleSet[,3] <- sample.int(directionPaths, size = particleNumber, replace = TRUE, prob = NULL)}

} # end of stepNumber iterations 

# ParticleSet
# pulling all the distances, x & y, for all particles into distX and distY then calculating the mean change in distance from center when 
# particles are initially positioned at the center.

if (particleDistribution == 1){dataSetD[cntr] <- mean(twoDDistanceCalc(ParticleSet[,1], ParticleSet[,2], centerX, centerY))}
cntr <- cntr+1
}

plot(ParticleSet[,1],ParticleSet[,2], xlim = c(0,Matrix_X*1.04), ylim = c(0,Matrix_Y*1.04),main = "Particle Distribution", xlab = "Grid Position X", ylab ="Grid Position Y", asp = Matrix_Y/Matrix_X)
if (particleDistribution == 1){
  write(timeSeriesLst, file = "dataA", ncolumns = 1)
  write(dataSetD, file = "dataB", ncolumns = 1)
  plot(timeSeriesLst, dataSetD, main = "Avg Distance Traveled (vs) Step Number", xlab = "Step Number (X)", ylab = "Average Particle Distance from Center (Y)")
  if (length(dataSetD) > 3){
    fit = nls(dataSetD ~ p1 + p2*timeSeriesLst^p3, start = list(p1 = 0.1, p2 = 1.0, p3 = 0.5))
    newt = data.frame(timeSeriesLst = seq(min(timeSeriesLst),max(timeSeriesLst),len=200))
    lines(newt$timeSeriesLst,predict(fit,newdata=newt))
    dataFitS <- matrix(summary(fit)[[10]],2,3)
    
    xLpos <- (max(timeSeriesLst) + min(timeSeriesLst))/1.5
    yLpos <- (max(dataSetD) + min(dataSetD))/4
    text(xLpos,yLpos, labels = sprintf("power (p) = %5.3f",dataFitS[3]))
    text(xLpos,yLpos*.8, labels = sprintf("amplitude (a) = %5.3f",dataFitS[2]))
    text(xLpos,yLpos*.6, labels = sprintf("offset (o) = %5.3f",dataFitS[1]))
    text(xLpos,yLpos*1.25, labels = sprintf("equation: Y = o + a*X^p"))}
}