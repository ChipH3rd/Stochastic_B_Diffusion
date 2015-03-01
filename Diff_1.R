# Setting up a 2D grid of diffusing particles moving one step per cycle and changing direction only after collision with a fixed barrier

# fundamental variables that will be changed depending on the experiment
particleNumber <- 4
barrierPercentage <- 1
stepNumber <- 5
Matrix_X <- 100
Matrix_Y <- 100
particleDistribution <- 1  # 0 for random particle placement within the grid, 1, for particles all placed in the center of the grid

# Misc. fixed values that are defined by the model or calculated from the experimental variables above
directionPaths <- 8
directionLst <- 1:8
stepLst <- list(X = c(0,1,1,1,0,-1,-1,-1), Y = c( 1,1,0,-1,-1,-1,0,1)) # movment in x and y for movment direction value
gridFociLst <- as.list(1:(Matrix_X*Matrix_Y))
barrierNumber <- round(barrierPercentage/100*Matrix_X*Matrix_Y, digits = 0)
particle <- data.frame(X=1, Y=1, D=1, location=1) ## define a particle located at x, y coordinates,  with a movement direction D. defaults 1 for all

# Misc functions
locationCalc <- function(part){(part[1]-1)*Matrix_Y+part[2]}
coordCalc <- function(part){xPos <- ceiling(part[4])
                           yPos <- part[4]-(xPos-1)*Matrix_Y
                           list(xPos,yPos)}

# create a set of diffusing particles
ParticleSet <- list(particle)
for (i in 2:particleNumber) {ParticleSet<-c(ParticleSet,list(particle))} 

# set values for the particle coordinates within the grid (0, for random; 1 for grid centered), and the randomly set particle movement direction
if (particleDistribution == 0){
for (i in 1:particleNumber) {ParticleSet[[i]][1] <- sample.int(Matrix_X, size=1, replace = TRUE, prob = NULL)
                             ParticleSet[[i]][2] <- sample.int(Matrix_Y, size=1, replace = TRUE, prob = NULL)
                             ParticleSet[[i]][3] <- sample.int(directionPaths, size = 1, replace = TRUE, prob = NULL)}}

if (particleDistribution == 1){
for (i in 1:particleNumber) {ParticleSet[[i]][1] <- round(Matrix_X/2, digits = 0)
                             ParticleSet[[i]][2] <- round(Matrix_Y/2, digits = 0)
                             ParticleSet[[i]][3] <- sample.int(directionPaths, size = 1, replace = TRUE, prob = NULL)}}

for (i in 1:particleNumber) {ParticleSet[[i]][4] <- locationCalc(ParticleSet[[i]])}

# set up obstacles in the grid space listed in barrierLst. Nonmoving, single point barriers that change the direction (must change direction) of the moving particle.
barrierLst <- sample(gridFociLst, barrierNumber, replace = FALSE, prob = NULL )

ParticleSet # print the particles' starting positions and movement direction

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
                      if (any(barrierLst == ParticleSet[[i]]$location)) {tmpLst <- directionLst[-unlist(ParticleSet[[i]]$D)]
                                                                          ParticleSet[[i]]$D <- sample(tmpLst, size = 1, replace = TRUE, prob = NULL)}

}
}

ParticleSet