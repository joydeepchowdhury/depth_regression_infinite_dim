# The function conditionaltrimmedmean returns the conditional 100r% trimmed mean given a fixed covariate
# value based on one of the halfspace depth, the spatial depth and the projection depth. Its arguments are
# the following:
#   (1) Response: A n-by-p matrix, where n is the sample size and p (> 1) is the dimension of the response.
#                 For projection depth, dimnesion of the response must be either 2 or 3.
#   (2) Covariate: A n-by-q matrix, where n is the sample size and q is either the length of the equispaced
#                 grid where the functional covariate is recorded, or the dimension of the finite dimensional
#                 covariate. If it is a vector, it is treated as a n-by-1 matrix.
#   (3) fixed.covariate: A vector whose length is equal to q, the number of columns of Covariate.
#   (4) trimming.proportion: A number between 0 and 1, equal to r in the definition of conditional 100r%
#                 trimmed mean. Default value is 0.5.
#   (5) depth.type: 1 for halfspace depth, 2 for spatial depth and 3 for projection depth. Default value
#                 is 1.

conditionaltrimmedmean = function(Response, Covariate, fixed.covariate, trimming.proportion = 0.5,
                                  depth.type = 1)
{
  if (is.vector(Covariate))
    dim(Covariate) = c(length(Covariate), 1)
  
  if (is.vector(Response))
    stop('Response must be multivariate.')
  
  if (trimming.proportion <= 0 || trimming.proportion >= 1)
    stop('trimming.proportion must be within 0 and 1.')
  
  X = Covariate
  Y = Response
  
  sample.size = nrow(Y)
  response.dim = ncol(Y)
  
  nearest.neighbour.size = ceiling((log(sample.size))^2)
  
  Distance.X = sqrt( rowSums((X - matrix(fixed.covariate, nrow = nrow(X),
                                         ncol = length(fixed.covariate), byrow = TRUE))^2) )
  
  h = sort(Distance.X)[nearest.neighbour.size]
  
  local.Y.values = Y[Distance.X <= h,]
  if (is.vector(local.Y.values))
    stop('ERROR: Neighbourhood size is 1.')
  
  nbhd.size = nrow(local.Y.values)
  
  after.trim.prop = 1 - trimming.proportion
  num.points.centralregion.trim.prop = ceiling( nbhd.size * after.trim.prop )
  
  ########### Defining functions depth.allpoints for halfspace depth, spatial depth and projection depth
  
  require(MASS)
  require(mrfDepth)
  require(ddalpha)
  require(projmed)
  
  depth.allpoints = function(Data, depth.type)
  {
    if (depth.type == 1){
      require(mrfDepth)
      depthsallpoints = hdepth(Data)$depthZ
      set.seed(NULL, kind = 'default', normal.kind = 'default')
      return(depthsallpoints)
    }else if (depth.type == 2){
      require(ddalpha)
      depthsallpoints = depth.spatial(Data, Data)
      return(depthsallpoints)
    }else if (depth.type == 3){
      require(projmed)
      if (ncol(Data) == 2){
        depthsallpoints = projdepthrcpp(Data, Data)
      }else if (ncol(Data) == 3){
        depthsallpoints = projdepth3drcpp(Data, Data)
      }else{
        stop('For projection depth, response dimension must be either 2 or 3.')
      }
      return(depthsallpoints)
    }else{
      stop('Enter correct depth.type.')
    }
  }
  
  ################################## Function definition complete ######################################
  
  depth.points = depth.allpoints(local.Y.values, depth.type)
  
  depth.points.sorted = sort(depth.points, decreasing = TRUE)
  depth.cutoff.trim.prop = depth.points.sorted[num.points.centralregion.trim.prop]
  
  maxdepthset.trim.prop = local.Y.values[(depth.points >= depth.cutoff.trim.prop),]
  if (is.vector(maxdepthset.trim.prop))
    stop('PROBLEM!! maxdepthset is a vector.')
  
  trimmed.mean = colMeans(maxdepthset.trim.prop)
  
  return(trimmed.mean)
}
