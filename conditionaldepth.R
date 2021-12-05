# The function conditionaldepth returns the conditional depth of a number of response value given a
# fixed covariate value based on one of the halfspace depth, the spatial depth and the projection depth.
# Its arguments are the following:
#   (1) Response: A n-by-p matrix, where n is the sample size and p (> 1) is the dimension of the response.
#                 For projection depth, dimnesion of the response must be either 2 or 3.
#   (2) Covariate: A n-by-q matrix, where n is the sample size and q is either the length of the equispaced
#                 grid where the functional covariate is recorded, or the dimension of the finite dimensional
#                 covariate. If it is a vector, it is treated as a n-by-1 matrix.
#   (3) fixed.covariate: A vector whose length is equal to q, the number of columns of Covariate.
#   (4) response.value: A m-by-p matrix, whose each row is a response value. If it a vector of length p,
#                 m is taken to be 1.
#   (5) depth.type: 1 for halfspace depth, 2 for spatial depth and 3 for projection depth. Default value
#                 is 1.

conditionaldepth = function(Response, Covariate, fixed.covariate, response.value, depth.type = 1)
{
  if (is.vector(Covariate))
    dim(Covariate) = c(length(Covariate), 1)
  
  if (is.vector(Response))
    stop('Response must be multivariate.')
  
  if (is.vector(response.value))
    dim(response.value) = c(1, length(response.value))
  
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
  
  if (depth.type == 1){
    require(mrfDepth)
    depthvalue = hdepth(local.Y.values, response.value)$depthZ
    set.seed(NULL, kind = 'default', normal.kind = 'default')
  }else if (depth.type == 2){
    require(ddalpha)
    depthvalue = depth.spatial(response.value, local.Y.values)
  }else if (depth.type == 3){
    require(projmed)
    if (ncol(Data) == 2){
      depthvalue = projdepthrcpp(local.Y.values, response.value)
    }else if (ncol(Data) == 3){
      depthvalue = projdepth3drcpp(local.Y.values, response.value)
    }else{
      stop('For projection depth, response dimension must be either 2 or 3.')
    }
  }else{
    stop('Enter correct depth.type.')
  }
  
  return(depthvalue)
}
