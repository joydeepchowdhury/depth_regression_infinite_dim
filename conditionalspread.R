# The function conditionalspread returns the value of the conditional spread measure Delta( r | x ) given
# a fixed covariate value 'x' and a number 'r' with 0 < r < 1, based on one of the halfspace depth, the
# spatial depth and the projection depth. Its arguments are the following:
#   (1) Response: A n-by-p matrix, where n is the sample size and p (> 1) is the dimension of the response.
#                 For projection depth, dimnesion of the response must be either 2 or 3.
#   (2) Covariate: A n-by-q matrix, where n is the sample size and q is either the length of the equispaced
#                 grid where the functional covariate is recorded, or the dimension of the finite dimensional
#                 covariate. If it is a vector, it is treated as a n-by-1 matrix.
#   (3) fixed.covariate: A vector whose length is equal to q, the number of columns of Covariate.
#   (4) r: A number between 0 and 1, equal to r in the definition of Delta( r | x ). Default value is 0.5.
#   (5) depth.type: 1 for halfspace depth, 2 for spatial depth and 3 for projection depth. Default value
#                 is 1.

conditionalspread = function(Response, Covariate, fixed.covariate, r = 0.5, depth.type = 1)
{
  if (is.vector(Covariate))
    dim(Covariate) = c(length(Covariate), 1)
  
  if (is.vector(Response))
    stop('Response must be multivariate.')
  
  if (r <= 0 || r >= 1)
    stop('r must be within 0 and 1.')
  
  X = Covariate
  Y = Response
  
  sample.size = nrow(Y)
  response.dim = ncol(Y)
  
  nearest.neighbour.size = ceiling((log(sample.size))^2)
  
  Distance.X = sqrt( rowSums((X - matrix(fixed.covariate, nrow = nrow(X),
                                         ncol = length(fixed.covariate), byrow = TRUE))^2) )
  
  h = sort(Distance.X)[nearest.neighbour.size]
  
  #### Defining function depths for halfspace depth, spatial depth and projection depth
  
  require(MASS)
  require(mrfDepth)
  require(ddalpha)
  require(projmed)
  
  depths = function(Y, depth.type)
  {
    if (depth.type == 1){
      halfspacedepth = function(Y)
      {
        require(mrfDepth)
        depths = hdepth(Y)$depthZ
        set.seed(NULL, kind = 'default', normal.kind = 'default')
        return(depths)
      }
      return(halfspacedepth(Y))
    }else if (depth.type == 2){
      require(ddalpha)
      return(depth.spatial(Y, Y))
    }else if (depth.type == 3){
      require(projmed)
      if (ncol(Y) == 2){
        return(projdepthrcpp(Y, Y))
      }else if (ncol(Y) == 3){
        return(projdepth3drcpp(Y, Y))
      }else{
        stop('For projection depth, response dimension must be either 2 or 3.')
      }
    }else{
      stop('Enter correct depth.type')
    }
  }
  
  ############################## Function definition complete ##############################
  
  local.Y.values = Y[Distance.X <= h,]
  
  nbhd.size = nrow(local.Y.values)
  
  num.points.boxplot.r = ceiling( nbhd.size * r )
  
  depth.points = depths(local.Y.values, depth.type)
  depth.points.sorted = sort(depth.points, decreasing = TRUE)
  depth.cutoff.r = depth.points.sorted[num.points.boxplot.r]
  
  maxdepthset.r = local.Y.values[(depth.points >= depth.cutoff.r),]
  
  Spread.r = max(dist(maxdepthset.r))
  
  return(Spread.r)
}
