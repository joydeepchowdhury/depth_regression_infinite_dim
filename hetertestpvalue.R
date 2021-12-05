# The function hetertestpvalue returns the p value for the test of heteroscedasticity based on random
# permutations and one of the halfspace depth, the spatial depth and the projection depth. The test is
# based on Delta( 0.5 | x ). See the chapter (section 4 and section 6) for the description of
# Delta( 0.5 | x ). Its arguments are the following:
#   (1) Response: A n-by-p matrix, where n is the sample size and p (> 1) is the dimension of the response.
#                 For projection depth, dimnesion of the response must be either 2 or 3.
#   (2) Covariate: A n-by-q matrix, where n is the sample size and is either the length of the equispaced
#                 grid where the functional covariate is recorded, or the dimension of the finite dimensional
#                 covariate. If it is a vector, it is treated as a n-by-1 matrix.
#   (3) depth.type: 1 for halfspace depth, 2 for spatial depth and 3 for projection depth. Default value
#                 is 1.
#   (4) randpermnum: Number of random permutations on which the test of heteroscedasticity is based on.
#                 Default value is 1000.

hetertestpvalue = function(Response, Covariate, depth.type = 1, randpermnum = 1000)
{
  if (is.vector(Covariate))
    dim(Covariate) = c(length(Covariate), 1)
  
  if (is.vector(Response))
    stop('Response must be multivariate.')
  
  randperm.indx.record = mat.or.vec(randpermnum, sample.size)
  for (i.rand in 1:randpermnum)
    randperm.indx.record[i.rand,] = sample(sample.size, sample.size)
  
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
  
  X.static = Covariate
  Y.static = Response
  
  sample.size = nrow(Y.static)
  
  X.distance = as.matrix( dist(X.static) )
  
  nearest.neighbour.size = ceiling((log(sample.size))^2)
  
  Spread.50pc = mat.or.vec(sample.size, 1)
  for (i in 1:sample.size)
  {
    Y = Y.static
    Distance.X = X.distance[i,]
    
    h = sort(Distance.X)[nearest.neighbour.size]
    
    local.Y.values = Y[Distance.X <= h,]
    if (is.vector(local.Y.values))
      stop('ERROR: Neighbourhood size is 1.')
    
    nbhd.size = nrow(local.Y.values)
    
    num.points.boxplot.50pc = ceiling( (nbhd.size * 50) / 100 )
    
    depth.points = depths(local.Y.values, depth.type)
    depth.points.sorted = sort(depth.points, decreasing = TRUE)
    depth.cutoff.50pc = depth.points.sorted[num.points.boxplot.50pc]
    
    maxdepthset.50pc = local.Y.values[(depth.points >= depth.cutoff.50pc),]
    
    Spread.50pc[i] = max(dist(maxdepthset.50pc))
  }
  
  var.Spread.50pc = var(Spread.50pc)
  
  # Computation for the random permutations
  
  var.Spread.50pc.perm = mat.or.vec(randpermnum, 1)
  for (i.rand in 1:randpermnum)
  {
    randperm.indx = randperm.indx.record[i.rand,]
    Y.perm = Y.static[randperm.indx,]
    
    Spread.50pc.perm = mat.or.vec(sample.size, 1)
    for (i in 1:sample.size)
    {
      Y = Y.perm
      Distance.X = X.distance[i,]
      
      h = sort(Distance.X)[nearest.neighbour.size]
      
      local.Y.values = Y[Distance.X <= h,]
      
      nbhd.size = nrow(local.Y.values)
      
      num.points.boxplot.50pc = ceiling( (nbhd.size * 50) / 100 )
      
      depth.points = depths(local.Y.values, depth.type)
      depth.points.sorted = sort(depth.points, decreasing = TRUE)
      depth.cutoff.50pc = depth.points.sorted[num.points.boxplot.50pc]
      
      maxdepthset.50pc = local.Y.values[(depth.points >= depth.cutoff.50pc),]
      
      Spread.50pc.perm[i] = max(dist(maxdepthset.50pc))
    }
    
    var.Spread.50pc.perm[i.rand] = var(Spread.50pc.perm)
  }
  
  pvalue.50pc.var = mean(var.Spread.50pc.perm >= var.Spread.50pc)
  
  return(pvalue.50pc.var)
}
