# The function conditionalskewness returns the value of one of the conditional skewness measures
# Psi_1( r_1, r_2 | x ) and Psi_2( r_2 | x ) given a fixed covariate value 'x' and numbers 'r_1' and 'r_2'
# with 0 < r_1, r_2 < 1, based on one of the halfspace depth, the spatial depth and the projection depth.
# Its arguments are the following:
#   (1) Response: A n-by-p matrix, where n is the sample size and p (> 1) is the dimension of the response.
#                 For projection depth, dimnesion of the response must be either 2 or 3.
#   (2) Covariate: A n-by-q matrix, where n is the sample size and q is either the length of the equispaced
#                 grid where the functional covariate is recorded, or the dimension of the finite dimensional
#                 covariate. If it is a vector, it is treated as a n-by-1 matrix.
#   (3) fixed.covariate: A vector whose length is equal to q, the number of columns of Covariate.
#   (4) measure.type: 1 for Psi_1( r_1, r_2 | x ) and 2 for Psi_2( r_1 | x ). Default is 1.
#   (5) r_1: A number between 0 and 1, equal to r_1 in the definition of Psi_1( r_1, r_2 | x ).
#                 Default value is 0.1.
#   (6) r_2: A number between 0 and 1, equal to r_2 in the definitions of Psi_1( r_1, r_2 | x ) and
#                 Psi_2( r_2 | x ). Default value is 0.5.
#   (7) depth.type: 1 for halfspace depth, 2 for spatial depth and 3 for projection depth. Default value
#                 is 1.

conditionalskewness = function(Response, Covariate, fixed.covariate, measure.type = 1,
                               r_1 = 0.1, r_2 = 0.5, depth.type = 1)
{
  if (is.vector(Covariate))
    dim(Covariate) = c(length(Covariate), 1)
  
  if (is.vector(Response))
    stop('Response must be multivariate.')
  
  if (r_1 <= 0 || r_1 >= 1)
    stop('r_1 must be within 0 and 1.')
  
  if (r_2 <= 0 || r_2 >= 1)
    stop('r_2 must be within 0 and 1.')
  
  X = Covariate
  Y = Response
  
  sample.size = nrow(Y)
  response.dim = ncol(Y)
  
  nearest.neighbour.size = ceiling((log(sample.size))^2)
  
  Distance.X = sqrt( rowSums((X - matrix(fixed.covariate, nrow = nrow(X),
                                         ncol = length(fixed.covariate), byrow = TRUE))^2) )
  
  h = sort(Distance.X)[nearest.neighbour.size]
  
  ########### Defining function depths for halfspace depth, spatial depth and projection depth
  
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
  
  ########### Defining functions depthmedian for halfspace depth, spatial depth and projection depth
  
  require(MASS)
  require(mrfDepth)
  require(ddalpha)
  require(projmed)
  require(lpSolveAPI)
  
  spatialmedian = function(Data, Weights)
  {
    if (sum(Weights) != 1)
      stop('Sum of weights must be 1.')
    
    ## Checking whether there is only one distinct observation in Data
    
    z = Data[1,]
    Difference = matrix(z, nrow = nrow(Data), ncol = length(z), byrow = TRUE) - Data
    norm_sq_Difference = rowSums(Difference^2)
    if (sum(norm_sq_Difference) < 1e-10)
    {
      median = z
      return(median)
    }
    
    ## Checking whether the weighted median is present in the data itself
    
    n = nrow(Data)
    
    for (i in 1:n)
    {
      X = Data
      x = X[i,]
      U = X - matrix(x, nrow = nrow(X), ncol = ncol(X), byrow = TRUE)
      weights_for_i = Weights
      
      weighted_norms_U = weights_for_i * sqrt(rowSums(U^2))
      all_indices_for_i = 1:n
      J_i = all_indices_for_i[weighted_norms_U < 1e-10]
      J_i_complement = setdiff(all_indices_for_i, J_i)
      J_i = setdiff(J_i, i)
      
      U_new = U[J_i_complement,]
      U_new = U_new / matrix(sqrt(rowSums(U_new^2)), nrow = nrow(U_new), ncol = ncol(U_new), byrow = FALSE)
      weights_proper = weights_for_i[J_i_complement]
      V = colSums(matrix(weights_proper, nrow = nrow(U_new), ncol = ncol(U_new), byrow = FALSE) * U_new)
      
      if (sqrt(sum(V^2)) <= sum(weights_for_i[J_i]))
      {
        median = x
        return(median)
      }
    }
    
    ## Checking whether the data lie on a straight line, and computing the median in that case
    
    x = Data[1,]
    for (i in 2:n)
    {
      y = Data[i,]
      direction_vector = y - x
      if (sum(direction_vector^2) > 0)
        break
    }
    Check_if_linear = mat.or.vec(n, 1)
    s = mat.or.vec(n, 1)
    for (i in 1:n)
    {
      z = Data[i,]
      s_vector = (z - x) / direction_vector
      s_vector = s_vector[-which((z - x) == 0 & direction_vector == 0)]
      if (length(unique(s_vector)) == 1)
      {
        Check_if_linear[i] = 1
        s[i] = unique(s_vector)
      }
    }
    if (sum(Check_if_linear) == n)
    {
      alpha = 0.5
      
      s_sorted_index = order(s, decreasing = FALSE)
      s_sorted = s[s_sorted_index]
      weights_sorted_index = Weights[s_sorted_index]
      cumulative_weights_sorted_index = cumsum(weights_sorted_index)
      index_weighted_median = which(cumulative_weights_sorted_index >= alpha)[1]
      s_weighted_median = s_sorted[index_weighted_median]
      
      median = x + (s_weighted_median * direction_vector)
      return(median)
    }
    
    ## Iteration procedure when the weighted median is not present in the data, or the data is not linear
    
    g_function_weighted = function(X_local, Q_local, weights_local)
    {
      g = sum(weights_local * sqrt(rowSums((matrix(Q_local, nrow = nrow(X_local), ncol = ncol(X_local),
                                                   byrow = TRUE) - X_local )^2))) / sum(weights_local)
      return(g)
    }
    
    X = Data
    
    Check = 0
    
    Q_1 = mat.or.vec(1, ncol(X))
    for (i in 1:ncol(X))
    {
      vector_concerned = X[,i]
      vector_concerned_sorted_index = order(vector_concerned)
      vector_concerned_sorted = vector_concerned[vector_concerned_sorted_index]
      weights_sorted_index = Weights[vector_concerned_sorted_index]
      cumulative_weights_sorted_index = cumsum(weights_sorted_index)
      index_weighted_median = which(cumulative_weights_sorted_index >= 0.5)[1]
      Q_1[i] = vector_concerned_sorted[index_weighted_median]
    }
    Q_best_till_now = Q_1
    g_best_till_now = g_function_weighted(Data, Q_best_till_now, Weights)
    
    Phi = mat.or.vec(ncol(X), ncol(X))
    for (i in 1:n)
    {
      t1 = X[i,] - Q_1
      if (sum(t1^2) > 0)
        Phi = Phi + Weights[i] * (( diag(ncol(X)) - ((t(t1) %*% t1) / sum(t1^2)) ) / sqrt(sum(t1^2)))
    }
    
    if (kappa(Phi) >= 10)
      warning('Bad initial value, final output may not be good estimate.')
    
    Threshold = 0.001
    iteration_number = 1
    maximum_iteration_number = 100
    while(iteration_number <= maximum_iteration_number)
    {
      X = Data
      Delta = mat.or.vec(ncol(X), 1)
      Phi = mat.or.vec(ncol(X), ncol(X))
      for (i in 1:n)
      {
        t1 = X[i,] - Q_1
        if (sum(t1^2) > 0)
        {
          Delta = Delta + Weights[i] * (t1 / sqrt(sum(t1^2)))
          Phi = Phi + Weights[i] * (( diag(ncol(X)) - ((t(t1) %*% t1) / sum(t1^2)) ) / sqrt(sum(t1^2)))
        }
      }
      Q_2 = Q_1 + t(solve(Phi) %*% t(Delta))
      
      difference_relative = sqrt(sum((Q_2 - Q_1)^2)) / max(sqrt(sum(Q_1^2)), sqrt(sum(Q_2^2)))
      
      if (difference_relative < Threshold){
        median = Q_2
        Check = 1
        break
      }else{
        g_at_Q_1 = g_function_weighted(Data, Q_1, Weights)
        g_at_Q_2 = g_function_weighted(Data, Q_2, Weights)
        if (g_at_Q_2 <= g_at_Q_1){
          if (g_at_Q_2 <= g_best_till_now)
          {
            Q_best_till_now = Q_2
            g_best_till_now = g_function_weighted(Data, Q_best_till_now, Weights)
          }
        }else{
          Q_2 = (g_at_Q_2 * Q_1 + g_at_Q_1 * Q_2) / (g_at_Q_1 + g_at_Q_2)
          g_at_Q_2 = g_function_weighted(Data, Q_2, Weights)
          if (g_at_Q_2 <= g_best_till_now)
          {
            Q_best_till_now = Q_2
            g_best_till_now = g_function_weighted(Data, Q_best_till_now, Weights)
          }
        }
        
        Q_1 = Q_2
      }
      
      iteration_number = iteration_number + 1
    }
    
    if (Check == 0)
      median = Q_best_till_now
    
    return(median)
  }
  
  projectionmedian.optimized = function(Data)
  {
    require(projmed)
    
    if (ncol(Data) == 2){
      Ab = projmedintermediatercpp(Data)
      
      A = Ab[, 1:(ncol(Ab)-1)]
      b = Ab[,ncol(Ab)]
      
      require(lpSolveAPI)
      
      lprec = make.lp(nrow = 0, ncol = 3, verbose = "neutral")
      set.objfn(lprec, c(1, 0, 0))
      for (i in 1:nrow(A))
        add.constraint(lprec, c(1, - A[i,]), ">=", - b[i])
      set.bounds(lprec, lower = c(-Inf, -Inf, -Inf), columns = 1:3)
      checksolve = solve(lprec)
      if (checksolve == 0){
        t = get.objective(lprec)
        t.x1.x2 = get.variables(lprec)
        
        projection.median = t.x1.x2[2:3]
      }else{
        stop('ERROR!!! checksolve nonzero!!!')
      }
      
      return(projection.median)
    }else{
      Ab = projmedintermediate3drcpp(Data)
      
      A = Ab[, 1:(ncol(Ab)-1)]
      b = Ab[,ncol(Ab)]
      
      require(lpSolveAPI)
      
      lprec = make.lp(nrow = 0, ncol = 4, verbose = "neutral")
      set.objfn(lprec, c(1, 0, 0, 0))
      for (i in 1:nrow(A))
        add.constraint(lprec, c(1, - A[i,]), ">=", - b[i])
      set.bounds(lprec, lower = c(-Inf, -Inf, -Inf, -Inf), columns = 1:4)
      checksolve = solve(lprec)
      if (checksolve == 0){
        t = get.objective(lprec)
        t.x1.x2.x3 = get.variables(lprec)
        
        projection.median = t.x1.x2.x3[2:4]
      }else{
        stop('ERROR!!! checksolve nonzero!!!')
      }
      
      return(projection.median)
    }
  }
  
  depthmedian = function(Data, depth.type)
  {
    if (depth.type == 1){
      require(mrfDepth)
      depth.median = hdepthmedian(Data)$median
      return(depth.median)
    }else if (depth.type == 2){
      Weights = rep(1, nrow(Data)) / nrow(Data)
      depth.median = spatialmedian(Data, Weights)
      return(depth.median)
    }else if (depth.type == 3){
      require(projmed)
      if (ncol(Data) %in% c(2,3)){
        depth.median = projectionmedian.optimized(Data)
      }else{
        stop('For projection depth, response dimension must be either 2 or 3.')
      }
      return(depth.median)
    }else{
      stop('Enter correct depth.type.')
    }
  }
  
  ################################## Function definition complete ##################################
  
  local.Y.values = Y[Distance.X <= h,]
  
  nbhd.size = nrow(local.Y.values)
  
  r1 = 1 - r_1
  num.points.boxplot.r1 = ceiling( nbhd.size * r_1 )
  num.points.boxplot.r2 = ceiling( nbhd.size * r_2 )
  
  depth.points = depths(local.Y.values, depth.type)
  depth.points.sorted = sort(depth.points, decreasing = TRUE)
  depth.cutoff.r1 = depth.points.sorted[num.points.boxplot.r1]
  depth.cutoff.r2 = depth.points.sorted[num.points.boxplot.r2]
  
  maxdepthset.r1 = local.Y.values[(depth.points >= depth.cutoff.r1),]
  maxdepthset.r2 = local.Y.values[(depth.points >= depth.cutoff.r2),]
  
  trimmed.mean.r1 = colMeans(maxdepthset.r1)
  Spread.r2 = max(dist(maxdepthset.r2))
  
  local.mean = colMeans(local.Y.values)
  
  if (depth.type == 1){
    depth.median = hdepthmedian(local.Y.values)$median
  }else if (depth.type == 2){
    Weights = rep(1, nrow(local.Y.values)) / nrow(local.Y.values)
    depth.median = spatialmedian(local.Y.values, Weights)
  }else{
    depth.median = projectionmedian.optimized(local.Y.values)
  }
  
  Psi1.r1.r2 = sqrt(sum((trimmed.mean.r1 - depth.median)^2)) / Spread.r2
  Psi2.r2 = sqrt(sum((local.mean - depth.median)^2)) / Spread.r2
  
  if (measure.type == 1){
    return(Psi1.r1.r2)
  }else{
    return(Psi2.r2)
  }
}
