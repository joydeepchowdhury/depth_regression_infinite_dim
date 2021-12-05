# The function skewtestpvalue returns the p value for the test of conditional skewness based on random
# samples and one of the halfspace depth, the spatial depth and the projection depth. Its arguments
# are the following:
#   (1) Response: A n-by-p matrix, where n is the sample size and p (> 1) is the dimension of the response.
#                 For projection depth, dimnesion of the response must be either 2 or 3.
#   (2) Covariate: A n-by-q matrix, where n is the sample size and q is either the length of the equispaced
#                 grid where the functional covariate is recorded, or the dimension of the finite dimensional
#                 covariate. If it is a vector, it is treated as a n-by-1 matrix.
#   (3) fixed.covariate: A vector whose length is equal to q, the number of columns of Covariate.
#   (4) depth.type: 1 for halfspace depth, 2 for spatial depth and 3 for projection depth. Default value
#                 is 1.
#   (5) pvalue.type: 1 for the p value of the test based on only Psi_1(0.1, 0.9 | x),
#                    2 for the p value of the test based on only Psi_2(0.5 | x),
#                    3 for both the above p values in a bivariate vector in the above order. See the chapter
#                 (section 5 and section 6) for the descriptions of Psi_1(0.1, 0.9 | x) and Psi_2(0.5 | x).
#   (6) bootstrap.repl.num: Number of random samples on which the test of conditional skewness is based on.
#                 Default value is 1000.

skewtestpvalue = function(Response, Covariate, fixed.covariate, depth.type = 1, pvalue.type = 3,
                           bootstrap.repl.num = 1000)
{
  if (is.vector(Covariate))
    dim(Covariate) = c(length(Covariate), 1)
  
  if (is.vector(Response))
    stop('Response must be multivariate.')
  
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
  
  num.points.boxplot.50pc = ceiling( (nbhd.size * 50) / 100 )
  num.points.boxplot.90pc = ceiling( (nbhd.size * 90) / 100 )
  
  ########### Defining functions depth.allpoints and depthmedian for halfspace depth, spatial depth
  ########### and projection depth
  
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
  
  depth.median = depthmedian(local.Y.values, depth.type)
  depth.points = depth.allpoints(local.Y.values, depth.type)
  
  depth.points.sorted = sort(depth.points, decreasing = TRUE)
  depth.cutoff.50pc = depth.points.sorted[num.points.boxplot.50pc]
  depth.cutoff.90pc = depth.points.sorted[num.points.boxplot.90pc]
  
  maxdepthset.50pc = local.Y.values[(depth.points >= depth.cutoff.50pc),]
  maxdepthset.90pc = local.Y.values[(depth.points >= depth.cutoff.90pc),]
  if (is.vector(maxdepthset.50pc) || is.vector(maxdepthset.90pc))
    stop('PROBLEM!! maxdepthset is a vector.')
  
  trimmed.mean.90pc = colMeans(maxdepthset.90pc)
  Spread.50pc = max(dist(maxdepthset.50pc))
  Spread.90pc = max(dist(maxdepthset.90pc))
  
  local.mean = colMeans(local.Y.values)
  
  statistic.value.90pc = sqrt(sum((trimmed.mean.90pc - depth.median)^2)) / Spread.90pc
  statistic.value.local = sqrt(sum((local.mean - depth.median)^2)) / Spread.50pc
  
  ############################## bootstrap calculation starts ##################################
  
  statistic.value.90pc.bootstrap = mat.or.vec(bootstrap.repl.num, 1)
  statistic.value.local.bootstrap = mat.or.vec(bootstrap.repl.num, 1)
  for (index.num in 1:bootstrap.repl.num)
  {
    local.Y.values.shifted = local.Y.values - matrix(depth.median, nrow = nrow(local.Y.values),
                                                     ncol = ncol(local.Y.values), byrow = TRUE)
    
    double.data = rbind(local.Y.values.shifted, -local.Y.values.shifted)
    sample.indices = sample(1:nrow(double.data), nrow(local.Y.values.shifted))
    
    local.Y.values.shifted.bootstrap.sample = double.data[sample.indices,]
    
    local.Y.values.new = local.Y.values.shifted.bootstrap.sample + 
      matrix(depth.median, nrow = nrow(local.Y.values), ncol = ncol(local.Y.values), byrow = TRUE)
    
    depth.median.new = depth.median
    depth.points.new = depth.allpoints(local.Y.values, depth.type)
    
    depth.points.sorted.new = sort(depth.points.new, decreasing = TRUE)
    depth.cutoff.50pc.new = depth.points.sorted.new[num.points.boxplot.50pc]
    depth.cutoff.90pc.new = depth.points.sorted.new[num.points.boxplot.90pc]
    
    maxdepthset.50pc.new = local.Y.values.new[(depth.points.new >= depth.cutoff.50pc.new),]
    maxdepthset.90pc.new = local.Y.values.new[(depth.points.new >= depth.cutoff.90pc.new),]
    
    trimmed.mean.90pc.new = colMeans(maxdepthset.90pc.new)
    Spread.50pc.new = max(dist(maxdepthset.50pc.new))
    Spread.90pc.new = max(dist(maxdepthset.90pc.new))
    
    local.mean.new = colMeans(local.Y.values.new)
    
    statistic.value.90pc.new = sqrt(sum((trimmed.mean.90pc.new - depth.median.new)^2)) / Spread.90pc.new
    statistic.value.local.new = sqrt(sum((local.mean.new - depth.median.new)^2)) / Spread.50pc.new
    
    statistic.value.90pc.bootstrap[index.num] = statistic.value.90pc.new
    statistic.value.local.bootstrap[index.num] = statistic.value.local.new
  }
  
  ############################### bootstrap calculation ends ###################################
  
  pvalue.skew.90pc = sum(statistic.value.90pc.bootstrap >= statistic.value.90pc) / bootstrap.repl.num
  pvalue.skew.local = sum(statistic.value.local.bootstrap >= statistic.value.local) / bootstrap.repl.num
  
  if (pvalue.type == 1){
    return(pvalue.skew.90pc)
  }else if (pvalue.type == 2){
    return(pvalue.skew.local)
  }else{
    return(c(pvalue.skew.90pc, pvalue.skew.local))
  }
}
