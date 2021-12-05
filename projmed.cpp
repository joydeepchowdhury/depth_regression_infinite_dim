#include <Rcpp.h>

using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericMatrix projmedintermediatercpp(NumericMatrix Data) {
  
  int n = Data.nrow(), p = Data.ncol();
  double normsqcurrent, normcurrent;
  NumericVector projectionmed(p);
  
  int i, j, k, l, counter = 0, counter_1, counter_2, counter_3;
  NumericMatrix subtractions((n * (n-1)) / 2, p), directionsubtractions((n * (n-1)) / 2, p);
  NumericMatrix directions((n * (n-1)) / 2, p);
  
  for (i = 0; i < (n - 1); ++i){
    for (k = (i + 1); k < n; ++k){
      normsqcurrent = 0;
      for (j = 0; j < p; ++j){
        subtractions(counter, j) = Data(i, j) - Data(k, j);
        normsqcurrent = normsqcurrent + pow(subtractions(counter, j), 2);
      }
      normcurrent = pow(normsqcurrent, 0.5);
      for (j = 0; j < p; ++j){
        if (normcurrent > 0){
          directionsubtractions(counter, j) = subtractions(counter, j) / normcurrent;
        }else{
          directionsubtractions(counter, j) = subtractions(counter, j);
        }
      }
      counter = counter + 1;
    }
  }
  for (i = 0; i < ((n * (n-1)) / 2); ++i){
    directions(i, 0) = directionsubtractions(i, 1);
    directions(i, 1) = - directionsubtractions(i, 0);
    if (directions(i, 1) < 0){
      for (j = 0; j < p; ++j){
        directions(i, j) = - directions(i, j);
      }
    }
  }
  
  NumericVector angles(((n * (n-1)) / 2) + 2);
  angles[0] = 0;
  for (i = 1; i < ((n * (n-1)) / 2) + 1; ++i){
    angles[i] = acos(directions(i-1, 0));
  }
  angles[((n * (n-1)) / 2) + 1] = M_PI; // M_PI is pi in rcpp
  
  NumericVector anglesorted = sort_unique(angles);
  int num_anglesorted = anglesorted.size();
  
  NumericVector midangles(num_anglesorted - 1);
  NumericMatrix middirections(num_anglesorted - 1, 2);
  for (i = 0; i < num_anglesorted - 1; ++i){
    midangles[i] = (anglesorted[i] + anglesorted[i+1]) / 2;
    
    middirections(i, 0) = cos(midangles[i]);
    middirections(i, 1) = sin(midangles[i]);
  }
  
  int num_midangles = midangles.size();
  
  NumericMatrix projectedData(num_midangles, n);
  IntegerMatrix orderprojectedData(num_midangles, n);
  for (i = 0; i < num_midangles; ++i){
    for (j = 0; j < n; ++j){
      projectedData(i, j) = (middirections(i, 0) * Data(j, 0)) + (middirections(i, 1) * Data(j, 1));
    }
  }
  NumericVector projectedData_i(n);
  int uniqueprojectedData_i;
  for (i = 0; i < num_midangles; ++i){
    for (j = 0; j < n; ++j){
      projectedData_i(j) = projectedData(i, j);
    }
    NumericVector sortedprojectedData_i = sort_unique(projectedData_i);
    uniqueprojectedData_i = sortedprojectedData_i.size();
    counter = 0;
    for (k = 0; k < uniqueprojectedData_i; ++k){
      for (j = 0; j < n; ++j){
        if (sortedprojectedData_i[k] == projectedData_i[j]){
          orderprojectedData(i, counter) = j;
          counter = counter + 1;
        }
      }
    }
  }
  
  IntegerVector indexconsequenceprelim(num_midangles);
  NumericMatrix medianDataprelim(num_midangles, p);
  if (n - 2*(n / 2) == 1){
    int medindproj[num_midangles];
    counter = 0;
    medindproj[counter] = orderprojectedData(0, (n+1)/2 - 1);
    indexconsequenceprelim[counter] = 0;
    for (i = 1; i < num_midangles; ++i){
      if (medindproj[counter] != orderprojectedData(i, (n+1)/2 - 1)){
        counter = counter + 1;
        medindproj[counter] = orderprojectedData(i, (n+1)/2 - 1);
        indexconsequenceprelim[counter] = i;
      }
    }
    for (i = 0; i <= counter; ++i){
      for (j = 0; j < p; ++j){
        medianDataprelim(i, j) = Data(medindproj[i], j);
      }
    }
  }else{
    IntegerMatrix medindproj(num_midangles, 2);
    counter = 0;
    medindproj(counter, 0) = orderprojectedData(0, n/2 - 1);
    medindproj(counter, 1) = orderprojectedData(0, n/2);
    indexconsequenceprelim[counter] = 0;
    for (i = 1; i < num_midangles; ++i){
      if (medindproj(counter, 0) != orderprojectedData(i, n/2 - 1) || medindproj(counter, 1) != orderprojectedData(i, n/2)){
        counter = counter + 1;
        medindproj(counter, 0) = orderprojectedData(i, n/2 - 1);
        medindproj(counter, 1) = orderprojectedData(i, n/2);
        indexconsequenceprelim[counter] = i;
      }
    }
    for (i = 0; i <= counter; ++i){
      for (j = 0; j < p; ++j){
        medianDataprelim(i, j) = (Data(medindproj(i, 0), j) + Data(medindproj(i, 1), j)) / 2;
      }
    }
  }
  
  IntegerVector indexconsequence(counter + 1);
  NumericMatrix medianData(counter + 1, p);
  for (i = 0; i < counter + 1; ++i){
    indexconsequence[i] = indexconsequenceprelim[i];
    for (j = 0; j < p; ++j){
      medianData(i, j) = medianDataprelim(i, j);
    }
  }
  
  NumericVector angles1(counter + 2);
  for (i = 0; i < counter + 1; ++i){
    angles1[i] = anglesorted(indexconsequence[i]);
  }
  angles1[counter + 1] = M_PI; // M_PI is pi in rcpp
  
  //return angles1;
  //////////////////////////////////
  
  int num_angles1 = counter + 2;
  
  NumericVector angles2prelim((int)(pow(n, 4) / 4) + 1);
  counter_1 = 0;
  NumericVector mediancurrent(p), t1(p), t2(p), t3(p);
  NumericMatrix Dataminusmediancurrent(n, p);
  double t2normsq, t3norm;
  for (i = 0; i < num_angles1 - 1; ++i){
    for (j = 0; j < p; ++j){
      mediancurrent[j] = medianData(i, j);
    }
    
    for (k = 0; k < n; ++k){
      for (j = 0; j < p; ++j){
        Dataminusmediancurrent(k, j) = Data(k, j) - mediancurrent[j];
      }
    }
    
    NumericMatrix newdirections((n * (n-1)) / 2, p);
    counter_2 = 0;
    for (j = 0; j < (n-1); ++j){
      for (l = 0; l < p; ++l){
        t1[l] = - Dataminusmediancurrent(j, l);
      }
      for (k = (j+1); k < n; ++k){
        t2normsq = 0;
        for (l = 0; l < p; ++l){
          t2[l] = t1[l] - Dataminusmediancurrent(k, l);
          t2normsq = t2normsq + pow(t2[l], 2);
        }
        if (t2normsq > 0){
          t3[0] = t2[1];
          t3[1] = - t2[0];
        }else{
          t3[0] = t1[0];
          t3[1] = t1[1];
        }
        t3norm = pow((pow(t3[0], 2) + pow(t3[1], 2)), 0.5);
        for (l = 0; l < p; ++l){
          if (t3norm > 0){
            if (t3[1] >= 0){
              newdirections(counter_2, l) = t3[l] / t3norm;
            }else{
              newdirections(counter_2, l) = - t3[l] / t3norm;
            }
          }else{
            newdirections(counter_2, l) = t3[l];
          }
        }
        counter_2 = counter_2 + 1;
      }
    }
    
    NumericVector newanglesprelim((n * (n-1)) / 2);
    for (j = 0; j < ((n * (n-1)) / 2); ++j){
      newanglesprelim[j] = acos(newdirections(j, 0));
    }
    IntegerVector newanglesindices((n * (n-1)) / 2);
    counter_2 = 0;
    for (j = 0; j < ((n * (n-1)) / 2); ++j){
      if (newanglesprelim[j] >= angles1[i] && newanglesprelim[j] < angles1[i+1]){
        newanglesindices[counter_2] = j;
        counter_2 = counter_2 + 1;
      }
    }
    
    if (counter_2 > 0){
      NumericVector newangles(counter_2);
      for (j = 0; j < counter_2; ++j){
        newangles[j] = newanglesprelim[newanglesindices[j]];
      }
      
      NumericVector newanglesortedprelim = sort_unique(newangles);
      int num_newanglesortedprelim = newanglesortedprelim.size();
      NumericVector newanglesorted(num_newanglesortedprelim + 2);
      newanglesorted[0] = angles1[i];
      for (j = 0; j < num_newanglesortedprelim; ++j){
        newanglesorted[j + 1] = newanglesortedprelim[j];
      }
      newanglesorted[num_newanglesortedprelim + 1] = angles1[i+1];
      
      int num_newanglesorted = newanglesorted.size();
      if (num_newanglesorted > 2){
        NumericVector midnewangles(num_newanglesorted - 1);
        NumericMatrix midnewdirections(num_newanglesorted - 1, 2);
        for (j = 0; j < num_newanglesorted - 1; ++j){
          midnewangles[j] = (newanglesorted[j] + newanglesorted[j + 1]) / 2;
          midnewdirections(j, 0) = cos(midnewangles[j]);
          midnewdirections(j, 1) = sin(midnewangles[j]);
        }
        
        NumericMatrix absZ(num_newanglesorted - 1, n);
        for (j = 0; j < num_newanglesorted - 1; ++j){
          for (k = 0; k < n; ++k){
            absZ(j, k) = 0;
            for (l = 0; l < p; ++l){
              absZ(j, k) = absZ(j, k) + midnewdirections(j, l) * Dataminusmediancurrent(k, l);
            }
            absZ(j, k) = fabs(absZ(j, k));
          }
        }
        
        IntegerMatrix orderabsZ(num_newanglesorted - 1, n);
        NumericVector absZ_j(n);
        for (j = 0; j < num_newanglesorted - 1; ++j){
          for (k = 0; k < n; ++k){
            absZ_j[k] = absZ(j, k);
          }
          NumericVector sortedZ_j = sort_unique(absZ_j);
          int uniqueZ_j = sortedZ_j.size();
          counter_3 = 0;
          for (k = 0; k < uniqueZ_j; ++k){
            for (l = 0; l < n; ++l){
              if (sortedZ_j[k] == absZ_j[l]){
                orderabsZ(j, counter_3) = l;
                counter_3 = counter_3 + 1;
              }
            }
          }
        }
        
        IntegerVector indexconsequence1(num_newanglesorted - 1);
        if (n - 2*(n / 2) == 1){
          int medindZ;
          counter_3 = 0;
          medindZ = orderabsZ(0, (n+1)/2 - 1);
          indexconsequence1[counter_3] = 0;
          for (j = 1; j < num_newanglesorted - 1; ++j){
            if (medindZ != orderabsZ(j, (n+1)/2 - 1)){
              counter_3 = counter_3 + 1;
              medindZ = orderabsZ(j, (n+1)/2 - 1);
              indexconsequence1[counter_3] = j;
            }
          }
        }else{
          int medindZ1, medindZ2;
          counter_3 = 0;
          medindZ1 = orderabsZ(0, n/2 - 1);
          medindZ2 = orderabsZ(0, n/2);
          indexconsequence1[counter_3] = 0;
          for (j = 1; j < num_newanglesorted - 1; ++j){
            if (medindZ1 != orderabsZ(j, n/2 - 1) || medindZ2 != orderabsZ(j, n/2)){
              counter_3 = counter_3 + 1;
              medindZ1 = orderabsZ(j, n/2 - 1);
              medindZ2 = orderabsZ(j, n/2);
              indexconsequence1[counter_3] = j;
            }
          }
        }
        
        NumericVector newanglesorted_toadd(counter_3 + 1);
        for (j = 0; j < counter_3 + 1; ++j){
          newanglesorted_toadd[j] = newanglesorted[indexconsequence1[j]];
        }
        
        for (j = 0; j < counter_3 + 1; ++j){
          angles2prelim[counter_1] = newanglesorted_toadd[j];
          counter_1 = counter_1 + 1;
        }
      }
    }
  }
  
  NumericVector angles2(counter_1);
  for (i = 0; i < counter_1; ++i){
    angles2[i] = angles2prelim[i];
  }
  
  int num_angles2 = angles2.size();
  NumericVector allanglesprelim(num_angles1 + num_angles2);
  for (i = 0; i < num_angles1; ++i){
    allanglesprelim[i] = angles1[i];
  }
  for (i = num_angles1; i < num_angles1 + num_angles2; ++i){
    allanglesprelim[i] = angles2[i - num_angles1];
  }
  
  NumericVector allangles = sort_unique(allanglesprelim);
  
  //return allangles;
  
  int num_allangles = allangles.size();
  NumericMatrix alldirections(2 * num_allangles, 2);
  for (i = 0; i < num_allangles; ++i){
    alldirections(i, 0) = cos(allangles[i]);
    alldirections(i, 1) = sin(allangles[i]);
    alldirections(num_allangles + i, 0) = - alldirections(i, 0);
    alldirections(num_allangles + i, 1) = - alldirections(i, 1);
  }
  
  NumericMatrix projectedDataall(2 * num_allangles, n);
  NumericMatrix absprojectedDataallminusmedians(2 * num_allangles, n);
  NumericVector projectedDataallmedians(2 * num_allangles);
  NumericVector MADprojectedDataall(2 * num_allangles);
  NumericVector t4(n);
  for (i = 0; i < 2 * num_allangles; ++i){
    for (j = 0; j < n; ++j){
      projectedDataall(i, j) = 0;
      for (k = 0; k < p; ++k){
        projectedDataall(i, j) = projectedDataall(i, j) + alldirections(i, k) * Data(j, k);
      }
    }
  }
  for (i = 0; i < 2 * num_allangles; ++i){
    for (j = 0; j < n; ++j){
      t4[j] = projectedDataall(i, j);
    }
    projectedDataallmedians[i] = median(t4);
  }
  for (i = 0; i < 2 * num_allangles; ++i){
    for (j = 0; j < n; ++j){
      absprojectedDataallminusmedians(i, j) = fabs(projectedDataall(i, j) - projectedDataallmedians[i]);
    }
  }
  for (i = 0; i < 2 * num_allangles; ++i){
    for (j = 0; j < n; ++j){
      t4[j] = absprojectedDataallminusmedians(i, j);
    }
    MADprojectedDataall[i] = median(t4);
  }
  
  ////////////////////////////////////////////////////////////////////////////////////////
  
  NumericMatrix A(2 * num_allangles, p);
  NumericVector b(2 * num_allangles);
  for (i = 0; i < 2 * num_allangles; ++i){
    if (MADprojectedDataall[i] > 0){
      for (j = 0; j < p; ++j){
        A(i, j) = alldirections(i, j) / MADprojectedDataall[i];
      }
      b[i] = projectedDataallmedians[i] / MADprojectedDataall[i];
    }else{
      for (j = 0; j < p; ++j){
        A(i, j) = 0;
      }
      b[i] = 0;
    }
  }
  NumericMatrix Abreturn(2 * num_allangles, p+1);
  for (i = 0; i < 2 * num_allangles; ++i){
    for (j = 0; j < p; ++j){
      Abreturn(i, j) = A(i, j);
    }
    Abreturn(i, p) = b[i];
  }
  
  return Abreturn;
}

// [[Rcpp::export]]
NumericVector projdepthrcpp(NumericMatrix Data, NumericMatrix Data_1) {
  
  int n = Data.nrow(), p = Data.ncol();
  double normsqcurrent, normcurrent;
  NumericVector projectionmed(p);
  
  int i, j, k, l, counter = 0, counter_1, counter_2, counter_3;
  NumericMatrix subtractions((n * (n-1)) / 2, p), directionsubtractions((n * (n-1)) / 2, p);
  NumericMatrix directions((n * (n-1)) / 2, p);
  
  for (i = 0; i < (n - 1); ++i){
    for (k = (i + 1); k < n; ++k){
      normsqcurrent = 0;
      for (j = 0; j < p; ++j){
        subtractions(counter, j) = Data(i, j) - Data(k, j);
        normsqcurrent = normsqcurrent + pow(subtractions(counter, j), 2);
      }
      normcurrent = pow(normsqcurrent, 0.5);
      for (j = 0; j < p; ++j){
        if (normcurrent > 0){
          directionsubtractions(counter, j) = subtractions(counter, j) / normcurrent;
        }else{
          directionsubtractions(counter, j) = subtractions(counter, j);
        }
      }
      counter = counter + 1;
    }
  }
  for (i = 0; i < ((n * (n-1)) / 2); ++i){
    directions(i, 0) = directionsubtractions(i, 1);
    directions(i, 1) = - directionsubtractions(i, 0);
    if (directions(i, 1) < 0){
      for (j = 0; j < p; ++j){
        directions(i, j) = - directions(i, j);
      }
    }
  }
  
  NumericVector angles(((n * (n-1)) / 2) + 2);
  angles[0] = 0;
  for (i = 1; i < ((n * (n-1)) / 2) + 1; ++i){
    angles[i] = acos(directions(i-1, 0));
  }
  angles[((n * (n-1)) / 2) + 1] = M_PI; // M_PI is pi in rcpp
  
  NumericVector anglesorted = sort_unique(angles);
  int num_anglesorted = anglesorted.size();
  
  NumericVector midangles(num_anglesorted - 1);
  NumericMatrix middirections(num_anglesorted - 1, 2);
  for (i = 0; i < num_anglesorted - 1; ++i){
    midangles[i] = (anglesorted[i] + anglesorted[i+1]) / 2;
    
    middirections(i, 0) = cos(midangles[i]);
    middirections(i, 1) = sin(midangles[i]);
  }
  
  int num_midangles = midangles.size();
  
  NumericMatrix projectedData(num_midangles, n);
  IntegerMatrix orderprojectedData(num_midangles, n);
  for (i = 0; i < num_midangles; ++i){
    for (j = 0; j < n; ++j){
      projectedData(i, j) = (middirections(i, 0) * Data(j, 0)) + (middirections(i, 1) * Data(j, 1));
    }
  }
  NumericVector projectedData_i(n);
  int uniqueprojectedData_i;
  for (i = 0; i < num_midangles; ++i){
    for (j = 0; j < n; ++j){
      projectedData_i(j) = projectedData(i, j);
    }
    NumericVector sortedprojectedData_i = sort_unique(projectedData_i);
    uniqueprojectedData_i = sortedprojectedData_i.size();
    counter = 0;
    for (k = 0; k < uniqueprojectedData_i; ++k){
      for (j = 0; j < n; ++j){
        if (sortedprojectedData_i[k] == projectedData_i[j]){
          orderprojectedData(i, counter) = j;
          counter = counter + 1;
        }
      }
    }
  }
  
  IntegerVector indexconsequenceprelim(num_midangles);
  NumericMatrix medianDataprelim(num_midangles, p);
  if (n - 2*(n / 2) == 1){
    int medindproj[num_midangles];
    counter = 0;
    medindproj[counter] = orderprojectedData(0, (n+1)/2 - 1);
    indexconsequenceprelim[counter] = 0;
    for (i = 1; i < num_midangles; ++i){
      if (medindproj[counter] != orderprojectedData(i, (n+1)/2 - 1)){
        counter = counter + 1;
        medindproj[counter] = orderprojectedData(i, (n+1)/2 - 1);
        indexconsequenceprelim[counter] = i;
      }
    }
    for (i = 0; i <= counter; ++i){
      for (j = 0; j < p; ++j){
        medianDataprelim(i, j) = Data(medindproj[i], j);
      }
    }
  }else{
    IntegerMatrix medindproj(num_midangles, 2);
    counter = 0;
    medindproj(counter, 0) = orderprojectedData(0, n/2 - 1);
    medindproj(counter, 1) = orderprojectedData(0, n/2);
    indexconsequenceprelim[counter] = 0;
    for (i = 1; i < num_midangles; ++i){
      if (medindproj(counter, 0) != orderprojectedData(i, n/2 - 1) || medindproj(counter, 1) != orderprojectedData(i, n/2)){
        counter = counter + 1;
        medindproj(counter, 0) = orderprojectedData(i, n/2 - 1);
        medindproj(counter, 1) = orderprojectedData(i, n/2);
        indexconsequenceprelim[counter] = i;
      }
    }
    for (i = 0; i <= counter; ++i){
      for (j = 0; j < p; ++j){
        medianDataprelim(i, j) = (Data(medindproj(i, 0), j) + Data(medindproj(i, 1), j)) / 2;
      }
    }
  }
  
  IntegerVector indexconsequence(counter + 1);
  NumericMatrix medianData(counter + 1, p);
  for (i = 0; i < counter + 1; ++i){
    indexconsequence[i] = indexconsequenceprelim[i];
    for (j = 0; j < p; ++j){
      medianData(i, j) = medianDataprelim(i, j);
    }
  }
  
  NumericVector angles1(counter + 2);
  for (i = 0; i < counter + 1; ++i){
    angles1[i] = anglesorted(indexconsequence[i]);
  }
  angles1[counter + 1] = M_PI; // M_PI is pi in rcpp
  
  //return angles1;
  //////////////////////////////////
  
  int num_angles1 = counter + 2;
  
  NumericVector angles2prelim((int)(pow(n, 4) / 4) + 1);
  counter_1 = 0;
  NumericVector mediancurrent(p), t1(p), t2(p), t3(p);
  NumericMatrix Dataminusmediancurrent(n, p);
  double t2normsq, t3norm;
  for (i = 0; i < num_angles1 - 1; ++i){
    for (j = 0; j < p; ++j){
      mediancurrent[j] = medianData(i, j);
    }
    
    for (k = 0; k < n; ++k){
      for (j = 0; j < p; ++j){
        Dataminusmediancurrent(k, j) = Data(k, j) - mediancurrent[j];
      }
    }
    
    NumericMatrix newdirections((n * (n-1)) / 2, p);
    counter_2 = 0;
    for (j = 0; j < (n-1); ++j){
      for (l = 0; l < p; ++l){
        t1[l] = - Dataminusmediancurrent(j, l);
      }
      for (k = (j+1); k < n; ++k){
        t2normsq = 0;
        for (l = 0; l < p; ++l){
          t2[l] = t1[l] - Dataminusmediancurrent(k, l);
          t2normsq = t2normsq + pow(t2[l], 2);
        }
        if (t2normsq > 0){
          t3[0] = t2[1];
          t3[1] = - t2[0];
        }else{
          t3[0] = t1[0];
          t3[1] = t1[1];
        }
        t3norm = pow((pow(t3[0], 2) + pow(t3[1], 2)), 0.5);
        for (l = 0; l < p; ++l){
          if (t3norm > 0){
            if (t3[1] >= 0){
              newdirections(counter_2, l) = t3[l] / t3norm;
            }else{
              newdirections(counter_2, l) = - t3[l] / t3norm;
            }
          }else{
            newdirections(counter_2, l) = t3[l];
          }
        }
        counter_2 = counter_2 + 1;
      }
    }
    
    NumericVector newanglesprelim((n * (n-1)) / 2);
    for (j = 0; j < ((n * (n-1)) / 2); ++j){
      newanglesprelim[j] = acos(newdirections(j, 0));
    }
    IntegerVector newanglesindices((n * (n-1)) / 2);
    counter_2 = 0;
    for (j = 0; j < ((n * (n-1)) / 2); ++j){
      if (newanglesprelim[j] >= angles1[i] && newanglesprelim[j] < angles1[i+1]){
        newanglesindices[counter_2] = j;
        counter_2 = counter_2 + 1;
      }
    }
    
    if (counter_2 > 0){
      NumericVector newangles(counter_2);
      for (j = 0; j < counter_2; ++j){
        newangles[j] = newanglesprelim[newanglesindices[j]];
      }
      
      NumericVector newanglesortedprelim = sort_unique(newangles);
      int num_newanglesortedprelim = newanglesortedprelim.size();
      NumericVector newanglesorted(num_newanglesortedprelim + 2);
      newanglesorted[0] = angles1[i];
      for (j = 0; j < num_newanglesortedprelim; ++j){
        newanglesorted[j + 1] = newanglesortedprelim[j];
      }
      newanglesorted[num_newanglesortedprelim + 1] = angles1[i+1];
      
      int num_newanglesorted = newanglesorted.size();
      if (num_newanglesorted > 2){
        NumericVector midnewangles(num_newanglesorted - 1);
        NumericMatrix midnewdirections(num_newanglesorted - 1, 2);
        for (j = 0; j < num_newanglesorted - 1; ++j){
          midnewangles[j] = (newanglesorted[j] + newanglesorted[j + 1]) / 2;
          midnewdirections(j, 0) = cos(midnewangles[j]);
          midnewdirections(j, 1) = sin(midnewangles[j]);
        }
        
        NumericMatrix absZ(num_newanglesorted - 1, n);
        for (j = 0; j < num_newanglesorted - 1; ++j){
          for (k = 0; k < n; ++k){
            absZ(j, k) = 0;
            for (l = 0; l < p; ++l){
              absZ(j, k) = absZ(j, k) + midnewdirections(j, l) * Dataminusmediancurrent(k, l);
            }
            absZ(j, k) = fabs(absZ(j, k));
          }
        }
        
        IntegerMatrix orderabsZ(num_newanglesorted - 1, n);
        NumericVector absZ_j(n);
        for (j = 0; j < num_newanglesorted - 1; ++j){
          for (k = 0; k < n; ++k){
            absZ_j[k] = absZ(j, k);
          }
          NumericVector sortedZ_j = sort_unique(absZ_j);
          int uniqueZ_j = sortedZ_j.size();
          counter_3 = 0;
          for (k = 0; k < uniqueZ_j; ++k){
            for (l = 0; l < n; ++l){
              if (sortedZ_j[k] == absZ_j[l]){
                orderabsZ(j, counter_3) = l;
                counter_3 = counter_3 + 1;
              }
            }
          }
        }
        
        IntegerVector indexconsequence1(num_newanglesorted - 1);
        if (n - 2*(n / 2) == 1){
          int medindZ;
          counter_3 = 0;
          medindZ = orderabsZ(0, (n+1)/2 - 1);
          indexconsequence1[counter_3] = 0;
          for (j = 1; j < num_newanglesorted - 1; ++j){
            if (medindZ != orderabsZ(j, (n+1)/2 - 1)){
              counter_3 = counter_3 + 1;
              medindZ = orderabsZ(j, (n+1)/2 - 1);
              indexconsequence1[counter_3] = j;
            }
          }
        }else{
          int medindZ1, medindZ2;
          counter_3 = 0;
          medindZ1 = orderabsZ(0, n/2 - 1);
          medindZ2 = orderabsZ(0, n/2);
          indexconsequence1[counter_3] = 0;
          for (j = 1; j < num_newanglesorted - 1; ++j){
            if (medindZ1 != orderabsZ(j, n/2 - 1) || medindZ2 != orderabsZ(j, n/2)){
              counter_3 = counter_3 + 1;
              medindZ1 = orderabsZ(j, n/2 - 1);
              medindZ2 = orderabsZ(j, n/2);
              indexconsequence1[counter_3] = j;
            }
          }
        }
        
        NumericVector newanglesorted_toadd(counter_3 + 1);
        for (j = 0; j < counter_3 + 1; ++j){
          newanglesorted_toadd[j] = newanglesorted[indexconsequence1[j]];
        }
        
        for (j = 0; j < counter_3 + 1; ++j){
          angles2prelim[counter_1] = newanglesorted_toadd[j];
          counter_1 = counter_1 + 1;
        }
      }
    }
  }
  
  NumericVector angles2(counter_1);
  for (i = 0; i < counter_1; ++i){
    angles2[i] = angles2prelim[i];
  }
  
  int num_angles2 = angles2.size();
  NumericVector allanglesprelim(num_angles1 + num_angles2);
  for (i = 0; i < num_angles1; ++i){
    allanglesprelim[i] = angles1[i];
  }
  for (i = num_angles1; i < num_angles1 + num_angles2; ++i){
    allanglesprelim[i] = angles2[i - num_angles1];
  }
  
  NumericVector allangles = sort_unique(allanglesprelim);
  
  //return allangles;
  
  int num_allangles = allangles.size();
  NumericMatrix alldirections(2 * num_allangles, 2);
  for (i = 0; i < num_allangles; ++i){
    alldirections(i, 0) = cos(allangles[i]);
    alldirections(i, 1) = sin(allangles[i]);
    alldirections(num_allangles + i, 0) = - alldirections(i, 0);
    alldirections(num_allangles + i, 1) = - alldirections(i, 1);
  }
  
  NumericMatrix projectedDataall(2 * num_allangles, n);
  NumericMatrix absprojectedDataallminusmedians(2 * num_allangles, n);
  NumericVector projectedDataallmedians(2 * num_allangles);
  NumericVector MADprojectedDataall(2 * num_allangles);
  NumericVector t4(n);
  for (i = 0; i < 2 * num_allangles; ++i){
    for (j = 0; j < n; ++j){
      projectedDataall(i, j) = 0;
      for (k = 0; k < p; ++k){
        projectedDataall(i, j) = projectedDataall(i, j) + alldirections(i, k) * Data(j, k);
      }
    }
  }
  for (i = 0; i < 2 * num_allangles; ++i){
    for (j = 0; j < n; ++j){
      t4[j] = projectedDataall(i, j);
    }
    projectedDataallmedians[i] = median(t4);
  }
  for (i = 0; i < 2 * num_allangles; ++i){
    for (j = 0; j < n; ++j){
      absprojectedDataallminusmedians(i, j) = fabs(projectedDataall(i, j) - projectedDataallmedians[i]);
    }
  }
  for (i = 0; i < 2 * num_allangles; ++i){
    for (j = 0; j < n; ++j){
      t4[j] = absprojectedDataallminusmedians(i, j);
    }
    MADprojectedDataall[i] = median(t4);
  }
  
  NumericVector x(p), projected_x(2 * num_allangles);
  double maxoutlyingness_x;
  NumericVector alloutlyingness_x(2 * num_allangles), projectiondepthvalues(Data_1.nrow());
  for (i = 0; i < Data_1.nrow(); ++i){
    for (j = 0; j < p; ++j){
      x[j] = Data_1(i, j);
    }
    for (j = 0; j < 2 * num_allangles; ++j){
      projected_x[j] = 0;
      for (k = 0; k < p; ++k){
        projected_x[j] = projected_x[j] + alldirections(j, k) * x[k];
      }
      if (MADprojectedDataall[j] > 0){
        alloutlyingness_x[j] = fabs(projected_x[j] - projectedDataallmedians[j]) / MADprojectedDataall[j];
      }else{
        alloutlyingness_x[j] = 0;
      }
    }
    maxoutlyingness_x = max(alloutlyingness_x);
    
    projectiondepthvalues[i] = 1 / (1 + maxoutlyingness_x);
  }
  
  return projectiondepthvalues;
}

// [[Rcpp::export]]
NumericMatrix projmedintermediate3drcpp(NumericMatrix Data) {
  
  int n = Data.nrow(), p = Data.ncol();
  
  int i, j, k, counter;
  NumericVector allangles_phi = runif(100, 0.0, 2 * M_PI);
  NumericVector allangles_theta = runif(50, 0.0, M_PI);
  
  NumericMatrix allangles(allangles_phi.size() * allangles_theta.size(), 2);
  counter = 0;
  for (i = 0; i < allangles_phi.size(); ++i){
    for (j = 0; j < allangles_theta.size(); ++j){
      allangles(counter, 0) = allangles_phi[i];
      allangles(counter, 1) = allangles_theta[j];
      counter = counter + 1;
    }
  }
  
  int num_allangles = allangles.nrow();
  int num_alldirections = 2 * num_allangles;
  NumericMatrix alldirections(num_alldirections, 3);
  for (i = 0; i < num_allangles; ++i){
    alldirections(i, 0) = sin(allangles(i, 1)) * cos(allangles(i, 0));
    alldirections(i, 1) = sin(allangles(i, 1)) * sin(allangles(i, 0));
    alldirections(i, 2) = cos(allangles(i, 1));
    alldirections(num_allangles + i, 0) = - alldirections(i, 0);
    alldirections(num_allangles + i, 1) = - alldirections(i, 1);
    alldirections(num_allangles + i, 2) = - alldirections(i, 2);
  }
  
  NumericMatrix projectedDataall(num_alldirections, n);
  NumericMatrix absprojectedDataallminusmedians(num_alldirections, n);
  NumericVector projectedDataallmedians(num_alldirections);
  NumericVector MADprojectedDataall(num_alldirections);
  NumericVector t4(n);
  for (i = 0; i < num_alldirections; ++i){
    for (j = 0; j < n; ++j){
      projectedDataall(i, j) = 0;
      for (k = 0; k < p; ++k){
        projectedDataall(i, j) = projectedDataall(i, j) + (alldirections(i, k) * Data(j, k));
      }
    }
  }
  for (i = 0; i < num_alldirections; ++i){
    for (j = 0; j < n; ++j){
      t4[j] = projectedDataall(i, j);
    }
    projectedDataallmedians[i] = median(t4);
  }
  for (i = 0; i < num_alldirections; ++i){
    for (j = 0; j < n; ++j){
      absprojectedDataallminusmedians(i, j) = fabs(projectedDataall(i, j) - projectedDataallmedians[i]);
    }
  }
  for (i = 0; i < num_alldirections; ++i){
    for (j = 0; j < n; ++j){
      t4[j] = absprojectedDataallminusmedians(i, j);
    }
    MADprojectedDataall[i] = median(t4);
  }
  
  ////////////////////////////////////////////////////////////////////////////////////////
  
  NumericMatrix A(num_alldirections, p);
  NumericVector b(num_alldirections);
  for (i = 0; i < num_alldirections; ++i){
    if (MADprojectedDataall[i] > 0){
      for (j = 0; j < p; ++j){
        A(i, j) = alldirections(i, j) / MADprojectedDataall[i];
      }
      b[i] = projectedDataallmedians[i] / MADprojectedDataall[i];
    }else{
      for (j = 0; j < p; ++j){
        A(i, j) = 0;
      }
      b[i] = 0;
    }
  }
  NumericMatrix Abreturn(num_alldirections, p+1);
  for (i = 0; i < num_alldirections; ++i){
    for (j = 0; j < p; ++j){
      Abreturn(i, j) = A(i, j);
    }
    Abreturn(i, p) = b[i];
  }
  
  return Abreturn;
}

// [[Rcpp::export]]
NumericVector projdepth3drcpp(NumericMatrix Data, NumericMatrix Data_1) {
  
  int n = Data.nrow(), p = Data.ncol();
  
  int i, j, k, counter;
  NumericVector allangles_phi = runif(100, 0.0, 2 * M_PI);
  NumericVector allangles_theta = runif(50, 0.0, M_PI);
  
  NumericMatrix allangles(allangles_phi.size() * allangles_theta.size(), 2);
  counter = 0;
  for (i = 0; i < allangles_phi.size(); ++i){
    for (j = 0; j < allangles_theta.size(); ++j){
      allangles(counter, 0) = allangles_phi[i];
      allangles(counter, 1) = allangles_theta[j];
      counter = counter + 1;
    }
  }
  
  int num_allangles = allangles.nrow();
  int num_alldirections = 2 * num_allangles;
  NumericMatrix alldirections(num_alldirections, 3);
  for (i = 0; i < num_allangles; ++i){
    alldirections(i, 0) = sin(allangles(i, 1)) * cos(allangles(i, 0));
    alldirections(i, 1) = sin(allangles(i, 1)) * sin(allangles(i, 0));
    alldirections(i, 2) = cos(allangles(i, 1));
    alldirections(num_allangles + i, 0) = - alldirections(i, 0);
    alldirections(num_allangles + i, 1) = - alldirections(i, 1);
    alldirections(num_allangles + i, 2) = - alldirections(i, 2);
  }
  
  NumericMatrix projectedDataall(num_alldirections, n);
  NumericMatrix absprojectedDataallminusmedians(num_alldirections, n);
  NumericVector projectedDataallmedians(num_alldirections);
  NumericVector MADprojectedDataall(num_alldirections);
  NumericVector t4(n);
  for (i = 0; i < num_alldirections; ++i){
    for (j = 0; j < n; ++j){
      projectedDataall(i, j) = 0;
      for (k = 0; k < p; ++k){
        projectedDataall(i, j) = projectedDataall(i, j) + (alldirections(i, k) * Data(j, k));
      }
    }
  }
  for (i = 0; i < num_alldirections; ++i){
    for (j = 0; j < n; ++j){
      t4[j] = projectedDataall(i, j);
    }
    projectedDataallmedians[i] = median(t4);
  }
  for (i = 0; i < num_alldirections; ++i){
    for (j = 0; j < n; ++j){
      absprojectedDataallminusmedians(i, j) = fabs(projectedDataall(i, j) - projectedDataallmedians[i]);
    }
  }
  for (i = 0; i < num_alldirections; ++i){
    for (j = 0; j < n; ++j){
      t4[j] = absprojectedDataallminusmedians(i, j);
    }
    MADprojectedDataall[i] = median(t4);
  }
  
  NumericVector x(p), projected_x(num_alldirections);
  double maxoutlyingness_x;
  NumericVector alloutlyingness_x(num_alldirections), projectiondepthvalues(Data_1.nrow());
  for (i = 0; i < Data_1.nrow(); ++i){
    for (j = 0; j < p; ++j){
      x[j] = Data_1(i, j);
    }
    for (j = 0; j < num_alldirections; ++j){
      projected_x[j] = 0;
      for (k = 0; k < p; ++k){
        projected_x[j] = projected_x[j] + (alldirections(j, k) * x[k]);
      }
      if (MADprojectedDataall[j] > 0){
        alloutlyingness_x[j] = fabs(projected_x[j] - projectedDataallmedians[j]) / MADprojectedDataall[j];
      }else{
        alloutlyingness_x[j] = 0;
      }
    }
    maxoutlyingness_x = max(alloutlyingness_x);
    
    projectiondepthvalues[i] = 1 / (1 + maxoutlyingness_x);
  }
  
  return projectiondepthvalues;
}
