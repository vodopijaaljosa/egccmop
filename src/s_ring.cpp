#include <Rcpp.h>
#include <cstdlib>
using namespace Rcpp;

// Utility - Perceptron: Return a number > 0.0 if the customer is taken
double perceptron(NumericVector weights, 
                NumericVector elevators, 
                NumericVector customers,
                int state) {
  
  int n = elevators.size();
  double x = 0;
  double y = 0;
  
  for (int i = 0; i < n; i++) {
    int j = (state + i) % n; 
    x += weights[i] * elevators[j];
    y += weights[i + n] * customers[j];
  }
  
  return x + y;
}

// Utility - Random: Return TRUE with a probability equal to prob
bool new_customer(double prob) {
  return rand() < prob * RAND_MAX;
}

// [[Rcpp::export]]
NumericVector s_ring(NumericVector weights,
                     int no_elevators,
                     double prob,
                     double alpha,
                     double beta,
                     int no_cycles) {
  
  // initialize state, elevators, customers and output vector
  int no_states = weights.size() / 2;
  int no_iter = no_cycles * no_states;
  int state = 0;
  int next_state = 1;
  int total_stops = 0;
  int r = no_states;
  int x = 0;
  int s = 0;
  
  NumericVector skips(no_states);
  NumericVector max_skips(no_states);
  NumericVector elevators(no_states);
  NumericVector customers(no_states);
  NumericVector output(3);
  
  for (int i = 0; i < no_elevators; i++) {
    elevators[i] = 1;
    customers[i] = 0;
    skips[i] = 0;
    max_skips[i] = 0; 
  }
  
  for (int i = no_elevators; i < no_states; i++) {
    elevators[i] = 0;
    customers[i] = 0;
    skips[i] = 0;
    max_skips[i] = 0;
  }
  
  // main loop
  for (int i = 0; i < no_iter; i++) {
    
    if (!customers[state]) {
      if (new_customer(prob)) {
        customers[state] = 1;
        r--;
      }
    }
    
    if (elevators[state]) {
      if (elevators[next_state]) {
        total_stops++;
        if (customers[state]) {
          max_skips[state] = std::max(max_skips[state], skips[state]);
          skips[state] = 0;
          customers[state] = 0;
          r++;
        } 
      } else {
        if (customers[state]) {
          if (perceptron(weights, elevators, customers, state) > 0.0) {
            total_stops++;
            max_skips[state] = std::max(max_skips[state], skips[state]);
            skips[state] = 0;
            customers[state] = 0;
            r++;
          } else {
            skips[state]++;
            elevators[state] = 0;
            elevators[next_state] = 1;
          }
        } else {
          elevators[state] = 0;
          elevators[next_state] = 1;
        } 
      }
    }
    
    x += r; 
    next_state = state;
    state = (state + no_states - 1) % no_states; 
  
  }
  
  for (int i = 0; i < max_skips.size(); i++) {
    max_skips[i] = std::max(max_skips[i], skips[i]);
    if (max_skips[i] > s) {
      s = max_skips[i];
    } 
  }
  
  output[0] = (double)no_states - (double)x / (double)no_iter;
  output[1] = (double)total_stops;
  output[2] = (double)s;
  
  // objective normalization
  output[0] = pow(output[0] / no_states, alpha);
  output[1] = pow(output[1]  * no_states / (no_elevators * no_iter), beta);
  return output;
  
}
