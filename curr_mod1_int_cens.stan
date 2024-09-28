// Stan file for modeling, intakes a data vector (processed in the modeling Rmd file) and computes samples of the joint posterior distribution
// using MCMC sampling (specifically NUTS sampling)

data {
  int<lower=0> N; // Number of rows in covariate matrices, equal to the maximum time of first observation over site-years
  int<lower=0> M; // Number of site-year combinations in the data
  int cens[M]; // Indicator variable for censoring for each site-year: 0 for uncensored, 1 for censored. This determines the form of the likelihood
  int obs_time[M]; // Time in day of year that first occurrence or censoring happened
  int int_cens_time[M]; // Previous sampling time before observation time for interval censored cases (i.e. non-left-censored)
  
  // Covariate matrices, where each column is a unique site-year
  matrix[N, M] x_heat;
  matrix[N, M] x_frost;
  matrix[N, M] x_min;
  matrix[N, M] x_rain;
  matrix[N, M] x_spring;
  matrix[N, M] x_photo;
  matrix[N, M] x_phen;
}

parameters {
  real alpha; // intercept
  real Beta_heat;
  real Beta_frost;
  real Beta_min;
  real Beta_photo;
  real Beta_rain;
  real Beta_spring;
  real Beta_phen;
}

model {
  alpha ~ normal(-20, 10); // negative intercept apriori
  Beta_heat ~ normal(0, 0.25);
  Beta_frost ~ normal(0, 0.25);
  Beta_min ~ normal(0, 0.25);
  Beta_rain ~ normal(0, 0.1);
  Beta_spring ~ normal(0, 1);
  Beta_phen ~ normal(0, 0.25);

  // Defining some variables
  real prob;
  real sum_prob;
  real eta;
  real p;
  real prob_complement;
  real before_cens_absences;
  int counter;
  
  for (j in 1:M) {
    prob = 0;
    if (cens[j] == 0) {                  // uncensored case
      
      if (obs_time[j] == 365) {         // absences case
        for (i in 1:365) {
          eta = alpha;                      // calculate the linear predictor
          eta = eta + Beta_heat*x_heat[i,j];
          eta = eta + Beta_rain*x_rain[i,j];
          eta = eta + Beta_min*x_min[i,j];
          eta = eta + Beta_frost*x_frost[i,j];
          eta = eta + Beta_spring*x_spring[i,j]; 
          eta = eta + Beta_phen*x_phen[i,j];
          eta = eta + Beta_photo*x_photo[i,j];

          p = (exp(eta))/(1 + exp(eta));       // apply the inverse link function
          
          if (i == 1) {
            prob = (1-p); // initialization on day 1
          } else {
            prob = prob*(1-p); // cumulative product of absences
          }
        }
        target += log(prob); // add this site-year's contribution to log likelihood
      } else { // now for the presence case
        for (i in 1:int_cens_time[j]) {    // absences up until and including previous sampling time (lower endpoint of censoring interval)
          eta = alpha;                      // calculate the linear predictor
          eta = eta + Beta_heat*x_heat[i,j];
          eta = eta + Beta_rain*x_rain[i,j];
          eta = eta + Beta_min*x_min[i,j];
          eta = eta + Beta_frost*x_frost[i,j];
          eta = eta + Beta_spring*x_spring[i,j]; 
          eta = eta + Beta_phen*x_phen[i,j];
          eta = eta + Beta_photo*x_photo[i,j];

          p = (exp(eta))/(1 + exp(eta));       // apply the inverse link function
          
          if (i == 1) {
            prob = (1-p); // initialize on day 1
          } else {
            prob = prob*(1-p); // cumulative product of absences
          }
        }
        
        before_cens_absences = prob; // store the cumulative product which represents the probability of absences before the censoring interval
        
        // now for the censored days
        counter = int_cens_time[j] + 1; // set up a counter to track the endpoint of the sums over interval censored days
        sum_prob = 0; // set the sum over the days to be zero initially
        
        while (counter <= obs_time[j]) { // contribute probabilities up until the first observation time (the upper endpoint of the censoring interval)
            for (m in int_cens_time[j]:counter) { // sum up the probabilities from the start of the interval to the counter
              eta = alpha;                      // calculate the linear predictor for each day m
              eta = eta + Beta_heat*x_heat[m,j];
              eta = eta + Beta_rain*x_rain[m,j];
              eta = eta + Beta_min*x_min[m,j];
              eta = eta + Beta_frost*x_frost[m,j];
              eta = eta + Beta_spring*x_spring[m,j]; 
              eta = eta + Beta_phen*x_phen[m,j];
              eta = eta + Beta_photo*x_photo[m,j];

              p = (exp(eta))/(1 + exp(eta));       // apply the inverse link function
            
              if (m == counter) { // possible presence contributes p_t
                if (m == int_cens_time[j]) {
                  prob = (p);  // initialization of prob on the first day of censoring
                } else {
                  prob = prob*(p); // cumulative product with a final presence when the time equals the counter variable
                }
              } else { // possible absence contributes (1 - p_t)
                if (m == int_cens_time[j]) {
                  prob = (1 - p); // initialization
                } else {
                  prob = prob*(1 - p); // cumulative product of absences until day "counter" is reached, whereupon a presence is multiplied as above
                }
              }
            }
            sum_prob = sum_prob + prob; // add the probability that the first occurrence is at day "counter" to the sum
            counter = counter + 1; // increment the counter over the censored days
          }
        sum_prob = sum_prob*(before_cens_absences); // multiply the prob. of absences before censoring and all possible presences during censoring interval
        target += log(sum_prob); // add log likelihood of the censored interval
      } // end interval censored presence case (i.e. not left censored)
      
    } else { // else do left-censored case
      for (i in 1:obs_time[j]) {
        eta = alpha;                      // calculate the linear predictor
        eta = eta + Beta_heat*x_heat[i,j];
        eta = eta + Beta_rain*x_rain[i,j];
        eta = eta + Beta_min*x_min[i,j];
        eta = eta + Beta_frost*x_frost[i,j];
        eta = eta + Beta_spring*x_spring[i,j]; 
        eta = eta + Beta_phen*x_phen[i,j];
        eta = eta + Beta_photo*x_photo[i,j];

        p = (exp(eta))/(1 + exp(eta));       // apply the inverse link function
        
        if (i == 1) {
          prob = (1 - p); // initialization on day 1
        } else {
          prob = prob*(1 - p); // cumulative product of absences
        }
      }
      prob_complement = 1 - prob; // probability that the occurrence time is greater than the censoring time
      target += log(prob_complement);
    }
  }
}
