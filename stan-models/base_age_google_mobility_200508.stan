functions {
    matrix[] country_model(matrix cntct_weekends_mean_local,
                         matrix cntct_weekdays_mean_local,
                         real avg_cntct_local,
                         int[] wkend_idx_local,
                         int WKEND_IDX_N_local, 
                         matrix covariates_res_local,
                         matrix covariates_nonres_local,
                         int N0,
                         int N2,
                         int A,
                         int init_A,
                         row_vector rev_serial_interval,
                         row_vector rev_ifr_daysSinceInfection,
                         int SI_CUT,
                         matrix ifr_age,
                         real ifr_country_scale_local,
                         real R0_local,
                         vector beta_res_wkend,
                         vector beta_res_wkday,
                         vector beta_nonres_wkday,
                         real ifr_noise_local,
                         real e_cases_N0_local
                         ) {
   // useful vectors for matrix multiplication
    vector[A] ones_vector_A= rep_vector(1.,A);
    row_vector[A] ones_row_vector_A= rep_row_vector(1.,A);
    // probability of infection given contact in location m
    real rho0;
    // expected deaths by calendar day (1st dim) age (2nd dim) and country (3rd dim), under self-renewal model
    // and expected deaths by calendar day (rows) and country (cols), under self-renewal model
    matrix[N2,A] E_deathsByAge;
    // expected new cases by calendar day, age, and location under self-renewal model
    // and a container to store the precomputed cases by age
    matrix[N2,A] E_casesByAge;
    row_vector[A] tmp_row_vector_A;
    // scaling of contacts after intervention effect on day t in location m
    vector[N2] impact_intv_nonres;
    vector[N2] impact_intv_res;
    
    // Rt for each age band and each location
    matrix[N2,A] RtByAge;
    matrix[N2,A] tmp;
    int comp_time[2] = rep_array(1, 2);
    int comp_time2[2] = rep_array(1, 2);

    // define probability of infection given contact in location m
    rho0 = R0_local /  avg_cntct_local;

    // define multipliers for residential contacts in each location for both weekdays and weekends
    // define multipliers for non-residential contacts in each location for both weekdays and weekends
    // multiply the multipliers with rho0 in each location
    impact_intv_res = exp( covariates_res_local * beta_res_wkday );
    impact_intv_nonres = exp( covariates_nonres_local * beta_nonres_wkday );
    impact_intv_res[ wkend_idx_local[1:WKEND_IDX_N_local]] = exp( covariates_res_local[ wkend_idx_local[1:WKEND_IDX_N_local], :] * beta_res_wkend );
    impact_intv_nonres[ wkend_idx_local[1:WKEND_IDX_N_local]] = rep_vector(0., WKEND_IDX_N_local);
    impact_intv_res *= rho0;
    impact_intv_nonres *= rho0;
    
    // define Rt by age and lcoation
    RtByAge = rep_matrix( ones_row_vector_A*cntct_weekends_mean_local , N2);
    RtByAge .*= rep_matrix( impact_intv_res, A);
    tmp = rep_matrix( ones_row_vector_A*cntct_weekdays_mean_local , N2);
    tmp .*= rep_matrix( impact_intv_nonres, A);
    RtByAge += tmp;
    
    // init expected cases by age and location
    // init expected cases by age and location in first N0 days
    E_casesByAge = rep_matrix( 0., N2, A );
    E_casesByAge[1:N0,init_A] = rep_vector( e_cases_N0_local, N0 );

    // calculate expected cases by age and country under self-renewal model after first N0 days
    // and adjusted for saturation
    for (t in (N0+1):N2)
    {
      comp_time[1] = SI_CUT - t+2;
      comp_time2[1] =t-SI_CUT;
      tmp_row_vector_A = rev_serial_interval[max(comp_time):SI_CUT] * E_casesByAge[max(comp_time2):(t-1),:];
      E_casesByAge[t,:] = tmp_row_vector_A * ( impact_intv_res[t] * cntct_weekends_mean_local + impact_intv_nonres[t] * cntct_weekdays_mean_local );
    }
  
    // calculate expected deaths by age and country  
    E_deathsByAge = rep_matrix( 0., N2, A );
    E_deathsByAge[1,:] = 1e-15 * E_casesByAge[1,:];
    for (t in 2:N2)
    {
      E_deathsByAge[t,:] = rev_ifr_daysSinceInfection[(N2-(t-1)+1):N2 ] * E_casesByAge[1:(t-1),:];
    }
    E_deathsByAge .*= ifr_age;
    E_deathsByAge *= (ifr_country_scale_local * ifr_noise_local);

    return({E_deathsByAge, E_casesByAge, RtByAge });
}

  
  // this is the partial sum function which calculates for
  // a subset of countries the log-lik contribution
  // so it computes for the countries m= start...end the
  // log lik
  real country_lpdf(real[] ifr_noise,
                    int start, int end,
                    int[,] deaths,
                    real phi,
                    int[] EpidemicStart,
                    matrix[] cntct_weekends_mean,
                    matrix[] cntct_weekdays_mean,
                    vector avg_cntct,
                    int[,] wkend_idx,
                    int[] WKEND_IDX_N, 
                    matrix[] covariates_res,
                    matrix[] covariates_nonres,
                    int[] N,
                    int N0,
                    int N2,
                    int A,
                    int init_A,
                    row_vector rev_serial_interval,
                    row_vector rev_ifr_daysSinceInfection,
                    int SI_CUT,
                    matrix ifr_age,
                    real[] ifr_country_scale,
                    vector R0,
                    vector beta_res_wkend,
                    vector beta_res_wkday,
                    vector beta_nonres_wkday,
                    real[] e_cases_N0
                    ) {
                      
    // log-lik of this subset
    vector[A] ones_vector_A= rep_vector(1.,A);
    vector[N2] E_deaths= rep_vector(0., N2);
    real log_lik = 0.0;
    
    for (m in start:end) {
      matrix[N2,A] E_deathsByAge = country_model(
          cntct_weekends_mean[m],
          cntct_weekdays_mean[m],
          avg_cntct[m],
          wkend_idx[:,m],
          WKEND_IDX_N[m], 
          covariates_res[m],
          covariates_nonres[m],
          N0,
          N2,
          A,
          init_A,
          rev_serial_interval,
          rev_ifr_daysSinceInfection,
          SI_CUT,
          ifr_age,
          ifr_country_scale[m],
          R0[m],
          beta_res_wkend,
          beta_res_wkday,
          beta_nonres_wkday,
          ifr_noise[m],
          e_cases_N0[m])[1];
      
     // calculate expected deaths by country
     E_deaths = E_deathsByAge * ones_vector_A;
    
      log_lik += neg_binomial_2_lpmf(deaths[EpidemicStart[m]:N[m], m] |
                                     E_deaths[EpidemicStart[m]:N[m]], phi );
    }

    return(log_lik);
  }
}

data {
  int<lower=1> M; // number of countries
  int<lower=1> N0; // number of initial days for which to estimate infections
  int<lower=1> N[M]; // days of observed data for country m. each entry must be <= N2
  int<lower=1> N2; // days of observed data + # of days to forecast
  int<lower=1> A; // number of age bands
  int<lower=1> SI_CUT; // number of days in serial interval to consider
  int<lower=1> P_RES; // number of predictors for residential contacts
  int<lower=1> P_NONRES; // number of predictors for non-residential contacts
  int WKEND_IDX_N[M]; // number of weekend indices in each location
  //	data
  real pop[M];
  matrix<lower=0, upper=1>[A,M] popByAge; // proportion of age bracket in population in location
  int epidemicStart[M];
  int deaths[N2, M]; // reported deaths -- the rows with i > N contain -1 and should be ignored
  int<lower=0> wkend_idx[N2,M]; //indices of 1:N2 that correspond to weekends in location m
  matrix[N2,P_RES] covariates_res[M]; // predictors for residential contacts
  matrix[N2,P_NONRES]  covariates_nonres[M]; // predictors for non-residential contacts
  //	priors
  matrix[A,A] cntct_weekdays_mean[M]; // mean of prior contact rates between age groups on weekdays
  matrix[A,A] cntct_weekends_mean[M]; // mean of prior contact rates between age groups on weekends
  real<lower=0> ifr_country_scale[M]; // relative probability of death for location, s days after infection, for age band a
  matrix<lower=0>[N2,A] ifr_age; // probability of death for age band a, stored N2 times
  row_vector[N2] rev_ifr_daysSinceInfection; // probability of death s days after infection in reverse order
  row_vector[SI_CUT] rev_serial_interval; // fixed pre-calculated serial interval using empirical data from Neil in reverse order
  int<lower=1, upper=A> init_A; // age band in which initial cases occur in the first N0 days
}

transformed data{
  vector<lower=0>[M] avg_cntct;
  vector[A] ones_vector_A;
  row_vector[A] ones_row_vector_A;
  ones_vector_A= rep_vector(1.,A);
  ones_row_vector_A= rep_row_vector(1.,A);
  
  for( m in 1:M )
  {
    avg_cntct[m] = popByAge[:,m]' * ( cntct_weekdays_mean[m] * ones_vector_A ) * 5./7.;
    avg_cntct[m] += popByAge[:,m]' * ( cntct_weekends_mean[m] * ones_vector_A ) * 2./7.;
  }
}

parameters {
  vector<lower=0>[M] R0; // R0
  //real<lower=0> kappa; // variance parameter for country-specific R0  
  real<lower=0> tau; // prior rate of expected number of cases per day in the first N0 days, for each country
  real<lower=0> e_cases_N0[M]; // expected number of cases per day in the first N0 days, for each country
  vector<lower=0>[P_RES] beta_res_wkend; // regression coefficients for time varying multipliers on residential contacts on weekends
  vector<lower=0>[P_RES] beta_res_wkday; // regression coefficients for time varying multipliers on residential contacts on weekdays
  vector<lower=0>[P_NONRES] beta_nonres_wkday; // regression coefficients for time varying multipliers on non-residential contacts on weekdays
  real<lower=0> phi; // overdispersion parameter for likelihood model
  real<lower=0> ifr_noise[M];
  real<lower=0> kappa;
}

model {
  // priors
  tau ~ exponential(0.03);
  e_cases_N0 ~ exponential(1/tau);
  phi ~ normal(0,5);
  beta_res_wkend ~ normal(0,2);
  beta_res_wkday ~ normal(0,2);
  beta_nonres_wkday ~ normal(0,2);
  ifr_noise ~ normal(1,0.1);
  kappa ~ normal(0,0.5);
  R0 ~ normal(3.28, kappa); // citation: https://academic.oup.com/jtm/article/27/2/taaa021/5735319
  
  // likelihood
  // reduce_sum requires CmdStan >= 2.23 and samples parallel when
  // STAN_THREADS=true is set in make/local
      
  target += reduce_sum_static(
    country_lpdf, ifr_noise, 1,
    deaths,
    phi,
    epidemicStart,
    cntct_weekends_mean,
    cntct_weekdays_mean,
    avg_cntct,
    wkend_idx,
    WKEND_IDX_N, 
    covariates_res,
    covariates_nonres,
    N,
    N0,
    N2,
    A,
    init_A,
    rev_serial_interval,
    rev_ifr_daysSinceInfection,
    SI_CUT,
    ifr_age,
    ifr_country_scale,
    R0,
    beta_res_wkend,
    beta_res_wkday,
    beta_nonres_wkday,
    e_cases_N0);
     
  /*
  // Compatible with Rstan
  target += country_lpdf(
    ifr_noise |
    1, M,
    deaths,
    phi,
    epidemicStart,
    cntct_weekends_mean,
    cntct_weekdays_mean,
    avg_cntct,
    wkend_idx,
    WKEND_IDX_N, 
    covariates_res,
    covariates_nonres,
    N,
    N0,
    N2,
    A,
    init_A,
    rev_serial_interval,
    rev_ifr_daysSinceInfection,
    SI_CUT,
    ifr_age,
    ifr_country_scale,
    R0,
    beta_res_wkend,
    beta_res_wkday,
    beta_nonres_wkday,
    e_cases_N0);
   */
}

generated quantities {
  matrix<lower=0>[N2,M] E_death;
  matrix<lower=0>[N2,A] E_deathsByAge[M];
  matrix<lower=0>[N2,A] E_casesByAge[M];
  matrix<lower=0>[N2,M] Rt;
  matrix<lower=0>[N2,A] RtByAge[M]; 

  for( m in 1:M )
  {
    matrix[N2,A] localA[3] 
      = country_model(
          cntct_weekends_mean[m],
          cntct_weekdays_mean[m],
          avg_cntct[m],
          wkend_idx[:,m],
          WKEND_IDX_N[m], 
          covariates_res[m],
          covariates_nonres[m],
          N0,
          N2,
          A,
          init_A,
          rev_serial_interval,
          rev_ifr_daysSinceInfection,
          SI_CUT,
          ifr_age,
          ifr_country_scale[m],
          R0[m],
          beta_res_wkend,
          beta_res_wkday,
          beta_nonres_wkday,
          ifr_noise[m],
          e_cases_N0[m]);
          
    E_deathsByAge[m] = localA[1];
    E_casesByAge[m] = localA[2];
    RtByAge[m] = localA[3];
    E_death[:,m] = E_deathsByAge[m] * ones_vector_A;
    Rt[:,m] = RtByAge[m] * popByAge[:,m];

  }
}
