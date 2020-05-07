library(rstan)
library(data.table)
library(lubridate)
library(gdata)
library(dplyr)
library(tidyr)
library(EnvStats)
library(scales)
library(stringr)
library(abind)

process_covariates <- function(countries, interventions, d, ifr.by.country,N2){
  serial.interval = read.csv("data/serial_interval.csv")
  # Pads serial interval with 0 if N2 is greater than the length of the serial
  # interval array
  if (N2 > length(serial.interval$fit)) {
    pad_serial.interval <- data.frame(
      "X"=(length(serial.interval$fit)+1):N2,
      "fit"=rep(1e-17, max(N2-length(serial.interval$fit), 0 ))
    )
    serial.interval = rbind(serial.interval, pad_serial.interval)
  }
  # various distributions required for modeling
  mean1 <- 5.1; cv1 <- 0.86; # infection to onset
  mean2 <- 17.8; cv2 <- 0.45 # onset to death
  x1 <- rgammaAlt(1e6,mean1,cv1) # infection-to-onset distribution
  x2 <- rgammaAlt(1e6,mean2,cv2) # onset-to-death distribution
  
  ecdf.saved <- ecdf(x1+x2)
  forecast <- 0
  dates <- list()
  reported_cases <- list()
  stan_data <- list(M=length(countries$Regions),N=NULL,deaths=NULL,f=NULL,
                   N0=6,cases=NULL,SI=serial.interval$fit[1:N2],features=NULL,
                   EpidemicStart = NULL, pop = NULL)
  deaths_by_country <- list()
  covariate_list <- list()
  k=1
  # going over each region
  for (Country in countries$Regions){
    IFR <- ifr.by.country$ifr[ifr.by.country$country == Country]
    region_pop <- ifr.by.country[ifr.by.country$country==Country,]
    region <- d[d$Country==Country,]
    region$DateRep <-region$DateRep
    region <-region[order(as.Date(region$DateRep)),]  # ensure date ordering
    
    # padding in raw data backwards ex. portugal
    date_min <- dmy('31/12/2019') 
    if (region$DateRep[1] > date_min){
      print(paste(Country,'In padding ECDC data'))
      pad_days <-region$DateRep[1] - date_min
      pad_dates <- date_min + days(1:pad_days[[1]]-1)
      padded_data <- data.frame("Country" = rep(Country, pad_days),
                                "DateRep" = pad_dates,
                                "Cases" = as.integer(rep(0, pad_days)),
                                "Deaths" = as.integer(rep(0, pad_days)),
                                stringsAsFactors=F)
      
     region <- bind_rows(padded_data,region)
    }
    index = which(region$Cases>0)[1]
    index1 = which(cumsum(region$Deaths)>=10)[1] # also 5
    index2 = index1-30
    
    print(sprintf("First non-zero cases is on day %d, and 30 days before 10 deaths is day %d",index,index2))
    region=region[index2:nrow(region),]
    stan_data$EpidemicStart = c(stan_data$EpidemicStart,index1+1-index2)
    stan_data$pop = c(stan_data$pop,region_pop$popt)
    # NPI interventionss are being used
    interventions_region <- interventions[interventions$Country == Country, c(2,3,4,5,6)] # school, self-isolation, public, lockdown, social-distancing
    for (ii in 1:ncol(interventions_region)) {
      covariate = names(interventions_region)[ii]
      region[covariate] <- (region$DateRep >= interventions_region[1,covariate])*1  # should this be > or >=?
    }
    
    dates[[Country]] =region$DateRep
    # hazard estimation
    N = length(region$Cases)
    print(sprintf("%s has %d days of data",Country,N))
    forecast = N2 - N
    if(forecast < 0) {
      print(sprintf("%s: %d", Country, N))
      print("ERROR!!!! increasing N2")
      N2 = N
      forecast = N2 - N
    }
    
    # IFR is the overall probability of dying given infection
    convolution = function(u) (IFR * ecdf.saved(u))
    
    f = rep(0,N2) # f is the probability of dying on day i given infection
    f[1] = (convolution(1.5) - convolution(0))
    for(i in 2:N2) {
      f[i] = (convolution(i+.5) - convolution(i-.5)) 
    }
    reported_cases[[Country]] = as.vector(as.numeric(region$Cases))
    deaths=c(as.vector(as.numeric(region$Deaths)),rep(-1,forecast))
    cases=c(as.vector(as.numeric(region$Cases)),rep(-1,forecast))
    deaths_by_country[[Country]] = as.vector(as.numeric(region$Deaths))
    region_intervention <- as.data.frame(region[, colnames(interventions_region)])
    region_intervention[N:(N+forecast),] <- region_intervention[N,]
    school =region_intervention[,1]
    selfIsolation =region_intervention[,2]
    publicEvents = region_intervention[,3]
    lockdown = region_intervention[,4]
    socialDistancing = region_intervention[,5]
    firstIntervention = 1*((school+ selfIsolation+ publicEvents+ lockdown + socialDistancing) >= 1)
    ## append data
    stan_data$N = c(stan_data$N,N)
    # stan_data$x = cbind(stan_data$x,x)
    stan_data$f = cbind(stan_data$f,f)
    stan_data$deaths = cbind(stan_data$deaths,deaths)
    stan_data$cases = cbind(stan_data$cases,cases)
    
    stan_data$N2=N2
    stan_data$x=1:N2
    if(length(stan_data$N) == 1) {
      stan_data$N = as.array(stan_data$N)
    }
    df_features = data.frame('school' = school, 'selfIsolation' = selfIsolation, 'publicEvents' = publicEvents,
                             'firstIntervention' = firstIntervention, 'lockdown' = lockdown, 'socialDistancing' = socialDistancing)
    features <- as.matrix(df_features)
    covariate_list[[k]] <- features
    k <- k+1
  }
  stan_data$P = dim(features)[2]
  stan_data$X = array(NA, dim = c(stan_data$M , stan_data$N2 ,stan_data$P ))
  for (i in 1:stan_data$M){
    stan_data$X[i,,] = covariate_list[[i]] 
  }
  return(list("stan_data" = stan_data, "dates" = dates, "reported_cases"=reported_cases, "deaths_by_country" = deaths_by_country))
}

process_covariates_age <- function(countries, google_mobility, d, ifr.by.country, ifr.by.age, serial.interval, pop_by_age, contact_tab, seedAge, N2, forecast)
{
  # various distributions required for modeling
  set.seed(5678) # make reproducible 
  mean1 <- 5.1; cv1 <- 0.86; # infection to onset
  mean2 <- 18.8; cv2 <- 0.45 # onset to death
  x1 <- rgammaAlt(1e6,mean1,cv1) # infection-to-onset distribution
  x2 <- rgammaAlt(1e6,mean2,cv2) # onset-to-death distribution
  ecdf_death <- ecdf(x1+x2)
  
  #	disrectise to days since infection
  ifr_by_dsi <- vector('numeric',N2)
  ifr_by_dsi[1] <- ecdf_death(1.5) - ecdf_death(0)
  for(s in 2:N2)
  {
    ifr_by_dsi[s] <- ecdf_death(s+.5) - ecdf_death(s-.5)
  }
  
  # serial intervention cut
  serial.interval.cut <- length( which(serial.interval$fit * 5e6 >= 1 ) )
  
  # contact bands
  cntct_bands <- nrow(contact_tab[[1]][[1]]) 
  
  # going over each region
  dates <- list()
  reported_cases <- list()
  deaths_by_country <- list()
  stan_data <- list(
    #	constants
    M = length(countries), 									# int<lower=1> M		number of countries
    N0 = 6L, 												# int<lower=1> N0		N0 = 6 to make it consistent with Rayleigh
    N = vector('integer',length(countries)),				# int<lower=1> N[M]		days of observed data for country m. each entry must be <= N2
    N2 = as.integer(N2), 									# int<lower=1> N2		days of observed data + # of days to forecast
    A = as.integer(cntct_bands), 							# int<lower=1> A		number of age bands
    SI_CUT = serial.interval.cut, 							# int<lower=1> SI_CUT	number of days in serial interval to consider
    P_RES = 1L,												#int<lower=1> P_RES; 	number of predictors for residential contacts		
    P_NONRES = 4L,											# int<lower=1> P_NONRES number of predictors for non-residential contacts
    WKEND_IDX_N= vector('integer',length(countries)),		# int WKEND_IDX_N[M];	number of weekend indices in each location
    #	data
    pop = vector('numeric',length(countries)),				# real pop[M];									population counts
    popByAge = matrix(NA_integer_, as.integer(cntct_bands), length(countries)),		# matrix<lower=0, upper=1>[A,M] popByAge	proportion of age bracket in population in location
    epidemicStart = vector('integer',length(countries)), 	# int epidemicStart[M];		
    deaths = matrix(NA_integer_, N2, length(countries)), 	# int deaths[N2, M];	 						reported deaths
    wkend_idx = matrix(NA_integer_, N2, length(countries)),	# int<lower=1> wkend_idx[N2,M]; 				indices of 1:N2 that correspond to weekends in location m
    covariates_nonres = vector('list', length(countries)),	# matrix[N2,P_NONRES]  covariates_nonres[M]		predictors for non-residential contacts
    covariates_res = vector('list', length(countries)),		# matrix[N2,P_RES] covariates_res[M]			predictors for residential contacts
    #	priors
    cntct_weekdays_mean = vector('list', length(countries)), 								# matrix[A,A] cntct_mean[M]						mean of prior contact rates between age groups
    cntct_weekends_mean = vector('list', length(countries)), 
    ifr_country_scale = vector('numeric', length(countries)),						# real<lower=0> ifr_country_scale[M]; 			relative probability of death for location, s days after infection, for age band a 
    ifr_age = matrix(data=ifr.by.age$ifr_mean, ncol=cntct_bands, nrow=N2, byrow=TRUE), 	# matrix<lower=0>[N2,A] ifr_age; 				probability of death for age band a, stored N2 times 
    hyperpara_ifr_age = ifr.by.age[, c("alpha", "beta")],
    rev_ifr_daysSinceInfection = rev( ifr_by_dsi ), 									# row_vector[N2] rev_ifr_daysSinceInfection;	probability of death s days after infection in reverse order
    rev_serial_interval = rev( serial.interval$fit[1:serial.interval.cut] ), 			# row_vector[SI_CUT] rev_serial_interval;		fixed pre-calculated serial interval using empirical data from Neil in reverse order
    init_A = ifelse(cntct_bands==101, 30L, ifelse(cntct_bands==20, as.numeric(seedAge), NA_integer_))	# int<lower=1, upper=A> init_A					age band in which initial cases occur in the first N0 days 
  ) 		
  
  cat('\nProcessing country data... \n')
  for(m in seq_along(countries)) 
  {
    #m <- 1	
    Country <- countries[m]
    
    #	create padded data of cases and deaths
    # start counting 30 days prior to reaching 10 deaths
    d1 <- d[d$Country==Country,]
    d1$date <- as.Date(d1$DateRep,format='%Y-%m-%d')
    d1$t <- decimal_date(d1$date) 
    d1 <- d1[order(d1$t),]
    date_min <- dmy('31/12/2019') 
    if (as.Date(d1$DateRep[1], format='%Y-%m-%d') > as.Date(date_min, format='%Y-%m-%d')){
      print(paste(Country,'In padding'))
      pad_days <- as.Date(d1$DateRep[1], format='%Y-%m-%d') - date_min
      pad_dates <- date_min + days(1:pad_days[[1]]-1)
      padded_data <- data.frame("Country" = rep(Country, pad_days),
                                "DateRep" = format(pad_dates, '%Y-%m-%d'),
                                "t" = decimal_date(as.Date(pad_dates,format='%Y-%m-%d')),
                                "date" = as.Date(pad_dates,format='%Y-%m-%d'),
                                "Cases" = as.integer(rep(0, pad_days)),
                                "Deaths" = as.integer(rep(0, pad_days)),
                                stringsAsFactors=F)
      padded_data$DateRep = as.Date(padded_data$DateRep, format='%Y-%m-%d')
      d1 <- bind_rows(padded_data, d1)
    }
    index <- which(d1$Cases>0)[1]
    index1 <- which(cumsum(d1$Deaths)>=10)[1] # also 5
    index2 <- index1-30L
    print(sprintf("First non-zero cases is on day %d, and 30 days before 10 deaths is day %d",index,index2))
    d1 <- d1[index2:nrow(d1),]
    
    #	add data on mobility indices
    google_mobility$date = as.Date(google_mobility$date,format='%Y-%m-%d')
    google_mobility$ISO_name = as.character(google_mobility$ISO_name)
    tmp <- google_mobility %>% 
      select(-c(ISO_code,county)) %>% 
      rename(Country:=ISO_name)
    d1 <- d1 %>% 
      group_by(Country,date) %>% 
      left_join(tmp) %>% 
      ungroup()
    for(x in c("grocery.pharmacy", "parks", "residential", "retail.recreation","transitstations","workplace"))
    {
      d1[[x]]
      z <- which(is.na(d1[[x]]))			
      z2 <- which(diff(z)!=1)
      if(!length(z2)) 
        z2 <- length(z)
      d1[[x]][1:z2] <- 0
      z <- which(is.na(d1[[x]]))
      if(length(z)) 
        d1[[x]][z] <- d1[[x]][z[1]-1] 
    }
    
    #	determine the index of weekend dates
    wkend_idx <- which(as.integer(format(d1$date, "%u"))>5)
    
    # store dates for data from this country
    dates[[Country]] <- d1$date  
    
    # hazard estimation
    N <- nrow(d1)
    print(sprintf("%s has %d days of data",Country,N))
    forecast <- N2 - N
    if(forecast < 0) {
      print(sprintf("%s: %d", Country, N))
      print("ERROR!!!! increasing N2")
      N2 <- N
      forecast <- N2 - N
    }
    
    #	calculate probability of deaths in age a after s days of infection for this country 
    # to obtain ifr.by.age.country, I scale here ifr.by.age with the ratio ifr.by.country/ifr.china
    # where ifr.china is 0.657% as in Verity 2020
    # TODO this will likely need to be revised
    ifr_china <- 0.00657
    ifr_country_scale <- ifr.by.country$ifr[ifr.by.country$country == Country] / ifr_china  
    
    reported_cases[[Country]] <- as.vector(as.numeric(d1$Cases))
    deaths <- c(as.vector(as.numeric(d1$Deaths)),rep(-1,forecast))
    cases <- c(as.vector(as.numeric(d1$Cases)),rep(-1,forecast))
    deaths_by_country[[Country]] <- as.vector(as.numeric(d1$Deaths))
    covariates <- d1 %>% select(-c(DateRep, Cases, Deaths, Country, date, t))
    covariates[N:(N+forecast),] <- covariates[N,]
    
    #	add constants
    stan_data$N[m] <- N
    stan_data$WKEND_IDX_N[m] <- length(wkend_idx)
    stan_data$wkend_idx[,m] <- c(wkend_idx, rep(0,N2-length(wkend_idx)))
    
    # add data  
    stan_data$pop[m] <- ifr.by.country$popt[ifr.by.country$country==Country]
    stan_data$popByAge[,m] <- pop_by_age[,Country]
    stan_data$epidemicStart[m] <- index1+1L-index2	#TODO confirm, this is always 31?
    stan_data$deaths[,m] <- deaths
    stan_data$covariates_res[[m]] <- as.matrix( covariates %>% select(residential) )
    stan_data$covariates_nonres[[m]] <- as.matrix( covariates %>% select(parks, retail.recreation, transitstations, workplace) )		
    
    # add priors				
    stan_data$cntct_weekends_mean[[m]] <- as.matrix(contact_tab[["United_Kingdom"]][["weekend"]])
    stan_data$cntct_weekdays_mean[[m]] <- as.matrix(contact_tab[["United_Kingdom"]][["weekday"]])
    
    stan_data$ifr_country_scale[m] <- ifr_country_scale
    
    if(length(stan_data$N) == 1) 
    {
      stan_data$N <- as.array(stan_data$N)
    }
  }
  
  # transform in array single vector 
  if(stan_data$M == 1){
    stan_data$pop = as.array(stan_data$pop)
    stan_data$epidemicStart = as.array(stan_data$epidemicStart)
    stan_data$ifr_country_scale = as.array(stan_data$ifr_country_scale)
  }
  
  return(list("stan_data" = stan_data, "dates" = dates, "reported_cases"=reported_cases, "deaths_by_country" = deaths_by_country))
} 
