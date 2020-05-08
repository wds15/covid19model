library(rstan)
library(data.table)
library(gdata)
library(EnvStats)

setwd("~/git/covid19model")

## command line parsing if any
args <- list(
  stanModelFile= 'base_age_google_mobility_200427',
  seed= 42,
  chain= 1,
  job_tag= 'myfirstrun',
  cntct_by = 5,
  seedAge = 7L
)

## set other args
args$file_stanModel <- file.path('stan-models',paste0(args$stanModelFile,'.stan'))
args$file_serial_interval <- file.path('data','serial_interval.csv')
args$file_covariates <- file.path('data','interventions.csv')
args$file_contact <- file.path('data','contact_UK.rda') # with weekday weekend
args$file_pop_age <- file.path('data','popByAge.csv')
args$file_age_ifr <- file.path('data','ifr_age.csv')
args$file_google_mobility <- file.path('data','google-mobility.csv')
tmp <- Sys.getenv("PBS_JOBID")
args$job_id <- ifelse(tmp!='', tmp, as.character(abs(round(rnorm(1) * 1e6))) )
args$job_dir <- file.path("results",paste0(args$stanModelFile,'-',args$job_tag,'-',args$job_id)) 
args$DEBUG <- TRUE
args$path_to_cmdstan <- filepath("~/git/cmdstan")

## start script
cat(sprintf("Running\n"))
str(args)
set.seed(args$seed)

source(file.path("utils", "read-data.r"))
source(file.path("utils", "process-covariates.r"))

## make job dir
dir.create( args$job_dir )

## save input args
saveRDS( args, file=file.path(args$job_dir, paste0(basename(args$job_dir), '_args.RDS')))

# countries to investigate
countries <- read.csv('data/regions.csv', stringsAsFactors = FALSE)

## Reading all data
d <- read_obs_data(countries)

# Read ifr 
ifr.by.country <- read_ifr_data()

#	read age-specific ifr
ifr.by.age = read.csv(args$file_age_ifr)

#	read contact rates
load(args$file_contact)

#	read age-specific pop counts
pop_by_age = read.csv(args$file_pop_age)

# Read interventions
covariates <- read_interventions(args$file_covariates)

# Read serial interval
serial.interval <- read.csv(args$file_serial_interval)

# read google mobility
google_mobility <- read.csv(args$file_google_mobility)

# Process covariates
forecast <- 0
N2 <- 90 # increase if you need more forecast

processed_data <- process_covariates_age( countries$Regions, 
                                         google_mobility, 
                                         d, 
                                         ifr.by.country, 
                                         ifr.by.age, 
                                         serial.interval, 
                                         pop_by_age, 		 
                                         contact_tab, 
                                         args$seedAge, 
                                         N2, 
                                         forecast)
stan_data = processed_data$stan_data
dates = processed_data$dates
deaths_by_country = processed_data$deaths_by_country
reported_cases = processed_data$reported_cases

## save image before running Stan
save.image( file=file.path(args$job_dir, paste0(basename(args$job_dir), '_stanin.RData')) )

# adapt format for cmdstan
stan_data$covariates_res = array(as.numeric(unlist(stan_data$covariates_res)), dim=c(stan_data$M, stan_data$N2, stan_data$P_RES))
stan_data$covariates_nonres = array(as.numeric(unlist(stan_data$covariates_nonres)), dim=c(stan_data$M, stan_data$N2, stan_data$P_NONRES))
stan_data$cntct_weekdays_mean = array(as.numeric(unlist(stan_data$cntct_weekdays_mean)), dim=c(stan_data$M, stan_data$A, stan_data$A))
stan_data$cntct_weekends_mean = array(as.numeric(unlist(stan_data$cntct_weekdays_mean)), dim=c(stan_data$M, stan_data$A, stan_data$A))

## run Stan
cat('\nRunning Stan... \n')
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
model <- stan_model(args$file_stanModel)

if(rstan){
  if(args$DEBUG) {
    fit <- sampling(model,data=stan_data,iter=10,warmup=5,chains=1,seed=args$seed,verbose=TRUE)
  } else { 
    # uncomment the line below for a full run to replicate results and comment the second line below 
    fit <- sampling(model,data=stan_data,iter=3000, warmup=2000, chains=1,seed=args$seed, thin=1, control = list(adapt_delta = 0.95, max_treedepth = 15))
    #fit = sampling(model,data=stan_data,iter=10,warmup=5,chains=1,seed=args$seed,thin=1,control = list(adapt_delta = 0.95, max_treedepth = 10))
  }
}
if(cmdstan){
  set_cmdstan_path(args$path_to_cmdstan)
  
  if(args$DEBUG) {
    fit <- sampling(model,data=stan_data,iter=10,warmup=5,chains=1,seed=args$seed,verbose=TRUE)
  } else { 
    # uncomment the line below for a full run to replicate results and comment the second line below 
    fit <- sampling(model,data=stan_data,iter=3000, warmup=2000, chains=1,seed=args$seed, thin=1, control = list(adapt_delta = 0.95, max_treedepth = 15))
    #fit = sampling(model,data=stan_data,iter=10,warmup=5,chains=1,seed=args$seed,thin=1,control = list(adapt_delta = 0.95, max_treedepth = 10))
  }
}

save(fit, file=file.path(args$job_dir, paste0(basename(args$job_dir), '_stanout.RData')) )
