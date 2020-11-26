# set the working directory
setwd("C:\\Users/TO/Desktop")

# load the dataset
data=read.csv("ten-roberto-wessels.csv",sep=";",header=TRUE) 

# Install all necessary packages
install.packages("ggplot2")
install.packages("dygraphs")
install.packages("tidyverse")
install.packages("lubridate")
install.packages("zoo")
install.packages("xts")
install.packages("reshape2")
install.packages("plotly")
install.packages("PerformanceAnalytics")
install.packages("corrplot")
set.seed(123)
# Load the packages
library(ggplot2)
library(dygraphs)
library(tidyverse)
library(lubridate)
library(zoo)
library(xts)
library(reshape2)
library(plotly)
library(PerformanceAnalytics)
library(corrplot)

###################################################################################################################
# STRATEGY-INDEPENDENT FUNCTIONS

# Function for calculating the means of every column/asset in the dataset:
get_means_vector=function(data){
  means_vector=apply(data,2,mean)
  return(means_vector)
}

# Function for calculating the SD of each column/asset in the dataset:
get_sd_vector=function(data){
  sd_vector=apply(data,2,sd)
  return(sd_vector)
}

# Function for calculating the relative weights of the assets in the dataset:
get_rel_weights_vector=function(weights_vector){
  
  # Calculate the absolute sum of the weights vector:
  abs_sum=abs(sum(weights_vector))

  # Divide each of the values in the absolute weights vector by the absolute sum:
  rel_weights=weights_vector/abs_sum

  #return the vector with the relative asset weights
  return (rel_weights)
}

# Function for getting the in-sample sharpe ratio OR the in-sample portfolio returns
get_insample_sharperatio_returns = function (data, absolute_weights_vector, sharperatio){
  
  # Create a new vector to be filled with the portfolio returns in sample in a loop rowise:
  pf_rtr_in_sample=c(length=length(data[,1]))
  for(i in 1:(length(data[,1]))){
    pf_rtr_in_sample[i]=unlist(data[i,])%*%(get_rel_weights_vector(absolute_weights_vector))
  }
  
  # an if-else construct to determine whether the in-sample sharpe ratio or the in-sample portfolio returns must be returned
  if (sharperatio == TRUE) { # compute the sharpe ratio in-sample applying the formula from DeMiguel page 1928
    return (mean(pf_rtr_in_sample)/sd(pf_rtr_in_sample)) 
  } else { # return the portfolio returns in-sample
    return (pf_rtr_in_sample)
  }
}

# Function for getting the out-of-sample sharpe ratio OR the out-of-sample portfolio returns
get_outofsample_sharperatio_returns = function (M, data, sharperatio, strategy) {
  
  # set the rolling window:
  rolling_window=M
  
  # Calculate length of the new vector with the portfolio returns:
  len_portfolio_returns=length(data[,1])-rolling_window
  
  # Create the vector with the respective length:
  portfolio_returns_outofsample=c(length=len_portfolio_returns)
  
  # Calculate the (in this case 264 - 120) excess returns and add each value in the portfolio_returns vector:
  for(i in 1:len_portfolio_returns){
    # Set the start index for each iteration:
    start_window=i
    
    # Set the start index for each iteration:
    end_window=rolling_window+i-1
    
    # Create a new "time"-matrix, which contains the data for a certain 120-days period
    # with the start index and the end index (rowise):
    time_matrix=data[start_window:end_window,]
    
    # Create the covariance matrix of the "time"-matrix:
    cov_time_matrix=cov(time_matrix)
    
    # Calculate the absolute weights of the assets in row end_window + 1 based on the last 120 rows:
    # an if-else construct to differentiate between both strategies, since these require different functions 
    if (strategy == "mv") { # mean-variance strategy
      weights_vct=get_weights_vector(get_means_vector(time_matrix),cov_time_matrix)
    } else { # Kan and Zhou three-fund portfolio strategy
      weights_vct=get_weights(time_matrix, length(time_matrix[,1]), length(time_matrix[1,]))
    }
    
    # Calculate the portfolio return using the excess returns in row 
    # end_window + 1 of the initial data and the computed relative weights:
    single_pf_return=unlist(data[end_window+1,])%*%get_rel_weights_vector(weights_vct)
    
    # Add each value in the vector portfolio returns:
    portfolio_returns_outofsample[start_window]=single_pf_return
  }
  
  portfolio_returns_outofsample = c(t(portfolio_returns_outofsample))
  
  # an if-else construct to determine whether the out-of-sample sharpe ratio or the out-of-sample portfolio returns must be returned
  if (sharperatio == TRUE) { # compute the out-of-sample sharperatio
    sharpe_ratio_out_of_sample=mean(portfolio_returns_outofsample)/sd(portfolio_returns_outofsample)
    return (sharpe_ratio_out_of_sample) 
  } else { # return the out-of-sample portfolio returns
    return (portfolio_returns_outofsample)
  }
}

# Function for calculating the dynamics of weights all of the T-M periods (T = total observations, M = rolling window)  

# description of the function's parameters
# 1. M: a numeric variable ->length of the rolling window
# 2. data: a data frame -> the relevant data 
# 3. assets: a numeric variable -> amount of assets considered
# 4. cov_matrix: a boolean variable -> the cov_matrix needed to compute the absolute weights vector for the mean-variance strategy;
#    the Khan and Zou strategy doesn't require a cov-matrix for the computation of the absolute weights vector
#    => if assigned TRUE, the function computes the absolute weights vector according to the mean-variance strategy
#       if assigned FALSE, the function computes the absolute weights vector according to the Kan and Zhou strategy


get_weights_dynamics = function(M, data, assets, cov_matrix) {
  
  # create an empty matrix which should collect the relative weights vector for each period considered
  collector = matrix(, nrow = assets, ncol = (length(data[,1])-M))
  colnames(collector) = c(1:(length(data[,1])-M))
  rownames(collector) = c(colnames(data))
  
  for(i in 1:(length(data[,1])-M)) {
    # Set the start index for each iteration:
    start_window=i
    
    # Set the start index for each iteration:
    end_window=M+i-1
    
    # Create a new "time"-matrix with the start index and the end index (rowise) to collect the absolute r.w. difference in t:
    time_matrix=data[start_window:end_window,]
    
    # an if-else construct to distinguish between the mean-variance and the Khan and Zou strategy, since the computation of
    # the abssolute weights vectors are different for both strategies
    
    if (cov_matrix == TRUE) { # mean-variance strategy
      # Create the covariance matrix of the time matrix
      cov_time_matrix=cov(time_matrix)
      # Calculate the absolute weights of the assets over time (always based on the previous 120 rows)
      weights_vct=get_weights_vector(get_means_vector(time_matrix), cov_time_matrix)
    } else { # Khan and Zou three-fund portfolio strategy
      weights_vct=get_weights(time_matrix, M, assets)
    }
    # Calculate the relative weights (the function for the computation of the relative weights is not strategy-dependent)
    rel_weights_vct=get_rel_weights_vector(weights_vct)
    
    # Collect the relative weights vectors for each period in the collector-matrix
    collector[,i] = rel_weights_vct
  }
  
  # transpose the collector matrix
  weight_matrix = t(collector)
  
  # plot a basic ggplot
  p1= ggplot(melt(weight_matrix), aes(x=Var1, y=value, col=Var2))+geom_point()+ggtitle("Dynamics of Weights")+ylab("Relative weights")+xlab("Periods")
  # make the basic ggplot interactive
  dynamics_return = ggplotly(p1)
  
  return(dynamics_return)
}

# Function for plotting the SD dynamics OR the development of a portfolio's returns over time

# description of the function's parameters
# 1. portfolio_returns: a data frame or a matrix, containing the in-sample or out of sample portfolio returns over time
# 2. width: the interval (2 months, 6 months, 1 year etc.), over which the function should be applied -> in months!
# 3. SD: a boolean variable 
#    => if assigned TRUE, the function plots the SD dynamics
#    => if assigned FALSE, the function plots the portfolio's returns dynamics

get_sd_dynamics_or_devreturn = function(portfolio_returns, width, sd, dates_seq) {
  
  pf_returns = portfolio_returns
  pf_returns = data.frame(pf_returns)
  
  # assign the sequence the rows of the return's data frame
  rownames(pf_returns)=dates_seq
  
  # format the data frame as time series
  pf_returns=as.xts(pf_returns,dateFormat="Date")
  
  #an if-else construct to differentiate which function should be applied
  if (sd == TRUE) { #plot the SD dynamics
    b = chart.RollingPerformance(R = pf_returns, width = width, FUN = "StdDev.annualized")
  } else { # plot the portfolio returns' dynamics
    b = chart.RollingPerformance(R = pf_returns, width = width, FUN = "Return.annualized")
  }
  return(b)
}

# Function for making a benchmark comparison: plot the development of the weighted portfolio's returns versus the returns
# of the benchmark portfolio

# hand over a benchmark matrix or data frame, which contains the weighted portfolio's returns and the returnss of the 
# benchmark portfolio over time
bm_comparison = function (benchmark) {
  
  bm_df = data.frame(benchmark)
  
  #turn the benchmark object into a matrix with appropriate length, width and names
  bm_mat = matrix(benchmark, nrow = (length(bm_df[,1])), ncol = (length(bm_df[1,])))
  rownames(bm_mat) = c(1:(length(bm_df[,1])))
  colnames(bm_mat) = c(colnames(benchmark))
  
  #plot a basic ggplot
  bm= ggplot(melt(bm_mat), aes(x=Var1, y=value, col=Var2))+geom_line(alpha = 0.7)+ggtitle("Benchmark comparison")+ylab("Returns")+xlab("Periods")
  
  
  #make the basic ggplot interactive
  plot = ggplotly(bm)
  
  return(plot)
}

# Function for getting the means and sd vector of the data frame from the function above
# hand over the same parameter as in the function above
get_means_and_sd_vector = function(benchmark_c) {
  
  df = data.frame(benchmark_c)
  # compute the means vector of the weighted portfolio's returns and the benchmark returns 
  bm_means = get_means_vector(df)
  
  # compute the SD vector of the weighted portfolio's returns and the benchmark returns
  bm_sd = get_sd_vector(df)
  
  # create an empty matrix with an appropriate length, width and names
  sum_matrix = matrix(, nrow = length(bm_means), ncol = 2)
  rownames(sum_matrix) = c(colnames(df))
  colnames(sum_matrix) = c("means_vector", "sd_vector")
  
  # fill the matrix with the computed vectors
  sum_matrix[,1] = bm_means
  sum_matrix[,2] = bm_sd
  
  # return the transposed matrix
  return(t(sum_matrix))
  
}

#####################################################################################################################
# MEAN VARIANCE STRATEGY SPECIFIC FUNCTIONS

# Function for calculating the absolute weights of the assets in the dataset:
get_weights_vector= function(means_vector,cov_matrix){
  
  # Create the inverse of the passed covariance matrix:
  inverse_cov=solve(cov_matrix)
  
  # Implement the formula for Xt (DeMiguel, Garlappi, Uppal; page 1922):
  x_t=inverse_cov%*%means_vector
  x_t_vector=c(x_t)
  
  # return the vector with the absolute asset weights
  return (x_t_vector)
}


####################################################################################################################
# KAN AND ZHOU THREE FUND STRATEGY PORTFOLIO SPECIFIC FUNCTIONS

# Function for getting the mu-g-parameter from the Kan and Zhou paper page 643
# by dimensions is meant the amout of assets in the dataset
get_mu_g = function(data, dimensions){
  
  # compute the inverse matrix of the covariance matrix of the data
  inverse_matrix = solve(cov(data))
  
  # create a n-dimensional 1-vector
  vec1 = numeric(0)
  i_vector = c(vec1, 1:dimensions)
  i_vector[1:dimensions] = 1
  
  # compute the nominator and the denominator needed for the final division
  nominator =  (t(get_means_vector(data))) %*%inverse_matrix%*%i_vector
  denominator = (t(i_vector)) %*%inverse_matrix%*%i_vector
  
  # get mu_g and make it a vector
  mu_g = nominator/denominator
  mu_g = as.vector(mu_g)
  
  return (mu_g)
}

# Function for computing the psi^2-parameter from the Kan and Zhou paper page 643
get_psi_square = function(data, mu_g_object, dimensions ){
  
  # compute the inverse matrix of the covariance matrix of the data
  inverse_matrix = solve(cov(data))
  
  # create a n-dimensional 1-vector
  vec1 = numeric(0)
  i_vector = c(vec1, 1:dimensions)
  i_vector[1:dimensions] = 1
  
  # get psi^2 and make it a vector
  psi_square = (t(get_means_vector(data)) - mu_g_object*i_vector) %*% inverse_matrix %*% (get_means_vector(data) - mu_g_object*i_vector)
  psi_square = as.vector(psi_square)
  
  return (psi_square)
}

# Function for computing the c3-parameter from the Kan and Zhou paper page 636
get_c3 = function(observations, dimensions){
  
  # get the first and the second term for the final multiplication
  first_term = (observations-dimensions-4)/observations
  second_term = (observations-dimensions-1)/(observations-2)
  
  # get the c3-parameter
  c3 = first_term*second_term
  
  return(c3)
}

# Function for computing the absolute weights vector of a Kan and Zhou three fund portfolio from the Kan and Zhou paper page 642
get_weights = function(data, observations, dimensions) {
  
  # get the c3 paremeter
  c3_object = get_c3(observations, dimensions)
  
  # get the mu_g parameter
  mu_g_object =get_mu_g(data, dimensions)
  
  # get the psi^2 parameter
  psi_square_object = get_psi_square(data, mu_g_object, dimensions)
  
  # compute the inverse matrix of the covariance matrix of the data
  inverse_matrix = solve(cov(data))
  
  # create a n-dimensional 1-vector
  ma = matrix(1, 1, dimensions) 
  i_vector = c(ma)
  
  # get the means vector of the data
  means_vector = get_means_vector(data)
  
  # compute the first and the second term for the final multiplication
  first_term = (psi_square_object/(psi_square_object+(dimensions/observations))) * inverse_matrix %*%means_vector
  second_term = ((dimensions/observations) / (psi_square_object + (dimensions/observations))) * mu_g_object * inverse_matrix %*% i_vector
  
  # get the absolute weights vector
  weights = c3_object * (first_term + second_term)
  
  return(weights)
}

#####################################################################################################################
#0) Subsetting and predefinition of constants

# a subset without the date-column
data.red=data[,-1]

# a subset, containing the excess returns (after subtracting the risk-free rate in column for each period -> column 13 in the original dataset)
# the return in the initial dataset still contain a risk free rate 
#data.new=data.red[,-12]-data.red[,12]
data.new=matrix(,nrow=264, ncol=11)
colnames(data.new)=c(colnames(data[,2:12]))
for (i in 1:11){
  data.new[,i]=data.red[,i]-data.red[,12]
}
data.new=data.frame(data.new)


# create a subset, containing the non-excessive returns (the original data)
data.probe = data.red[,-12]

# create a subset, containing the benchmark SP500-portfolio
bm_SP500 = data.new$S.P500

# determine the amount of assets considered
assets = length(data.new[1,])

# determine the amount of observations considered
observations = length(data[,1])

#####################################################################################################################
# APPLYING THE FUNCTIONS TO THE MEAN VARIANCE STRATEGY

#1) Calculate the Sharpe ratio -> excess returns needed => use the data.new subset

#1.1) out-of-sample

sharpe_ratio_out_of_sample_mv=get_outofsample_sharperatio_returns(120, data.new, TRUE, "mv")
round(sharpe_ratio_out_of_sample_mv, digits=4)
#1.2) in-sample

# 1. alternative
sharpe_ratio_in_sample_mv=get_insample_sharperatio_returns(data.new, get_weights_vector((get_means_vector(data.new)), cov(data.new)), TRUE)
round(sharpe_ratio_in_sample_mv, digits=4)
# 2. alternative: apply the integrated SharpeRatio-function in R as a check
dates=seq(as.Date("1981/01/30"), as.Date("2002/12/31"), by = "1 month",tzone="GMT")-1
rownames(data.new)=dates
time_series=as.xts(data.new,dateFormat="Date")

# compute the relative weights vector: needed for the weights-parameter of the integrated SharpeRatio function in R
aw = get_weights_vector((get_means_vector(data.new)), cov(data.new)) # absolute
rw = get_rel_weights_vector(aw) # relative
SharpeRatio(R = time_series, Rf = 0, p = 0.95, FUN = c("StdDev"),weights = rw, annualize = FALSE)

#3) calculate the certainty-equivalent -> unadjusted returns needed => use the data.probe subset

#3.1) out-of-sample 

# compute the out-of-sample portfolio returns
returns_out= get_outofsample_sharperatio_returns(120, data.probe, FALSE, "mv")

# compute the out-of-sample certainty-equivalent:
certainty_equivalent_out_of_sample_mv=mean(returns_out) - ((1/2)*(var(returns_out)))
round(certainty_equivalent_out_of_sample_mv,digits=4)
#3.2) in-sample

# compute the absolute weights vector of the returns: needed as a parameter for the following function
absolute_weights_in = get_weights_vector((get_means_vector(data.probe)), cov(data.probe))

# compute the in-sample portfolio returns
returns_in = get_insample_sharperatio_returns(data.probe, absolute_weights_in, FALSE )

certainty_equivalent_in_sample_mv=mean(returns_in) - (1/2)*(var(returns_in))
round(certainty_equivalent_in_sample_mv,digits=4)
#4) Plotting the dynamics of a portfolio returns' SD: the plot is based on a 3-month interval 
#   adjusted returns => data.new subset

#4.1) in-sample portfolio returns

# compute the absolute weights vector of the adjusted returns: needed as a parameter for the following function
absolute_weights_adj = get_weights_vector((get_means_vector(data.new)), cov(data.new))

# compute the adjusted in-sample adjusted portfolio returns 
adjusted_returns_in = get_insample_sharperatio_returns(data.new, absolute_weights_adj, FALSE )

# create a sequence of the dates of the observations: needed as a parameter for the following function
dates_in=seq(as.Date("1981/01/30"), as.Date("2002/12/31"), by = "1 month",tzone="GMT")-1

# plot SD dynamics of the in-sample portfolio returns
get_sd_dynamics_or_devreturn(adjusted_returns_in, 3, TRUE, dates_in)

#4.2) out-of-sample portfolio returns

# compute the adjusted out-of-sample adjusted portfolio returns
adjusted_returns_out = get_outofsample_sharperatio_returns(120, data.new, FALSE, "mv")

# create a sequence of the dates of the observations: needed as a parameter for the following function
data.new[121,]
data.new[264,]
dates_out=seq(as.Date("1991/01/29"), as.Date("2002/12/29"), by = "1 month",tzone="GMT")-1

# plot SD dynamics of the out-of-sample portfolio returns
get_sd_dynamics_or_devreturn(adjusted_returns_out, 3, TRUE, dates_out)

#6) Plotting the dynamics of weights of all assets

# based on unadjusted portfolio returns
get_weights_dynamics(120, data.probe, 11, TRUE)

# based on adjusted portfolio returns
get_weights_dynamics(120, data.new, 11, TRUE)

#7) Plotting the dynamics of a portfolio's returns: the plot is based on a 6-month interval

# based on unadjusted in-sample portfolio returns
get_sd_dynamics_or_devreturn(returns_in, 6, FALSE, dates_in)

# based on unadjusted out-of-sample portfolio returns
get_sd_dynamics_or_devreturn(returns_out, 6, FALSE, dates_out)

# based on adjusted in-sample portfolio returns
get_sd_dynamics_or_devreturn(adjusted_returns_in, 6, FALSE, dates_in)

# based on adjusted out-of-sample portfolio returns
get_sd_dynamics_or_devreturn(adjusted_returns_out, 6, FALSE, dates_out)

#8) Benchmark comparison in sample 

benchmark_in_sample = cbind(adjusted_returns_in, bm_SP500)

bm_comparison(benchmark_in_sample)
get_means_and_sd_vector(benchmark_in_sample)

#9) benchmark comparison out of sample 

W = 120 + 1
L = length(data[,1])
benchmark_out = cbind(adjusted_returns_out, bm_SP500[W:L])

bm_comparison(benchmark_out)
get_means_and_sd_vector(benchmark_out)

#####################################################################################################################
# APPLYING THE FUNCTIONS TO THE KAN AND ZHOU THREE FUND PORTFOLIO STRATEGY


#1) Calculate the Sharpe ratio -> excessive returns needed => use the data.new subset

#1.1) out-of-sample

sharpe_ratio_out_of_sample_kz=get_outofsample_sharperatio_returns(120, data.new, TRUE, "kz")
round(sharpe_ratio_out_of_sample_kz, digits=4)
#1.2) in-sample

# 1. alternative
sharpe_ratio_in_sample_kz=get_insample_sharperatio_returns(data.new, get_weights(data.new, observations, assets), TRUE)
round(sharpe_ratio_in_sample_kz, digits=4)

# 2. alternative: apply the integrated SharpeRatio-function in R as a check
dates_kz=seq(as.Date("1981/01/30"), as.Date("2002/12/31"), by = "1 month",tzone="GMT")-1
rownames(data.new)=dates_kz
time_series_kz=as.xts(data.new,dateFormat="Date")

# compute the relative weights vector: needed for the weights-parameter of the integrated SharpeRatio function in R
aw_kz = as.vector(get_weights(data.new, observations, assets)) # absolute
rw_kz = get_rel_weights_vector(aw_kz) # relative
SharpeRatio(R = time_series_kz, Rf = 0, p = 0.95, FUN = c("StdDev"),weights = rw_kz, annualize = FALSE)

#3) calculate the certainty-equivalent -> unadjusted returns needed => use the data.probe subset

#3.1) out-of-sample (passt)

# compute the out-of-sample portfolio returns
returns_out_kz= get_outofsample_sharperatio_returns(120, data.probe, FALSE, "kz")

# compute the out-of-sample sharpe ratio
certainty_equivalent_out_of_sample_kz=mean(returns_out_kz) - ((1/2)*(var(returns_out_kz)))
round(certainty_equivalent_out_of_sample_kz,digits=4)
#3.2) in-sample

# compute the absolute weights vector of the returns: needed as a parameter for the following function
absolute_weights_kz = get_weights(data.probe, observations, assets)

# compute the in-sample portfolio returns
returns_in_kz = get_insample_sharperatio_returns(data.probe, absolute_weights_kz, FALSE )

certainty_equivalent_in_sample_kz=mean(returns_in_kz) - (1/2)*(var(returns_in_kz))
round(certainty_equivalent_in_sample_kz,digits=4)
#4) Plotting the dynamics of a portfolio returns' SD: the plot is based on a 3-month interval 
#   adjusted returns => data.new subset

#4.1) in-sample portfolio returns

# compute the absolute weights vector of the adjusted returns: needed as a parameter for the following function
absolute_weights_adj_kz = get_weights(data.new, observations, assets)

# compute the adjusted in-sample adjusted portfolio returns 
adjusted_returns_in_kz = get_insample_sharperatio_returns(data.new, absolute_weights_adj_kz, FALSE )

# create a sequence of the dates of the observations: needed as a parameter for the following function
dates_in_kz=seq(as.Date("1981/01/30"), as.Date("2002/12/31"), by = "1 month",tzone="GMT")-1

# plot SD dynamics of the in-sample portfolio returns
get_sd_dynamics_or_devreturn(adjusted_returns_in_kz, 3, TRUE, dates_in_kz)

#4.2) out-of-sample portfolio returns

# compute the adjusted out-of-sample adjusted portfolio returns
adjusted_returns_out_kz = get_outofsample_sharperatio_returns(120, data.new, FALSE, "kz")

# create a sequence of the dates of the observations: needed as a parameter for the following function
data.new[121,]
data.new[264,]
dates_out_kz=seq(as.Date("1991/01/29"), as.Date("2002/12/29"), by = "1 month",tzone="GMT")-1

# plot SD dynamics of the out-of-sample portfolio returns
get_sd_dynamics_or_devreturn(adjusted_returns_out_kz, 3, TRUE, dates_out_kz)

#6) Plotting the dynamics of weights of all assets

# based on unadjusted portfolio returns
get_weights_dynamics(120, data.probe, 11, FALSE)

# based on adjusted portfolio returns
get_weights_dynamics(120, data.new, 11, FALSE)

#7) Plotting the dynamics of a portfolio's returns: the plot is based on a 6-month interval

# based on unadjusted in-sample portfolio returns
get_sd_dynamics_or_devreturn(returns_in_kz, 6, FALSE, dates_in_kz)

# based on unadjusted out-of-sample portfolio returns
get_sd_dynamics_or_devreturn(returns_out_kz, 6, FALSE, dates_out_kz)

# based on adjusted in-sample portfolio returns
get_sd_dynamics_or_devreturn(adjusted_returns_in_kz, 6, FALSE, dates_in_kz)

# based on adjusted out-of-sample portfolio returns
get_sd_dynamics_or_devreturn(adjusted_returns_out_kz, 6, FALSE, dates_out_kz)

#8) Benchmark comparison in sample 

benchmark_in_sample_kz = cbind(adjusted_returns_in_kz, bm_SP500)

bm_comparison(benchmark_in_sample_kz)
get_means_and_sd_vector(benchmark_in_sample_kz)

#9) benchmark comparison out of sample 

W = 120 + 1
L = length(data[,1])
benchmark_out_kz = cbind(adjusted_returns_out_kz, bm_SP500[W:L])

bm_comparison(benchmark_out_kz)
get_means_and_sd_vector(benchmark_out_kz)

#10) correlation matrix of the 11 assets + the one benchmark 

corr_matrix = cor(data.new)
corrplot(corr_matrix, method="number", number.cex = 0.5)
