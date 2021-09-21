library(HMM)
options(digits=2)
library(ggplot2)
library("dplyr")
library("tidyr")
# The possible states that the robot might be in. The states goes in a circle.
states <- c('S0', 'S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8', 'S9')

# Symbols, i.e in what state does our signal say we are? 
observations <- c('S0', 'S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8', 'S9')

# The robot starts in a state with equal probability 
start_probability <- rep(0.1, 10)

# Function for creating the transition matrix. There is a 50% chance that 
# the robot stays in the current state and a 50% chance that the robot moves 
# on to the next state. If the state is 10, the robot moves to state 1. 
create_trans <- function(){
  matrix <- matrix(rep(0,100), nrow=10, ncol=10)
  diag(matrix) <- 0.5
  diag(matrix[1:9,2:10]) <- 0.5
  matrix[10,1] <- 0.5
  return(matrix)
}

transProbs <- create_trans()

# Constructing a matrix that contains the probabilities that the robot is in a 
# state SX when the signal comes from state SY. 
emissionProbs <- matrix(rep(0,100), nrow=10, ncol=10)
diag(emissionProbs) <- 0.2
diag(emissionProbs[1:9,2:10]) <- 0.2
diag(emissionProbs[1:8,3:10]) <- 0.2
diag(emissionProbs[2:10,1:9]) <- 0.2
diag(emissionProbs[3:10,1:8]) <- 0.2
emissionProbs[1,9:10] <- 0.2
emissionProbs[2,10] <- 0.2
emissionProbs[9:10,1] <- 0.2
emissionProbs[10,2] <- 0.2


# The hidden markov model
hmm_model <- initHMM(states, observations, start_probability, transProbs, emissionProbs)

# Computing the position where the robot IS most likely to be at time t.
# Filtered algorithm
filtered_algo <- function(alpha, forward){
  filtered <- matrix(rep(0,1000), ncol=100, nrow=10)
  for (i in 1:dim(alpha)[2]){
    filtered[,i] <- alpha[,i] / sum(alpha[,i])
  }
  return(filtered)
}

# Computing the position where the robot WAS most likely to be at time t.
# Smoothed algorithm
smoothed_algo <- function(alpha, beta){
  smoothed <- matrix(rep(0,1000), ncol=100, nrow=10)
  for (i in 1:dim(alpha)[2]){
    smoothed[, i] <- alpha[, i] * beta[, i] / sum(alpha[, i] * beta[, i]) 
  }
  return(smoothed)
}

# Finding the most probable state that the robot was in at time t with 
# regards to the observations 
state_func <- function(algo){
  state <- c()
  for (i in 1:dim(algo)[2]){
    # taking -1 to get the states in the form of 0-9
    state[i] <- which.max(algo[, i] )-1
  }
  return(state)
}

# Misclassification function
missclass <- function(true, pred){
  1 - sum(diag(table(true, pred)))/sum(table(true, pred))
}

# Function for printing the confusion matrices and misclassification rates
genHeatmap <- function (true_states, algo_state, name){
  conf_matr <- table(true_states, algo_state)
  prop.table(table(true_states, algo_state), margin=1)
  missclass(true_states, algo_state)
  conf_matr <- as.data.frame(conf_matr)
  names(conf_matr) <- c("TrueHiddenState","Prediction", "Freq")
  heatmap <- ggplot(conf_matr, aes(Prediction, TrueHiddenState, fill= Freq)) + 
    geom_tile(aes(fill=Freq)) +
    geom_text(aes(label = round(Freq, 1))) + 
    ggtitle(name) + 
    theme(plot.title = element_text(size = 20, face = "bold"))+
    geom_abline(intercept = 0, slope = 1, color="white", size=0.1)  
  return(heatmap)
}


simulation <- function(hmm){
  misclass_vec <- c()
  # 100 time steps simulated from our hidden markov model 
  simulated <- simHMM(hmm, 100)
  # Our observations from our simulation, e.g our signal
  obs <- simulated$observation
  # Forward algo
  forward <- forward(hmm, obs)
  # Backward algo
  backward <- backward(hmm, obs)
  # Computing alpha as e^forward since forward is in log form
  alpha <- exp(forward)
  # Computing beta as e^backward since backward is in log form
  beta <- exp(backward)
  # calculating the probabilities the robot was in state x at a certain time point
  filtered <- filtered_algo(alpha)
  smoothed <- smoothed_algo(alpha, beta)
  # And most likely path
  viterbi <- viterbi(hmm, obs)
  # The predicted states
  filtered_state <- state_func(filtered)
  smoothed_state <- state_func(smoothed)
  # The true states
  true_states <- simulated$states
  
  # Misclassification rates
  misclass_vec[1] <- missclass(true_states, filtered_state)
  misclass_vec[2] <- missclass(true_states, smoothed_state)
  misclass_vec[3] <- missclass(true_states, viterbi)
  return(list("misclass_vec" = misclass_vec,
              "filtered_state" = filtered_state,
              "smoothed_state" = smoothed_state,
              "true_states" = true_states,
              "viterbi" =viterbi,
              "filtered"= filtered))
}

simulate <- simulation(hmm_model)
# What confusion matrices could look like for one simulation
filtered_state <- simulate$filtered_state
smoothed_state <- simulate$smoothed_state
viterbi <- simulate$viterbi
true_states <- simulate$true_states
print("filtered")
filter_confmatr <- genHeatmap(true_states, filtered_state, "Filtered Confusion Matrix")

print("smoothed")
smooth_confmatr <- genHeatmap(true_states, smoothed_state, "Smoothed Confusion Matrix")

print("viterbi")
viterbi_confmatr <- genHeatmap(true_states, viterbi, "Viterbi Confusion Matrix")

misclass_matrix <- matrix(0, nrow=100, ncol=3)
for (i in 1:dim(misclass_matrix)[1]){
  misclass_matrix[i,] <- simulation(hmm_model)$misclass_vec
}

# Plotting the misclassification distribution for the different algorithms
df <- data.frame(misclass_matrix)

misclass_plot <- ggplot(data=df)+
  ggtitle("Missclassification Error") +
  geom_density(aes(x=X2, colour="Smoothed"), fill = "palegreen3", alpha=0.6 )+
  geom_density(aes(x=X3, color="Viterbi"), fill = "light blue", alpha=0.5)+
  geom_density(aes(x=X1, colour="Filtered"), fill="tomato", alpha=0.3)+
  
  labs(colour="Algorithm") + 
  xlab('Misclass error') + ylab('Frequency') + 
  xlim(0,1)

# Empirical Estimators of Entropy and Mutual Information and Related Quantities
library(entropy)

# Calculating the entropy at different time steps for 100 simulation and 
# averaging the results
entropy_matrix <- matrix(0, ncol=100, nrow=100)
for (j in 1:100){
  entropy <- c()
  filtered <- simulation(hmm_model)$filtered
  for (i in 1:dim(filtered)[2]){
    
    entropy[i] <- entropy.empirical(filtered[, i], unit="log2")  
  }
  entropy_matrix[,j] <- entropy
}

entropy <- c()
for (i in 1:100){
  entropy[i] <- sum(entropy_matrix[i,])/dim(entropy_matrix)[2]
  
}


entropy <- as.data.frame(entropy)
entropy$x <- seq(1,100,1)
entropy_plot <- ggplot() +
  geom_line(data=entropy, aes(y=entropy, x=x)) +
  xlab("Time step") + 
  ggtitle("Entropy of the filtered distribution at different time points")

# Calculating the probability that robot is in a certain hidden state at time 
# step 101 given the probabilities at time 100.
state_101 <- rep(0, 10)
for (i in 1:10){
  state_101[i] <- state_101[i] + 0.5 * filtered[i, 100]
  if (i==10){
    state_101[1] <- state_101[1] + 0.5 * filtered[i, 100]
  }else{
    state_101[i+1] <- state_101[i+1] + 0.5 * filtered[i, 100]
    
  }
  
}
state_101
filtered[,100]


