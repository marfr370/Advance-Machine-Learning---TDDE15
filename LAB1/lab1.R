setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(bnlearn)
library(RBGL)
library(Rgraphviz)
library(gRain)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("RBGL")
#BiocManager::install("Rgraphviz")
#BiocManager::install("gRain")

##############################################################################
##########################  ASSIGNMENT 1 #####################################
##############################################################################
# Show that multiple runs of the hill-climbing algorithm can return non-equivalent Bayesian network (BN) structures. 
# Two DAGs represent the same independencies according to the separation criterion (i.e. they are equivalent) 
# if and only if they have the same adjacencies and unshielded colliders.
# The causal variables influencing the collider are themselves not necessarily associated. 
# If the parents are not adjacent, the collider is unshielded.

data("asia")
nodes = names(asia)
e = empty.graph(nodes)


hc_restart0 = hc(asia, restart=0, score="bde", iss=15) 
hc_restart100 = hc(asia, restart=100, score="bde", iss=15) 

information <- function(hc1, hc2){
  print(arcs(hc1))
  print(arcs(hc2))
  print('VSTRUCTS')
  #moral = false is the default value and it makes vstructs return unshielded colliders
  print(vstructs(hc1, moral=FALSE))
  print(vstructs(hc2, moral=FALSE))
  print('CPDAG:')
  cpdag(hc1)
  cpdag(hc2)
  all.equal(hc1,hc2)
}

information(hc_restart0, hc_restart100)
library(Rgraphviz)


#graphviz.compare(hc_restart0, hc_restart100, layout = "dot", shape = "circle", main = NULL,
#                sub = NULL, diff = "from-first", diff.args = list())




# Explain why this happens.
# Since the hill climbing algorithm uses a randomized search approach and is not 
# asymptotically correct under faithfulness
# it may get stuck in a local maxima where the solution can not be approved upon by
# neighboring states since the cost function is higher in those. 
# To attempt to avoid getting stuck in local optima, one could use restarts (i.e. repeated local search)
# which is what we define in the hc algorithm. 

# iss is the user-defined imaginary sample size (the higher the less regularization). 

# Bayesian score favours models that trade off fit of data and model complexity.

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##############################################################################
##########################  ASSIGNMENT 2 #####################################
##############################################################################
library(bnlearn)
library(RBGL)
library(Rgraphviz)
library(gRain)

data('asia')
n=dim(asia)[1]
set.seed(12345)
id=sample(1:n, floor(n*0.8))
data=asia[id,]
test = asia[-id,]


information <- function(hc1, hc2){
  #print(arcs(hc1))
  #print(arcs(hc2))
  print('VSTRUCTS')
  #moral = false is the default value and it makes vstructs return unshielded colliders
  print(vstructs(hc1, moral=FALSE))
  print(vstructs(hc2, moral=FALSE))
  print('CPDAG:')
  #print(cpdag(hc1))
  #print(cpdag(hc2))
  all.equal(hc1,hc2)
}


hc_restart100 = hc(data, restart=100, score="bde", iss=5)

true_dag = model2network("[A][S][T|A][L|S][B|S][D|B:E][E|T:L][X|E]")

information(hc_restart100, true_dag)

#graphviz.compare(hc_restart100, true_dag, layout = "dot", shape = "circle", main = NULL,
#                 sub = NULL, diff = "from-first", diff.args = list())

#graphviz.plot(true_dag)

bn_fitted = bn.fit(hc_restart100, data = data)
true_bn_fitted = bn.fit(true_dag, data= data)

grain_obj = as.grain(bn_fitted)
true_grain_obj = as.grain(true_bn_fitted)

juncTree = compile(grain_obj)
true_juncTree = compile(true_grain_obj)
summary(juncTree)

#plot(juncTree)

predictfunc <- function(testdata, tree, xvars){
  prediction = c()
  for (i in 1:dim(testdata)[1]){
    xvals <- testdata[i, xvars]
    xvals <- as.vector(unlist(xvals, use.names = FALSE))
    states <- setEvidence(tree, nodes = xvars, states = xvals)
    probs <- querygrain(states, nodes="S")
    prediction[i] <- names(which.max(probs$S))
    
  }
  return(prediction)
}

missclass <- function(true, pred){
  1 - sum(diag(table(true, pred)))/sum(table(true, pred))
}

xvars <- colnames(test[, -which(names(test) == "S")])
pred <- predictfunc(test, juncTree, xvars)
asia_Pred <- predictfunc(test, true_juncTree, xvars)
true_s <- test$S
table(true_s, pred)
table(true_s, asia_Pred)
table(pred, asia_Pred)

missclass(true_s, pred)
missclass(true_s, asia_Pred)
##############################################################################
##########################  ASSIGNMENT 3 #####################################
##############################################################################
#graphviz.plot(bn_fitted)
# A subset that contains all the useful information is called a Markov blanket. 
# Since S has no children and therefore no parents of it's children,
# the Markov blanket only consists of S's parents themselves. 
# In an undirected graph (markov network) the markov blanket is all of the variables
# that are connected to the variable in question via an edge. 
mb_nodes = mb(x=bn_fitted, node='S')
xvars <- colnames(test[, mb_nodes])
xvars

mb_pred <- predictfunc(test, juncTree, xvars)
table(true_s, pred)
missclass(true_s, pred)


##############################################################################
##########################  ASSIGNMENT 4 #####################################
##############################################################################
# Classification is then done by applying Bayes rule to compute the probability of C given
# the particular instance of A1,...,An, and then predicting the class with the highest posterior probability. 
# This computation is rendered feasible by making a strong independence
# assumption: all the attributes Ai are conditionally independent given the value of the class C
# For the Naive Bayesian Network which is more well-known nowadays, 
# all features are considered as attributes and are independent given the class.
naive_dag = model2network("[S][A|S][B|S][T|S][L|S][E|S][X|S][D|S]")
#graphviz.plot(naive_dag, main="Naïve Bayes Classifier")
naive_fitted = bn.fit(naive_dag, data = data)
naiveGrain_obj = as.grain(naive_fitted)
naiveJuncTree = compile(naiveGrain_obj)
xvars <- colnames(test[, -which(names(test) == "S")])
naive_prediction <- predictfunc(test, naiveJuncTree, xvars)
table(true_s, naive_prediction)
missclass(true_s, naive_prediction)



