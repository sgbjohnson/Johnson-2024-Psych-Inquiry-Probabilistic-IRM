### FUNCTIONS ###
# Function used in simulation to adjust 'global' levels of trait means to create correlated trait means across groups
tau.calc <- function(global,change,scaling){
  for(i in (1:length(global))){
    if(change[i] < 0){global[i] <- global[i]+global[i]*change[i]*scaling}
    else{global[i] <- global[i] + (1-global[i])*change[i]*scaling}
  }
  return(global)
}
# Function used in simulation to adjust self-sample to account for over/undersampling of desirable/undesirable traits
wishful.adjust <- function(tau,gamma,delta){
  if(delta > 0.5){tau + gamma*(1-tau)*(2*delta-1)}
  else{tau - gamma*tau*(1-2*delta)}
}

### SIMULATION PROCEDURE ###
# Perform a single run of the simulation
simulate <- function(
    ntraits=8, # Number of traits (set to 10 in paper)
    change.amount=0.65, # Higher values of this parameter makes ingroup and outgroup more distinct (set to 0.65 in paper)
    homogeneity=5, # Higher values indicate more homogeneous groups (set to 5 in paper)
    group.size=1, # Number of individuals within each group (set to 20 in paper)
    sample.size=5, # Number of observations of each trait for each individual (set to 5 in paper)
    gamma=0.56, # Parameter controlling extent of "wishful thinking" in self-sampling (set to 0.56 in paper)
    sample.size.self=10, # Parameter controlling precision of self-knowledge (set to 10 in paper)
    pi.ingroup=0.5, # Parameter controlling projection to the in-group (set to 0.5 in paper)
    pi.outgroup=0,
    noise.size=9 # Parameter controlling amount of noise
     ){ 
  # Set trait desirability (uniform from 0 to 1)
    deltas <- runif(ntraits,0,1)
     
  # Set trait means for the ingroup and outgroup
  # Parameters used in paper yield a correlation of about 0.5 between ingroup and outgroup means
    tau.means.global <- runif(ntraits,0,1)
    tau.means.ingroup.change <- runif(ntraits,-1,1)
    tau.means.outgroup.change <- runif(ntraits,-1,1)
    tau.means.ingroup <- tau.calc(tau.means.global,tau.means.ingroup.change,change.amount)
    tau.means.outgroup <- tau.calc(tau.means.global,tau.means.outgroup.change,change.amount)
  
  # Re-parameterize group-level trait distributions as beta distributions
    alphas.ingroup <- tau.means.ingroup*homogeneity
    betas.ingroup <- homogeneity - alphas.ingroup
    alphas.outgroup <- tau.means.outgroup*homogeneity
    betas.outgroup <- homogeneity - alphas.outgroup
  
  # Create individual-level trait levels (taus) from group-level distributions
    taus.ingroup <- t(replicate(group.size,apply(as.matrix(c(1:ntraits)),1,function(X){rbeta(1,alphas.ingroup[X],betas.ingroup[X])})))
    taus.outgroup <- t(replicate(group.size,apply(as.matrix(c(1:ntraits)),1,function(X){rbeta(1,alphas.outgroup[X],betas.outgroup[X])})))
  
  # Create observational data for each individual
    data.ingroup <- apply(taus.ingroup,2,function(X){apply(as.matrix(X),1,function(Y){rbinom(1,sample.size,Y)})})
    data.outgroup <- apply(taus.outgroup,2,function(X){apply(as.matrix(X),1,function(Y){rbinom(1,sample.size,Y)})})
    if(group.size==1){
      data.ingroup <- t(matrix(data.ingroup))
      data.outgroup <- t(matrix(data.outgroup))
    }
    
  # Create sample for the self
    # Set true taus, drawn from ingroup distribution
    tau.self <- apply(as.matrix(c(1:ntraits)),1,function(X){rbeta(1,alphas.ingroup[X],betas.ingroup[X])})
    # Adjust sample for wishful thinking
    prob.self <- apply(as.matrix(c(1:ntraits)),1,function(X){wishful.adjust(tau.self[X],gamma,deltas[X])})
    sample.self <- apply(as.matrix(c(1:ntraits)),1,function(X){rbinom(1,sample.size.self,prob.self[X])})
    judgments.self <- sample.self/sample.size.self
    
  # Create noise
    noise.ingroup <- matrix(runif(group.size*ntraits,0,noise.size),ncol=ntraits)
    noise.outgroup <- matrix(runif(group.size*ntraits,0,noise.size),ncol=ntraits)
    
## CALCULATE MEAN OF POSTERIOR DISTRIBUTIONS ##
    # Follows Table 1 in paper, plus noise term
      # Case 1: Stereotypes (No), Observation (No), Projection (No)
        ingroup.N.N.N <- t(apply(as.matrix(c(1:group.size)),1, # Iterating over each group member
                             function(X){apply(as.matrix(c(1:ntraits)),1, # Iterating over each trait
                                               function(Y){(1 + noise.ingroup[Y])/(2 + noise.size)})})) # Mean of posterior on Line 1 of Table 1
        outgroup.N.N.N <- t(apply(as.matrix(c(1:group.size)),1, # Iterating over each group member
                                 function(X){apply(as.matrix(c(1:ntraits)),1, # Iterating over each trait
                                                   function(Y){(1 + noise.outgroup[Y])/(2 + noise.size)})})) # Mean of posterior on Line 1 of Table 1
        
      # Case 2: Stereotypes (Yes), Observation (No), Projection (No)
        ingroup.Y.N.N <- t(apply(as.matrix(c(1:group.size)),1,
                             function(X){apply(as.matrix(c(1:ntraits)),1,
                                               function(Y){(alphas.ingroup[Y] + noise.ingroup[Y])/(homogeneity + noise.size)})})) # Mean of posterior on Line 2 of Table 1
        outgroup.Y.N.N <- t(apply(as.matrix(c(1:group.size)),1,
                                 function(X){apply(as.matrix(c(1:ntraits)),1,
                                                   function(Y){(alphas.outgroup[Y] + noise.outgroup[Y])/(homogeneity + noise.size)})})) # Mean of posterior on Line 2 of Table 1    
        
      # Case 3: Stereotypes (No), Observation (Yes), Projection (No)
        ingroup.N.Y.N <- t(apply(as.matrix(c(1:group.size)),1,
                             function(X){apply(as.matrix(c(1:ntraits)),1,
                                               function(Y){(1 + data.ingroup[X,Y] + noise.ingroup[Y])/(2 + sample.size + noise.size)})})) # Mean of posterior on Line 3 of Table 1
        outgroup.N.Y.N <- t(apply(as.matrix(c(1:group.size)),1,
                                 function(X){apply(as.matrix(c(1:ntraits)),1,
                                                   function(Y){(1 + data.outgroup[X,Y] + noise.outgroup[Y])/(2 + sample.size + noise.size)})})) # Mean of posterior on Line 3 of Table 1
        
      # Case 4: Stereotypes (Yes), Observation (Yes), Projection (No)
        ingroup.Y.Y.N <- t(apply(as.matrix(c(1:group.size)),1,
                             function(X){apply(as.matrix(c(1:ntraits)),1,
                                               function(Y){(alphas.ingroup[Y] + data.ingroup[X,Y] + noise.ingroup[Y])/(homogeneity + sample.size + noise.size)})})) # Mean of posterior on Line 4 of Table 1
        outgroup.Y.Y.N <- t(apply(as.matrix(c(1:group.size)),1,
                                 function(X){apply(as.matrix(c(1:ntraits)),1,
                                                   function(Y){(alphas.outgroup[Y] + data.outgroup[X,Y] + noise.outgroup[Y])/(homogeneity + sample.size + noise.size)})})) # Mean of posterior on Line 4 of Table 1
        
      # Case 5: Stereotypes (No), Observation (No), Projection (Yes) [Situation of original IRM]
        ingroup.N.N.Y <- t(apply(as.matrix(c(1:group.size)),1, # Iterating over each group member
                              function(X){apply(as.matrix(c(1:ntraits)),1, # Iterating over each trait
                                                function(Y){(1 + pi.ingroup*sample.self[Y] + noise.ingroup[Y])/(2 + pi.ingroup*sample.size.self + noise.size)})})) # Mean of posterior on Line 5 of Table 1
        outgroup.N.N.Y <- t(apply(as.matrix(c(1:group.size)),1, # Iterating over each group member
                                 function(X){apply(as.matrix(c(1:ntraits)),1, # Iterating over each trait
                                                   function(Y){(1 + pi.outgroup*sample.self[Y] + noise.outgroup[Y])/(2 + pi.outgroup*sample.size.self + noise.size)})})) # Mean of posterior on Line 5 of Table 1 (adds random perturbation to account for no outgroup projection to ensure variance)
        
      # Case 6: Stereotypes (Yes), Observation (No), Projection (Yes)
        ingroup.Y.N.Y <- t(apply(as.matrix(c(1:group.size)),1,
                              function(X){apply(as.matrix(c(1:ntraits)),1,
                                                function(Y){(alphas.ingroup[Y] + pi.ingroup*sample.self[Y] + noise.ingroup[Y])/(homogeneity + pi.ingroup*sample.size.self + noise.size)})})) # Mean of posterior on Line 6 of Table 1
        outgroup.Y.N.Y <- t(apply(as.matrix(c(1:group.size)),1,
                                 function(X){apply(as.matrix(c(1:ntraits)),1,
                                                   function(Y){(alphas.outgroup[Y] + pi.outgroup*sample.self[Y] + noise.outgroup[Y])/(homogeneity + pi.outgroup*sample.size.self + noise.size)})})) # Mean of posterior on Line 6 of Table 1
        
      # Case 7: Stereotypes (No), Observation (Yes), Projection (Yes)
        ingroup.N.Y.Y <- t(apply(as.matrix(c(1:group.size)),1,
                              function(X){apply(as.matrix(c(1:ntraits)),1,
                                                function(Y){(1 + data.ingroup[X,Y] + pi.ingroup*sample.self[Y] + noise.ingroup[Y])/(2 + sample.size + pi.ingroup*sample.size.self + noise.size)})})) # Mean of posterior on Line 7 of Table 1
        outgroup.N.Y.Y <- t(apply(as.matrix(c(1:group.size)),1,
                                 function(X){apply(as.matrix(c(1:ntraits)),1,
                                                   function(Y){(1 + data.outgroup[X,Y] + pi.outgroup*sample.self[Y] + noise.outgroup[Y])/(2 + sample.size + pi.outgroup*sample.size.self + noise.size)})})) # Mean of posterior on Line 7 of Table 1
        
      # Case 8: Stereotypes (Yes), Observation (Yes), Projection (Yes)
        ingroup.Y.Y.Y <- t(apply(as.matrix(c(1:group.size)),1,
                              function(X){apply(as.matrix(c(1:ntraits)),1,
                                                function(Y){(alphas.ingroup[Y] + data.ingroup[X,Y] + pi.ingroup*sample.self[Y] + noise.ingroup[Y])/(homogeneity + sample.size + pi.ingroup*sample.size.self + noise.size)})})) # Mean of posterior on Line 8 of Table 1
        outgroup.Y.Y.Y <- t(apply(as.matrix(c(1:group.size)),1,
                                 function(X){apply(as.matrix(c(1:ntraits)),1,
                                                   function(Y){(alphas.outgroup[Y] + data.outgroup[X,Y] + pi.outgroup*sample.self[Y] + noise.outgroup[Y])/(homogeneity + sample.size + pi.outgroup*sample.size.self + noise.size)})})) # Mean of posterior on Line 8 of Table 1

## CALCULATING CORRELATIONS ##
  # Calculate r(s,i) for the average group member [to match original IRM, should be ~0.5]
        si.mean.cors.N.N.N <- cor(colMeans(ingroup.N.N.N),judgments.self)
        si.mean.cors.Y.N.N <- cor(colMeans(ingroup.Y.N.N),judgments.self)
        si.mean.cors.N.Y.N <- cor(colMeans(ingroup.N.Y.N),judgments.self)
        si.mean.cors.Y.Y.N <- cor(colMeans(ingroup.Y.Y.N),judgments.self)
        si.mean.cors.N.N.Y <- cor(colMeans(ingroup.N.N.Y),judgments.self)
        si.mean.cors.Y.N.Y <- cor(colMeans(ingroup.Y.N.Y),judgments.self)
        si.mean.cors.N.Y.Y <- cor(colMeans(ingroup.N.Y.Y),judgments.self)
        si.mean.cors.Y.Y.Y <- cor(colMeans(ingroup.Y.Y.Y),judgments.self)
  # Calculate r(s,o) for the average group member [to match original IRM, should be ~0]
        so.mean.cors.N.N.N <- cor(colMeans(outgroup.N.N.N),judgments.self)
        so.mean.cors.Y.N.N <- cor(colMeans(outgroup.Y.N.N),judgments.self)
        so.mean.cors.N.Y.N <- cor(colMeans(outgroup.N.Y.N),judgments.self)
        so.mean.cors.Y.Y.N <- cor(colMeans(outgroup.Y.Y.N),judgments.self)
        so.mean.cors.N.N.Y <- cor(colMeans(outgroup.N.N.Y),judgments.self)
        so.mean.cors.Y.N.Y <- cor(colMeans(outgroup.Y.N.Y),judgments.self)
        so.mean.cors.N.Y.Y <- cor(colMeans(outgroup.N.Y.Y),judgments.self)
        so.mean.cors.Y.Y.Y <- cor(colMeans(outgroup.Y.Y.Y),judgments.self)
  # Calculate r(i,o) for the *average* ingroup member and the *average* outgroup member
        io.mean.cors.N.N.N <- cor(colMeans(ingroup.N.N.N),colMeans(outgroup.N.N.N))
        io.mean.cors.Y.N.N <- cor(colMeans(ingroup.Y.N.N),colMeans(outgroup.Y.N.N))
        io.mean.cors.N.Y.N <- cor(colMeans(ingroup.N.Y.N),colMeans(outgroup.N.Y.N))
        io.mean.cors.Y.Y.N <- cor(colMeans(ingroup.Y.Y.N),colMeans(outgroup.Y.Y.N))
        io.mean.cors.N.N.Y <- cor(colMeans(ingroup.N.N.Y),colMeans(outgroup.N.N.Y))
        io.mean.cors.Y.N.Y <- cor(colMeans(ingroup.Y.N.Y),colMeans(outgroup.Y.N.Y))
        io.mean.cors.N.Y.Y <- cor(colMeans(ingroup.N.Y.Y),colMeans(outgroup.N.Y.Y))
        io.mean.cors.Y.Y.Y <- cor(colMeans(ingroup.Y.Y.Y),colMeans(outgroup.Y.Y.Y))
  # Calculate r(s,d) [to match original IRM, should be ~0.5]
        sd.cor <- cor(judgments.self,deltas)
  # Calculate r(i,d) for the *average* ingroup member
        id.mean.cors.N.N.N <- cor(colMeans(ingroup.N.N.N),deltas)
        id.mean.cors.Y.N.N <- cor(colMeans(ingroup.Y.N.N),deltas)
        id.mean.cors.N.Y.N <- cor(colMeans(ingroup.N.Y.N),deltas)
        id.mean.cors.Y.Y.N <- cor(colMeans(ingroup.Y.Y.N),deltas)
        id.mean.cors.N.N.Y <- cor(colMeans(ingroup.N.N.Y),deltas)
        id.mean.cors.Y.N.Y <- cor(colMeans(ingroup.Y.N.Y),deltas)
        id.mean.cors.N.Y.Y <- cor(colMeans(ingroup.N.Y.Y),deltas)
        id.mean.cors.Y.Y.Y <- cor(colMeans(ingroup.Y.Y.Y),deltas)
  # Calculate r(o,d) for the *average* ingroup member
        od.mean.cors.N.N.N <- cor(colMeans(outgroup.N.N.N),deltas)
        od.mean.cors.Y.N.N <- cor(colMeans(outgroup.Y.N.N),deltas)
        od.mean.cors.N.Y.N <- cor(colMeans(outgroup.N.Y.N),deltas)
        od.mean.cors.Y.Y.N <- cor(colMeans(outgroup.Y.Y.N),deltas)
        od.mean.cors.N.N.Y <- cor(colMeans(outgroup.N.N.Y),deltas)
        od.mean.cors.Y.N.Y <- cor(colMeans(outgroup.Y.N.Y),deltas)
        od.mean.cors.N.Y.Y <- cor(colMeans(outgroup.N.Y.Y),deltas)
        od.mean.cors.Y.Y.Y <- cor(colMeans(outgroup.Y.Y.Y),deltas)
  # Calculate validity of ingroup judgments for the *average* ingroup member
        i.val.N.N.N <- cor(colMeans(ingroup.N.N.N),colMeans(taus.ingroup))
        i.val.Y.N.N <- cor(colMeans(ingroup.Y.N.N),colMeans(taus.ingroup))
        i.val.N.Y.N <- cor(colMeans(ingroup.N.Y.N),colMeans(taus.ingroup))
        i.val.Y.Y.N <- cor(colMeans(ingroup.Y.Y.N),colMeans(taus.ingroup))
        i.val.N.N.Y <- cor(colMeans(ingroup.N.N.Y),colMeans(taus.ingroup))
        i.val.Y.N.Y <- cor(colMeans(ingroup.Y.N.Y),colMeans(taus.ingroup))
        i.val.N.Y.Y <- cor(colMeans(ingroup.N.Y.Y),colMeans(taus.ingroup))
        i.val.Y.Y.Y <- cor(colMeans(ingroup.Y.Y.Y),colMeans(taus.ingroup))
  # Calculate validity of outgroup judgments for the *average* outgroup member
        o.val.N.N.N <- cor(colMeans(outgroup.N.N.N),colMeans(taus.outgroup))
        o.val.Y.N.N <- cor(colMeans(outgroup.Y.N.N),colMeans(taus.outgroup))
        o.val.N.Y.N <- cor(colMeans(outgroup.N.Y.N),colMeans(taus.outgroup))
        o.val.Y.Y.N <- cor(colMeans(outgroup.Y.Y.N),colMeans(taus.outgroup))
        o.val.N.N.Y <- cor(colMeans(outgroup.N.N.Y),colMeans(taus.outgroup))
        o.val.Y.N.Y <- cor(colMeans(outgroup.Y.N.Y),colMeans(taus.outgroup))
        o.val.N.Y.Y <- cor(colMeans(outgroup.N.Y.Y),colMeans(taus.outgroup))
        o.val.Y.Y.Y <- cor(colMeans(outgroup.Y.Y.Y),colMeans(taus.outgroup))   
  
## SUMMARIZING DATA FOR OUTPUT ##           
  # Create correlation vectors
  si.cors <- c(si.mean.cors.N.N.N,si.mean.cors.Y.N.N,si.mean.cors.N.Y.N,si.mean.cors.Y.Y.N,si.mean.cors.N.N.Y,si.mean.cors.Y.N.Y,si.mean.cors.N.Y.Y,si.mean.cors.Y.Y.Y)
  so.cors <- c(so.mean.cors.N.N.N,so.mean.cors.Y.N.N,so.mean.cors.N.Y.N,so.mean.cors.Y.Y.N,so.mean.cors.N.N.Y,so.mean.cors.Y.N.Y,so.mean.cors.N.Y.Y,so.mean.cors.Y.Y.Y)
  io.cors <- c(io.mean.cors.N.N.N,io.mean.cors.Y.N.N,io.mean.cors.N.Y.N,io.mean.cors.Y.Y.N,io.mean.cors.N.N.Y,io.mean.cors.Y.N.Y,io.mean.cors.N.Y.Y,io.mean.cors.Y.Y.Y)
  sd.cors <- c(rep(sd.cor,8))
  id.cors <- c(id.mean.cors.N.N.N,id.mean.cors.Y.N.N,id.mean.cors.N.Y.N,id.mean.cors.Y.Y.N,id.mean.cors.N.N.Y,id.mean.cors.Y.N.Y,id.mean.cors.N.Y.Y,id.mean.cors.Y.Y.Y)
  od.cors <- c(od.mean.cors.N.N.N,od.mean.cors.Y.N.N,od.mean.cors.N.Y.N,od.mean.cors.Y.Y.N,od.mean.cors.N.N.Y,od.mean.cors.Y.N.Y,od.mean.cors.N.Y.Y,od.mean.cors.Y.Y.Y)
  i.val.cors <- c(i.val.N.N.N,i.val.Y.N.N,i.val.N.Y.N,i.val.Y.Y.N,i.val.N.N.Y,i.val.Y.N.Y,i.val.N.Y.Y,i.val.Y.Y.Y)
  o.val.cors <- c(o.val.N.N.N,o.val.Y.N.N,o.val.N.Y.N,o.val.Y.Y.N,o.val.N.N.Y,o.val.Y.N.Y,o.val.N.Y.Y,o.val.Y.Y.Y)
  
  # Outputs
  model.checks <- c(cor(tau.means.ingroup,tau.means.outgroup),# IO true group-level trait distribution correlation
                    cor(colMeans(taus.ingroup),colMeans(taus.outgroup)), # IO true individual-level trait distribution correlation
                    cor(judgments.self,deltas), # SD correlation (Krueger et al. assume 0.5)
                    cor(judgments.self,colMeans(ingroup.N.N.Y)), # SI correlation with social projection (Krueger et al. assume 0.5)
                    cor(judgments.self,colMeans(outgroup.N.N.Y)), # SO correlation with social projection (Krueger et al. assume 0)
                    cor(judgments.self,colMeans(taus.ingroup)), # SI validity (Kreuger et al. assume about 0.5)
                    cor(judgments.self,colMeans(taus.outgroup))) # SO validity (Krueger et al. assume 0, but we set this around 0.25)
  model.outputs <- c(si.cors,so.cors,io.cors,sd.cors,id.cors,od.cors,i.val.cors,o.val.cors)
  return(c(model.checks,model.outputs))
} 

### GENERATE SIMULATION RESULTS ###

# Overall model results
set.seed(101011)
dat <- t(replicate(50000,simulate()))

summ.checks <- rbind(round(apply(dat[,1:7],2,mean),2),round(apply(dat[,1:7],2,sd),2))
colnames(summ.checks) <- c("IO mean true","IO tau true","SD","SI proj","SO proj","SI valid","SO valid")
row.names(summ.checks) <- c("M","SD")

summ.outputs <- matrix(round(colMeans(dat[,8:dim(dat)[2]]),2),nrow=8,byrow=F)
row.names(summ.outputs) <- c("N.N.N","Y.N.N","N.Y.N","Y.Y.N","N.N.Y","Y.N.Y","N.Y.Y","Y.Y.Y") 
colnames(summ.outputs) <- c("SI","SO","IO","SD","ID","OD","Val I","Val O")

summ.checks;summ.outputs
