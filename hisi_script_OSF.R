# Barnby, J.M. & Moutoussis, M., 2020 #

# -------------------------------------------------------------------------

#rm(list=ls()) #to clear environment

## Define Functions

noisyBino <- function(pSucc, U, binN ){
  # this is the whole mdf for applying extra uncertainty U to binomial 
  # distro of binN-1 draws w. success param pSucc as per
  #  MDF = binopdf(0:(binN-1),binN-1,pSuccs); MDF = MDF ^ 1/U ... etc.
  n = binN-1;
  MDF = dbinom(0:n,n,pSucc);
  MDF = MDF ^ (1/U);
  return( MDF / sum(MDF));
}

# -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
# beta distribution with scaled x values, so that instead of bet. 0 and 1,
# the random var obtains values between lo and hi. Can have lo > hi, in
# which case they are mirrored.
dbetasc <- function(x, shape1, shape2, lo=0, hi=1, ncp=0, log=FALSE){
  
  xtr <- (x-lo)/(hi-lo); # will work even if hi<lo  
  if (log==FALSE) {
    return( dbeta( xtr, shape1, shape2, ncp, log)/abs(hi-lo) );
  }
  else {
    return( dbeta( xtr, shape1, shape2, ncp, log) - log(abs(hi-lo)) );
  }
}
# -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
rbetasc <- function(n, shape1, shape2, lo=0, hi=1, ncp=0){
  auxv <- rbeta(n, shape1, shape2, ncp);
  return(  lo+ auxv*(hi-lo) );
}

rmdf <- function(k, mdf,check=FALSE){
  cdf <- cumsum(mdf);
  n <- length(mdf);
  if (check){
    Tol=1e-8;  # for rough checking ...
    if (abs(cdf[n] - 1) > Tol)
    { stop('mdf provided does not add up to 1 within 1e-16'); };
  }
  cdf[n] <- 1.0;  # force it to 1, for good measure ...
  x <- runif(k);
  CDF <- repRowMat(cdf,k);
  return( 1+n - rowSums(CDF >= x));
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
repRowMat <- function(x, renN) {
  #
  if (!(is.vector(x))){
    stop('non - vector 1st argument in repRowMat')
  }
  return(t(matrix(rep(x,renN),nrow=length(x)))); 
}
###libraries
library(patchwork)
library(dplyr)
library(readr)
library(plyr)
library(tidyr)
library(psych)
library(car)
library(reshape)
library(reshape2)

#Visualisation
library(sjPlot)
library(sjlabelled)
library(sjmisc)
library(gtable)
library(ggplot2)
library(ggpubr)
library(ggforce)
library(grid)
library(gridExtra)
library(qgraph)
library(ggvis)
library("see")
library(purrr)
library(pacman)
library(RColorBrewer)
library(lattice)
library(corrplot)
library(dabestr)
library(ggstatsplot)
library(jtools)
library(ggstance)

######### GENERATIVE MODEL FUNCTIONS (LOG-LIKELIHOOD) ################

# 0. Toy (basic demo) model Log-likelihood ~~~~~~~~~~~~~~~~~~~~~~~~~~
#
L0 <- function(resp,stim,param){
  # Assume that the pt. has been presented with three 'partners'
  # which are the rows of 'stim'. Their params here are: 
  # param[1] = alpha and param[2] = beta of a Beta distribution.
  # so that the Beta distribution just encodes their belief that the 
  # Other is 'nice'=1 (vs. 'nasty'=0)
  # param[3] is a learning rate that shifts expectations depending on experience
  # param[4] is a decision noise
 
  
  sessN = 3;   maxTr = 6;                  # hard-coded as same for all.
  sll = 0;                                 # Initialize main output, sum log lik  
  a0 = param[1]; b0=param[2]; e=param[3];  tau=param[4];   # copies for clarity
  evo <- matrix(NA,3*6,5) # a record of how belief params will evolve
  colnames(evo) <- c('a','b','expectation','response','lnLikDensity')
  pNice <- 0;
  for (sess in 1:3){  
    # rough way to update beliefs about Others that can be expected in task
    # after each person has been encountered.
    if (sess > 1){
       a0 <- a0 + e*pNice;  b0 <- b0 + e*(1-pNice);  # update starting params for this 'Other'
    }
    a <- a0;     b <- b0;                         
    for (trial in 1:maxTr){
     # Record belief indicators of start of trial:
     evo[(sess-1)*maxTr+trial,1:3] <- c(a,b,a/(a+b))
     # here the participant holds beliefs about the Other beta(a,b) and has 
     # a response function beta(a/tau,b/tau)
     lp <- dbeta(resp[sess,trial],a/tau,b/tau,log=TRUE)
     sll <- sll + lp;
     a <- a + stim[sess, trial]
     b <- b + (1-stim[sess,trial])
     # For easy refernce, also record choice and likelihood density of that choice
     evo[(sess-1)*maxTr+trial,4] <- resp[sess,trial];
     evo[(sess-1)*maxTr+trial,5] <- lp; 
    }
    pNice <- a/(a+b);  # expectation for this Other
  }
  
  outp <- list(); outp[[1]] <- sll; outp[[2]] <- evo;
  
  return(outp)
}

dotest =0;
if (dotest){
# Now let's give it some data and some parameters:
rsp1 <- t(matrix(c(rep(0.5,6),rep(0.1,6),rep(0.2,6)),6,3)); rsp1
st1  <- t(matrix(c(rep(0,10),rep(1,8)),6,3)); st1
p1 <- c(2,2,0.66,0.5)

test <- L0(rsp1,st1,p1)
plot(test[[2]][,3],ylim=c(0,1),main='evolving expectation about Other',xlab='trials',ylab='expectation'); 
abline(0.5,0,col='gray'); abline(v=6.5); abline(v=12.5)
print(test[[2]])
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1.a Inference Model A Log-likelihood with full range of possible returns ~~~~~~
#     binned in 9 intervals from 0% to 100%. Used July - Sept 2019

# Define likelihood function for inference model. This gives
# the probability density of the responses given the participant's 
# parameters. Offers are col 1 of d, responses are col 2:3 . d contains
# all 3 partners x 6 trials each, in order of presentation, down the rows.
infHISIll9 <- function(ptp,d,details=0,plots=0) {
  # rem ptp is of the form c(pHI0,uHI0,pSI0,uSI0,upi,eta)
  tn = 6;  # As per Joe's experiment, 6 trials
  on = 3;  # three Others (partners)
  
  # Parameters to help describe personality and action. These are
  # not free parameters to be fitted, i.e. don't vary across pts.
  Nb = 9;      # number of bins for each personality dim.
  Na = Nb;     # number of action bins (within which density is flat)
               # i.e. proportions that may be received by the Other.
               # In general will be set to Nb. 
  
  #   Prior beliefs of the pt ** about the population **
  #   Prior beliefs about greed on its own, and malice on its own:
  PSI0 = noisyBino(ptp[3],ptp[4],Nb); PHI0 = noisyBino(ptp[1],ptp[2],Nb);
  # In this version, these baseline beliefs are considered as independent,
  # so a simple matrix multiplication gives the joint:
  PSIHI0 = PSI0 %*% t(PHI0); 
  
  if (plots){
    # Action names:
    anames = c('selfish','v.tight','tight','fair','kind','v.kind','altruistic')
  }
  
  # Now to formulate the policies. There is one individual parameter
  # here, which determines the uncertainty with which say a high-HI person will
  # indeed choose the action most fitting to their type (i.e., keep everything),
  # or will show choices more variable over options:
  upi = ptp[5];  # convenience copy of policy variability (inverse precision)
  eta = ptp[6];
  if (length(ptp) > 6){piu=ptp[7]} else {piu = 1} # arbitrary policy-inverse-uncertainty 
      # param; could set to e.g. 1/upi to reference it to self, roughly. 
  if (length(ptp) > 7){err=ptp[8]} else {err =0.02/(Nb*Nb)} # arbitrary lapse-rate-like 
      # param; Note scaling by the number of attribution states considered. 
  # Set up the map between 'attributes' and actions :
  pi = array(NA,c(Nb,Nb,Na));   # This will hold all possible policies
  # pinit, pstep, bu, mu fine-tune the map ... see below.
  pinit = 0.05; pstep= (1-2*pinit)/(Nb-1); 
  bu = 2.5; mu= -2*pstep;     # more handwavey constants to modulae u. Should have
  # bu+mu*Nb and bu+mu both valid values for u in the noisyBino. Negative values of
  # mu over-weigh high values of SI and HI, so that SI=1, HI=9 is skewed towards small
  # returns, while m=0 would be symmetric around the middle. 
  for (SI in 1:Nb){
    for (HI in 1:Nb) {
      x = noisyBino(pinit+(SI-1)*pstep, bu+mu*SI,Na) * 
        noisyBino(pinit+(HI-1)*pstep, bu+mu*HI,Na);
      pi[SI,HI,] = fliplr(x^(1/piu) / sum(x^(1/piu)))
    }
  }
  # The above probably not worth adding a free param at the indiv. level.
  
  
  
  # Run the inference. Here the posterior of one
  # trial will form the prior for the next trial, staring from the population prior beliefs PSIHI0.
  #
  
  sll = 0; 
  pri0 = PSIHI0; # prior at the very start of encounter with 'partner'.
  # Construct an output/results object, which is just the sum log lik 
  # if details==0. If not, record trial-by-trial data, attributions,
  # prob. of data given params (likelihood), log-lik, most likely attrib given
  # the parameters and data, and also a set of simulated data (incl. decision variability)
  if (details){
    llout = list(); 
    llout[[1]]=0; 
    hd <- c('pHI0','uHI0','pSI0','uSI0','upi','eta','piu','err')
    llout[[2]] = c(ptp[1:6],piu,err);  names(llout[[2]]) <- hd;
    llout[[3]] = matrix(NA,tn*on+1,10); 
    llout[[3]][,1] = c(0,1:(tn*on))
    colnames(llout[[3]]) = c('trial','ret','hi','si','lik','ll','HImode','SImode','HIsim','SIsim')
    llout[[3]][2:(1+tn*on),2:4] <- d
    # Hypothetical policy before any data seen:
    pol = pri0^(1/upi); pol = pol / sum(as.vector(pol));
    pol = (pol+err)/(1+err*length(pol))
    llout[[4]] <- array(dim=c(dim(pri0),(1+tn*on)))
    llout[[4]][,,1] <- pri0;
    names(llout) <- c('sll','param', 'evo','policy')
  }
    
  for (other in 1:3){    # counter of 3 different 'partners'
    if (other > 1){
      pri0 <- pri0*(1-eta) + post*eta;  # learning over partners
      # modifies the starting prior if this isn't the first partner.
    }
    post <- pri0; # this is the beief that will be updated with each trial
    
    # rows to be processed and convenience copies of other's actions,
    # malice and greed attributions:
    ro = ((other-1)*tn+1):(other*tn); # rows of data matrix
    as = d[ro,1];  aind = round((Na-1)*as+1) 
    hi = d[ro,2];  hind = round((Na-1)*hi+1)
    si = d[ro,3];  sind = round((Na-1)*si+1)
    
    for (t in 1:tn){  # loop
      if(plots){
        heatmap(log(post),Rowv=NA, Colv=NA, col = heat.colors(128), scale="none", 
                margins=c(4,8),asp=1,labRow=0:(Nb-1),labCol=0:(Nb-1), 
                main = paste('\n lnPost. at start of trial ',t),xlab='HI',ylab='SI')
      }
      pri = post;              # new prior is last posterior
      # In the next line, the pt. uses the pi entry as likelihood, pri as prior,
      # over the character of the partner. This post is their post. beliefs 
      post = pi[,,aind[t]] * pri; post = post / sum(as.vector(post))  # Bayes
      # Now the probability of the response, incl. the lapse-rate-like err:
      pol = post^(1/upi);  pol= pol/sum(as.vector(pol));
      pol = (pol+err)/(1+err*length(pol))
      lik = pol[sind[t],hind[t]];
      sll = sll + log(lik);         # accumulate sum log lik
      if (details){
        llout$evo[((other-1)*tn+t+1),'lik'] <- lik
        llout$evo[((other-1)*tn+t+1),'ll'] <- log(lik)
        # find mode of pol
        c = max.col(pol);  # this finds the max. col. pos. for each row
        m=c(); for (r in 1:Nb){m[r]=pol[r,c[r]]};  # retrieve the values ...
        r=which.max(m);  # ... so was to find the row with the mode of the distro.
        llout$evo[((other-1)*tn+t+1),c('HImode','SImode')] <- c(c[r]-0.5,r-0.5)/Nb
        # Now sample randomly 
        hisim <- rmdf(1,colSums(pol)) # joint=marginal*conditional, so sample first dim from the marginal ...
        sisim <- rmdf(1,pol[,hisim]/sum(pol[,hisim]))  # ... and the other from the corresponding conditional.
        # debug: barplot(colSums(pol));       barplot(pol[,hisim]); 
        llout$evo[((other-1)*tn+t+1),c('HIsim','SIsim')] <- c((hisim-0.5),(sisim-0.5))/Nb
        llout$policy[,,(other-1)*tn+t+1] <- pol
      }
    }
    if (plots) {
      heatmap(log(post),Rowv=NA, Colv=NA, col = heat.colors(128), scale="none", 
            margins=c(4,8),asp=1,labRow=0:(Nb-1),labCol=0:(Nb-1), 
            main = paste('\n lnPost. after triall ',t),xlab='HI',ylab='SI')
    }
  }
  if (details){llout$sll <- sll} else {llout <- sll}
  return(llout)
} 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1.b Inference Model A Log-likelihood ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     Only 2 possible returns, .5 and 0, as per pt instructions
#     Nov 2019 on

# Define likelihood function for inference model. This gives
# the probability density of the responses given the participant's 
# parameters. Offers are col 1 of d, responses are col 2:3 . d contains
# all 3 partners x 6 trials each, in order of presentation, down the rows.
infHISIll <- function(ptp,d,details=0,plots=0) {
  # rem ptp is of the form c(pHI0,uHI0,pSI0,uSI0,upi,eta)
  tn = 6;  # As per Joe's experiment, 6 trials
  on = 3;  # three Others (partners)
  
  # Parameters to help describe personality and action. These are
  # not free parameters to be fitted, i.e. don't vary across pts.
  Nb = 9;     # number of bins for each personality dim.
  Nf = 9;     # 'full' number of action bins (within which density is flat) for detailed sim work
  Na = 2;     # actual number of actions in Joe's expt.
  
  # i.e. proportions that may be received by the Other.
  # In general will be set to Nb. 
  
  #   Prior beliefs of the pt ** about the population **
  #   Prior beliefs about greed on its own, and malice on its own:
  PSI0 = noisyBino(ptp[3],ptp[4],Nb); PHI0 = noisyBino(ptp[1],ptp[2],Nb);
  # In this version, these baseline beliefs are considered as independent,
  # so a simple matrix multiplication gives the joint:
  PSIHI0 = PSI0 %*% t(PHI0); 
  
  if (plots){
    # Action names (was:  anames = c('selfish','v.tight','tight','fair','kind','v.kind','altruistic') )
    anames = c('selfish','fair')
  }
  
  # Now to formulate the policies. There is one individual parameter
  # here, which determines the uncertainty with which say a high-HI person will
  # indeed choose the action most fitting to their type (i.e., keep everything),
  # or will show choices more variable over options:
  upi = ptp[5];  # convenience copy of policy variability (inverse precision)
  eta = ptp[6];
  if (length(ptp) > 6){piu=ptp[7]} else {piu = 1} # arbitrary policy-inverse-uncertainty 
  # param; could set to e.g. 1/upi to reference it to self, roughly. 
  if (length(ptp) > 7){err=ptp[8]} else {err =0.02/(Nb*Nb)} # arbitrary lapse-rate-like 
  # param; Note scaling by the number of attribution states considered. 
  # Set up the map between 'attributes' and actions :
  pif = array(NA,c(Nb,Nb,Nf));      # This will hold all possible policies for full range of actions
  pi  = array(NA,c(Nb,Nb,Na));    # Possible policies for actual range of actions
  fp = c(1,1,1,0,0,0,0,0,0);      # auxiliary vector to further bin down pi, here to a two-bin vector
  f2p = t(matrix(c(fp,1-fp),Nb,Na)); 
  # To plot average policy / returns : 
  if (plots) {
    piav = array(NA,c(Nb,Nb)); # This will hold average policies
    pifav = piav;
  }
  # pinit, pstep, bu, mu fine-tune the map ... see below.
  pinit = 0.1; pstep= (1-2*pinit)/(Nb-1); 
  bu = 2.5; mu= -2*pstep;     # more handwavey constants to modulate u. Should have
  # bu+mu*Nb and bu+mu both valid values for u in the noisyBino. Negative values of
  # mu over-weigh high values of SI and HI, so that SI=1, HI=9 is skewed towards small
  # returns, while m=0 would be symmetric around the middle. 
  for (HI in 1:Nb){
    for (SI in 1:Nb) {
      x = noisyBino(pinit+(HI-1)*pstep, bu+mu*HI,Nf) * 
        noisyBino(pinit+(SI-1)*pstep, bu+mu*SI,Nf);
      pif[HI,SI,] = fliplr(x^(1/piu) / sum(x^(1/piu)))
      pi[HI,SI,]  = as.vector(f2p %*% pif[HI,SI,]);      # further bin down!
      if (plots){
        piav[HI,SI] = sum(pi[HI,SI,]*c(0,0.5))
        pifav[HI,SI] = sum(pif[HI,SI,]*(1:Nb))
      }
    }
  }
  if (plots){           # Display the average policy / returns as a heatmap
    heatmap( piav,
             Rowv=NA, Colv=NA, col = topo.colors(512), 
             scale="none", margins=c(5,8),asp=1,
             labRow=((0:(Nb-1))+0.5)*0.1,
             labCol=((0:(Nb-1))+0.5)*0.1, 
             main = paste('\n attributes vs. mean policy'),
             xlab='HI',ylab='SI\n')
  }
  
  #              ******** Run the inference ********
  # Here the posterior of one trial will form the prior for the next trial, 
  # staring from the population prior beliefs PSIHI0.
  #
  
  sll = 0; 
  pri0 = PSIHI0; # prior at the very start of encounter with 'partner'.
  # Construct an output/results object, which is just the sum log lik 
  # if details==0. If not, record trial-by-trial data, attributions,
  # prob. of data given params (likelihood), log-lik, most likely attrib given
  # the parameters and data, and also a set of simulated data (incl. decision variability)
  if (details){
    llout = list(); 
    llout[[1]]=0; 
    hd <- c('pHI0','uHI0','pSI0','uSI0','upi','eta','piu','err')
    llout[[2]] = c(ptp[1:6],piu,err);  names(llout[[2]]) <- hd;
    llout[[3]] = matrix(NA,tn*on+1,10); 
    llout[[3]][,1] = c(0,1:(tn*on))
    colnames(llout[[3]]) = c('trial','ret','hi','si','lik','ll','HImode','SImode','HIsim','SIsim')
    llout[[3]][2:(1+tn*on),2:4] <- d
    # Hypothetical (attribn. reporting) policy before any data seen:
    pol = pri0^(1/upi); pol = pol / sum(as.vector(pol));
    pol = (pol+err)/(1+err*length(pol))
    llout[[4]] <- array(dim=c(dim(pri0),(1+tn*on)))
    llout[[4]][,,1] <- pri0;
    names(llout) <- c('sll','param', 'evo','policy')
  }
  
  for (other in 1:3){    # counter of 3 different 'partners'
    if (other > 1){
      pri0 <- pri0*(1-eta) + post*eta;  # learning over partners
      # modifies the starting prior if this isn't the first partner.
    }
    post <- pri0; # this is the beief that will be updated with each trial
    
    # rows to be processed and convenience copies of other's actions,
    # malice and greed attributions:
    ro = ((other-1)*tn+1):(other*tn); # rows of data matrix
    as = d[ro,1];  aind = round((Na-1)*as+1) 
    hi = d[ro,2];  hind = round((Nb-1)*hi+1)
    si = d[ro,3];  sind = round((Nb-1)*si+1)
    
    for (t in 1:tn){  # loop
      if(plots){
        heatmap(log(post),Rowv=NA, Colv=NA, col = heat.colors(128), scale="none", 
                margins=c(4,8),asp=1,labRow=0:(Nb-1),labCol=0:(Nb-1), 
                main = paste('\n lnPost. at start of trial ',t),xlab='HI',ylab='SI')
      }
      pri = post;              # new prior is last posterior
      # In the next line, the pt. uses the pi entry as likelihood, pri as prior,
      # over the character of the partner. This post is their post. beliefs 
      post = pi[,,aind[t]] * pri; post = post / sum(as.vector(post))  # Bayes
      # Now the probability of the response, incl. the lapse-rate-like err:
      pol = post^(1/upi);  pol= pol/sum(as.vector(pol));
      pol = (pol+err)/(1+err*length(pol))
      lik = pol[sind[t],hind[t]];
      sll = sll + log(lik);         # accumulate sum log lik
      if (details){
        llout$evo[((other-1)*tn+t+1),'lik'] <- lik
        llout$evo[((other-1)*tn+t+1),'ll'] <- log(lik)
        # find mode of pol
        c = max.col(pol);  # this finds the max. col. pos. for each row
        m=c(); for (r in 1:Nb){m[r]=pol[r,c[r]]};  # retrieve the values ...
        r=which.max(m);  # ... so was to find the row with the mode of the distro.
        llout$evo[((other-1)*tn+t+1),c('HImode','SImode')] <- c(c[r]-0.5,r-0.5)/Nb
        # Now sample randomly 
        hisim <- rmdf(1,colSums(pol)) # joint=marginal*conditional, so sample first dim from the marginal ...
        sisim <- rmdf(1,pol[,hisim]/sum(pol[,hisim]))  # ... and the other from the corresponding conditional.
        # debug: barplot(colSums(pol));       barplot(pol[,hisim]); 
        llout$evo[((other-1)*tn+t+1),c('HIsim','SIsim')] <- c((hisim-0.5),(sisim-0.5))/Nb
        llout$policy[,,(other-1)*tn+t+1] <- pol
      }
    }
    if (plots) {
      heatmap(log(post),Rowv=NA, Colv=NA, col = heat.colors(128), scale="none", 
              margins=c(4,8),asp=1,labRow=0:(Nb-1),labCol=0:(Nb-1), 
              main = paste('\n lnPost. after triall ',t),xlab='HI',ylab='SI')
    }
  }
  if (details){llout$sll <- sll} else {llout <- sll}
  return(llout)
} 


#######  WRAPPERS TO FACILITATE USE OF REGULARIZATION ###############
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2.a Warpper for Inference Model A Log-likelihood ~~~~~~~~~~~~~~~~~~~

# REM               SEb, acc0max, acc0min,  eta,  gam,  wexp,  wrpe, Tresp,  sig
#      parMat <- c( 0.8,  0.8,    0.6,      0.1,  0.8,   0.3,   0.4,  0.2,   0.1)
#      if fixAcc (fixed acceptance proportions) are given, then reset eta to zero,
#      acc0min and max to NA and use fixAcc (for all ptN) in hapPol[1,expInd,ptN].
#      usu. leave fixAcc at default c(0.85,0.70,0.30,0.15)

# We'll try to use beta distribution with scaled x values, so that instead of 
# bet. 0 and 1, x obtains values between lo and hi, which is the interval of
# psychobiological interest. Can have lo > hi, in which case they are mirrored.
# dbetasc <- function(x, shape1, shape2, lo=0, hi=1, ncp=0, log=FALSE)
# e.g. for uncertaint param, may use: 
# x <- 0.02*(0:1000);  plot(x, dbetasc(x,1.2,3.6,0,25),t='l')
# for a pretty uniform probability or learning rate prior:
# x <- 0.002*(0:1000)-0.5;  y <-  dbetasc(x,1.05,1.05,0,1); plot(x,y,t='l')
# for a normal looking prior with mean 10 and sd about 108:
# x <- (0:1000)-500;  y <-  dbetasc(x,10.41,10,-500,500); plot(x,y,t='l')
# So each col. of ParM has the form ashape, bshape, lowlim, hilim., 
# This example has weakly informative priors for all except piu, which is
# essentially fixed to the value of 2, and err, which is infomative:
#          pHI0   uHI0 pSI0  uSI0   upi   eta    piu    err
# ashape   1.05   1.2  1.05   1.2   1.2   1.05   100   1.05
# bshape   1.05   3.6  1.05   3.6   3.6   1.05   100    4
# lowlim    0      0    0      0     0     0      0     0
# hilim     1     25    1     25    25     1      4     1

msLPhisi1a <- function(ParM, datAr, scbeta0=NA,check=0){

  parM <- as.vector(ParM); # in case it's inputed in another format
  parn <- length(parM)
  
  if ((scbeta0[1] < 0) && !is.na(scbeta0)){ 
    # i.e. a number, but not a valid scaled distr. param.,
    # which means 'use default, weak regularizing priors'
    scbeta0 <- matrix(c(1.05,1.05,0, 1,
                      1.2, 3.6, 0, 25,
                      1.05,1.05,0, 1,
                      1.2, 3.6, 0, 25,
                      1.2, 3.6, 0, 25,
                      1.05,1.05,0, 1,
                      1.2, 3.6, 0, 25,
                      1.05,1.05,0, 1  ), 4, 8)
    if(check){
      colnames(scbeta0) <- c('pHI0','uHI0','pSI0','uSI0','upi','eta','piu','err')
    }
  }
  
  # Cacl. the log prior for MAP purposes etc, all calc'd in short form:
  mSLPrior <- 0;
  if (length(scbeta0)>1){  # legit prior must have 24 elements or so!
      mSLPrior <- mSLPrior - sum(dbetasc( parM, 
                    scbeta0[1,1:parn],scbeta0[2,1:parn],
                    scbeta0[3,1:parn],scbeta0[4,1:parn], log=TRUE)); 
  }

  if (mSLPrior == Inf){  # If we are in an a priori prohibited parameter region
    # do not attempt to calculate the likelihood - it will be nonsense anyway.
    return(Inf); 
  } else {
      return(mSLPrior - infHISIll(ParM,datAr))
  }


} # end of msLPhisi1a

                  #######  FIGURE 8 PLOT   ###############
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

plots <- 1;
# opar <- par();

uHI0 = 1
uSI0 = 0.2
pHI0 = 0.33
pSI0 = 0.75
upi  = 2  # decision variability for participant (not [imagined] for Other)
eta  = 0.25 # learning rate to cater for order effects
# optional params:
piu  = 2  # how noisly Other's attributes map to choices
xi   = 0.01  # generic lapse rate like param
testpar = c(pHI0,uHI0,pSI0,uSI0,upi,eta)   # ,piu,xi)
parhd = c("pHI0","uHI0","pSI0","uSI0","upi","eta","piu","xi")

# Set up some toy data:
as <- c(rep(c(0.5,0),3),rep(0,6),rep(0.5,6))
toyD = matrix(c(as,
                rep(0,6),rep(0.5,6),rep(0.1,6),
                rep(0.2,6),rep(0.4,6),rep(0.67,6)),
              18,3); 
colnames(toyD) = c('ret','hi','si'); toyD

# demo
# Parameters to help describe personality and action. These are
# not free parameters to be fitted, i.e. don't vary across pts.
Nb = 9;     # number of bins for each personality dim.
Nf = 9;     # 'full' number of action bins (within which density is flat) for detailed sim work
Na = 2;     # actual number of actions in Joe's expt.

# i.e. proportions that may be received by the Other.
# In general will be set to Nb. 

#   Prior beliefs of the pt ** about the population **
#   Prior beliefs about greed on its own, and malice on its own:
PSI0 = noisyBino(testpar[3],testpar[4],Nb); PHI0 = noisyBino(testpar[1],testpar[2],Nb);
# In this version, these baseline beliefs are considered as independent,
# so a simple matrix multiplication gives the joint:
PSIHI0 = PSI0 %*% t(PHI0); 

if (plots){
  # Action names (was:  anames = c('selfish','v.tight','tight','fair','kind','v.kind','altruistic') )
  anames = c('selfish','fair')
}

# Now to formulate the policies. There is one individual parameter
# here, which determines the uncertainty with which say a high-HI person will
# indeed choose the action most fitting to their type (i.e., keep everything),
# or will show choices more variable over options:
upi = testpar[5];  # convenience copy of policy variability (inverse precision)
eta = testpar[6];
if (length(testpar) > 6){piu=testpar[7]} else {piu = 1} # arbitrary policy-inverse-uncertainty 
# param; could set to e.g. 1/upi to reference it to self, roughly. 
if (length(testpar) > 7){err=testpar[8]} else {err =0.02/(Nb*Nb)} # arbitrary lapse-rate-like 
# param; Note scaling by the number of attribution states considered. 
# Set up the map between 'attributes' and actions :
pif = array(NA,c(Nb,Nb,Nf));      # This will hold all possible policies for full range of actions
pi  = array(NA,c(Nb,Nb,Na));    # Possible policies for actual range of actions
fp = c(1,1,1,0,0,0,0,0,0);      # auxiliary vector to further bin down pi, here to a two-bin vector
f2p = t(matrix(c(fp,1-fp),Nb,Na)); 
# To plot average policy / returns : 
if (plots) {
  piav = array(NA,c(Nb,Nb)); # This will hold average policies
  pifav = piav;
}
# pinit, pstep, bu, mu fine-tune the map ... see below.
pinit = 0.1; pstep= (1-2*pinit)/(Nb-1); 
bu = 2.5; mu= -2*pstep;     # more handwavey constants to modulate u. Should have
# bu+mu*Nb and bu+mu both valid values for u in the noisyBino. Negative values of
# mu over-weigh high values of SI and HI, so that SI=1, HI=9 is skewed towards small
# returns, while m=0 would be symmetric around the middle. 
for (HI in 1:Nb){
  for (SI in 1:Nb) {
    x = noisyBino(pinit+(HI-1)*pstep, bu+mu*HI,Nf) * 
      noisyBino(pinit+(SI-1)*pstep, bu+mu*SI,Nf);
    pif[HI,SI,] = fliplr(x^(1/piu) / sum(x^(1/piu)))
    pi[HI,SI,]  = as.vector(f2p %*% pif[HI,SI,]);      # further bin down!
    piav[HI,SI] = sum(pi[HI,SI,]*c(0,0.5))
    pifav[HI,SI] = sum(pif[HI,SI,]*(1:Nb))
  }
}
if (plots){           # Display the average policy / returns as a heatmap
  heatmap( piav,
           Rowv=NA, Colv=NA, col = topo.colors(512), 
           scale="none", margins=c(5,8),asp=1,
           labRow=((0:(Nb-1))+0.5)*0.1,
           labCol=((0:(Nb-1))+0.5)*0.1, 
           main = paste('\n attributes vs. mean policy'),
           xlab='HI',ylab='SI\n')
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

                #####TESTING and MODELLING SIMULATED DATA #######

#install.packages("R.matlab")
library(matlab)    # This provides some convenient matlab/octave functions
library(R.matlab)  # Functions for handling mat files etc.
library(Hmisc)
library(ppcor)
library(tidyr); library(dplyr)

setwd("/Users/josephbarnby/Downloads/")
pfits <- read.csv("hisi04b.csv")
as <- c(rep(c(0.5,0),3),rep(0,6),rep(0.5,6))
toyd = matrix(c(as,rep(0,6),rep(0.5,6),rep(0.1,6),rep(0.2,6),rep(0.4,6),rep(0.67,6)),18,3); 
colnames(toyd) <- c("ret", "HI", "SI")
load("alldetc.RData")
alld[[1]][1]

testfull <- list()
simulated <- list()
testfullparams <- list()
simulatedfull <- list()
test.full.fits <- list()

set.seed(1)

names(pfits)

for (i in 1:1763) {
  
  test                <- infHISIll(as.numeric(pfits[i,2:7]),  
                                  matrix(unlist(alld[[i]][1]), 
                                         ncol = 3, byrow = F), 
                                  details=1,plots=0)
  
  testfull[[i]]       <- test[[3]][2:19,]
  testfullparams[[i]] <- test[[2]]
  
  simulated[[i]]      <- cbind(testfull[[i]][1:18,c(2,9:10)], 
                               alld[[i]]$nontask[1,c(2:4,6:8)],
                               matrix(unlist(testfullparams[[i]]), ncol = 8, byrow = F))
  
  test.full.fits[[i]] <- cbind(testfull[[i]][1:18,c(1, 5:6)], 
                               alld[[i]]$nontask[1,c(2:4,6:8)])
  
}

  #modelled fits

set.seed(2)

simulated[[1]] # one participant of simulated data

simulated.df <- do.call(rbind, simulated)
simulated.df <- as.data.frame(simulated.df)
colnames(simulated.df)[10:17] <- c("pHI0", "uHI0", "pSI0", "uSI0", "upi", "eta", "piu", "err")

test.full.fits.df <- as.data.frame(do.call(rbind, test.full.fits))
test.full.fits.df$Dictator <- c(rep("Fair", 6), rep("Partially Fair", 6), rep("Unfair", 6))
test.full.fits.df2 <- cbind(test.full.fits.df[c(1:3, 10)], simulated.df)
test.full.fits.df2$Level <- car::recode(test.full.fits.df2$GPTSTotal,
                                        "32:36 = 'Low';
                                37:44 = 'Medium';
                                45:61 = 'High';
                                62:101.9 = 'Very High';
                                102:170 = 'Clinical'")
test.full.fits.df2$Level <- factor(test.full.fits.df2$Level,levels = c("Low", "Medium", "High", "Very High", "Clinical"), ordered = T)
names(test.full.fits.df2)

summary.fits <- as.data.frame(format(Rmisc::summarySE(test.full.fits.df2, measurevar = "ll",
                                                      groupvars=c("Dictator", "trial")), digits =3 ))
summary.fits2 <- as.data.frame(format(Rmisc::summarySE(test.full.fits.df2, measurevar = "ll",
                                                       groupvars=c("Dictator", "Level")), digits =3 ))
summary.fits[4:7] <- sapply(summary.fits[4:7], as.numeric)

#script for figure 1:

ll1 <- ggplot(test.full.fits.df2, aes(trial, ll, color = Dictator))+
  geom_jitter(alpha = 0.01)+
  stat_summary(fun.data = "mean_cl_boot")+
  scale_x_continuous(breaks = c(1:18))+
  geom_text(data = Rmisc::summarySE(test.full.fits.df2, measurevar = "ll",
                                    groupvars=c("Dictator", "trial")),
            aes(label= format(ll, digits = 3), y = ll + 1, color = Dictator))+
  xlab("Trial")+
  ylab("Log Likelihood")+
  ggtitle("A | Loglikelihood for each trial")+
  geom_hline(yintercept = -4.394, alpha = 0.5)+
  scale_color_brewer(palette = "Dark2")+
  sjPlot::theme_blank()+
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.title = element_blank())

ll2 <- ggplot(test.full.fits.df2, aes(ll, col = Dictator))+
  geom_density()+
  xlab("Log Likelihood")+
  ylab("Density")+
  ggtitle("B | Density plot of loglikelihood for each trial")+
  scale_color_brewer(palette = "Dark2")+
  geom_vline(xintercept = -4.394, alpha = 0.5)+
  geom_vline(aes(xintercept = mean(ll), color = Dictator), alpha = 0.5)+
  facet_wrap(.~trial, nrow = 3, ncol = 6)+
  sjPlot::theme_blank()+
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.title = element_blank())

ll3 <- ggplot(test.full.fits.df2, aes(Level, ll, color = Level))+
  #geom_smooth(method = "lm", show.legend = FALSE)+
  stat_summary(fun.data = "mean_cl_boot")+
  geom_jitter(alpha = 0.01, show.legend = FALSE)+
  ylab("Log Likelihood")+
  geom_hline(yintercept = -4.394, alpha = 0.5)+
  scale_color_brewer(palette = "Reds")+
  geom_text(data = Rmisc::summarySE(test.full.fits.df2, measurevar = "ll",
                                    groupvars=c("Level")),
            aes(label= format(ll, digits = 3), y = ll + 1, color = Level))+
  stat_cor(label.y.npc = 'bottom')+
  sjPlot::theme_blank()+
  ggtitle("C | Loglikelihood for each division of GPTS score")+
  theme(legend.position = c(0.7, 0.15),
        legend.direction = "horizontal")+
  theme(axis.title.x = element_blank())+
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.title = element_blank())

ll4 <-  ggplot(test.full.fits.df2, aes(ll, col = Level))+
  geom_density()+
  xlab("Log Likelihood")+
  ylab("Density")+
  scale_color_brewer(palette = "Reds")+
  geom_vline(xintercept = -4.394, alpha = 0.5)+
  ggtitle("D | Density plot of loglikelihood for each trial")+
  facet_wrap(.~trial, nrow = 3, ncol = 6)+
  sjPlot::theme_blank()+
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.title = element_blank())

ll5 <- ggplot(test.full.fits.df2, aes(GPTSTotal, ll))+
  geom_smooth(method = "lm", show.legend = FALSE)+
  stat_summary(fun.data = "mean_cl_boot")+
  geom_jitter(alpha = 0.01, show.legend = FALSE)+
  xlab("GPTSTotal")+
  ylab("Log Likelihood")+
  geom_hline(yintercept = -4.394, alpha = 0.5)+
  scale_x_continuous(breaks = c(32, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160))+
  scale_color_brewer(palette = "Dark2")+
  stat_cor(label.y.npc = 'bottom')+
  ggtitle("E | Loglikelihood value for each gradation of GPTS score")+
  sjPlot::theme_blank()+
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.title = element_blank())

summary(test.full.fits.df2$ll)
mean(test.full.fits.df2$lik)

lllist1 <- ggpubr::ggarrange( ll1, ll2, labels = c("", ""), common.legend = T, nrow = 1) 
lllist2 <- ggpubr::ggarrange( ll3, ll4, labels = c("", ""), common.legend = T, nrow = 1) 
ggarrange(lllist1, lllist2, ll5, nrow = 3, labels = c("", "", ""))

# log likelihood by trials and gpts divisions - regression models

summary(lme4::lmer(ll ~ trial + Dictator + GPTSTotal + (1|ID), data = test.full.fits.df2))
confint(lme4::lmer(ll ~ trial + Dictator + GPTSTotal + (1|ID), data = test.full.fits.df2))
effectsize::effectsize(lme4::lmer(ll ~ trial + Dictator + GPTSTotal + (1|ID), data = test.full.fits.df2))

#wrangle main simmed data for use in analyses 

set.seed(2)

simulated[[i]] # one participant of simulated data

simulated.df <- do.call(rbind, simulated)
simulated.df <- as.data.frame(simulated.df)
colnames(simulated.df)[10:17] <- c("pHI0", "uHI0", "pSI0", "uSI0", "upi", "eta", "piu", "err")
simulated.df.unfair <- simulated.df[which(simulated.df$Receiver=="PlayerStingy"),]
simulated.df.pfair <- simulated.df[which(simulated.df$Receiver=="PlayerEqual"),]
simulated.df.fair <- simulated.df[which(simulated.df$Receiver=="PlayerFair"),]

simulated.df.unfair$Dictator <- c(rep("Unfair", 6), rep("PFair", 6), rep("Fair", 6))
simulated.df.pfair$Dictator <- c(rep("PFair", 6), rep("Fair", 6), rep("Unfair", 6))
simulated.df.fair$Dictator <- c(rep("Fair", 6), rep("Unfair", 6), rep("PFair", 6))

unfairsims <- rbind(simulated.df.unfair[which(simulated.df.unfair$Dictator == "Unfair"),],
                    simulated.df.pfair[which(simulated.df.pfair$Dictator == "Unfair"),],
                    simulated.df.fair[which(simulated.df.fair$Dictator == "Unfair"),])
fairsims <- rbind(simulated.df.unfair[which(simulated.df.unfair$Dictator == "Fair"),],
                  simulated.df.pfair[which(simulated.df.pfair$Dictator == "Fair"),],
                  simulated.df.fair[which(simulated.df.fair$Dictator == "Fair"),])
pfairsims <- rbind(simulated.df.unfair[which(simulated.df.unfair$Dictator == "PFair"),],
                   simulated.df.pfair[which(simulated.df.pfair$Dictator == "PFair"),],
                   simulated.df.fair[which(simulated.df.fair$Dictator == "PFair"),])

unfairsims$TrialHI <- rep(c(1,2,3,4,5,6))
unfairsims$TrialSI <- rep(c(1,2,3,4,5,6))
fairsims$TrialHI <- rep(c(1,2,3,4,5,6))
fairsims$TrialSI <- rep(c(1,2,3,4,5,6))
pfairsims$TrialHI <- rep(c(1,2,3,4,5,6))
pfairsims$TrialSI <- rep(c(1,2,3,4,5,6))

#tag by dictator 

unfairsims %>%
  group_by(ID) %>%
  mutate(meanHI = mean(HIsim),
         meanSI = mean(SIsim)) -> unfairsims.calc
quantile(simulated.df$GPTSTotal)
unfairsims.calc$Level <- as.numeric(unfairsims.calc$GPTSTotal)
unfairsims.calc$Level <- car::recode(unfairsims.calc$Level,
                                     "32:36 = 'Low';
                                37:44 = 'Medium';
                                45:61 = 'High';
                                62:101.9 = 'Very High';
                                102:159 = 'Clinical'")
fairsims %>%
  group_by(ID) %>%
  mutate(meanHI = mean(HIsim),
         meanSI = mean(SIsim)) -> fairsims.calc
fairsims.calc$Level <- as.numeric(fairsims.calc$GPTSTotal)
fairsims.calc$Level <- car::recode(fairsims.calc$Level,
                                   "32:36 = 'Low';
                                37:44 = 'Medium';
                                45:61 = 'High';
                                62:101.9 = 'Very High';
                                102:159 = 'Clinical'")
pfairsims %>%
  group_by(ID) %>%
  mutate(meanHI = mean(HIsim),
         meanSI = mean(SIsim)) -> pfairsims.calc
pfairsims.calc$Level <- as.numeric(pfairsims.calc$GPTSTotal)
pfairsims.calc$Level <- car::recode(pfairsims.calc$Level,
                                    "32:36 = 'Low';
                                37:44 = 'Medium';
                                45:61 = 'High';
                                62:101.9 = 'Very High';
                                102:159 = 'Clinical'")
#combined data

combsims <- rbind(fairsims.calc, pfairsims.calc, unfairsims.calc)
combsims$Level <- factor(combsims$Level,levels = c("Low", "Medium", "High", "Very High", "Clinical"), ordered = T)

# simulated data analysis

names(combsims)
combsims$OrdinalValueHI <- factor(car::recode(combsims$meanHI,
                                              "0.00 : 0.2 = '1';
                            0.2000000001 : 0.4 = '2';
                            0.4000000001 : 0.6 = '3';
                            0.6000000001 : 0.8 = '4'; 
                            0.8000000001 : 1.00 = '5'"), levels = c("1", "2", "3", "4", "5"), ordered = T)
combsims$OrdinalValueSI <- factor(car::recode(combsims$meanSI,
                                              "0.00 : 0.2 = '1';
                            0.2000000001 : 0.4 = '2';
                            0.4000000001 : 0.6 = '3';
                            0.6000000001 : 0.8 = '4'; 
                            0.8000000001 : 1.00 = '5'"), levels = c("1", "2", "3", "4", "5"), ordered = T)
combsims$Age <- as.numeric(combsims$Age)
ordered(combsims$Dictator)

#generated models

model.HI.sims <- combsims[c(4, 5, 7, 8, 18, 24)]
model.SI.sims <- combsims[c(4, 5, 7, 8, 18, 25)]

model.HI.sims$Age <- as.numeric(model.HI.sims$Age)
model.HI.sims$GPTSTotal <- scale(model.HI.sims$GPTSTotal)
model.HI.sims$Dictator <- factor(ordered(model.HI.sims$Dictator))
model.HI.sims$Receiver <- factor(model.HI.sims$Receiver, levels = c("PlayerFair", "PlayerEqual", "PlayerStingy"), ordered = T)
model.SI.sims$Age <- as.numeric (model.SI.sims$Age)
model.SI.sims$GPTSTotal <- scale(model.SI.sims$GPTSTotal)
model.SI.sims$Dictator <- factor(ordered(model.SI.sims$Dictator))
model.SI.sims$Receiver <- factor(model.SI.sims$Receiver, levels = c("PlayerFair", "PlayerEqual", "PlayerStingy"), ordered = T)

ordinalUF.gen.HI <- ordinal::clm(OrdinalValueHI ~ GPTSTotal + Age + Receiver,
                                 data = model.HI.sims[model.HI.sims$Dictator=="Unfair",], na.action =  na.fail)
ordinalF.gen.HI  <- ordinal::clm(OrdinalValueHI ~ GPTSTotal + Age + Receiver,
                                 data = model.HI.sims[model.HI.sims$Dictator=="Fair",], na.action =  na.fail)
ordinalPF.gen.HI <- ordinal::clm(OrdinalValueHI ~ GPTSTotal + Age + Receiver,
                                 data = model.HI.sims[model.HI.sims$Dictator=="PFair",], na.action =  na.fail)

effectsize::effectsize(ordinalUF.gen.HI)

effectsize::effectsize(ordinalF.gen.HI)

effectsize::effectsize(ordinalPF.gen.HI)

ordinalUF.gen.SI <- ordinal::clm(OrdinalValueSI ~ GPTSTotal + Age + Receiver,
                                 data = model.SI.sims[model.SI.sims$Dictator=="Unfair",], na.action =  na.fail)
ordinalF.gen.SI  <- ordinal::clm(OrdinalValueSI ~ GPTSTotal + Age + Receiver,
                                 data = model.SI.sims[model.SI.sims$Dictator=="Fair",], na.action =  na.fail)
ordinalPF.gen.SI <- ordinal::clm(OrdinalValueSI ~ GPTSTotal + Age + Receiver,
                                 data = model.SI.sims[model.SI.sims$Dictator=="PFair",], na.action =  na.fail)

effectsize::effectsize(ordinalUF.gen.SI)

effectsize::effectsize(ordinalF.gen.SI)

effectsize::effectsize(ordinalPF.gen.SI)

combsims$Receiver <- factor(combsims$Receiver, levels = c("PlayerFair", "PlayerEqual", "PlayerStingy"), ordered = T)

###Attributional global models

#globalmodel.gen.HI <- ordinal::clmm(OrdinalValueHI ~ scale(GPTSTotal) + scale(Age) + ordered(Dictator) + Receiver + (1|ID),
#                                 data = combsims, na.action =  na.fail)
#globalmodel.gen.SI <- ordinal::clmm(OrdinalValueSI ~ scale(GPTSTotal) + scale(Age) + ordered(Dictator) + Receiver + (1|ID),
#                                    data = combsims, na.action =  na.fail)

globalmodel.gen.HI <- lme4::lmer(scale(meanHI) ~ scale(GPTSTotal) + scale(Age) + ordered(Dictator) + Receiver + (1|ID),
                                 data = combsims, na.action =  na.fail)
globalmodel.gen.SI <- lme4::lmer(scale(meanSI) ~ scale(GPTSTotal) + scale(Age) + ordered(Dictator) + Receiver + (1|ID),
                                 data = combsims, na.action =  na.fail)

globalmodelHI.dredge <- MuMIn::dredge(globalmodel.gen.HI, REML = F, trace = 2)
globalmodelHI.models<-MuMIn::get.models(globalmodelHI.dredge, subset=delta<2)
globalmodelHIa<-MuMIn::model.avg(globalmodelHI.models, adjusted=FALSE, revised.var=TRUE)

summary(globalmodel.gen.HI) #only one model supplied
confint(globalmodel.gen.HI)

globalmodelSI.dredge <- MuMIn::dredge(globalmodel.gen.SI, REML = F, trace = 2)
globalmodelSI.models<-MuMIn::get.models(globalmodelSI.dredge, subset=delta<2)
globalmodelSIa<-MuMIn::model.avg(globalmodelSI.models, adjusted=FALSE, revised.var=TRUE)

summary(globalmodelSIa)
confint(globalmodelSIa)

# tests for generating new data

# script for figure 2:

sumrepReg <- ggplot(gather(combsims, "Type", "Value", 21:22), aes(GPTSTotal, Value, color = Type))+ #clinical groups (from 101.9 - 159) dip in their HI estimates as GPTS value increases
  geom_smooth(method = "lm")+
  coord_cartesian(ylim = c(0,1))+
  ylab("Attribution")+
  facet_wrap(~Dictator)+
  scale_x_continuous(breaks = c(32, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160))+
  stat_cor(label.y.npc = 'bottom')+
  cowplot::theme_cowplot()+
  labs(title = "E | Mean Simulated Attributions by GPTS Score")+
  theme(legend.position = c(0.04, 0.86), legend.direction = 'horizontal')+
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.title = element_blank())

density.HI <- ggplot(combsims[which(combsims$Dictator!="PFair"),], aes(HIsim, colour = Level))+
  geom_density()+
  xlab("Harmful Intent Simulations")+
  facet_wrap(Dictator~TrialHI, nrow = 2, ncol = 6)+
  scale_color_brewer(palette = "Reds")+
  sjPlot::theme_blank()+
  labs(title = "B | Density of Simulated Harmful Intent Attributions")+
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.title = element_blank())
density.SI <- ggplot(combsims[which(combsims$Dictator!="PFair"),], aes(SIsim, colour = Level))+
  geom_density()+
  facet_wrap(Dictator~TrialHI, nrow = 2, ncol = 6)+
  xlab("Self Interest Simulations")+
  scale_color_brewer(palette = "Blues")+
  sjPlot::theme_blank()+
  labs(title = "D | Density of Simulated Self Interest Attributions")+
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.title = element_blank())

sumrepdatHI <- Rmisc::summarySE(combsims, measurevar = "HIsim",
                                groupvars=c("TrialHI", "Dictator","Level"))
sumrepdatSI <- Rmisc::summarySE(combsims, measurevar = "SIsim",
                                groupvars=c("TrialSI", "Dictator","Level"))

sumrepplotHI <- 
  ggplot(              combsims   [combsims$Dictator!="PFair",],    aes(x = TrialHI, y = HIsim, group = Level)) +
  geom_line(data =     sumrepdatHI[sumrepdatHI$Dictator!="PFair",], aes(x = TrialHI, y = HIsim, colour =Level, group =  Level), linetype = 3)+
  geom_point(data =    sumrepdatHI[sumrepdatHI$Dictator!="PFair",], aes(x = TrialHI, y = HIsim, group = Level, colour = Level), shape = 18) +
  geom_errorbar(data = sumrepdatHI[sumrepdatHI$Dictator!="PFair",], aes(x = TrialHI, y = HIsim, group = Level, colour = Level, ymin = HIsim-se, ymax = HIsim+se), width = .05)+
  ylab("HI Attributions (scale 0-1)")+
  labs(color = "Level")+
  scale_x_continuous(breaks = c(1,2,3,4,5,6))+
  scale_color_brewer(palette = "Reds")+
  facet_wrap(.~Dictator)+
  cowplot::theme_cowplot()+
  labs(title = "A | Simulated Harmful Intent Attributions for Each Trial")+
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.title = element_blank())
sumrepplotSI <- 
  ggplot(              combsims   [combsims$Dictator!="PFair",],    aes(x = TrialSI, y = SIsim, group = Level)) +
  geom_line(data =     sumrepdatSI[sumrepdatSI$Dictator!="PFair",], aes(x = TrialSI, y = SIsim, colour =Level, group =  Level), linetype = 3)+
  geom_point(data =    sumrepdatSI[sumrepdatSI$Dictator!="PFair",], aes(x = TrialSI, y = SIsim, group = Level, colour = Level), shape = 18) +
  geom_errorbar(data = sumrepdatSI[sumrepdatSI$Dictator!="PFair",], aes(x = TrialSI, y = SIsim, group = Level, colour = Level, ymin = SIsim-se, ymax = SIsim+se), width = .05)+
  ylab("SI Attributions (scale 0-1)")+
  labs(color = "Level")+
  scale_color_brewer(palette = "Blues")+
  scale_x_continuous(breaks = c(1,2,3,4,5,6))+
  facet_wrap(.~Dictator)+
  cowplot::theme_cowplot()+
  labs(title = "C | Simulated Self Interest Attributions for Each Trial")+
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.title = element_blank())

HIcomb <- ggpubr::ggarrange(sumrepplotHI, density.HI, ncol = 2, common.legend = T, labels = c("", ""))
SIcomb <- ggpubr::ggarrange(sumrepplotSI, density.SI, ncol = 2, common.legend = T, labels = c("", ""))
ggpubr::ggarrange(HIcomb, SIcomb, sumrepReg, nrow = 3, labels = c("", "", ""))

# tests for parameters with paranoia

names(combsims)
quantile(combsims$GPTSTotal, prob = seq(0, 1, length = 11), type = 5)
combsims$POrdinal <- car::recode(combsims$GPTSTotal,
                                 "0:33 = '1';
                                33:34 = '2';
                                35:37 = '3';
                                38:40 = '4';
                                41:44 = '5';
                                45:49 = '6';
                                50:56 = '7';
                                57:67 = '8';
                                68:83 = '9';
                                84:170 = '10'")

colnames(combsims)[c(11, 13, 14)] <- c("Uncertainty of HI (uHI0)", "Uncertainty of SI (uSI0)", "Uncertainty Overall (uΠ)")

# figure 3: 

g1 <- ggplot(gather(combsims, "Parameter", "Value", c(11,13, 14)), aes(meanHI, Value, color = Parameter))+ #clinical groups (from 101.9 - 159) dip in their HI estimates as GPTS value increases
    geom_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1)+
    geom_point(alpha = 1/250)+
    ylab("Uncertainty of Partner Policies")+
    xlab("Mean Harmful Intent Attribution Over 18 Trials")+
    coord_cartesian(ylim = c(0,10))+
    #facet_wrap(.~Dictator)+
    ggpubr::stat_cor(aes(color = Parameter), label.y.npc = "middle", method = "spearman")+
    scale_color_brewer(palette = "Set1")+
    cowplot::theme_cowplot()+
    theme(legend.position = c(0.55, 0.75))

g2 <-  ggplot(gather(combsims, "Parameter", "Value", c(11,13, 14)), aes(meanSI, Value, color = Parameter))+ #clinical groups (from 101.9 - 159) dip in their HI estimates as GPTS value increases
    geom_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1)+
    geom_point(alpha = 1/250)+
    ylab("Uncertainty of Partner Policies")+
    xlab("Mean Self-Interest Attributions Over 18 Trials")+
    ggpubr::stat_cor(aes(color = Parameter), label.y.npc = "middle", method = "spearman")+
    coord_cartesian(ylim = c(0,10))+
    #facet_wrap(.~Dictator)+
    scale_color_brewer(palette = "Set1")+
    cowplot::theme_cowplot()+
    theme(legend.position = c(0.55, 0.75))

g3 <-  ggplot(gather(combsims, "Parameter", "Value", c(11,13, 14)), aes(GPTSTotal, Value, color = Parameter))+ #clinical groups (from 101.9 - 159) dip in their HI estimates as GPTS value increases
    #geom_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1)+
    geom_smooth(method = "lm", size = 1, alpha = 0.1)+
    geom_point(alpha = 1/250)+
    ylab("Uncertainty of Partner Policies")+
    xlab("Paranoia (GPTS Total)")+
    ggpubr::stat_cor(aes(color = Parameter), label.y.npc = "middle", method = "spearman")+
    coord_cartesian(ylim = c(0,10))+
    #facet_wrap(.~Dictator)+
    scale_color_brewer(palette = "Set1")+
    cowplot::theme_cowplot()+
  theme(legend.position = c(0.55, 0.75))

g4 <- ggplot(gather(combsims, "Attribution", "Value", c(21,22)), aes(Value, eta, color = Attribution))+ #clinical groups (from 101.9 - 159) dip in their HI estimates as GPTS value increases
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1)+ ##   ETA VALUES BY MEAN HI AND SI SCORES
  geom_point(alpha = 1/250)+
  ylab("Learning Rate (η)")+
  xlab("Mean Attribution Over 18 Trials")+
  coord_cartesian(ylim = c(0,1))+
  #facet_wrap(.~Dictator)+
  ggpubr::stat_cor(aes(color = Attribution), label.y.npc = "top", method = "spearman")+
  scale_color_brewer(palette = "Set1")+
  cowplot::theme_cowplot()+
  theme(legend.position = c(0.75, 0.75))


(g1 + g2) / (g3 + g4)


# script for Figure 9:

names(combsims)

melt.combsims <- melt(combsims, measure.vars = c("pHI0","Uncertainty of HI (uHI0)", "pSI0",                    
                                            "Uncertainty of SI (uSI0)", "Uncertainty Overall (uΠ)", "eta"))

ggarrange(
  ggplot(melt.combsims[melt.combsims$variable==c("pHI0", "pSI0"),], aes(value, colour = variable))+
    geom_density()+
    sjPlot::theme_blank()+
    theme(legend.position = c(0.75, 0.5))+
    scale_colour_brewer(palette = "Dark2"),
  ggplot(melt.combsims[melt.combsims$variable==c("Uncertainty of HI (uHI0)", "Uncertainty of SI (uSI0)", "Uncertainty Overall (uΠ)"),], aes(value, colour = variable))+
    geom_density()+
    sjPlot::theme_blank()+
    theme(legend.position = c(0.75, 0.5))+
    scale_colour_brewer(palette = "Dark2"),
  ggplot(melt.combsims[melt.combsims$variable==c("eta"),], aes(value, colour = variable))+
    geom_density()+
    sjPlot::theme_blank()+
    theme(legend.position = c(0.75, 0.5))+
    scale_colour_brewer(palette = "Dark2"))

# create df for latent statistical analysis

combsims.modelling <- combsims
combsims.modelling[c(7, 10:14, 21, 22)] <- sapply(combsims.modelling[c(7, 10:14, 21, 22)], scale)
combsims.modelling <- as.data.frame(combsims.modelling)
str(combsims.modelling)
combsims.modelling$Age <- as.numeric(combsims.modelling$Age)
combsims.modelling$Dictator <- ordered(combsims.modelling$Dictator)
combsims.modelling[combsims.modelling$Sex=="Male",]$Sex <- 1
combsims.modelling[combsims.modelling$Sex=="Female",]$Sex <- 2
head(combsims.modelling, 10)
combsims.modelling$HIsimO <- ordered(car::recode(combsims.modelling$HIsim, 
                                                 "0:0.200 =      1;
                                   0.2001: 0.400 = 2;
                                   0.4001: 0.600 = 3;
                                   0.6001: 0.800 = 4;
                                   0.8001: 1.00 = 5", as.factor = T))
combsims.modelling$SIsimO <- ordered(car::recode(combsims.modelling$SIsim, 
                                                 "0:0.200 =      1;
                                   0.2001: 0.400 = 2;
                                   0.4001: 0.600 = 3;
                                   0.6001: 0.800 = 4;
                                   0.8001: 1.00 = 5", as.factor = T))

# TO CONFIRM PARAMETRIC MIXED MODEL RESULTS WITH ORDINAL MODELS

#combsims.modelling %>% 
#  mutate(uHIOrd=cut(`Uncertainty of HI`, breaks=c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf), 
#                      labels=c('1', '2', '3', '4', '5')),
#         uSIOrd=cut(`Uncertainty of SI`, breaks=c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf), 
#                    labels=c('1', '2', '3', '4', '5')),
#         etaOrd=cut(`eta`, breaks=c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf), 
#                    labels=c('1', '2', '3', '4', '5')),
#         upiOrd=cut(`upi`, breaks=c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf), 
#                    labels=c('1', '2', '3', '4', '5')),) -> combsims.o

# modelling latent parameters using glms

colnames(combsims.modelling)[c(11, 13, 14)] <- c("uHI0", "uSI0", "uΠ")

uncertHI0.model <- glm(scale(uHI0) ~ GPTSTotal + Age + Sex #+ scale(meanHI)
                       , 
            data = combsims.modelling, na.action =  na.fail)
uncertHI0.global <- MuMIn::model.avg(MuMIn::get.models(MuMIn::dredge(uncertHI0.model, REML = F), subset = delta<2))
summary(uncertHI0.global)
confint(uncertHI0.global)
effectsize::effectsize(uncertHI0.global)

uncertSI0.model <- glm(scale(uSI0) ~ GPTSTotal + Age + Sex #+ scale(meanSI)
                       , 
            data = combsims.modelling, na.action =  na.fail)
uncertSI0.global <- MuMIn::model.avg(MuMIn::get.models(MuMIn::dredge(uncertSI0.model, REML = F), subset = delta<2))
summary(uncertSI0.global)
confint(uncertSI0.global)
effectsize::effectsize((uncertSI0.global))

upi.model1 <- glm(scale(uΠ) ~ GPTSTotal + Age + Sex #+ scale(meanHI)
                  ,
                       data = combsims.modelling, na.action =  na.fail)
upi.global1 <- MuMIn::model.avg(MuMIn::get.models(MuMIn::dredge(upi.model1, REML = F), subset = delta<2))
summary(upi.global1)
confint(upi.global1)
effectsize::effectsize(upi.global1)

eta.model1 <- glm(scale(eta) ~ GPTSTotal + Age + Sex #+ scale(meanHI)
                  , 
            data = combsims.modelling, na.action =  na.fail)
eta.global1 <- MuMIn::model.avg(MuMIn::get.models(MuMIn::dredge(eta.model1, REML = F), subset = delta<2))
summary(eta.global1)
confint(eta.global1)
effectsize::effectsize(eta.global1)

#pHI0.model <- glm(scale(pHI0) ~ GPTSTotal + Age + Sex #+ scale(meanHI)
#                       , 
#                       data = combsims.modelling, na.action =  na.fail)
#pHI0.global <- MuMIn::model.avg(MuMIn::get.models(MuMIn::dredge(pHI0.model, REML = F), subset = delta<2))
#summary(pHI0.model)
#confint(pHI0.model)
#
#pSI0.model <- glm(scale(pSI0) ~ GPTSTotal + Age + Sex #+ scale(meanHI)
#                  , 
#                  data = combsims.modelling, na.action =  na.fail)
#pSI0.global <- MuMIn::model.avg(MuMIn::get.models(MuMIn::dredge(pSI0.model, REML = F), subset = delta<2)) #only one model supplied
#summary(pSI0.model) 
#confint(pSI0.model)

#modelling simulated HI and SI attributions with parameters

uHI.model.trial <- lme4::lmer(scale(meanHI) ~ GPTSTotal + Age + Sex + uHI0 + eta + uΠ + (1|ID),
                                data = combsims.modelling, na.action =  na.fail)
#effectsize::effectsize(uHI.model.trial)

uSI.model.trial <- lme4::lmer(scale(meanSI) ~ GPTSTotal + Age + Sex + uSI0 + eta + uΠ + (1|ID),
                                data = combsims.modelling, na.action =  na.fail)
#effectsize::effectsize(uSI.model.trial)

uHI0.global <- MuMIn::model.avg(MuMIn::get.models(MuMIn::dredge(uHI.model.trial, REML = F), subset = delta<2))
summary(uHI0.global)
confint(uHI0.global)
MuMIn::importance(uHI0.global)
effectsize::effectsize(uHI0.global)

uSI0.global <- MuMIn::model.avg(MuMIn::get.models(MuMIn::dredge(uSI.model.trial, REML = F), subset = delta<2))
summary(uSI0.global)
confint(uSI0.global)
MuMIn::importance(uSI0.global)
effectsize::effectsize(uSI0.global)

#correlations and densities of parameters#

library(ggcorrplot)
library(GGally)

ggplot(gather(combsims, "Parameter", "Value", 10:15), aes(Value, group = Level, color = Level)) + 
  geom_density() + 
  facet_wrap(~Parameter, scales="free") +
  theme_blank()+
  scale_color_brewer(palette = "Reds")

names(combsims)

colnames(combsims)[c(11, 13)] <- c("uHI0", "uSI0")


s.corr <- round(cor(combsims[c(7,10:15)], method = "spearman"), 2)
head(s.corr)
p.corr <- round(cor(combsims[c(7,10:15)]), 2)
head(p.corr)

s.mat <- cor_pmat(combsims[c(7, 10:15)], method = "spearman")
head(s.mat[, 1:4])
p.mat <- cor_pmat(combsims[c(7, 10:15)])
head(p.mat[, 1:4])

p.corrplot <- ggcorrplot(p.corr, hc.order = TRUE, type = "lower",
           outline.col = "white",
           lab = T,
           ggtheme = sjPlot::theme_blank,
           colors = c("#6D9EC1", "white", "#E46726"),
           p.mat = p.mat, insig = "pch") + theme(legend.position = c(0.15, 0.70))
s.corrplot <- ggcorrplot(s.corr, hc.order = TRUE, type = "lower",
                         outline.col = "white",
                         lab = T,
                         ggtheme = sjPlot::theme_blank,
                         colors = c("#6D9EC1", "white", "#E46726"),
                         p.mat = p.mat)

s.corrplot + theme(legend.position = c(0.2, 0.75))

## Parameter clusters for simulations

den.hisi04b_1 <- ggplot(simulated.df, aes(pHI0, pSI0))+
              geom_point()+
              geom_density_2d()+
              theme_blank()
den.hisi04b_2 <- ggplot(simulated.df, aes(uHI0, uSI0))+
              ylim(0, 10)+
              xlim(0, 10)+
              geom_point()+
              geom_density_2d()+
  theme_blank()
den.hisi04b_3 <- ggplot(simulated.df, aes(pHI0, uHI0))+
              ylim(0, 10)+
              geom_point()+
              geom_density_2d()+
              theme_blank()
den.hisi04b_4 <- ggplot(simulated.df, aes(pSI0, uSI0))+
              ylim(0, 10)+
              geom_point()+
              geom_density_2d()+
              theme_blank()
den.hisi04b_5 <- ggplot(simulated.df, aes(pHI0, upi))+
  ylim(0, 10)+
  geom_point()+
  geom_density_2d()+
  theme_blank()
den.hisi04b_6 <- ggplot(simulated.df, aes(pHI0, eta))+
  geom_point()+
  geom_density_2d()+
  theme_blank()
(den.hisi04b_1|den.hisi04b_2)/
  (den.hisi04b_3|den.hisi04b_4)/
  (den.hisi04b_5|den.hisi04b_6)

combsims[10:15] <- sapply(combsims[10:15], scale)

##bayes bits
#setwd('/Users/josephbarnby/Dropbox/PhD/Helsinki Summit/Data/')
#p1 <- read.csv("GPTSscores.csv")
#p2 <- read.csv("dataBASELINE.csv")
#p  <- rbind(p1[1:33], p2[1:33])
#p$socref <- rowSums(p[2:17])
#p$persec <- rowSums(p[18:33])
#p$convic <- rowSums(p[c(5, 7, 10, 11, 18:24, 32)])
#p$preocc <- rowSums(p[c(4, 9, 14, 15, 16, 25, 28, 29, 33)])
#p$distre <- rowSums(p[c(2, 3, 6, 8, 12, 13, 17, 26, 27, 30, 31)])
#
#colnames(p)[1] <- 'ID'
#
#combsims.modelling.x <- join(combsims.modelling, p[,c(1, 34:38)], by = "ID")
#
#names(combsims.modelling.x)[29:33] <- c("Social Reference", "Persecution", "Conviction", "Preoccupation", "Distress")
#subscale.cor <- pcor(combsims.modelling.x[c(10:15, 7,29:33)], method = "spearman")
#subscale.mat <- cor_pmat(combsims.modelling.x[c(10:15, 7,29:33)], method = "spearman")
#ggcorrplot(subscale.cor$estimate, hc.order = TRUE, type = "lower",
#           outline.col = "white",
#           lab = T,
#           ggtheme = sjPlot::theme_blank,
#           colors = c("#6D9EC1", "white", "#E46726"),
#           p.mat = subscale.mat)+ theme(legend.position = c(0.2, 0.75))
#
##correlations between simulated and real data
#
#names(combsims.modelling)
#
#realdat1 <- read.csv("AnxHelsinkiData_w_decisions.csv")
#realdat2 <- read.csv("HelsinkiData.csv")
#colnames(realdat1)[1] <- "ID"
#realdat  <- rbind(realdat1[c(1, 4:15, 18:29, 32:43)], realdat2[c(1, 4:15, 18:29, 32:43)])
#names(realdat)
#
#realdat$HImean_real <- rowSums(realdat[c(3, 5, 7, 9, 11, 13, 
#                                15, 17, 19, 21, 23, 25,
#                                27, 29, 31, 33, 35, 37)]/18)
#realdat$SImean_real <- rowSums(realdat[c(2, 4, 6, 8, 10, 12, 
#                                14, 16, 18, 20, 22, 24, 
#                                26, 28, 30, 32, 34, 36)]/18)
#
#combsims.modelling.y <- join(combsims.modelling, realdat[c(1, 38, 39)], by = "ID")
#names(combsims.modelling.y)
#
#ggarrange(
#  
#ggplot(na.omit(combsims.modelling.y), aes(meanHI, scale(HImean_real), color = Level))+
#  geom_point(alpha = 0.05, show.legend = F)+
#  geom_smooth(method = "lm", aes(fill = Level))+
#  stat_cor(method = "spearman", show.legend = F)+
#  
#  xlab("Simulated")+
#  ylab("Real")+
#  
#  scale_color_brewer(palette = "Reds")+
#  scale_fill_brewer(palette = "Reds")+
#  
#  bbplot::bbc_style()+
#  theme(legend.position = c(0.8, 0.2),
#        panel.grid.major.y = element_blank(),
#        axis.title = element_text(size = 14)),
#
#ggplot(na.omit(combsims.modelling.y), aes(meanSI, scale(SImean_real), color = Level))+
#  geom_point(alpha = 0.05, show.legend = F)+
#  geom_smooth(method = "lm", aes(fill = Level))+
#  stat_cor(method = "spearman", show.legend = F)+
#  
#  scale_color_brewer(palette = "Blues")+
#  scale_fill_brewer(palette = "Blues")+
#  
#  xlab("Simulated")+
#  ylab("Real")+
#  
#  bbplot::bbc_style()+
#  theme(legend.position = c(0.8, 0.2),
#        panel.grid.major.y = element_blank(),
#        axis.title = element_text(size = 14))
#)
#### simulations 

#simulated params to check for model fitting

#simp <- matrix(NA,nrow=200,ncol=6)
#
#simp[,1] <- runif(200,0.01,0.1) ; hist(simp[,2],xlim=c(0,1)) #pHI
#simp[,2] <- runif(200,0.9,1.1)  ; hist(simp[,3],xlim=c(0,10)) #uHI
#simp[,3] <- runif(200,0.87,0.93); hist(simp[,1],xlim=c(0,1)) #pSI
#simp[,4] <- runif(200,1,1.2)    ; hist(simp[,4],xlim=c(0,10)) #uSI
#simp[,5] <- runif(200,1.3,1.7)  ; hist(simp[,5],xlim=c(0,10)) #upi
#simp[,6] <- runif(200,0.03,0.07); hist(simp[,6],xlim=c(0,1)) #eta
#
#colnames(simp) <- c("pHI", "uHI", "pSI", "uSI", "upi", "eta")
#head(simp)
#
## Simple log-likelihood demos
#
#test0 <- infHISIll(simp,toyD,details=1,plots=1);
#
##simD <- test0$evo[2:(tn*on+1),c('ret','HImode','S.Imode')]
#simD <- test0$evo[2:(tn*on+1),c('ret','HIsim','SIsim')]
#print('Simulated data:')
#print(round(t(simD[1:6,]),3));   # Parter 1
#print(round(t(simD[7:12,]),3));  # Parter 2
#print(round(t(simD[13:18,]),3)); # Parter 3
#simll <- infHISIll(simp,simD,details = 1,plots = 0);
#
#plot(0:18, simll$evo[,'HImode'],t='b',ylim=c(0,1),col='green4',
#     ylab = 'green: HI     gold: SI',xlab='trial',
#     main=paste('Synthetic data with correct param.,\n sll=',round(simll[[1]],3)))
#lines(0:18,simll$evo[,'SImode'],t='b',ylim=c(0,1),col='gold4')
#lines(0:18,simll$evo[,'HIsim'],t='b',ylim=c(0,1),col='palegreen3')
#lines(0:18,simll$evo[,'SIsim'],t='b',ylim=c(0,1),col='yellow3')
#
#abline(v=6.5,col='gray'); abline(v=12.5,col='gray');
#
#testfull <- list()
#simulated <- list()
#testfullparams <- list()
#simulatedfull <- list()
#test.full.fits <- list()
#
#
#for (i in 1:200) {
#  
#  simll                <- infHISIll(as.numeric(simp[i,1:6]),  
#                                   matrix(unlist(alld[[i]][1]), 
#                                          ncol = 3, byrow = F), 
#                                   details=1,plots=0)
#  
#  testfull[[i]]       <- simll[[3]][2:19,]
#  testfullparams[[i]] <- simll[[2]]
#  
#  simulated[[i]]      <- cbind(testfull[[i]][1:18,c(2,9:10)], 
#                               alld[[i]]$nontask[1,c(2:4,6:8)],
#                               matrix(unlist(testfullparams[[i]]), ncol = 8, byrow = F))
#  
#  test.full.fits[[i]] <- cbind(testfull[[i]][1:18,c(1, 5:6)], 
#                               alld[[i]]$nontask[1,c(2:4,6:8)]
#                               )
#  
#}
#
#simulated.df <- do.call(rbind, simulated)
#simulated.df <- as.data.frame(simulated.df)
#colnames(simulated.df)[10:17] <- c("pHI0", "uHI0", "pSI0", "uSI0", "upi", "eta", "piu", "err")
#colnames(simulated.df)
#
#cor_auto(simulated.df[10:15])
