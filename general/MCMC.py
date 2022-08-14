import numpy as np
from numpy import random
from scipy import stats

class Proposal():
    def __init__(self,dist_f):
        self.dist_f = dist_f

    def draw_proposal(self,old_par):
        forward_dist = self.dist_f(old_par)
        new_par = forward_dist.rvs()
        backward_dist = self.dist_f(new_par)
        log_ratio = backward_dist.logpdf(old_par) - forward_dist.logpdf(new_par)
        return new_par,log_ratio

class MultiProposal(Proposal):
    def __init__(self,dist_fs,which=None):
        if which == None:
            self.which = [i for i in range(len(dist_fs))]
        else:
            self.which = which
        self.proposals = [Proposal(dist_f) for dist_f in dist_fs]

    def draw_proposal(self,old_pars):
        new_pars = list(old_pars)
        joint_log_ratio = 0
        for i in range(len(self.which)):
            par_ind = self.which[i]
            old_par = old_pars[par_ind]
            new_par,log_ratio = self.proposals[i].draw_proposal(old_par)
            new_pars[par_ind] = new_par
            joint_log_ratio += log_ratio
        return new_pars,joint_log_ratio
            
def norm_f(prop_std):
    def prop_dist(par):
        return stats.norm(loc=par,scale=prop_std)
    return prop_dist

def truncnorm_f(prop_std,a=-np.inf,b=np.inf):
    def prop_dist(par):
        return stats.truncnorm(a,b,loc=par,scale=prop_std)
    return prop_dist

def lognorm_f(prop_std):
    def prop_dist(par):
        return stats.lognorm(s=prop_std,scale=par)
    return prop_dist

def gibbs_update(updatef,log_likf,**kwargs):
    def update(priors,pars,*args):
        new_pars = updatef(priors,pars,**kwargs)
        return new_pars,1,log_likf(priors,new_pars)
    return update

def MH_update(proposal,log_likf,**kwargs):

    def update(priors,pars,loglik_old):
        prop,logratio = proposal.draw_proposal(pars,**kwargs)
        loglik_new = log_likf(priors,prop)
        print("~~~~~~")
        print("new, old, prop ratio")
        print(loglik_new)
        print(loglik_old)
        print(logratio)
        log_p = loglik_new - loglik_old + logratio
        if log_p > 0:
            accept_prob = 1
        else:
            accept_prob = np.exp(log_p)
        rand = random.uniform()
        if rand <= accept_prob:
            return prop,1,loglik_new
        else:
            return pars,0,loglik_old
    return update

def MCMC_cycle(updates):
    def cycle(priors,pars,loglik_old):
        successes = np.zeros(len(updates))
        for i in range(len(updates)):
            upd_f = updates[i]
            pars,success,loglik_old = upd_f(priors,pars,loglik_old)
            successes[i] += success
        return pars,successes,loglik_old
    return cycle

def MCMC_init(init_f,priors,loglik_f,**kwargs):
    def init():
        pars = init_f(priors,loglik_f,**kwargs)
        loglik = loglik_f(priors,pars)
        return pars,loglik
    return init
    
def MCMC_chain(num_warmup,num_sample,thin_rate,priors,init_f,cycle_f):
    pars,loglik_old = init_f()
    pars,success,loglik_old = cycle_f(priors,pars,loglik_old)
    for i in range(num_warmup-1):
        pars,success,loglik_old = cycle_f(priors,pars,loglik_old)

    samples = [pars for i in range(num_sample)]
    successes = np.zeros(len(success))
    for i in range(num_sample):
        for j in range(thin_rate-1):
            pars,success,loglik_old = cycle_f(priors,pars,loglik_old)
            successes += success
        pars,success,loglik_old = cycle_f(priors,pars,loglik_old)
        samples[i] = pars
        successes += success

    return samples,successes
