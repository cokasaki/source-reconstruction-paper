from fenics import *
import numpy as np
from numpy import random
from scipy import sparse
from scipy.special import kv,gamma,loggamma
from sksparse import cholmod as ch

from util import *
###
#
# The Matern Covariance
#
##

def d2(x1,x2):
    x1 = np.array(x1)
    x2 = np.array(x2)
    return np.sqrt(np.sum((x1-x2)**2))

def dists(pts):
    return np.array([[d2(x1,x2) for x1 in pts] for x2 in pts])

# Neumann BC is default for fenics
def make_K1K2LinvL(d,V):
    u = TrialFunction(V)
    v = TestFunction(V)
    stiff1 = dot(u,v)*dx
    stiff2 = inner(grad(u), grad(v))*dx
    mass = inner(u,v)*dx

    K1 = assemble_mat(stiff1)
    K2 = assemble_mat(stiff2)
    L = assemble_mat(mass)
    invL = sparse.spdiags(1/L.sum(axis=0),0,L.shape[0],L.shape[1])
    return K1,K2,L,invL
mat_dict = {}
def make_matern(alpha,kappa,d,V,sigma2,threshold=None,space_name=None,downscale=False):
    if space_name is not None and space_name in mat_dict.keys():
        K1,K2,L,invL = mat_dict[space_name]
    else:
        K1,K2,L,invL = make_K1K2LinvL(d,V)
        mat_dict[space_name] = K1,K2,L,invL

    K = kappa**2*K1 + K2
    if downscale:
        K = K1 + K2/kappa**2
    maternMat = K

    i = 1
    while alpha > i:
        maternMat = maternMat @ invL @ K
        i += 1
    sigma2_anal = analytical_matern_var(alpha-d/2,kappa,d)

    max_entry = np.max(maternMat[maternMat.nonzero()])
    min_entry = np.min(np.abs([maternMat[maternMat.nonzero()]]))

    maternMat *= sigma2_anal/sigma2
    max_entry = np.max(np.abs(maternMat[maternMat.nonzero()]))
    min_entry = np.min(np.abs([maternMat[maternMat.nonzero()]]))
    
    I,J,V = [],[],[]
    maternMat = maternMat.tocoo()
    for i,j,v in zip(maternMat.row,maternMat.col,maternMat.data):
        if threshold is None or np.abs(v/max_entry) >= threshold:
            I += [i]
            J += [j]
            V += [v]
    maternMat = sparse.coo_matrix((V,(I,J))).tocsr()

    return maternMat

def make_matern_factors(alpha,kappa,d,V,sigma2):
    u = TrialFunction(V)
    v = TestFunction(V)
    stiffness = kappa**2*u*v*dx + inner(grad(u), grad(v))*dx
    mass = inner(u,v)*dx

    # Assemble precision matrices
    K = assemble(stiffness)
    K = as_backend_type(K).mat()
    K = sparse.csr_matrix(K.getValuesCSR()[::-1],shape=K.size)

    L = assemble(mass)
    L = as_backend_type(L).mat()
    (m,n) = L.size
    L = sparse.csr_matrix(L.getValuesCSR()[::-1],shape=L.size)
    invL = sparse.spdiags(1/L.sum(axis=0),0,m,n) # the approx from Lindgren

    sigma2_anal = analytical_matern_var(alpha-d/2,kappa,d)
    if alpha % 2 == 1:
        return [K,invL]*(alpha//2) + [ch.cholesky(K.tocsc()).L()] + [sparse.eye(K.shape[0])*np.sqrt(sigma2_anal/sigma2)]
    if alpha % 2 == 0:
        sqrt_invL = sparse.diags(np.sqrt(np.diag(invL.todense())))
        return [K] + [invL,K]*(alpha//2 - 1) + [sqrt_invL] + [sparse.eye(K.shape[0])*np.sqrt(sigma2_anal/sigma2)]

def analytical_matern_cov(nu,kappa,d,sigma2,dist,threshold=1e-10):
    '''returns the analytical value of the matern covariance at distance dist
    and with parameters nu and kappa'''
    assert nu > 0, "nu must be greater than zero"
    assert kappa > 0, "kappa must be greater than zero"

    if np.isinf(kappa):
        return sigma2*np.eye(dist.shape[0])

    logsigma2 = np.log(sigma2)

    logconst = np.log(2)+logsigma2+nu*np.log(kappa/2.0)-loggamma(nu)
    const = np.exp(logconst)
    dist_exp = (dist**nu)*kv(nu,kappa*dist)
    S = const*dist_exp*(const*dist_exp*sigma2 > threshold)

    # if you're having troubles check if varratio (sigma2) is inf
    # print(sigma2)

    it = np.nditer(dist,flags=['multi_index'])
    for x in it:
        mi = it.multi_index
        if dist[mi] == 0:
            S[mi] = sigma2

    return S

    
def analytical_matern_var(nu,kappa,d):
    if nu==0:
        return 1/(4*np.pi)
    else:
        logsigma2 = loggamma(nu)-loggamma(nu+d/2)-(d/2)*np.log(4*np.pi)-(2*nu)*np.log(kappa)
        sigma2 = np.exp(logsigma2)
        return sigma2

def sim_matern(alpha,kappa,d,V,sigma2):
    u = TrialFunction(V)
    v = TestFunction(V)
    stiffness = kappa**2*u*v*dx + inner(grad(u), grad(v))*dx
    mass = inner(u,v)*dx

    # Assemble precision matrices
    K = assemble(stiffness)
    K = as_backend_type(K).mat()
    K = sparse.csr_matrix(K.getValuesCSR()[::-1],shape=K.size)

    L = assemble(mass)
    L = as_backend_type(L).mat()
    (m,n) = L.size
    L = sparse.csr_matrix(L.getValuesCSR()[::-1],shape=L.size)
    invL = sparse.spdiags(1/L.sum(axis=0),0,m,n) # the approx from Lindgren

    sigma2_anal = analytical_matern_var(alpha-d/2,kappa,d)
    sim = random.normal(size=m,scale=np.sqrt(sigma2/sigma2_anal))
    if alpha % 2 == 1:
        Kchol = ch.cholesky(K.tocsc())
        sim = Kchol.solve_A(sim)
    if alpha % 2 == 0:      
        sqrt_invL = sparse.diags(np.sqrt(np.diag(invL.todense())))
        sim = splinalg.spsolve(K,sqrt_invL @ sim)
    for i in range((alpha-1)//2):
        sim = splinalg.spsolve(K,invL @ sim)
    return sim

def sim_matern_time(alpha,kappa,tau,d,V,sigma2,T,dt=None,max_num=None,threshold=1e-10,tol=1e-5):
    assert alpha%2==0,"alpha must be even"

    u = TrialFunction(V)
    v = TestFunction(V)
    stiffness = kappa**2*u*v*dx + inner(grad(u), grad(v))*dx
    mass = inner(u,v)*dx

    # Assemble precision matrices
    K = assemble(stiffness)
    K = as_backend_type(K).mat()
    K = sparse.csr_matrix(K.getValuesCSR()[::-1],shape=K.size)

    L = assemble(mass)
    L = as_backend_type(L).mat()
    (m,n) = L.size
    L = sparse.csr_matrix(L.getValuesCSR()[::-1],shape=L.size)
    invL = sparse.spdiags(1/L.sum(axis=0),0,m,n) # the approx from Lindgren 
    sqrt_invL = sparse.diags(np.sqrt(np.diag(invL.todense())))

    invLK = invL @ K
    diag = invLK[0,0]
    if dt is not None:
        print(1/(2*diag))
        print(dt)
        assert 1/(2*diag) > dt*(1-tol),"fixed step size is too large"
    else:
        dt = 1/(2*diag)
    N = int(np.ceil(T/dt))

    print(N," time steps required")
    if max_num is not None:
        assert N <= max_num,"too many time steps required"

    M = sparse.eye(m) - dt*invLK
    sigma2_anal = analytical_matern_var(alpha-d/2,kappa,d)
    sim = random.normal(size=(N,m),scale=np.sqrt(sigma2/sigma2_anal))
    simold = np.zeros((N,m))
    simnew = np.zeros((N,m))
    for i in range(N):
        simold[i,:] = sqrt_invL @ sim[i,:]
    for i in range(alpha//2):
        simnew[0,:] = simold[0,:]
        for i in range(1,N):
            simnew[i,:] = M @ simnew[i-1,:] + simold[i,:]
        simold = dt*simnew
        simnew = np.zeros((N,m))
    return simold

###
#
# The ARMA Precisions
#
###

def make_ar1_ch(tau,dt,sigma,num_time_points,neumann=False):
    rho = 1-dt/tau
    Rinv = sparse.diags([np.ones(num_time_points),-np.ones(num_time_points-1)*rho],offsets=[0,-1]).tocsc()
    if neumann:
        Rinv[0,0] = dt/tau
    return Rinv/sigma

def make_ar1(phi,sigma,N):
    Q = sparse.diags([np.ones(N)*(1+phi**2),-np.ones(N-1)*phi,-np.ones(N-1)*phi],offsets=[0,1,-1])
    Q = Q - sparse.coo_matrix(([phi**2,phi**2],([0,N-1],[0,N-1])))
    return Q/sigma**2

def ar1logdet(phi,sigma,N):
    return np.log(1-phi**2)-2*N*np.log(sigma)
