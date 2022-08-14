import numpy as np
from numpy import random
from scipy import sparse
from scipy.sparse import linalg as splinalg
from sksparse import cholmod as ch
from fenics import *

# simulate from the cholesky decomposition of a precision matrix
def sim_from_Qchol(Qchol,n):
    sim = random.normal(size=n)
    sim = Qchol.apply_Pt(Qchol.solve_Lt(sim,use_LDLt_decomposition=False))
    return sim

# calculate a simulataneous confidence interval from simulations
def simultaneous_confint_from_sims(sims,p=0.05,mean=None,std=None):
    if mean is None:
        mean = np.mean(sims,axis=0)
    if std is None:
        std = np.std(sims,axis=0)
    scaled = (sims-mean)/std
    max_scaled = np.amax(np.abs(scaled),axis=0)
    q = np.quantile(max_scaled,1-p)
    return q*std

# calculate a simulataneous confidence interval from a covariance matrix
def simultaneous_confint(cov,num_sim,p=0.05):
    dim = cov.shape[0]
    sim = random.multivariate_normal(np.zeros(dim),cov,size=num_sim)
    return simultaneous_confint_from_sims(sim,p)

# calculate a simultaneous confidence interval from a sparse precision matrix
def simultaneous_confint_prec(prec,num_sim,L=None,p=0.05):
    dim = prec.shape[0]
    R = ch.cholesky(prec.tocsc())
    if L is None:
        L = sparse.eye(dim)
    std = np.sqrt(np.diag((L @ splinalg.inv(prec.tocsc()) @ L.transpose()).todense()))
    sim = random.multivariate_normal(np.zeros(dim),np.eye(dim),size=num_sim)
    max_scaled = []
    for i in range(num_sim):
        corr_sim = L @ R.solve_A(sim[i,:]) # shouldnt this be solve_Lt!!!
        corr_sim_scaled = corr_sim/std
        max_scaled += [np.max(np.abs(corr_sim_scaled))]
    q = np.quantile(max_scaled,1-p)
    return q*std

# calculate a simultaneous confidence interval from a sparse precision matrix
def simultaneous_confint_prec_2(fac,num_sim,p=0.05):
    dim = prec.shape[0]

    # how many extreme simulations do we need to keep
    num_to_keep = np.ceil(p*num_sim)

    sims = np.abs(random.multivariate_normal(size=(num_to_keep,dim**2)))
    extreme_sims = np.zeros((num_to_keep,dim**2))
    for i in range(num_to_keep):
        extreme_sims[i,:] = fac.solve_L(sims[i,:])

    per_sim_max = np.amax(extreme_sims,axis=1)
    least_extreme_sim = np.argmin(per_sim_max)
    least_extreme_dev = per_sim_max[least_extreme_sim]
    discard_max = np.zeros(dim**2)

    sims = np.abs(random.multivariate_normal(size=(num_sim-num_to_keep,dim**2)))
    for i in range(num_sim-num_to_keep):
        corr_sim = fac.solve_L(sims[i,:])
        if np.max(corr_sim) < least_extreme_dev:
            discard_max = np.amax([discard_max,corr_sim])
        else:
            discard_max = np.amax([discard_max,extreme_sims[least_extreme_sim]])

            per_sim_max[least_extreme_sim] = np.max(corr_sim)
            extreme_sims[least_extreme_sim] = corr_sim
            least_extreme_sim = np.argmin(per_sim_max)
            least_extreme_dev = per_sim_max[least_extreme_sim]

    return discard_max

# simulates a multivariate normal with sparse inverse cholesky factors [R_i]
def mvn_invchol(chol_L,dim,size=1):
    sim = random.multivariate_normal(np.zeros(dim),np.eye(dim),size=size)
    for i in range(len(chol_L)):
        for j in range(size):
            sim[j,:] = splinalg.spsolve(chol_L[-1-i].tocsc(), sim[j,:])
    return sim

class SparseMVN():
    def __init__(self,mean=None,prec=None):
        if mean is not None:
            self.mean=mean
        if prec is not None:
            self.prec=prec.tocsc()
            self.prec_fac = ch.cholesky(prec.tocsc())


# Define functions to transfer vector to function and vice versa
def vec_to_func(vec,space,tol=1e-5):
    f = Function(space)
    f.vector()[:] = vec.copy()
    return f

def func_to_vec(func,tol=1e-5):
    vec = func.vector()[:]
    return vec

def spatialvec_to_func_slow(pts,vec,space,tol=1e-5):
    f = Function(space)

    dofmap = space.dofmap()
    nvertices = space.ufl_cell().num_vertices()

    # Set up a vertex_2_dof list
    dim = len(vec[0])
    for i in range(dim):
        indices = [dofmap.tabulate_entity_dofs(0, j)[i] for j in range(nvertices)]

        vertex_2_dof = dict()
        [vertex_2_dof.update(dict(vd for vd in zip(cell.entities(0),
                                                dofmap.cell_dofs(cell.index())[indices])))
                                for cell in cells(space.mesh())]

        # Get the vertex coordinates
        X = space.mesh().coordinates()

        for j in range(len(pts)):
            pt = pts[j]

            # Find the matching vertex (if it exists)
            vertex_idx = np.where(np.apply_along_axis(lambda x: np.allclose(x,pt),1,X))
            if not vertex_idx:
                print('No matching vertex!')
            else:
                vertex_idx = vertex_idx[0][0]
                dof_idx = vertex_2_dof[vertex_idx]
            
            f.vector()[dof_idx] = vec[j][i]

    return f

def spatialvec_to_func(pts,vec,space,tol=1e-5):
    f = Function(space)

    dofmap = space.dofmap()
    nvertices = space.ufl_cell().num_vertices()

    # Set up a vertex_2_dof list
    dim = len(vec[0])
    for i in range(dim):
        indices = [dofmap.tabulate_entity_dofs(0, j)[i] for j in range(nvertices)]

        for j in range(len(pts)):
            pt = pts[j]

            # Find the matching vertex (if it exists)
            cell_id = space.mesh().bounding_box_tree().compute_first_entity_collision(Point(pt))
            cell = Cell(space.mesh(),cell_id)
            coords = space.element().tabulate_dof_coordinates(cell)
            dofs = dofmap.cell_dofs(cell.index())
            which_index = np.where([np.allclose(coords[i],pt) for i in indices])[0]
            if which_index.size != 1:
                print("Number of matching vertices: "+str(which_index.size))
                print("Coordinates of non-matching point: "+str(pt))
            else:
                index = indices[which_index[0]]
                dof = dofs[index]
            
            f.vector()[dof] = vec[j][i]

    return f


# if something breaks or we need to speed it up we can try this
# try this https://fenicsproject.discourse.group/t/pass-function-to-c-expression/1081 

# a class to represent a numerically estimated function that can be used as a PDE coefficient                 
# function is assumed to be constant over cells                                          
class DataExpression(UserExpression):
    def __init__(self,data=None,**kwargs):
        if data is not None:
            self._data = data
        super(DataExpression,self).__init__(**kwargs)
    def eval_cell(self,values,x,ufc_cell):
        values[0] = self._data[ufc_cell.index]
    def value_shape(self):
        return (1,)

# a class to represent a finite-element discretized function that can be used as a PDE coefficient                    
# function is assumed to be constant over cells, so each cell is assigned the value at its center of mass     
class FunctionExpression(UserExpression):
    def __init__(self,function=None,mesh=None,**kwargs):
        if function is not None:
            self._function = function
        if mesh is not None:
            self._mesh = mesh
        super(FunctionExpression,self).__init__(**kwargs)
    def eval_cell(self,values,x,ufc_cell):
        vertex_indices = self._mesh.cells()[ufc_cell.index]
        vertex_coordinates = map(lambda i: self._mesh.coordinates()[i],vertex_indices)
        center = np.mean(list(vertex_coordinates))
        values[0] = self._function(center)
    def value_shape(self):
        return (1,)

# a function to assemble a stiffness or mass matrix
def assemble_mat(equation):
    M = assemble(equation)
    M = as_backend_type(M).mat()
    M = sparse.csr_matrix(M.getValuesCSR()[::-1],shape=M.size)
    return M

# a function to assemble an offset vector
def assemble_vec(equation):
    v = assemble(equation)
    v = as_backend_type(v).vec()
    v = v.getArray()
    return v

###
#
# Function to apply dirichlet boundary conditions
# they need to be applied in different ways to
# mu, K, L, and Q
#
###

# a function to split up a symmetric matrix into three parts
#   - the part referring to the boundary of a PDE
#   - the part referring to the interior of the PDE
#   - the part coupling the boundary to the interior
def split_mat_internal_external(M,BC):
    
    BC_dict = BC.get_boundary_values()
    BC_indices = list(BC_dict.keys())
    
    num_pts = M.shape[0]
    internal_indices = [i for i in range(num_pts) if i not in BC_indices]

    M11 = M[internal_indices,:][:,internal_indices]
    M12 = M[internal_indices,:][:,BC_indices]
    M22 = M[BC_indices,:][:,BC_indices]
    return M11,M22,M12 # interior, exterior, cross

# a function to split up a rectangular matrix into two parts
#   - the part referring to the boundary of a PDE
#   - the part referring to the interior of the PDE
def split_mat_rows_internal_external(M,BC):
    
    BC_dict = BC.get_boundary_values()
    BC_indices = list(BC_dict.keys())
    
    num_pts = M.shape[0]
    internal_indices = [i for i in range(num_pts) if i not in BC_indices]

    M1 = M[internal_indices,:]
    M2 = M[BC_indices,:]
    return M1,M2 # interior, exterior

# a function to split up a rectangular matrix into two parts
#   - the part referring to the boundary of a PDE
#   - the part referring to the interior of the PDE
def split_mat_cols_internal_external(M,BC):
    
    BC_dict = BC.get_boundary_values()
    BC_indices = list(BC_dict.keys())
    
    num_pts = M.shape[0]
    internal_indices = [i for i in range(num_pts) if i not in BC_indices]

    M1 = M[internal_indices,:]
    M2 = M[BC_indices,:]
    return M1,M2 # interior, exterior

# a function to split up a vector into two parts
#  - the part referring to the boudary of the PDE
#  - the part referring to the interior of the PDE
def split_vec_internal_external(v,BC):

    BC_dict = BC.get_boundary_values()
    BC_indices = list(BC_dict.keys())

    num_pts = v.shape[0]
    internal_indices = [i for i in range(num_pts) if i not in BC_indices]

    v_int = v[internal_indices]
    v_ext = v[BC_indices]
    return v_int,v_ext

# given a source problem with u=soln and f=source
# and given a prior on the source with mean mu_source
# and given the FEM matrices K (stiffness) and L (mass)
# calculate the prior mean of u
def calc_prior_solninterior_mean_with_BC(K,L,mu_source,BC,const=None):
    K11,K22,K12 = split_mat_internal_external(K,BC)
    L11,L22,L12 = split_mat_internal_external(L,BC)
    mu_int,mu_ext = split_vec_internal_external(mu_source,BC)

    BC_dict = BC.get_boundary_values()
    BC_values = np.array(list(BC_dict.values()))

    term1 = L11 @ mu_int + L12 @ mu_ext
    term2 = K12 @ BC_values
    prior_mean = splinalg.spsolve(K11,term1-term2)

    if const is not None:
        c_int,c_ext = split_vec_internal_external(const,BC)
        prior_mean -= splinalg.spsolve(K11,c_int)
    return prior_mean

# given a source problem with u=soln and f=source
# and given a prior on the source with precision Q_source
# and given the FEM matrices K (stiffness) and L (mass)
# calculate the prior precision of u
def calc_prior_solninterior_precision_with_BC(K,Ltilde,Q_source,BC):
    K11,K22,K12 = split_mat_internal_external(K,BC)
    L11,L22,L12 = split_mat_internal_external(Ltilde,BC)
    Q11,Q22,Q12 = split_mat_internal_external(Q_source,BC)

    L11inv = sparse.diags(1/L11.diagonal())
    Qu = K11.transpose() @ L11inv @ Q11 @ L11inv @ K11
    return Qu

# given a vector for the interior of the solution and the BC
# calculate the whole solution
def calc_soln_vec(interior,BC):

    BC_dict = BC.get_boundary_values()
    BC_keys = list(BC_dict.keys())
    BC_keys.sort()
    solution = interior.copy()
    for ind in BC_keys:
        solution = np.insert(solution,ind,BC_dict[ind])
    return solution

# given a source problem with u=soln and f=source
# and given a solution vector u_1
# and given the FEM matrices K (stiffness) and L (mass)
# and given the mass approximations Ltilde and Lbar_inv
# and given the source mean and precision matrices
# calculate the posterior distribution of f_boundary 
def calc_sourceboundary_dist_with_BC(u1,K,L,Ltilde,Lbar_inv,mu_source,Q_source,BC):
    K11,K22,K12 = split_mat_internal_external(K,BC)
    L11,L22,L12 = split_mat_internal_external(L,BC)
    Ltilde11,Ltilde22,Ltilde12 = split_mat_internal_external(Ltilde,BC)
    Q11,Q22,Q12 = split_mat_internal_external(Q_source,BC)

    S22inv = Q22 - Q21 @ splinalg.spsolve(Q11,Q12)
    S11inv = Q11 - Q12 & splinalg.spsolve(Q22,Q21)

    # calculate posterior precision
    Qpost = S22inv + Ltilde12.transpose() @ Lbar_inv @ Q11 @ Lbar_inv @ Ltilde12

    # calculate posterior mean
    BC_dict = BC.get_boundary_values()
    BC_inds = np.array(list(BC_dict.keys()))
    BC_values = np.array(list(BC_dict.values()))
    term1a = L11 @ mu_source - K12 @ BC_values - K11 @ u1
    term1b = L12.transpose() @ Lbar_inv @ S11inv @ Lbar_inv @ term1a
    term2 = S22inv @ mu_source[BC_inds]
    mu_post = - splinalg.spsolve(Qpost,term1b - term2)

    return mu_post,Qpost

def calc_sourceinterior_with_BC(u1,f2,K,L,BC):
    K11,K22,K12 = split_mat_internal_external(K,BC)
    L11,L22,L12 = split_mat_internal_external(L,BC)
    u2 = np.array(list(BC.get_boundary_values().values()))
    
    f1 = splinalg.spsolve(L11,K11 @ u1 + K12 @ u2 - L12 @ f2)
    return f1

    
