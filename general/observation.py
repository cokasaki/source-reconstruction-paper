from fenics import *
#from mshr import *
from scipy import sparse
import numpy as np

def get_covered_cells(mesh_small,mesh_big):
    # the multimesh object allows us to calculate...
    intersector = MultiMesh()
    intersector.add(mesh_big)
    intersector.add(mesh_small)

    # which cells of the mesh are covered
    covered = intersector.covered_cells(0)
    return covered

def get_cut_cells(mesh_small,mesh_big):
    # the multimesh object allows us to calculate...
    intersector = MultiMesh()
    intersector.add(mesh_big)
    intersector.add(mesh_small)

    # which cells of the mesh are cut
    cut = intersector.cut_cells(0)
    return cut

def refine_cut_cells(mesh_small,mesh_big,fem):
    return None

def make_obs_mat_areas(polygons,fem,mesh,space,resolution=10,threshold=5e-3):
    V = []
    I = []
    J = []

    for i in range(len(polygons)):
        verts = polygons[i]
        pts = [Point(vert) for vert in verts]
        area = Polygon(pts)
        area_mesh = generate_mesh(domain,resolution)

        # we need to include all the cells that are covered...
        covered = get_covered_cells(area_mesh,mesh_big)
        v = [map(lambda cell_id: Cell(mesh,cell_id).volume(),covered)]
        
        for k in range(len(covered)):
            if v[k] >= threshold:
                I += [i]
                J += [covered[k]]
                V += [v[k]]

                # and which are partially covered by the area in question
                # however, I do not yet know how to calculate cell intersection volumes
                # cut = intersector.cut_cells(0)

    O = sparse.coo_matrix((V,(I,J)),shape=(len(pts),mesh.num_vertices())).tocsr()
    return O

def make_obs_mat_pts(pts,fem,mesh,space,threshold=5e-3):
    V,I,J = [],[],[]

    tree = mesh.bounding_box_tree()
    for j in range(len(pts)):
        pt = pts[j]
        mesh_pt = Point(pt)

        # find the cell containing this point
        cell_id = tree.compute_first_entity_collision(mesh_pt)
        cell = Cell(mesh,cell_id)

        # calculate the coordinates of the vertices of the cell
        coordinate_dofs = cell.get_vertex_coordinates()

        # array of vertex indices associates with this cell
        i = list(space.dofmap().cell_dofs(cell.index()))

        # basis function value for each vertex at this point
        v = [fem.evaluate_basis(k,pt,coordinate_dofs,cell.orientation())[0] for k in range(fem.space_dimension())]

        for k in range(fem.space_dimension()):
            if v[k] >= threshold:
                I += [i[k]]
                J += [j]
                V += [v[k]]

    O = sparse.coo_matrix((V,(I,J)),shape=(space.dim(),len(pts))).tocsr()
    return O

def make_obs_mat_points(points,fem,mesh,space,threshold=5e-3):
    V,I,J = [],[],[]

    tree = mesh.bounding_box_tree()
    for j in range(len(points)):
        mesh_pt = points[j]

        # find the cell containing this point
        cell_id = tree.compute_first_entity_collision(mesh_pt)
        cell = Cell(mesh,cell_id)

        # calculate the coordinates of the vertices of the cell
        coordinate_dofs = cell.get_vertex_coordinates()

        # array of vertex indices associates with this cell
        i = list(space.dofmap().cell_dofs(cell.index()))

        # basis function value for each vertex at this point
        v = [fem.evaluate_basis(k,mesh_pt.coordinates(),coordinate_dofs,cell.orientation())[0] for k in range(fem.space_dimension())]

        for k in range(fem.space_dimension()):
            if v[k] >= threshold:
                I += [i[k]]
                J += [j]
                V += [v[k]]

    O = sparse.coo_matrix((V,(I,J)),shape=(mesh.num_vertices(),len(pts))).tocsr()
    return O

def make_obs_mat_time(obs_times,tL,tR,num_pts):
    dt = (tR-tL)/(num_pts-1)
    pts = np.linspace(tL,tR,num_pts)
    V,I,J = [],[],[]
    for j in range(len(obs_times)):
        obs_t = obs_times[j]
        left_i = int(np.floor((obs_t-tL)/dt))
        right_i = int(np.ceil((obs_t-tL)/dt))
        left_t = pts[left_i]
        right_t = pts[right_i]
        V+=[(right_t-obs_t)/dt,(obs_t-left_t)/dt]
        I+=[left_i,right_i]
        J+=[j,j]
    
    O = sparse.coo_matrix((V,(I,J)),shape=(num_pts,len(obs_times))).tocsr()
    return O
