from firedrake import *
from netgen import *
from netgen.occ import *
from netgen.gui import *
from netgen.meshing import meshsize, MeshingParameters

geo = OCCGeometry('scaled_reactor.step')
mp = MeshingParameters(minh=0.1)  # or meshsize.very_coarse, or meshsize.moderate, or meshsize.fine, or meshsize.very_fine ngmesh = geo.GenerateMesh(mp) mesh = Mesh(Mesh(ngmesh).curve_field(3))  # if you want a high-order mesh
VTKFile("output/mesh.pvd").write(mesh)
