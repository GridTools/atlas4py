#!/usr/bin/env python3
import atlas4py as atlas
import math

atlas.initialize() # Required

mesh = atlas.Mesh( atlas.Grid("H16") )

fs_nodes = atlas.functionspace.NodeColumns(mesh)

field = fs_nodes.create_field(name="myfield")

lonlat = atlas.make_view(fs_nodes.lonlat)

view = atlas.make_view(field)
for n in range(field.shape[0]):
   lat = lonlat[n,1] * math.pi / 180.
   view[n] = math.cos(4.*lat)

gmsh = atlas.Gmsh("example1.msh", coordinates='xyz')
gmsh.write(mesh)
gmsh.write(field)

def plot_pyvista(mesh,field,title):
   """ Plot with pyvista only mpi-serial for now"""
   import pyvista
   ds = atlas.pyvista.dataset_create(mesh, coordinates='xyz', field=field)
   pyvista.set_plot_theme('paraview')
   pl = pyvista.Plotter()
   pl.add_mesh(ds, show_edges=True)
   pl.add_title(title)
   pl.show(interactive=True)

# plot_pyvista(mesh,field,'title')

atlas.finalize()

