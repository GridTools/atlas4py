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

gmsh = atlas.Gmsh("out.msh", coordinates='xyz')
gmsh.write(mesh)
gmsh.write(field)

atlas.finalize()

