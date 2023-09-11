#!/usr/bin/env python3

# Example program that initializes a spectral field with a
# spherical harmonic. Then the inverse transform is performed.
# Then the result can be visualized with Gmsh and pyvista
#    gmsh example2.msh

# ------------------------------------------------------------------------
#
# Limitations with Trans.backend('ectrans')
# -----------------------------------------
# 1) Only implemented for global Gaussian grids
#
# ------------------------------------------------------------------------
#
# Limitations with Trans.backend('local')
# ---------------------------------------
# 1) Only `invtrans` is implemented, and only MPI-serial
#
# ------------------------------------------------------------------------
#
# Known issues with Trans.backend('local'), to be fixed
# -----------------------------------------------------
# 1) We need to initialize trans with
#         trans = atlas.Trans(grid, truncation)
#    instead of
#         trans = atlas.Trans(fs_gp, fs_sp)
# 2) Using atlas Field API, only level=0 is supported (rank-1).
# 3) The Healpix grid is identified as unstructured
#    because the grid wrongly gets cropped to a "west", and its
#    structure gets lost
#
# ------------------------------------------------------------------------

import atlas4py as atlas

# ------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------
def set_spherical_harmonic(field_sp, n, m):
   """ Set field_sp to a spherical harmonic """
   spectral_functionspace = atlas.functionspace.Spectral(field_sp.functionspace)
   if not spectral_functionspace:
      raise TypeError("'field_sp.functionspace' could not be cast to "
                    + "'atlas.functionspace.Spectral'. It is of type "
                    + "'atlas.functionspace."+field_sp.functionspace.type+"'.")
   sp = atlas.make_view(field_sp)

   sp[:] = 0.
   def spherical_harmonic(re,im,_n,_m):
      if _n == n and _m == m:
         sp[re] = 1.
         sp[im] = 1.

   import numpy.random as random
   random.seed(1)
   rand1=random.rand(len(sp))
   random.seed(2)
   rand2=random.rand(len(sp))
   def spectrum(re,im,_n,_m):
      if _m==0:
         sp[re] = pow(_n+1,-5/6)*(rand1[re]-0.5)
      else:
         sp[re] = pow(_n+1,-5/6)*(rand1[re]-0.5)
         sp[im] = pow(_n+1,-5/6)*(rand2[im]-0.5) 


   #spectral_functionspace.parallel_for( spectrum )
   spectral_functionspace.parallel_for( spherical_harmonic )

# ------------------------------------------------------------------------

def plot_pyvista(mesh,field,title):
   """ Plot using pyviste, only works with mpi-serial for now """
   import pyvista
   ds = atlas.pyvista.dataset_create(mesh, coordinates='xyz', field=field)
   pyvista.set_plot_theme('paraview')
   pl = pyvista.Plotter()
   pl.add_mesh(ds, show_edges=False)
   pl.add_title(title)
   pl.show(interactive=True)

# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Start here
# ------------------------------------------------------------------------

atlas.initialize() # Required

init_total_wave_number=16
init_zonal_wave_number=8

grid = atlas.Grid("F400")
truncation = 127
atlas.Trans.backend("ectrans")
partitioner = atlas.Partitioner("ectrans")
  # 'ectrans' partitioner needs to be used for 'ectrans' Trans backend only.
  # Partitions are nearly identical to the equal_regions partitioner

if not atlas.GaussianGrid(grid) or not atlas.Trans.has_backend("ectrans"):
   # Fallback to "local" backend
   atlas.Trans.backend("local")
   partitioner = atlas.Partitioner("equal_regions")
atlas.Trans.backend("local")
partitioner = atlas.Partitioner("equal_regions")

# Create function spaces
fs_sp = atlas.functionspace.Spectral(truncation)
fs_gp = atlas.functionspace.StructuredColumns(grid, partitioner, halo=1)

# Create fields (already distributed)
field_sp = fs_sp.create_field(name="sp")
field_gp = fs_gp.create_field(name="gp")

# Initial condition for spectral field
set_spherical_harmonic(field_sp, n=init_total_wave_number, m=init_zonal_wave_number)

# Spectral transform
trans = atlas.Trans(grid, truncation)
trans.invtrans(field_sp, field_gp)

# Visualisation of gridpoint field
mesh = atlas.MeshGenerator(type='structured').generate(grid)
atlas.Gmsh("example2.msh", coordinates='xyz').write(mesh).write(field_gp)

# plot_pyvista(mesh, field_gp, title=grid.name + " - T" + str(truncation) + "; n="+str(init_total_wave_number)+" m="+str(init_zonal_wave_number))
#    --> only works with mpi-serial for now

atlas.finalize()
# ------------------------------------------------------------------------
