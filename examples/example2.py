#!/usr/bin/env python3

# Example program that initializes a spectral field with a
# spherical harmonic. Then the inverse transform is performed
# Then the result can be visualized with Gmsh:
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
# 2) Using atlas Field API, only level=1 is supported.
# 3) The Healpix grid is identified as unstructured
#    because the grid wrongly gets cropped to a "west", and its
#    structure gets lost
#
# ------------------------------------------------------------------------

import atlas4py as atlas

# ------------------------------------------------------------------------
# Functions
def set_spherical_harmonic(field_sp, n, m):
   """ Set field_sp to a spherical harmonic """
   nflds = field_sp.shape[1]
   spectral_functionspace = atlas.functionspace.Spectral(field_sp.functionspace)
   if not spectral_functionspace:
      raise TypeError("'field_sp.functionspace' could not be cast to "
                    + "'atlas.functionspace.Spectral'. It is of type "
                    + "'atlas.functionspace."+field_sp.functionspace.type+"'.")

   sp = atlas.make_view(field_sp)
   sp[:,:] = 0.
   def spherical_harmonic(re,im,_n,_m):
      for jfld in range(nflds):
         if n*(jfld+1) == _n and m*(jfld+1) == _m:
            sp[re,jfld] = 1.
            sp[im,jfld] = 1.

   spectral_functionspace.parallel_for( spherical_harmonic )
# ------------------------------------------------------------------------


# ------------------------------------------------------------------------
# Start here
# ------------------------------------------------------------------------

atlas.initialize() # Required

grid = atlas.Grid("O128")
levels = 10
truncation = 255
atlas.Trans.backend("ectrans")
partitioner = atlas.Partitioner("ectrans")
  # 'ectrans' partitioner needs to be used for 'ectrans' Trans backend only.
  # It is essentially nearly identical to the equal_regions partitioner

if not atlas.GaussianGrid(grid) or not atlas.Trans.has_backend("ectrans"):
   # Fallback to "local" backend
   atlas.Trans.backend("local")
   levels = 1
   partitioner = atlas.Partitioner("equal_regions")

# Create function spaces
fs_sp = atlas.functionspace.Spectral(truncation)
fs_gp = atlas.functionspace.StructuredColumns(grid, partitioner)

# Create fields (already distributed)
field_sp = fs_sp.create_field(name="sp", levels=levels)
field_gp = fs_gp.create_field(name="gp", levels=levels)

# Initial condition for spectral field
set_spherical_harmonic(field_sp, n=4, m=2)

# Spectral transform
trans = atlas.Trans(grid, truncation)
trans.invtrans(field_sp, field_gp)

# Visualisation of gridpoint field
meshgenerator = atlas.MeshGenerator(type='structured')
mesh = meshgenerator.generate(grid)
gmsh = atlas.Gmsh("example2.msh", coordinates='xy')
gmsh.write(mesh)
gmsh.write(field_gp)

atlas.finalize()
# ------------------------------------------------------------------------
