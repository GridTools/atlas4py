
try:
   import pyvista as pv
   _atlas_pyvista_support = True
except ImportError:
   _atlas_pyvista_support = False

import numpy as np
import atlas4py as atlas

def _atlas_pyvista_connectivity(conn):
   num_neighbors = np.add.reduce(np.where(conn != -1, 1, 0), axis=1)
   tmp = np.concatenate((num_neighbors[:, None], conn), axis=1)
   return tmp[np.concatenate((np.full((conn.shape[0], 1), True, dtype=bool), conn != -1), axis=1)]

def _atlas_connectivity_to_numpy(atlas_conn, *, out=None ):
   if isinstance(atlas_conn, atlas.BlockConnectivity):
      shape = (atlas_conn.rows, atlas_conn.cols)
      out = np.zeros(shape, dtype=np.int64) if out is None else out
      assert out.shape == shape

      for i in range(atlas_conn.rows):
         for nb in range(atlas_conn.cols):
            out[i, nb] = atlas_conn[i, nb]

      return out

   shape = (atlas_conn.rows, atlas_conn.maxcols)
   out = np.zeros(shape, dtype=np.int64) if out is None else out
   assert out.shape == shape

   for i in range(atlas_conn.rows):
      cols = atlas_conn.cols(i)
      for nb in range(cols):
         out[i, nb] = atlas_conn[i, nb]
      out[i, cols:]=-1

   return out

def _atlas_pyvista_cells(mesh):
   np_conn = _atlas_connectivity_to_numpy(mesh.cells.node_connectivity)
   return _atlas_pyvista_connectivity(np_conn)

def _atlas_pyvista_points(mesh, coordinates='lonlat', **kwargs):
   if coordinates=='lonlat':
      return np.insert(mesh.nodes.lonlat, 2, 0., axis=1)
   elif coordinates=='xy':
      return np.insert(mesh.nodes.field('xy'), 2, 0., axis=1)
   elif coordinates=='xyz':
      deg2rad = np.pi/180.
      xyrad = atlas.make_view(mesh.nodes.lonlat) * deg2rad
      phi, theta = xyrad[:, 1], xyrad[:, 0]
      return np.stack((np.cos(phi) * np.cos(theta), np.cos(phi) * np.sin(theta), np.sin(phi)), axis=1)
   else:
      raise IndexError("No coordinates "+coordinates+" supported")

def _atlas_pyvista_mesh(mesh, coordinates):
   import pyvista as pv
   points = _atlas_pyvista_points( mesh, coordinates) 
   cells = _atlas_pyvista_cells(mesh)
   return pv.PolyData(points, cells)

def dataset_update(dataset,field,level=0,variable=0):
    if atlas.functionspace.NodeColumns(field.functionspace):
        field.halo_exchange()
        if field.rank == 1:
            dataset.point_data[field.name] = atlas.make_view(field)[:]
        if field.rank == 2:
            dataset.point_data[field.name] = atlas.make_view(field)[:,level]
        if field.rank == 3:
            dataset.point_data[field.name] = atlas.make_view(field)[:,level,variable]
    else:
        if not hasattr(dataset, 'atlas_mesh'):
            raise Exception("no atlas_mesh stored. This is not a atlas-compatible pyvista dataset")
        if not hasattr(dataset, 'nodes_fs'):
            dataset.nodes_fs = atlas.functionspace.NodeColumns(dataset.atlas_mesh)
        tmp_field = dataset.nodes_fs.create_field(dtype=field.datatype)
        n_points = min(field.shape[0],dataset.number_of_points)
        if field.rank == 1:
            atlas.make_view(tmp_field)[:n_points] = atlas.make_view(field)[:n_points]
        if field.rank == 2:
            atlas.make_view(tmp_field)[:n_points] = atlas.make_view(field)[:n_points,level]
        if field.rank == 3:
             atlas.make_view(tmp_field)[:n_points] = atlas.make_view(field)[:n_points,level,variable]
        tmp_field.halo_exchange()
        dataset.point_data[field.name] = atlas.make_view(tmp_field)

def dataset_create(mesh,coordinates='lonlat',field=None,level=0,variable=0):
   if not _atlas_pyvista_support:
      raise ModuleNotFoundError("No module named 'pyvista'")

   dataset = _atlas_pyvista_mesh(mesh, coordinates)
   dataset.atlas_mesh = mesh
   if field:
      dataset_update(dataset,field=field,level=level,variable=variable)
   return dataset

