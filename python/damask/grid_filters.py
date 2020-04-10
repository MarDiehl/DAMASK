from scipy import spatial as _spatial
import numpy as _np

def _ks(size,grid,first_order=False):
    """
    Get wave numbers operator.

    Parameters
    ----------
    size : numpy.ndarray
        physical size of the periodic field.

    """
    k_sk = _np.where(_np.arange(grid[0])>grid[0]//2,_np.arange(grid[0])-grid[0],_np.arange(grid[0]))/size[0]
    if grid[0]%2 == 0 and first_order: k_sk[grid[0]//2] = 0                                         # Nyquist freq=0 for even grid (Johnson, MIT, 2011)

    k_sj = _np.where(_np.arange(grid[1])>grid[1]//2,_np.arange(grid[1])-grid[1],_np.arange(grid[1]))/size[1]
    if grid[1]%2 == 0 and first_order: k_sj[grid[1]//2] = 0                                         # Nyquist freq=0 for even grid (Johnson, MIT, 2011)

    k_si = _np.arange(grid[2]//2+1)/size[2]

    kk, kj, ki = _np.meshgrid(k_sk,k_sj,k_si,indexing = 'ij')
    return _np.concatenate((ki[:,:,:,None],kj[:,:,:,None],kk[:,:,:,None]),axis = 3)


def curl(size,field):
    """
    Calculate curl of a vector or tensor field in Fourier space.

    Parameters
    ----------
    size : numpy.ndarray
        physical size of the periodic field.

    """
    n = _np.prod(field.shape[3:])
    k_s = _ks(size,field.shape[:3],True)

    e = _np.zeros((3, 3, 3))
    e[0, 1, 2] = e[1, 2, 0] = e[2, 0, 1] = +1.0                                                     # Levi-Civita symbol
    e[0, 2, 1] = e[2, 1, 0] = e[1, 0, 2] = -1.0

    field_fourier = _np.fft.rfftn(field,axes=(0,1,2))
    curl_ = (_np.einsum('slm,ijkl,ijkm ->ijks', e,k_s,field_fourier)*2.0j*_np.pi if n == 3 else       # vector, 3   -> 3
             _np.einsum('slm,ijkl,ijknm->ijksn',e,k_s,field_fourier)*2.0j*_np.pi)                     # tensor, 3x3 -> 3x3

    return _np.fft.irfftn(curl_,axes=(0,1,2),s=field.shape[:3])


def divergence(size,field):
    """
    Calculate divergence of a vector or tensor field in Fourier space.

    Parameters
    ----------
    size : numpy.ndarray
        physical size of the periodic field.

    """
    n = _np.prod(field.shape[3:])
    k_s = _ks(size,field.shape[:3],True)

    field_fourier = _np.fft.rfftn(field,axes=(0,1,2))
    div_ = (_np.einsum('ijkl,ijkl ->ijk', k_s,field_fourier)*2.0j*_np.pi if n == 3 else               # vector, 3   -> 1
            _np.einsum('ijkm,ijklm->ijkl',k_s,field_fourier)*2.0j*_np.pi)                             # tensor, 3x3 -> 3

    return _np.fft.irfftn(div_,axes=(0,1,2),s=field.shape[:3])


def gradient(size,field):
    """
    Calculate gradient of a vector or scalar field in Fourier space.

    Parameters
    ----------
    size : numpy.ndarray
        physical size of the periodic field.

    """
    n = _np.prod(field.shape[3:])
    k_s = _ks(size,field.shape[:3],True)

    field_fourier = _np.fft.rfftn(field,axes=(0,1,2))
    grad_ = (_np.einsum('ijkl,ijkm->ijkm', field_fourier,k_s)*2.0j*_np.pi if n == 1 else              # scalar, 1 -> 3
             _np.einsum('ijkl,ijkm->ijklm',field_fourier,k_s)*2.0j*_np.pi)                            # vector, 3 -> 3x3

    return _np.fft.irfftn(grad_,axes=(0,1,2),s=field.shape[:3])


def cell_coord0(grid,size,origin=_np.zeros(3)):
    """
    Cell center positions (undeformed).

    Parameters
    ----------
    grid : numpy.ndarray
        number of grid points.
    size : numpy.ndarray
        physical size of the periodic field.
    origin : numpy.ndarray, optional
        physical origin of the periodic field. Defaults to [0.0,0.0,0.0].

    """
    start = origin        + size/grid*.5
    end   = origin + size - size/grid*.5
    return _np.mgrid[start[0]:end[0]:grid[0]*1j,start[1]:end[1]:grid[1]*1j,start[2]:end[2]:grid[2]*1j].T


def cell_displacement_fluct(size,F):
    """
    Cell center displacement field from fluctuation part of the deformation gradient field.

    Parameters
    ----------
    size : numpy.ndarray
        physical size of the periodic field.
    F : numpy.ndarray
        deformation gradient field.

    """
    integrator = 0.5j*size/_np.pi

    k_s = _ks(size,F.shape[:3],False)
    k_s_squared = _np.einsum('...l,...l',k_s,k_s)
    k_s_squared[0,0,0] = 1.0

    displacement = -_np.einsum('ijkml,ijkl,l->ijkm',
                              _np.fft.rfftn(F,axes=(0,1,2)),
                              k_s,
                              integrator,
                             ) / k_s_squared[...,_np.newaxis]

    return _np.fft.irfftn(displacement,axes=(0,1,2),s=F.shape[:3])


def cell_displacement_avg(size,F):
    """
    Cell center displacement field from average part of the deformation gradient field.

    Parameters
    ----------
    size : numpy.ndarray
        physical size of the periodic field.
    F : numpy.ndarray
        deformation gradient field.

    """
    F_avg = _np.average(F,axis=(0,1,2))
    return _np.einsum('ml,ijkl->ijkm',F_avg-_np.eye(3),cell_coord0(F.shape[:3][::-1],size))


def cell_displacement(size,F):
    """
    Cell center displacement field from deformation gradient field.

    Parameters
    ----------
    size : numpy.ndarray
        physical size of the periodic field.
    F : numpy.ndarray
        deformation gradient field.

    """
    return cell_displacement_avg(size,F) + cell_displacement_fluct(size,F)


def cell_coord(size,F,origin=_np.zeros(3)):
    """
    Cell center positions.

    Parameters
    ----------
    size : numpy.ndarray
        physical size of the periodic field.
    F : numpy.ndarray
        deformation gradient field.
    origin : numpy.ndarray, optional
        physical origin of the periodic field. Defaults to [0.0,0.0,0.0].

    """
    return cell_coord0(F.shape[:3][::-1],size,origin) + cell_displacement(size,F)


def cell_coord0_gridSizeOrigin(coord0,ordered=True):
    """
    Return grid 'DNA', i.e. grid, size, and origin from array of cell positions.

    Parameters
    ----------
    coord0 : numpy.ndarray
        array of undeformed cell coordinates.
    ordered : bool, optional
        expect coord0 data to be ordered (x fast, z slow).

    """
    coords    = [_np.unique(coord0[:,i]) for i in range(3)]
    mincorner = _np.array(list(map(min,coords)))
    maxcorner = _np.array(list(map(max,coords)))
    grid      = _np.array(list(map(len,coords)),'i')
    size      = grid/_np.maximum(grid-1,1) * (maxcorner-mincorner)
    delta     = size/grid
    origin    = mincorner - delta*.5

    # 1D/2D: size/origin combination undefined, set origin to 0.0
    size  [_np.where(grid==1)] = origin[_np.where(grid==1)]*2.
    origin[_np.where(grid==1)] = 0.0

    if grid.prod() != len(coord0):
        raise ValueError('Data count {} does not match grid {}.'.format(len(coord0),grid))

    start = origin + delta*.5
    end   = origin - delta*.5 + size

    if not _np.allclose(coords[0],_np.linspace(start[0],end[0],grid[0])) and \
           _np.allclose(coords[1],_np.linspace(start[1],end[1],grid[1])) and \
           _np.allclose(coords[2],_np.linspace(start[2],end[2],grid[2])):
        raise ValueError('Regular grid spacing violated.')

    if ordered and not _np.allclose(coord0.reshape(tuple(grid[::-1])+(3,)),cell_coord0(grid,size,origin)):
        raise ValueError('I_nput data is not a regular grid.')

    return (grid,size,origin)


def coord0_check(coord0):
    """
    Check whether coordinates lie on a regular grid.

    Parameters
    ----------
    coord0 : numpy.ndarray
        array of undeformed cell coordinates.

    """
    cell_coord0_gridSizeOrigin(coord0,ordered=True)


def node_coord0(grid,size,origin=_np.zeros(3)):
    """
    Nodal positions (undeformed).

    Parameters
    ----------
    grid : numpy.ndarray
        number of grid points.
    size : numpy.ndarray
        physical size of the periodic field.
    origin : numpy.ndarray, optional
        physical origin of the periodic field. Defaults to [0.0,0.0,0.0].

    """
    return _np.mgrid[origin[0]:size[0]+origin[0]:(grid[0]+1)*1j,
                    origin[1]:size[1]+origin[1]:(grid[1]+1)*1j,
                    origin[2]:size[2]+origin[2]:(grid[2]+1)*1j].T


def node_displacement_fluct(size,F):
    """
    Nodal displacement field from fluctuation part of the deformation gradient field.

    Parameters
    ----------
    size : numpy.ndarray
        physical size of the periodic field.
    F : numpy.ndarray
        deformation gradient field.

    """
    return cell_2_node(cell_displacement_fluct(size,F))


def node_displacement_avg(size,F):
    """
    Nodal displacement field from average part of the deformation gradient field.

    Parameters
    ----------
    size : numpy.ndarray
        physical size of the periodic field.
    F : numpy.ndarray
        deformation gradient field.

    """
    F_avg = _np.average(F,axis=(0,1,2))
    return _np.einsum('ml,ijkl->ijkm',F_avg-_np.eye(3),node_coord0(F.shape[:3][::-1],size))


def node_displacement(size,F):
    """
    Nodal displacement field from deformation gradient field.

    Parameters
    ----------
    size : numpy.ndarray
        physical size of the periodic field.
    F : numpy.ndarray
        deformation gradient field.

    """
    return node_displacement_avg(size,F) + node_displacement_fluct(size,F)


def node_coord(size,F,origin=_np.zeros(3)):
    """
    Nodal positions.

    Parameters
    ----------
    size : numpy.ndarray
        physical size of the periodic field.
    F : numpy.ndarray
        deformation gradient field.
    origin : numpy.ndarray, optional
        physical origin of the periodic field. Defaults to [0.0,0.0,0.0].

    """
    return node_coord0(F.shape[:3][::-1],size,origin) + node_displacement(size,F)


def cell_2_node(cell_data):
    """Interpolate periodic cell data to nodal data."""
    n = (  cell_data + _np.roll(cell_data,1,(0,1,2))
         + _np.roll(cell_data,1,(0,))  + _np.roll(cell_data,1,(1,))  + _np.roll(cell_data,1,(2,))
         + _np.roll(cell_data,1,(0,1)) + _np.roll(cell_data,1,(1,2)) + _np.roll(cell_data,1,(2,0)))*0.125

    return _np.pad(n,((0,1),(0,1),(0,1))+((0,0),)*len(cell_data.shape[3:]),mode='wrap')


def node_2_cell(node_data):
    """Interpolate periodic nodal data to cell data."""
    c = (  node_data + _np.roll(node_data,1,(0,1,2))
         + _np.roll(node_data,1,(0,))  + _np.roll(node_data,1,(1,))  + _np.roll(node_data,1,(2,))
         + _np.roll(node_data,1,(0,1)) + _np.roll(node_data,1,(1,2)) + _np.roll(node_data,1,(2,0)))*0.125

    return c[:-1,:-1,:-1]


def node_coord0_gridSizeOrigin(coord0,ordered=False):
    """
    Return grid 'DNA', i.e. grid, size, and origin from array of nodal positions.

    Parameters
    ----------
    coord0 : numpy.ndarray
        array of undeformed nodal coordinates.
    ordered : bool, optional
        expect coord0 data to be ordered (x fast, z slow).

    """
    coords    = [_np.unique(coord0[:,i]) for i in range(3)]
    mincorner = _np.array(list(map(min,coords)))
    maxcorner = _np.array(list(map(max,coords)))
    grid      = _np.array(list(map(len,coords)),'i') - 1
    size      = maxcorner-mincorner
    origin    = mincorner

    if (grid+1).prod() != len(coord0):
        raise ValueError('Data count {} does not match grid {}.'.format(len(coord0),grid))

    if not _np.allclose(coords[0],_np.linspace(mincorner[0],maxcorner[0],grid[0]+1)) and \
           _np.allclose(coords[1],_np.linspace(mincorner[1],maxcorner[1],grid[1]+1)) and \
           _np.allclose(coords[2],_np.linspace(mincorner[2],maxcorner[2],grid[2]+1)):
        raise ValueError('Regular grid spacing violated.')

    if ordered and not _np.allclose(coord0.reshape(tuple((grid+1)[::-1])+(3,)),node_coord0(grid,size,origin)):
        raise ValueError('I_nput data is not a regular grid.')

    return (grid,size,origin)


def regrid(size,F,new_grid):
    """
    Return mapping from coordinates in deformed configuration to a regular grid.

    Parameters
    ----------
    size : numpy.ndarray
        physical size
    F : numpy.ndarray
        deformation gradient field
    new_grid : numpy.ndarray
        new grid for undeformed coordinates

    """
    c = cell_coord0(F.shape[:3][::-1],size) \
      + cell_displacement_avg(size,F) \
      + cell_displacement_fluct(size,F)

    outer = _np.dot(_np.average(F,axis=(0,1,2)),size)
    for d in range(3):
        c[_np.where(c[:,:,:,d]<0)]        += outer[d]
        c[_np.where(c[:,:,:,d]>outer[d])] -= outer[d]

    tree = _spatial.cKDTree(c.reshape(-1,3),boxsize=outer)
    return tree.query(cell_coord0(new_grid,outer))[1].flatten()
