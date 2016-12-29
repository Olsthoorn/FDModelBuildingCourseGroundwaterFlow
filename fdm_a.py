# write the function in this cell to disk as file fdm.py

import numpy as np
#import pdb  # in case we need to debug this function

def fdm3(x, y, z, kx, ky, kz, FQ, HI, IBOUND):
    '''Returns computed heads of steady state 3D finite difference grid.
    
    Steady state 3D Finite Difference Model that computes the heads a 3D ndarray.
    
    Parameters
    ----------
    `x` : ndarray, shape: Nx+1, [L]
        `x` coordinates of grid lines perpendicular to rows, len is Nx+1
    `y` : ndarray, shape: Ny+1, [L]
        `y` coordinates of grid lines along perpendicular to columns, len is Ny+1
    `z` : ndarray, shape: Nz+1, [L]
        `z` coordinates of layers tops and bottoms, len = Nz+1
    `kx`, `ky`, `kz` : ndarray, shape: (Ny, Nx, Nz) [L/T]
        hydraulic conductivities along the three axes, 3D arrays.
    `FQ` : ndarray, shape: (Ny, Nx, Nz), [L3/T]
        prescrived cell flows (injection positive, zero of no inflow/outflow)
    `IH` : ndarray, shape: (Ny, Nx, Nz), [L]
        initial heads. `IH` has the prescribed heads for the cells with prescribed head.
    `IBOUND` : ndarray of int, shape: (Ny, Nx, Nz), dim: [-]
        boundary array like in MODFLOW with values denoting
        * IBOUND>0  the head in the corresponding cells will be computed
        * IBOUND=0  cells are inactive, will be given value NaN
        * IBOUND<0  coresponding cells have prescribed head
    
    Returns
    -------    
    `Phi` : ndarray, shape: (Ny, Nx, Nz), [L]
        the 3D array with the final heads with `NaN` at inactive cells.
    
    TO 160905
    '''

    import numpy as np
    import scipy.sparse as sp
    from scipy.sparse.linalg import spsolve # to use its short name
#    pdb.set_trace()
    x = np.sort(np.array(x))       # enforce ascending
    y = np.sort(np.array(y))[::-1] # enforce descending
    z = np.sort(np.array(z))[::-1] # enforce descending
    
    # as well as the number of cells along the three axes
    SHP = Ny, Nx, Nz = len(y)-1, len(x)-1, len(z)-1

    Nod = np.prod(SHP)
    
    if Nod == 0:
        raise AssertationError("Nx, Ny and Nz must be >= 1")

    # assert correct shape of input arrays
    if kx.shape != SHP:
        raise AssertionError("shape of kx {0} differs from that of model {1}".format(kx.shape,SHP))
    if ky.shape != SHP:
        raise AssertionError("shape of ky {0} differs from that of model {1}".format(ky.shape,SHP))
    if kz.shape != SHP:
        raise AssertionError("shape of kz {0} differs from that of model {1}".format(kz.shape,SHP))
    
    # from this we have the width of columns, rows and layers
    dx = np.diff(x).reshape(1,Nx,1)
    dy = np.abs(np.diff(y).reshape(Ny,1,1)) # enforce positive
    dz = np.abs(np.diff(z)).reshape(1,1,Nz) # enforce positive
    
    active = (IBOUND>0).reshape(Nod,)  # boolean vector denoting the active cells
    inact  = (IBOUND==0).reshape(Nod,) # boolean vector denoting inacive cells
    fxhd   = (IBOUND<0).reshape(Nod,)  # boolean vector denoting fixed-head cells

    # half cell flow resistances
    Rx = 0.5 * dx / (dy * dz) / kx
    Ry = 0.5 * dy / (dz * dx) / ky
    Rz = 0.5 * dz / (dx * dy) / kz
    
    # set flow resistance in inactive cells to infinite
    Rx = Rx.reshape(Nod,); Rx[inact] = np.Inf; Rx=Rx.reshape(SHP)
    Ry = Ry.reshape(Nod,); Ry[inact] = np.Inf; Ry=Ry.reshape(SHP)
    Rz = Rz.reshape(Nod,); Rz[inact] = np.Inf; Rz=Rz.reshape(SHP)
    
    # conductances between adjacent cells
    Cx = 1 / (Rx[:,:-1,:] + Rx[:,1:,:])
    Cy = 1 / (Ry[:-1,:,:] + Ry[1:,:,:])
    Cz = 1 / (Rz[:,:,:-1] + Rz[:,:,1:])
    
    NOD = np.arange(Nod).reshape(SHP)
    
    IE = NOD[:,1:,:]  # east neighbor cell numbers
    IW = NOD[:,:-1,:] # west neighbor cell numbers
    IN = NOD[:-1,:,:] # north neighbor cell numbers
    IS = NOD[1:,:,:]  # south neighbor cell numbers
    IT = NOD[:,:,:-1] # top neighbor cell numbers
    IB = NOD[:,:,1:]  # bottom neighbor cell numbers
    
    R = lambda x : x.ravel()  # generate anonymous function R(x) as shorthand for x.ravel()

    # notice the call  csc_matrix( (data, (rowind, coind) ), (M,N))  tuple within tupple
    A = sp.csc_matrix(( -np.concatenate(( R(Cx), R(Cx), R(Cy), R(Cy), R(Cz), R(Cz)) ),\
                        (np.concatenate(( R(IE), R(IW), R(IN), R(IS), R(IB), R(IT)) ),\
                         np.concatenate(( R(IW), R(IE), R(IS), R(IN), R(IT), R(IB)) ),\
                      )),(Nod,Nod))
    
    # to use the vector of diagonal values in a call of sp.diags() we need to have it aa a 
    # standard nondimensional numpy vector.
    # To get this:
    # - first turn the matrix obtained by A.sum(axis=1) into a np.array by np.array( .. )
    # - then take the whole column to loose the array orientation (to get a dimensionless numpy vector)
    adiag = np.array(-A.sum(axis=1))[:,0]
    
    Adiag = sp.diags(adiag)  # diagonal matrix with a[i,i]
    
    RHS = FQ.reshape(Nod,1) - A[:,fxhd].dot(HI.reshape(Nod,1)[fxhd]) # Right-hand side vector
    
    Phi = HI.flatten() # allocate space to store heads
    
    Phi[active] = spsolve( (A+Adiag)[active][:,active] ,RHS[active] ) # solve heads at active locations
    
    Phi[inact] = np.NaN # put NaN at inactive locations
    
    return Phi.reshape(SHP) # reshape vector to 3D size of original model