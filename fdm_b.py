
import numpy as np
import pdb
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve # to use its short name
from collections import namedtuple

class InputError(Exception):
    pass

def quivdata(Out, x, y, z=None, iz=0):
    """Returns coordinates and velocity components to show velocity with quiver
    
    Compute arrays to display velocity vectors using matplotlib's quiver.
    The quiver is always drawn in the xy-plane for a specific layer. Default iz=0
    
    Parameters
    ----------
    `Out` : namedtuple holding arrays `Qx`, `Qy`, `Qz` as defined in `fdm3`
        `Qx` : ndarray, shape: (Ny, Nx-1, Nz), [L3/T]
            Interfacial flows in finite difference model in x-direction from `fdm3'
        `Qy` : ndarray, shape: (Ny-1, Nx, Nz), [L3/T]
            Interfacial flows in finite difference model in y-direction from `fdm3`
        `Qz` : ndarray, shape: (Ny, Nx, Nz-1), [L3/T]
            Interfacial flows in finite difference model in z-direction from `fdm3`            
    `x` : ndarray, [m]
        Grid line coordinates of columns
    'y' : ndarray, [m]
        Grid line coordinates of rows
    `z` : ndaray [L] | int [-]
        If z == None, then iz must be given (default = 0)
        If z is an ndarray vector of floats
            z will be interpreted as the elvations of uniform layers.
            iz will be ignored
        If z is a full 3D ndarray of floats
            z will be interpreted as the elevations of the tops and bottoms of all cells.
            iz will be ignored
    `iz` : int [-]
            iz is ignored if z ~= None
            iz is the number of the layer for which the data are requested,
            and all output arrays will be 2D for that layer.

    Returns
    -------
    `Xm` : ndarray, shape: (Nz, Ny, Nx), [L]
        x-coordinates of cell centers
    `Ym` : ndarray, shape: (Nz, Ny, Nx), [L]
        y-coodinates of cell centers
    `ZM` : ndarray, shape: (Nz, Ny, Nx), [L]
        `z`-coordinates at cell centers
    `U` : ndarray, shape: (Nz, Ny, Nx), [L3/d]
        Flow in `x`-direction at cell centers
    `V` : ndarray, shape: (Nz, Ny, Nx), [L3/T]
        Flow in `y`-direction at cell centers
    `W` : ndarray, shape: (Nz, Ny, Nx), [L3/T]
        Flow in `z`-direction at cell centers.
    
    """
    Ny = len(y)-1
    Nx = len(x)-1
    
    xm = 0.5 * (x[:-1] + x[1:])
    ym = 0.5 * (y[:-1] + y[1:])
    
    X, Y = np.meshgrid(xm, ym) # coordinates of cell centers
    
    # Flows at cell centers
    U = np.concatenate((Out.Qx[iz, :, 0].reshape((1, Ny, 1)), \
                        0.5 * (Out.Qx[iz, :, :-1].reshape((1, Ny, Nx-2)) +\
                               Out.Qx[iz, :, 1: ].reshape((1, Ny, Nx-2))), \
                        Out.Qx[iz, :, -1].reshape((1, Ny, 1))), axis=2).reshape((Ny,Nx))
    V = np.concatenate((Out.Qy[iz, 0, :].reshape((1, 1, Nx)), \
                        0.5 * (Out.Qy[iz, :-1, :].reshape((1, Ny-2, Nx)) +\
                               Out.Qy[iz, 1:,  :].reshape((1, Ny-2, Nx))), \
                        Out.Qy[iz, -1, :].reshape((1, 1, Nx))), axis=1).reshape((Ny,Nx))
    return X, Y, U, V


def unique(x, tol=0.0001):
    """return sorted unique values of x, keeping ascending or descending direction"""
    if x[0]>x[-1]:  # vector is reversed
        x = np.sort(x)[::-1]  # sort and reverse
        return x[np.hstack((np.diff(x) < -tol, True))]
    else:
        x = np.sort(x)
        return x[np.hstack((np.diff(x) > +tol, True))]

    
def fdm3(x, y, z, kx, ky, kz, FQ, HI, IBOUND):
    '''Steady state 3D Finite Difference Model returning computed heads and flows.
        
    Heads and flows are returned as 3D arrays as specified under output parmeters.
    
    Parameters
    ----------
    `x` : ndarray,[L]
        `x` coordinates of grid lines perpendicular to rows, len is Nx+1
    `y` : ndarray, [L]
        `y` coordinates of grid lines along perpendicular to columns, len is Ny+1
    `z` : ndarray, [L]
        `z` coordinates of layers tops and bottoms, len = Nz+1
    `kx`, `ky`, `kz` : ndarray, shape: (Ny, Nx, Nz), [L/T]
        hydraulic conductivities along the three axes, 3D arrays.
    `FQ` : ndarray, shape: (Ny, Nx, Nz), [L3/T]
        prescrived cell flows (injection positive, zero of no inflow/outflow)
    `IH` : ndarray, shape: (Ny, Nx, Nz), [L]
        initial heads. `IH` has the prescribed heads for the cells with prescribed head.
    `IBOUND` : ndarray, shape: (Ny, Nx, Nz) of int
        boundary array like in MODFLOW with values denoting
        * IBOUND>0  the head in the corresponding cells will be computed
        * IBOUND=0  cells are inactive, will be given value NaN
        * IBOUND<0  coresponding cells have prescribed head
    
    outputs
    -------    
    `Out` : namedtuple containing heads and flows:
        `Out.Phi` : ndarray, shape: (Ny, Nx, Nz), [L3/T] 
            computed heads. Inactive cells will have NaNs
        `Out.Q`   : ndarray, shape: (Ny, Nx, Nz), [L3/T]
            net inflow in all cells, inactive cells have 0
        `Out.Qx   : ndarray, shape: (Ny, Nx-1, Nz), [L3/T] 
            intercell flows in x-direction (parallel to the rows)
        `Out.Qy`  : ndarray, shape: (Ny-1, Nx, Nz), [L3/T] 
            intercell flows in y-direction (parallel to the columns)
        `Out.Qz`  : ndarray, shape: (Ny, Nx, Nz-1), [L3/T] 
            intercell flows in z-direction (vertially upward postitive)
        the 3D array with the final heads with `NaN` at inactive cells.
    
    TO 160905
    '''

    # define the named tuple to hold all the output of the model fdm3
    Out = namedtuple('Out',['Phi', 'Q', 'Qx', 'Qy', 'Qz'])
    Out.__doc__ = """fdm3 output, <namedtuple>, containing fields Phi, Qx, Qy and Qz\n \
                    Use Out.Phi, Out.Q, Out.Qx, Out.Qy and Out.Qz"""                            

    x = unique(x)
    y = unique(y)[::-1]  # unique and descending
    z = unique(z)[::-1]  # unique and descending
        
    # as well as the number of cells along the three axes
    SHP = Nz, Ny, Nx = len(z)-1, len(y)-1, len(x)-1
    
    Nod = np.prod(SHP)
 
    if Nod == 0:
        raise AssetError(
            "Grid shape is (Ny, Nx, Nz) = {0}. Number of cells in all 3 direction must all be > 0".format(SHP))
                                
    if kx.shape != SHP:
        raise AssertionError("shape of kx {0} differs from that of model {1}".format(kx.shape,SHP))
    if ky.shape != SHP:
        raise AssertionError("shape of ky {0} differs from that of model {1}".format(ky.shape,SHP))
    if kz.shape != SHP:
        raise AssertionError("shape of kz {0} differs from that of model {1}".format(kz.shape,SHP))
        
    # from this we have the width of columns, rows and layers
    dx = np.abs(np.diff(x)).reshape(1, 1, Nx)
    dy = np.abs(np.diff(y)).reshape(1, Ny, 1)
    dz = np.abs(np.diff(z)).reshape(Nz, 1, 1)
    
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
    Cx = 1 / (Rx[:, :, :-1] + Rx[:, :, 1:])
    Cy = 1 / (Ry[:, :-1, :] + Ry[:, 1:, :])
    Cz = 1 / (Rz[:-1, :, :] + Rz[1:, :, :])
    
    NOD = np.arange(Nod).reshape(SHP)
    
    IE = NOD[:, :, 1:]  # east neighbor cell numbers
    IW = NOD[:, :, :-1] # west neighbor cell numbers
    IN = NOD[:, :-1, :] # north neighbor cell numbers
    IS = NOD[:, 1:, :]  # south neighbor cell numbers
    IT = NOD[:-1, :, :] # top neighbor cell numbers
    IB = NOD[1:, :, :]  # bottom neighbor cell numbers
    
    R = lambda x : x.ravel()  # generate anonymous function R(x) as shorthand for x.ravel()

    # notice the call  csc_matrix( (data, (rowind, coind) ), (M,N))  tuple within tupple
    # also notice that Cij = negative but that Cii will be postive, namely -sum(Cij)
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
    
    Out.Phi = HI.flatten() # allocate space to store heads
    
    Out.Phi[active] = spsolve( (A+Adiag)[active][:,active] ,RHS[active] ) # solve heads at active locations
    
    # net cell inflow
    Out.Q  = (A+Adiag).dot(Out.Phi).reshape(SHP)

    # set inactive cells to NaN
    Out.Phi[inact] = np.NaN # put NaN at inactive locations
    
    # reshape Phi to shape of grid
    Out.Phi = Out.Phi.reshape(SHP)
    
    #Flows across cell faces
    Out.Qx =  -np.diff( Out.Phi, axis=2) * Cx
    Out.Qy =  +np.diff( Out.Phi, axis=1) * Cy
    Out.Qz =  +np.diff( Out.Phi, axis=0) * Cz
    
    return Out # all outputs in a named tuple for easy access