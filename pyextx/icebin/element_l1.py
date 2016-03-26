import numpy as np
import scipy.sparse

"""Computes a regridding matrix between an L1 (piecewise linear)
finite elements and a L0 ("regular") grid with constant-value grid
cells.

These functions are meant to work with ibgrid"""


def eqn_plane_subelement(element, basis_vertex):
    """Provides the equation for the sub-element with vertices at three points:
        element: (p0, p1, p2) = ((p0x,p0y), (p1x,p1y), (p2x,p2y))
    Equation is z = Ax+By+C, returned as (A,B,C)

    If basis_vertex == 0, we solve for A,B,C in:
        Ap0x + Bp0y + C = p0z
        Ap1x + Bp1y + C = 0
        AP2x + Bp2y + C = 0

    basis_vertex: Index (0-2) of the vertex that is set to 1 for this
    sub-basis function; other vertices equal 0.
    """

    M = np.zeros((3,3))
    M[0,0] = element[0].x
    M[0,1] = element[0].y
    M[1,0] = element[1].x
    M[1,1] = element[1].y
    M[2,0] = element[2].x
    M[2,1] = element[2].y
    M[:,2] = 1.
    rhs = np.zeros((3,))
    rhs[basis_vertex] = 1.      # p0z
    return np.linalg.solve(M, rhs)



def integrate_subelement(element, basis_vertex, polygon):
    """Integrates a sub-element's plane around a polygon.
    A sub-element is defined as the plane z=Ax+By+C defined by one
        vertex in one element.

    element: [p0, p1, p2] list(Vertex)
        The points in the element's triangle.
    polygon: [q0, q1, ...]
        The polygon to integrate around.  It will often be the same as
        element, but will sometimes be smaller."""

    A,B,C = eqn_plane_subelement(element, basis_vertex)
    # print('ABC', element, A,B,C)
    total = 0.
    q0=polygon[-1]
    for q1 in polygon:
        # These formulae were obtained from element3.maxima.
        # Run:   maxima <element3.maxima
        A_coeff = (1./6.)  *  (0 \
            + 2 * q0.x * q0.x * q0.y \
            +   q0.x * q0.x * q1.y \
            -   q0.x * q0.y * q1.x \
            +   q0.x * q1.x * q1.y \
            -   q0.y * q1.x * q1.x \
            - 2 * q1.x * q1.x * q1.y )

        B_coeff = (1./6.)  *  (0 \
            - 2  * q0.x * q0.y * q0.y \
            +   q0.x * q0.y * q1.y \
            -   q0.y * q0.y * q1.x \
            +   q0.x * q1.y * q1.y \
            -   q0.y * q1.x * q1.y \
            + 2 * q1.x * q1.y * q1.y)

        C_coeff = .5  *  (q0.x * q1.y - q0.y * q1.x)

        total += A*A_coeff + B*B_coeff + C*C_coeff
        q0 = q1
    return total


def compute_AvI(gridX, gridI):
    """Computes A<-I regridding matrix.
    returns: AvI, weightsA, weightsI

    To regrid using this matrix:
        fA = (1/weightsA) AvI fI
        fI = (1/weightsI) IvA fA    where IvA = transpose(AvI)

    NOTE: weightsA is the area of A that OVERLAPS I.  Similar for
        weightsI.  If you wish to be weighting by the FULL area of A,
        compute that elsewhere.

    Returns: AvI, weightsA, weightsI"""


    # Iterate through polygons in exchange grid.  This is the efficient
    # way to find ice basis functions relevant to teach atmosphere grid
    # cell.  Each of these polygons might involve multiple ice basis
    # functions
    data = []
    ii = []
    jj = []
    for cellX in gridX.cells.values():
        # Get the polygon we wish to integrate around
        polygon = cellX.vertices
        ix_A = cellX.i      # Index in the atmosphere grid

        # Fish out the element definition from gridI
        element = gridI.cells[cellX.j].vertices

        # Loop over the partial basis functions for this element
        # Each partial basis function is 1 at one vertex and 0 at the other two
        for basis_vertex in range(0,len(element)):
            # Index of this basis function in the ice mesh
            ix_I = element[basis_vertex].index

            val = integrate_subelement(element, basis_vertex, polygon)

            #AvI.add((ix_A, ix_I), val) # Equiv. to AvI[poly_iA, poly_iI] += val
            data.append(val)
            ii.append(ix_A)
            jj.append(ix_I)

    AvI = scipy.sparse.coo_matrix((data,(ii,jj)),
        shape=(gridX.grid1_ncells_full, gridX.grid2_nvertices_full))

    weightsA = np.asarray(AvI.sum(1)).transpose()[0]
    weightsI = np.asarray(AvI.sum(0))[0]
    return (AvI, weightsA, weightsI)
