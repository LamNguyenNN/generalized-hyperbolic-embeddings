import numpy as np
from .Circles_Lines import Circle, Segment

def entry_points_isom_cir(isom_cir):

    intersection = isom_cir.get_intersection(Circle(0,1))

    endpt1_arg = np.longdouble(get_arg(intersection[0]))
    endpt2_arg = np.longdouble(get_arg(intersection[1]))

    if np.abs(endpt1_arg - endpt2_arg) < np.pi:
        if endpt1_arg > endpt2_arg:
            init_point = intersection[0]
            term_point = intersection[1]
        else:
            init_point = intersection[1]
            term_point = intersection[0]

    else:
        if endpt1_arg < endpt2_arg:
            init_point = intersection[0]
            term_point = intersection[1]
        else:
            init_point = intersection[1]
            term_point = intersection[0]

    return init_point, term_point

def isProperVertex(z):
    if isclose(np.absolute(z), 1):
        return False
    else:
        return True

def NormalizedBoundary(elts):

    G = list(elts)

    G.sort(key=lambda T: get_arg(entry_points_isom_cir(T.get_isometric_circle())[1])) # Sort G based on the argument of terminal point of isometric circles

    print(G)

    L = Segment(0., 1.)
    intersect_pts = []
    intersect_angles = []
    min_angle = 2*np.pi
    min_x = 1
    min_idx = 0

    intersect_segment = False

    # Find isometric circles that intersect line segment (0,0) to (1,0)
    for (j, T) in enumerate(G):
        isom_cir = T.get_isometric_circle()
        intersection = L.get_intersection(isom_cir)

    if len(intersection) != 0: # Find isometric circle with smallest intersection with segment. If multiple circles with same intersection, find the one with smallest angle.
        intersect_segment = True
        angle = L.intersection_angle_in_unit_circle(T.get_isometric_circle())

        if intersection[0] < min_x or isclose(intersection[0], min_x): # non-strict inequality to allow for angle checks as well.

            if intersection[0] < min_x or angle < min_angle: # strict inequality means we can update to this circle for sure
              min_idx = j
              min_angle = angle

            min_x = intersection[0] # Always update the intersection

    if intersect_segment: # Cyclically shift G so that element intersecting segment with smallest intersection and angle is first
        shift_indices = [i % len(G) for i in range(min_idx, min_idx + len(G))]
        G = [G[i] for i in shift_indices]
        V = [L.get_intersection(G[0].get_isometric_circle())[0]]
    else: # Otherwise we start with G as is, i.e. arguments of terminal points in [0,2pi) sorted least to greatest
        V = [entry_points_isom_cir(G[0].get_isometric_circle())[1]]

    U = [G[0]] # U will track elements in boundary; V tracks vertices of boundary
    idx = 1

    while idx != len(G) + 1: # This condition is step 5

        # Looping through step 6...

        #print('U =', U)
        #print(G[idx % len(G)])

        isom_cir_T = U[-1].get_isometric_circle() # Isometric circle of last element of U
        isom_cir_Ti = G[idx % len(G)].get_isometric_circle() # Isometric circle of index idx element of G
        intersection = isom_cir_T.get_intersection(isom_cir_Ti) # Examine intersection point of these two circles

        if len(intersection) == 0: # If no intersection, examine terminal point isom_cir_Ti relative to isom_cir_T
        _, term_pt = entry_points_isom_cir(isom_cir_Ti)
            if isom_cir_T.contains_interior_point(term_pt): # If terminal point in the interior, increment idx and continue
                idx += 1
                continue
            else: # Otherwise, isom_cir_Ti completely outside, so add it to the list
                U.append(G[idx % len(G)])
                V.append(entry_points_isom_cir(isom_cir_T)[0])
                V.append(entry_points_isom_cir(isom_cir_Ti)[1])
                idx += 1

        else: # If there is an intersection point, compare distance between intersection and initial point of isom_cir_T and distance between last vertex of V and same initial point.
            v = intersection[0] if np.absolute(intersection[0]) < 1 else intersection[1]
            init_pt, _ = entry_points_isom_cir(isom_cir_T)

            if idx == 1 or len(U) == 1 or np.absolute(v - init_pt) < np.absolute(V[-1] - init_pt):

                if np.isclose(np.absolute(v - init_pt) - np.absolute(V[-1] - init_pt), 0): # Check if strict inequality actually holds
                    U = U[:-1]
                    V = V[:-1]
                    continue

                U.append(G[idx % len(G)])
                V.append(v)
                idx += 1
                continue

            else:
                U = U[:-1]
                V = V[:-1]
                continue

    elements = U[:-1]
    vertices = np.array(V[1:], dtype=np.clongdouble)

    vertex_angles = [get_arg(v) for v in vertices]

    sides = []
    for i in range(len(vertices)):
        if not isProperVertex(vertices[i-1]) and not isProperVertex(vertices[i]):
            continue
        else:
            sides.append([vertices[i], vertices[i-1]])

    sides = np.array(sides)

    icircs = [T.get_isometric_circle() for T in elements]

    return elements, vertices, vertex_angles, sides, icircs



def NormalizedBoundary_Outside(U, z):
    """Check if z lies on the outside or inside of a normalized boundary U. Return -1 if z in the interior or boundary of exterior domain. Otherwise, return index of isometric circle containing it."""

    for idx, icirc in enumerate(U[4]):
        if icirc.contains_interior_point(z):
            return idx

    # If we complete loop without returning, then z must lie in the exterior domain.
    return -1


def NormalizedBoundary_ShiftedPoint(icirc, scale, endpoint):
    # endpoint = 0 for initial
    # endpoint = 1 for terminal

    init_point, term_point = entry_points_isom_cir(icirc)
    init_arg = get_arg(init_point - icirc.center)
    term_arg = get_arg(term_point - icirc.center)

    if endpoint == 0:
        new_angle = term_arg - (term_arg - init_arg)*scale
    else:
        new_angle = init_arg + (term_arg - init_arg)*scale

    shifted_point = icirc.center + icirc.radius * np.exp(1j*new_angle)

    return shifted_point

def NormalizedBoundary_Plot(elements, show=True):

    plt.figure(figsize=(5, 5))
    plt.axis('equal')
    plt.xlim(xmin=-2, xmax=2)
    plt.ylim(ymin=-2, ymax=2)

    ts = np.linspace(0, 2*pi, 100)
    xs = np.cos(ts)
    ys = np.sin(ts)
    plt.plot(xs, ys)

    for T in elements:
        T.plot_isometric_circle()

    if show:
        plt.show()


def Reduction_Full(G, S, z, plot=False):
    S_red = S
    decomp = []
    z_ = S(z)

    while True:

        dists = np.array([diskDistFromOrigin(T(z_)) for T in G], np.clongdouble)

        if plot:

            plt.figure(figsize=(3, 3))
            plt.axis('equal')
            plt.xlim(xmin=-2, xmax=2)
            plt.ylim(ymin=-2, ymax=2)

            ts = np.linspace(0, 2*pi, 100)
            xs = np.cos(ts)
            ys = np.sin(ts)
            plt.plot(xs, ys)

            xs = np.real([T(z_) for T in G])
            ys = np.imag([T(z_) for T in G])

            xs = np.insert(xs, 0, np.real(z_))
            ys = np.insert(ys, 0, np.imag(z_))

            plt.scatter(xs[0], ys[0], c='red')
            plt.scatter(xs[1:], ys[1:], c='blue')

            for i in range(len(xs)):
                if i ==0:
                    plt.text(xs[i], ys[i], str(i), c='red')
                else:
                    plt.text(xs[i], ys[i], str(i), c='blue')

            for T in G:
                T.plot_isometric_circle()

            plt.show()


        if np.all((dists > diskDistFromOrigin(z_)) | (isclose(dists, diskDistFromOrigin(z_)))):
            break

        min_dist_idx = np.argmin(dists)

        z_ = G[min_dist_idx](z_)
        S_red = Mobius.compose(G[min_dist_idx], S_red)
        decomp.append(G[min_dist_idx])

    delta = Mobius.compose(S_red, S.invert())

    return S_red, delta, decomp


def Reduction_FromNormBound(U, S, z):

    # Note: If U is normalized basis, then reducing identity map returns element in group that sends z to fundamental domain.

    z_ = S(z)

    S_red = S
    decomp = []

    T = Mobius(1,0,0,1)

    while(True):
        outside = NormalizedBoundary_Outside(U, z_)

        if outside == -1:
            z_ = S_red(z)
            outside = NormalizedBoundary_Outside(U, z_)
            if outside == -1:
                break

        z_ = U[0][outside](z_)
        S_red = Mobius.compose(U[0][outside], S_red)
        decomp.append(U[0][outside])

    delta = Mobius.compose(S_red, S.invert())

    return S_red, delta, decomp


def EdgePairing(U, return_pair_unpaired = False):
    unpair = []
    pair = []

    for i, side in enumerate(U[3]):

        v0 = side[0]
        v1 = side[1]

        v0_image = U[0][i](v0)
        v0_im_arg = get_arg(v0_image)

        i1 = int(np.min([j for (j, arg) in enumerate(U[2]) if isclose(v0_im_arg, arg)], initial=len(U[2])))

        if i1 != len(U[2]):
            if not isclose(v0_image, U[1][i1]):
                i1 = len(U[2])

        v1_image = U[0][i](v1)
        v1_im_arg = get_arg(v1_image)

        if i1 != len(U[2]):
            i2 = i1 + 1
            if i2 == len(U[1]):
                i2 = 0

            if not isclose(v1_image, U[1][i2]):
                i2 = len(U[2])


            if i2 != len(U[2]):
                paired_side_idx = isclose(U[3], np.array([v1_image, v0_image], dtype=np.clongdouble))
                paired_side_idx = np.logical_and(paired_side_idx[:, 0], paired_side_idx[:, 1])
                paired_side_idx = np.where(paired_side_idx)[0][0]

                pair.append([side, np.array([v0_image, v1_image], dtype=np.clongdouble), U[0][i], U[0][paired_side_idx]])

            else:
                unpair.append([side, (True,False), U[4][i], U[0][i]]) # (True,False) means first vertex is paired, but second is unpaired

        else:
            i2 = int(np.min([j for (j, arg) in enumerate(U[2]) if np.isclose(arg, v1_im_arg)], initial=len(U[2])))
            if i2 != len(U[2]):
                if not isclose(v1_image, U[1][i2]):
                    i2 = len(U[2])

            if i2 != len(U[2]):
                unpair.append([side, (False,True), U[4][i], U[0][i]]) # First vertex is unpaired, but second is paired
            else:
                unpair.append([side, (False,False), U[4][i], U[0][i]]) # Both are unpaired

    if return_pair_unpaired:
        return unpair, pair
    else:
        return unpair


def NormalizedBasis(G):
    U = NormalizedBoundary(G)
    scale = 0.9999

    print('Computed initial boundary. Starting loop....')

    while(True):
        print('----------------------------------------------------------------------')
        print('Adding reductions to boundary')
        G_ = list(U[0])

        for S in G:
            S_red, _, _ = Reduction_Full(U[0], S, 0.)

            S_red_inv = S_red.invert()

            if S_red_inv == Mobius(1,0,0,1):
                print(S_red_inv.mat())
                print('Reduced element is identity. Not appending.')
            if S_red_inv != Mobius(1,0,0,1) and S_red_inv not in G_:
                G_.append(S_red_inv)

        print('Completed adding reductions to boundary.')

        U_ = NormalizedBoundary(G_)
        G = list(G_)

        NormalizedBoundary_Plot(U_[0])
        NormalizedBoundary_Plot(U[0])

        print('Comparing U and U_ (plots above)...')

        if not Mobius.compare_list(U_[0], U[0]):
            U = list(U_)
            print('No Match. Restarting...')
            print('----------------------------------------------------------------------')
            continue
        else:
            U = list(U_)
            print('Match found! Now checking for side pairing')

        print(U[0])
        NormalizedBoundary_Plot(U[0])

        unpair, pair = EdgePairing(U, return_pair_unpaired=True)

        if len(unpair) == 0:
            print('All sides are paired!')
            break

        print('----------------------------------------------------------------------')

        print(len(unpair), 'unpaired sides.')
        print(len(pair), 'paired sides.')

        print(unpair)
        print(pair)

        print('No side pairing. Computing reductions...')

        for idx, side_pairing_icirc_elt in enumerate(unpair):

            side = side_pairing_icirc_elt[0]
            pairing = side_pairing_icirc_elt[1]
            icirc = side_pairing_icirc_elt[2]
            S = side_pairing_icirc_elt[3]

            for i in range(2):
                if pairing[i]:
                    print('Vertex already paired. Continuing...')
                    continue

                v = side[i]

                if not isProperVertex(v):
                    v = NormalizedBoundary_ShiftedPoint(icirc, scale, i)

                print('Reducing', S, 'w.r.t', v)

                S_red, _, _ = Reduction_Full(U[0], S, v)

                print('Reduction completed')
                print(S_red)

                if S_red == Mobius(1,0,0,1):
                    print('Reduced element is identity. Not appending.')
                    print(S_red.mat())
                elif S_red not in G:
                    print('Unique reduction. Adding', S_red, 'to G...')
                    G.append(S_red)
                else:
                    print('Redundant reduction. Continuing...')


        print('Completed adding reductions for unpaired sides.')

        #return U
        U_new = NormalizedBoundary(G)

        if len(U_new[0]) == len(U[0]):

            if Mobius.compare_list(U_new[0], U[0]):
            scale = (9 + scale)/10

        U = list(U_new)

        print('Computed new boundary. Restarting loop...')
        #U = NormalizedBoundary(G)
        NormalizedBoundary_Plot(U[0])

    return U, pair