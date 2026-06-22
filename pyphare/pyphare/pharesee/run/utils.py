#
#
#

import numpy as np


def _current1d(by, bz, xby, xbz):
    # jx = 0
    # jy = -dxBz
    # jz = dxBy
    # the following hard-codes yee layout
    # which is not true in general
    # we should at some point provide proper
    # derivation routines in the gridlayout
    dx = xbz[1] - xbz[0]
    jy = np.zeros(by.size + 1)
    jy[1:-1] = -(bz[1:] - bz[:-1]) / dx
    dx = xby[1] - xby[0]
    jz = np.zeros(bz.size + 1)
    jz[1:-1] = (by[1:] - by[:-1]) / dx
    jy[0] = jy[1]
    jy[-1] = jy[-2]
    jz[0] = jz[1]
    jz[-1] = jz[-2]
    return jy, jz


def _current2d(bx, by, bz, dx, dy):
    # jx = dyBz
    # jy = -dxBz
    # jz = dxBy - dyBx
    # the following hard-codes yee layout
    # which is not true in general
    # we should at some point provide proper
    # derivation routines in the gridlayout
    jx = np.zeros(by.shape)
    jy = np.zeros(bx.shape)
    jz = np.zeros((bx.shape[0], by.shape[1]))

    jx[:, 1:-1] = (bz[:, 1:] - bz[:, :-1]) / dy
    jy[1:-1, :] = -(bz[1:, :] - bz[:-1, :]) / dx
    jz[1:-1, 1:-1] = (by[1:, 1:-1] - by[:-1, 1:-1]) / dx - (
        bx[1:-1, 1:] - bx[1:-1, :-1]
    ) / dy

    jy[0, :] = jy[1, :]
    jy[:, 0] = jy[:, 1]
    jy[-1, :] = jy[-2, :]
    jy[:, -1] = jy[:, -2]

    jz[0, :] = jz[1, :]
    jz[:, 0] = jz[:, 1]
    jz[-1, :] = jz[-2, :]
    jz[:, -1] = jz[:, -2]

    jx[0, :] = jx[1, :]
    jx[:, 0] = jx[:, 1]
    jx[-1, :] = jx[-2, :]
    jx[:, -1] = jx[:, -2]

    return jx, jy, jz


def _compute_current(patch, **kwargs):
    reference_pd = patch["Bx"]  # take Bx as a reference, but could be any other

    ndim = reference_pd.box.ndim
    if ndim == 1:
        By = patch["By"]
        xby = By.x
        Bz = patch["Bz"]
        xbz = Bz.x
        Jy, Jz = _current1d(By[:], Bz[:], xby, xbz)
        return (
            # Jx is invariant in 1D, so this value is invalid and should not be used
            reference_pd.copy_as(0, name="Jx", centering="dual"),
            By.copy_as(Jy, name="Jy", centering="primal"),
            Bz.copy_as(Jz, name="Jz", centering="primal"),
        )

    elif ndim == 2:
        Bx = patch["Bx"]
        By = patch["By"]
        Bz = patch["Bz"]

        dx, dy = reference_pd.dl
        Jx, Jy, Jz = _current2d(Bx[:], By[:], Bz[:], dx, dy)

        centering = {
            component: [
                reference_pd.layout.centering[direction][component]
                for direction in ("X", "Y")
            ]
            for component in ("Jx", "Jy", "Jz")
        }

        return (
            Bx.copy_as(Jx, name="Jx", centering=centering["Jx"]),
            By.copy_as(Jy, name="Jy", centering=centering["Jy"]),
            Bz.copy_as(Jz, name="Jz", centering=centering["Jz"]),
        )

    raise RuntimeError("dimension not implemented")


def _divB2D(Bx, By, xBx, yBy):
    dxbx = (Bx[1:, :] - Bx[:-1, :]) / (xBx[1] - xBx[0])
    dyby = (By[:, 1:] - By[:, :-1]) / (yBy[1] - yBy[0])
    return dxbx + dyby


def _compute_divB(patch, **kwargs):
    reference_pd = patch["Bx"]  # take Bx as a reference, but could be any other
    ndim = reference_pd.box.ndim

    if ndim == 1:
        raise ValueError("divB is 0 by construction in 1D")

    centering = ["dual"] * ndim

    if ndim == 2:
        By = patch["By"].dataset[:]
        Bx = patch["Bx"].dataset[:]
        xBx = patch["Bx"].x
        yBy = patch["By"].y
        divB = reference_pd.copy_as(
            _divB2D(Bx, By, xBx, yBy), name="divB", centering=centering
        )
        return (divB,)

    raise RuntimeError("dimension not implemented")


from pyphare.pharesee.hierarchy.primal import (  # noqa: F401, E402
    _inner_slices,
    _ppp_to_ppp_domain_slicing,
    _pdd_to_ppp_domain_slicing,
    _dpd_to_ppp_domain_slicing,
    _ddp_to_ppp_domain_slicing,
    _dpp_to_ppp_domain_slicing,
    _pdp_to_ppp_domain_slicing,
    _ppd_to_ppp_domain_slicing,
    _ddd_to_ppp_domain_slicing,
    slices_to_primal_,
    merge_centerings,
    slices_to_primal,
    _compute_to_primal_patch_data,
    _compute_to_primal,
)


def _get_rank(patch, **kwargs):
    """
    make a field dataset cell centered coding the MPI rank
    rank is obtained from patch global id == "rank#local_patch_id"
    """
    from pyphare.core.box import grow

    reference_pd = patch["Bx"]  # Bx as a ref, but could be any other
    ndim = reference_pd.box.ndim

    layout = reference_pd.layout
    centering = ["dual"] * ndim
    nbrGhosts = reference_pd.ghosts_nbr[0]
    shape = grow(reference_pd.box, [nbrGhosts] * 2).shape

    if ndim == 1:
        pass

    elif ndim == 2:
        data = np.zeros(shape) + int(patch.id.strip("p").split("#")[0])
        pd = reference_pd.copy_as(data, name="rank", centering=centering)
        return (pd,)
    else:
        raise RuntimeError("Not Implemented yet")


def _compute_pressure(patch, **kwargs):
    Mxx = patch["Mxx"]
    Mxy = patch["Mxy"]
    Mxz = patch["Mxz"]
    Myy = patch["Myy"]
    Myz = patch["Myz"]
    Mzz = patch["Mzz"]
    massDensity = patch["rho"][:]
    Vix = patch["Vx"][:]
    Viy = patch["Vy"][:]
    Viz = patch["Vz"][:]

    Pxx = Mxx[:] - Vix * Vix * massDensity
    Pxy = Mxy[:] - Vix * Viy * massDensity
    Pxz = Mxz[:] - Vix * Viz * massDensity
    Pyy = Myy[:] - Viy * Viy * massDensity
    Pyz = Myz[:] - Viy * Viz * massDensity
    Pzz = Mzz[:] - Viz * Viz * massDensity

    return (
        Mxx.copy_as(Pxx, name="Pxx"),
        Mxy.copy_as(Pxy, name="Pxy"),
        Mxz.copy_as(Pxz, name="Pxz"),
        Myy.copy_as(Pyy, name="Pyy"),
        Myz.copy_as(Pyz, name="Pyz"),
        Mzz.copy_as(Pzz, name="Pzz"),
    )


def _compute_pop_pressure_xx(patch, mass, popname):
    Mxx = patch[popname + "_Mxx"]
    Fx = patch[popname + "_Fx"][:]
    N = patch[popname + "_rho"][:]
    return {popname + "_Pxx": Mxx.copy_as(Mxx[:] - Fx * Fx * mass / N)}


def _compute_pop_pressure_xy(patch, mass, popname):
    Mxy = patch[popname + "_Mxy"]
    Fx = patch[popname + "_Fx"][:]
    Fy = patch[popname + "_Fy"][:]
    N = patch[popname + "_rho"][:]
    return {popname + "_Pxy": Mxy.copy_as(Mxy[:] - Fx * Fy * mass / N)}


def _compute_pop_pressure_xz(patch, mass, popname):
    Mxz = patch[popname + "_Mxz"]
    Fx = patch[popname + "_Fx"][:]
    Fz = patch[popname + "_Fz"][:]
    N = patch[popname + "_rho"][:]
    return {popname + "_Pxz": Mxz.copy_as(Mxz[:] - Fx * Fz * mass / N)}


def _compute_pop_pressure_yy(patch, mass, popname):
    Myy = patch[popname + "_Myy"]
    Fy = patch[popname + "_Fy"][:]
    N = patch[popname + "_rho"][:]
    return {popname + "_Pyy": Myy.copy_as(Myy[:] - Fy * Fy * mass / N)}


def _compute_pop_pressure_yz(patch, mass, popname):
    Myz = patch[popname + "_Myz"]
    Fy = patch[popname + "_Fy"][:]
    Fz = patch[popname + "_Fz"][:]
    N = patch[popname + "_rho"][:]
    return {popname + "_Pyz": Myz.copy_as(Myz[:] - Fy * Fz * mass / N)}


def _compute_pop_pressure_zz(patch, mass, popname):
    Mzz = patch[popname + "_Mzz"]
    Fz = patch[popname + "_Fz"][:]
    N = patch[popname + "_rho"][:]
    return {popname + "_Pzz": Mzz.copy_as(Mzz[:] - Fz * Fz * mass / N)}


def _compute_pop_pressure(patch, **kwargs):
    """
    computes the pressure tensor for a given population
    this method is different from _compute_pressure in that:
        P = M - V*V*n*mass
    where V is the bulk velocity, n is the density and mass is the particle mass
    but for populations we don't have V we have the Flux.
    so the formula becomes:
        P = M - F*F/N * mass
    """
    popname = kwargs["popname"]
    Mxx = patch[popname + "_Mxx"]
    Mxy = patch[popname + "_Mxy"]
    Mxz = patch[popname + "_Mxz"]
    Myy = patch[popname + "_Myy"]
    Myz = patch[popname + "_Myz"]
    Mzz = patch[popname + "_Mzz"]
    Fx = patch[popname + "_Fx"][:]
    Fy = patch[popname + "_Fy"][:]
    Fz = patch[popname + "_Fz"][:]
    N = patch[popname + "_rho"][:]

    mass = kwargs["mass"]

    Pxx = Mxx[:] - Fx * Fx * mass / N
    Pxy = Mxy[:] - Fx * Fy * mass / N
    Pxz = Mxz[:] - Fx * Fz * mass / N
    Pyy = Myy[:] - Fy * Fy * mass / N
    Pyz = Myz[:] - Fy * Fz * mass / N
    Pzz = Mzz[:] - Fz * Fz * mass / N

    return (
        Mxx.copy_as(Pxx, name=popname + "_Pxx"),
        Mxy.copy_as(Pxy, name=popname + "_Pxy"),
        Mxz.copy_as(Pxz, name=popname + "_Pxz"),
        Myy.copy_as(Pyy, name=popname + "_Pyy"),
        Myz.copy_as(Pyz, name=popname + "_Pyz"),
        Mzz.copy_as(Pzz, name=popname + "_Pzz"),
    )


def make_interpolator(data, coords, interp):
    """Build a scipy interpolator from flat field data and coordinates."""
    dim = coords.ndim

    if dim == 1:
        from scipy.interpolate import interp1d

        return interp1d(
            coords, data, kind=interp, fill_value="extrapolate", assume_sorted=False
        )

    elif dim == 2:
        from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator

        if interp == "nearest":
            return NearestNDInterpolator(coords, data)
        elif interp == "bilinear":
            return LinearNDInterpolator(coords, data)
        else:
            raise ValueError("interp can only be 'nearest' or 'bilinear'")

    else:
        raise ValueError("make_interpolator is not yet 3d")


def finest_coords_for(domain, dl, qty, nbrGhosts, origin=None):
    """Return coordinate arrays for the finest-level output grid."""
    from pyphare.core.gridlayout import yeeCoordsFor

    dim = len(dl)
    if origin is None:
        origin = [0] * dim

    if dim == 1:
        nx = int(np.round(domain[0] / dl[0]))
        return (
            yeeCoordsFor(origin, nbrGhosts[0], dl, [nx], qty, "x", withGhosts=True),
        )

    elif dim == 2:
        nCells = [int(np.round(d / dl_i)) for d, dl_i in zip(domain, dl)]
        x = yeeCoordsFor(origin, nbrGhosts[0], dl, nCells, qty, "x", withGhosts=True)
        y = yeeCoordsFor(origin, nbrGhosts[1], dl, nCells, qty, "y", withGhosts=True)
        return (x, y)

    else:
        raise ValueError("finest_coords_for is not yet 3d")


def interpolate_hierarchy(hier, box=None, quantity=None, interp="nearest"):
    """
    Interpolate a hierarchy onto a uniform finest-level grid, returning UniformGrids.

    box: optional Box with physical coordinates to restrict the output domain.
         If None, uses hier.selection_box if one is set (and it differs from the
         full domain_box), otherwise the full domain is used.
    """

    from pyphare.core.box import Box
    from pyphare.pharesee.hierarchy import hierarchy_utils as hootils
    from pyphare.pharesee.hierarchy import uniformgrid as uniform
    from pyphare.pharesee.hierarchy.func import GetDl, GetDomainSize

    if isinstance(box, (list, tuple)):
        raise TypeError(
            "interpolate_hierarchy takes a single box; call it once per box"
        )

    times = hier.times()
    if len(times) > 1:
        raise ValueError("interpolate_hierarchy does not support multiple times")
    time = times[0]

    if 0 not in hier.levels(time):
        raise ValueError("interpolate_hierarchy only supports coarse timesteps")

    dl = GetDl(hier, time)

    if box is None and hier.selection_box is not None:
        candidate = hier.selection_box.get(time)
        if candidate is not None and candidate != hier.domain_box:
            box = candidate

    nlevels = len(hier.levels(time))

    if box is not None:
        origin = np.array(box.lower, dtype=float)
        domain = np.array(box.upper - box.lower, dtype=float)
        n_cells = np.round(domain / dl).astype(int)
        layout_box = Box(np.zeros(hier.ndim, dtype=int), n_cells - 1)
    else:
        origin = np.zeros(hier.ndim)
        domain = GetDomainSize(hier)
        layout_box = hier.level_domain_box(nlevels - 1)

    ref_patch = hier.levels(time)[0].patches[0]
    layout = ref_patch.layout.copy_as(box=layout_box, origin=list(origin), dl=dl)

    datas = {}
    for qty in [quantity] if quantity else hier.quantities():
        nbrGhosts = ref_patch[qty].ghosts_nbr
        data, coords = hootils.flat_finest_field(hier, qty, time=time)
        interpolator = make_interpolator(data, coords, interp)
        finest_coords = finest_coords_for(domain, dl, qty, nbrGhosts, origin)
        mesh = np.meshgrid(*finest_coords, indexing="ij")
        datas[qty] = uniform.UniformGrid(
            layout=layout,
            name=ref_patch[qty].name,
            data=interpolator(*mesh),
            ghosts_nbr=nbrGhosts,
            centering=ref_patch[qty].centerings,
        )

    return uniform.UniformGrids(datas)
