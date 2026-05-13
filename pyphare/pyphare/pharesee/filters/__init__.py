#
#
#

import numpy as np
from copy import deepcopy


def gaussian(hier, qty=None, sigma=2, **kwargs):
    from pyphare.pharesee import hierarchy as harch

    time = harch.func.GetTime(hier)
    finest = harch.func.GetFinest(hier, time, qty)
    grids = deepcopy(finest)
    for key, grid in finest.items():
        grids[key] = gaussian_filter_uniform_grid(grid, sigma=sigma)
    return grids if len(grids) > 1 else next(iter(grids.values()))


def gaussian_filter_uniform_grid(grid, sigma=2):
    from scipy.ndimage import gaussian_filter

    ndim = grid.box.ndim
    nb_ghosts = grid.ghosts_nbr[0]
    ds = np.asarray(grid[:], dtype=float)
    nan_mask = np.isnan(ds)
    filled = np.where(nan_mask, 0.0, ds)
    weights = np.where(nan_mask, 0.0, 1.0)
    val_filt = gaussian_filter(filled, sigma=sigma)
    w_filt = gaussian_filter(weights, sigma=sigma)
    result = np.where(w_filt > 1e-6, val_filt / w_filt, np.nan)
    select = tuple([slice(nb_ghosts or None, -nb_ghosts or None) for _ in range(ndim)])
    ds_ = np.full(ds.shape, np.nan)
    ds_[select] = result[select]
    copy = deepcopy(grid)
    copy.dataset = ds_
    return copy
