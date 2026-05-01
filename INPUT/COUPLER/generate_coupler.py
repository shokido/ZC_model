import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as ncdf
fname_ocn="../OCN/grid_eqpac_30.nc"
fname_atm="../ATM/grid_256.nc"
fname_coupler="coupler_eqpac_30_atm_256.nc"
nc_ocn=ncdf.Dataset(fname_ocn,"r")
lon_ocn=nc_ocn.variables["lon_p"][:]
lat_ocn=nc_ocn.variables["lat_p"][:]
mask_ocn=nc_ocn.variables["mask_sst"][:]
nc_ocn.close()

nc_atm=ncdf.Dataset(fname_atm,"r")
lon_atm=nc_atm.variables["lon"][:]
lat_atm=nc_atm.variables["lat"][:]
nc_atm.close()

def _wrap_lon_360(lon):
    return np.mod(np.asarray(lon, dtype=np.float64), 360.0)

def _bracket_indices_1d(grid_sorted, values):
    """
    For each value, find left index i0 such that grid[i0] <= v <= grid[i0+1],
    clipped to interior [0, n-2]. Returns i0 and fractional t in [0,1].
    grid_sorted must be increasing.
    """
    g = grid_sorted
    v = np.asarray(values, dtype=np.float64)

    # insertion point: first index where g[idx] >= v
    idx = np.searchsorted(g, v, side="left")
    idx = np.clip(idx, 1, len(g) - 1)
    i0 = idx - 1
    i1 = idx

    g0 = g[i0]
    g1 = g[i1]

    # avoid divide by zero (should not happen for proper grids)
    denom = (g1 - g0)
    denom = np.where(denom == 0.0, np.nan, denom)

    t = (v - g0) / denom
    t = np.clip(t, 0.0, 1.0)
    return i0.astype(np.int32), t.astype(np.float64)
def _bracket_indices_1d_increasing(grid, values):
    """
    grid: 1D increasing array
    values: array-like

    Returns
    -------
    i0 : int32 array
        left index such that grid[i0] <= v <= grid[i0+1] (clipped internally)
    t : float64 array
        fraction in [0,1] inside the bracket interval
    inside : bool array
        True if v is within [grid[0], grid[-1]] (inclusive)
    """
    g = np.asarray(grid, dtype=np.float64)
    v = np.asarray(values, dtype=np.float64)

    inside = (v >= g[0]) & (v <= g[-1])

    # compute bracket indices (safe even if outside; we'll mask by inside later)
    idx = np.searchsorted(g, v, side="left")
    idx = np.clip(idx, 1, len(g) - 1)
    i0 = (idx - 1).astype(np.int32)
    i1 = i0 + 1

    g0 = g[i0]
    g1 = g[i1]
    denom = g1 - g0
    denom = np.where(denom == 0.0, np.nan, denom)

    t = (v - g0) / denom
    t = np.clip(t, 0.0, 1.0)

    return i0, t, inside
def _prepare_monotonic_1d(grid):
    g = np.asarray(grid, dtype=np.float64)
    order = np.argsort(g)
    return g[order], order
def get_AtoO_bilinear(lon_ocn, lat_ocn, mask_ocn, lon_atm, lat_atm, *,
                      nearest_only=False,
                      periodic_lon=True,
                      fill_value=1):
    """
    Ocean (mask==1) -> Atmos (4 corners) bilinear mapping on rectilinear ATM grid.

    Assumptions:
      - lon_atm, lat_atm are 1D (rectilinear atmosphere grid)
      - mask_ocn is 2D (ny,nx)
      - lon_ocn/lat_ocn are either both 1D (nx),(ny) or both 2D (ny,nx)
      - weights are based on lon/lat linear fractions (classic bilinear)

    Returns (Fortran 1-based indices):
      ind_AtoO_x: (4,ny,nx) int32
      ind_AtoO_y: (4,ny,nx) int32
      wgt_AtoO  : (4,ny,nx) float64
    Slot convention (k = 0..3):
      0: (i0, j0)
      1: (i1, j0)
      2: (i0, j1)
      3: (i1, j1)
    """
    lon_ocn = np.asarray(lon_ocn, dtype=np.float64)
    lat_ocn = np.asarray(lat_ocn, dtype=np.float64)
    mask_ocn = np.asarray(mask_ocn)

    lon_atm = np.asarray(lon_atm, dtype=np.float64)
    lat_atm = np.asarray(lat_atm, dtype=np.float64)

    # ocean lon/lat -> 2D
    if lon_ocn.ndim == 1 and lat_ocn.ndim == 1:
        lon2_o, lat2_o = np.meshgrid(lon_ocn, lat_ocn)  # (ny,nx)
    else:
        lon2_o, lat2_o = lon_ocn, lat_ocn

    ny, nx = mask_ocn.shape

    ind_x = np.full((4, ny, nx), fill_value, dtype=np.int32)
    ind_y = np.full((4, ny, nx), fill_value, dtype=np.int32)
    wgt   = np.zeros((4, ny, nx), dtype=np.float64)

    tgt = (mask_ocn == 1)
    if not np.any(tgt):
        return ind_x, ind_y, wgt

    # --- prepare atm grids (sorted) ---
    # lat: just sort
    lat_sorted, lat_order = _prepare_monotonic_1d(lat_atm)

    # lon: optionally wrap to 0..360 and sort
    if periodic_lon:
        lon_wrapped = _wrap_lon_360(lon_atm)
        lon_sorted, lon_order = _prepare_monotonic_1d(lon_wrapped)
        lon_q = _wrap_lon_360(lon2_o[tgt])
    else:
        lon_sorted, lon_order = _prepare_monotonic_1d(lon_atm)
        lon_q = lon2_o[tgt]

    lat_q = lat2_o[tgt]

    # --- bracketing indices in sorted space ---
    i0_s, tx = _bracket_indices_1d(lon_sorted, lon_q)
    j0_s, ty = _bracket_indices_1d(lat_sorted, lat_q)

    i1_s = i0_s + 1
    j1_s = j0_s + 1

    # map sorted indices -> original indices (0-based)
    i0 = lon_order[i0_s]
    i1 = lon_order[i1_s]
    j0 = lat_order[j0_s]
    j1 = lat_order[j1_s]

    # bilinear weights
    w00 = (1.0 - tx) * (1.0 - ty)
    w10 = tx * (1.0 - ty)
    w01 = (1.0 - tx) * ty
    w11 = tx * ty

    # If nearest_only: pick the max weight corner and make it 1, others 0
    if nearest_only:
        W = np.stack([w00, w10, w01, w11], axis=0)  # (4, Npts)
        kmax = np.argmax(W, axis=0)                 # (Npts,)
        w00 = (kmax == 0).astype(np.float64)
        w10 = (kmax == 1).astype(np.float64)
        w01 = (kmax == 2).astype(np.float64)
        w11 = (kmax == 3).astype(np.float64)

    # write back into (4,ny,nx) arrays at tgt locations
    # slot 0: (i0,j0)
    ind_x[0][tgt] = i0 + 1
    ind_y[0][tgt] = j0 + 1
    wgt[0][tgt]   = w00

    # slot 1: (i1,j0)
    ind_x[1][tgt] = i1 + 1
    ind_y[1][tgt] = j0 + 1
    wgt[1][tgt]   = w10

    # slot 2: (i0,j1)
    ind_x[2][tgt] = i0 + 1
    ind_y[2][tgt] = j1 + 1
    wgt[2][tgt]   = w01

    # slot 3: (i1,j1)
    ind_x[3][tgt] = i1 + 1
    ind_y[3][tgt] = j1 + 1
    wgt[3][tgt]   = w11

    return ind_x, ind_y, wgt
def get_OtoA_bilinear(lon_ocn, lat_ocn, mask_ocn, lon_atm, lat_atm, *,
                      periodic_lon=True,
                      nearest_only=False,
                      fill_index=1):
    """
    Ocean -> Atmos mapping using 4-point bilinear interpolation on a rectilinear
    OCEAN grid. Land/undef is handled by mask_ocn (1=valid ocean, 0=invalid).

    IMPORTANT POLICY:
      - No extrapolation. Atmos points outside the ocean lon/lat ranges are invalid:
        indices remain fill_index and weights are all 0.

    Parameters
    ----------
    lon_ocn, lat_ocn : 1D arrays (increasing). Ocean grid coordinates.
    mask_ocn : 2D array (ny_o, nx_o). 1 = sea, 0 = land/invalid.
    lon_atm, lat_atm : 1D arrays (increasing). Atmos grid coordinates (targets).
    periodic_lon : bool
        If True, longitudes are wrapped to [0,360). Use this when your grids are in 0..360.
        NOTE: This assumes lon_ocn (after wrap) is still increasing and continuous (no dateline break inside array).
    nearest_only : bool
        If True, select only one corner among the valid ones (max renormalized weight):
          wgt[0]=1, wgt[1:]=0 for that chosen slot (slot index varies internally but we store it in 0..3 weights).
    fill_index : int
        Index value for invalid/unmapped points 

    Returns
    -------
    ind_OtoA_x : (4, ny_a, nx_a) int32, Fortran 1-based i index into ocean grid
    ind_OtoA_y : (4, ny_a, nx_a) int32, Fortran 1-based j index into ocean grid
    wgt_OtoA   : (4, ny_a, nx_a) float64, weights (sum to 1 where valid; 0 where invalid)
    Slot convention:
      0: (i0, j0)
      1: (i1, j0)
      2: (i0, j1)
      3: (i1, j1)
    """
    lon_ocn = np.asarray(lon_ocn, dtype=np.float64)
    lat_ocn = np.asarray(lat_ocn, dtype=np.float64)
    mask_ocn = np.asarray(mask_ocn)
    lon_atm = np.asarray(lon_atm, dtype=np.float64)
    lat_atm = np.asarray(lat_atm, dtype=np.float64)

    if lon_ocn.ndim != 1 or lat_ocn.ndim != 1:
        raise ValueError("lon_ocn and lat_ocn must be 1D for rectilinear ocean grid.")
    if lon_atm.ndim != 1 or lat_atm.ndim != 1:
        raise ValueError("lon_atm and lat_atm must be 1D for rectilinear atmosphere grid.")
    if mask_ocn.shape != (len(lat_ocn), len(lon_ocn)):
        raise ValueError("mask_ocn shape must be (len(lat_ocn), len(lon_ocn)).")

    ny_a = len(lat_atm)
    nx_a = len(lon_atm)

    ind_x = np.full((4, ny_a, nx_a), fill_index, dtype=np.int32)
    ind_y = np.full((4, ny_a, nx_a), fill_index, dtype=np.int32)
    wgt   = np.zeros((4, ny_a, nx_a), dtype=np.float64)

    # Atmos target mesh
    lon2_a, lat2_a = np.meshgrid(lon_atm, lat_atm)  # (ny_a, nx_a)

    # Longitude wrap (no extrapolation; only affects coordinate system)
    if periodic_lon:
        lon_ocn_w = _wrap_lon_360(lon_ocn)
        lon_q = _wrap_lon_360(lon2_a)
    else:
        lon_ocn_w = lon_ocn
        lon_q = lon2_a

    # Bracket indices in ocean grid
    i0, tx, inside_lon = _bracket_indices_1d_increasing(lon_ocn_w, lon_q)
    j0, ty, inside_lat = _bracket_indices_1d_increasing(lat_ocn,   lat2_a)

    i1 = i0 + 1
    j1 = j0 + 1

    # Only points inside both lon/lat ranges are eligible (no extrapolation)
    inside = inside_lon & inside_lat

    # Bilinear weights (computed everywhere, but will be masked by inside & mask_ocn)
    w00 = (1.0 - tx) * (1.0 - ty)
    w10 = tx * (1.0 - ty)
    w01 = (1.0 - tx) * ty
    w11 = tx * ty

    # Validity from ocean mask at the four corners
    v00 = (mask_ocn[j0, i0] == 1) & inside
    v10 = (mask_ocn[j0, i1] == 1) & inside
    v01 = (mask_ocn[j1, i0] == 1) & inside
    v11 = (mask_ocn[j1, i1] == 1) & inside

    # Masked weights
    W0 = np.where(v00, w00, 0.0)
    W1 = np.where(v10, w10, 0.0)
    W2 = np.where(v01, w01, 0.0)
    W3 = np.where(v11, w11, 0.0)

    Wsum = W0 + W1 + W2 + W3
    ok = (Wsum > 0.0)  # at least one valid ocean corner exists

    # Renormalize among valid corners
    W0n = np.zeros_like(W0); W1n = np.zeros_like(W1); W2n = np.zeros_like(W2); W3n = np.zeros_like(W3)
    W0n[ok] = W0[ok] / Wsum[ok]
    W1n[ok] = W1[ok] / Wsum[ok]
    W2n[ok] = W2[ok] / Wsum[ok]
    W3n[ok] = W3[ok] / Wsum[ok]

    # nearest_only: keep only the max-weight valid corner
    if nearest_only:
        stackW = np.stack([W0n, W1n, W2n, W3n], axis=0)  # (4, ny_a, nx_a)
        kmax = np.argmax(stackW, axis=0)                # (ny_a, nx_a)

        W0n = ((kmax == 0) & ok).astype(np.float64)
        W1n = ((kmax == 1) & ok).astype(np.float64)
        W2n = ((kmax == 2) & ok).astype(np.float64)
        W3n = ((kmax == 3) & ok).astype(np.float64)

    # Fill indices and weights where ok
    i0F = i0 ; i1F = i1 
    j0F = j0 ; j1F = j1 

    ind_x[0][ok] = i0F[ok]; ind_y[0][ok] = j0F[ok]; wgt[0][ok] = W0n[ok]
    ind_x[1][ok] = i1F[ok]; ind_y[1][ok] = j0F[ok]; wgt[1][ok] = W1n[ok]
    ind_x[2][ok] = i0F[ok]; ind_y[2][ok] = j1F[ok]; wgt[2][ok] = W2n[ok]
    ind_x[3][ok] = i1F[ok]; ind_y[3][ok] = j1F[ok]; wgt[3][ok] = W3n[ok]

    # outside or no-valid-corner points remain fill_index with wgt=0
    return ind_x, ind_y, wgt
ind_AtoO_x,ind_AtoO_y,wgt_AtoO=get_AtoO_bilinear(lon_ocn, lat_ocn, mask_ocn, lon_atm, lat_atm)
ind_OtoA_x,ind_OtoA_y,wgt_OtoA=get_OtoA_bilinear(lon_ocn, lat_ocn, mask_ocn, lon_atm, lat_atm)
kcoord = np.arange(4, dtype=np.int32)  # 0,1,2,3

with ncdf.Dataset(fname_coupler, "w", format="NETCDF4") as nc:
    # Dimensions
    nc.createDimension("k", 4)
    nc.createDimension("lon_atm", len(lon_atm))
    nc.createDimension("lat_atm", len(lat_atm))
    nc.createDimension("lon_ocn", len(lon_ocn))
    nc.createDimension("lat_ocn", len(lat_ocn))

    # Coordinates
    v_k = nc.createVariable("k", "i4", ("k",))
    v_k.long_name = "corner index (0..3)"
    v_k[:] = kcoord

    v_lon_atm = nc.createVariable("lon_atm", "f8", ("lon_atm",))
    v_lon_atm.units = "degrees_east"
    v_lon_atm[:] = lon_atm

    v_lat_atm = nc.createVariable("lat_atm", "f8", ("lat_atm",))
    v_lat_atm.units = "degrees_north"
    v_lat_atm[:] = lat_atm

    v_lon_ocn = nc.createVariable("lon_ocn", "f8", ("lon_ocn",))
    v_lon_ocn.units = "degrees_east"
    v_lon_ocn[:] = lon_ocn

    v_lat_ocn = nc.createVariable("lat_ocn", "f8", ("lat_ocn",))
    v_lat_ocn.units = "degrees_north"
    v_lat_ocn[:] = lat_ocn

    # Mapping variables (A->O): (k, lat_ocn, lon_ocn)
    v_ind_AtoO_x = nc.createVariable("ind_AtoO_x", "i4", ("k", "lat_ocn", "lon_ocn"), zlib=True)
    v_ind_AtoO_y = nc.createVariable("ind_AtoO_y", "i4", ("k", "lat_ocn", "lon_ocn"), zlib=True)
    v_wgt_AtoO   = nc.createVariable("wgt_AtoO",   "f4", ("k", "lat_ocn", "lon_ocn"), zlib=True)

    v_ind_AtoO_x.long_name = "ATM i-index for A->O (Fortran 1-based; fill=0 means invalid)"
    v_ind_AtoO_y.long_name = "ATM j-index for A->O (Fortran 1-based; fill=0 means invalid)"
    v_wgt_AtoO.long_name   = "weights for A->O bilinear interpolation (sum over k is 1 where valid)"

    v_ind_AtoO_x[:] = ind_AtoO_x
    v_ind_AtoO_y[:] = ind_AtoO_y
    v_wgt_AtoO[:]   = wgt_AtoO

    # Mapping variables (O->A): (k, lat_atm, lon_atm)
    v_ind_OtoA_x = nc.createVariable("ind_OtoA_x", "i4", ("k", "lat_atm", "lon_atm"), zlib=True)
    v_ind_OtoA_y = nc.createVariable("ind_OtoA_y", "i4", ("k", "lat_atm", "lon_atm"), zlib=True)
    v_wgt_OtoA   = nc.createVariable("wgt_OtoA",   "f4", ("k", "lat_atm", "lon_atm"), zlib=True)

    v_ind_OtoA_x.long_name = "OCN i-index for O->A (Fortran 1-based; fill=0 means invalid)"
    v_ind_OtoA_y.long_name = "OCN j-index for O->A (Fortran 1-based; fill=0 means invalid)"
    v_wgt_OtoA.long_name   = "weights for O->A bilinear interpolation (mask-aware; sum over k is 1 where valid)"
    v_ind_OtoA_x[:] = ind_OtoA_x
    v_ind_OtoA_y[:] = ind_OtoA_y
    v_wgt_OtoA[:]   = wgt_OtoA

    # Global attrs
    nc.title = "Zebiak-Cane coupler mapping indices and weights"
    nc.history = "Created by Python (netCDF4)"
    nc.comment = "A->O: Atmos->Ocean; O->A: Ocean->Atmos. Indices are Fortran 1-based for A-O, and 0-based for O-A. Invalid points have ind=0, wgt=0."



