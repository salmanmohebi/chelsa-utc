"""SAGA-exact numpy versions of the operators used by CHELSA V2.1 daily.

CHELSA V2.1 runs SAGA GIS 7.6 tools (``chelsa_main.py``, ``saga_functions.py``). To
rebuild its output exactly we copy what SAGA does, not what a "normal" interpolation
would do. Verified against the published tifs (to the uint16 step), including the
dateline, small islands and 84 N; see ``docs/method.md`` and ``docs/validation.md``.

Grids
-----
* CHELSA 30" grid (``dem_latlong.sgrd``): NX x NY = 43200 x 20880 pixel centres
  ``x_i = XMIN + i * CS``, ``y_r = Y_TOP - r * CS`` (r = tif row, north -> south).
  ``CS`` is SAGA's cell size 0.0083333333, not exactly 1/120.
* ERA5 0.25 deg node grid after ``change_latlong(PATCH=True)``: 721 x 1441 nodes,
  lat -90..90 ascending, lon -180..180 (column 1440 = copy of column 0). Use ``patch``.

Sampling (``NodeGrid.sample``) = SAGA ``CSG_Grid::Get_Value``:
* the target must lie inside the cell-edge extent, and its NEAREST node must be in
  the grid and not NoData; otherwise NoData;
* 4 x 4 stencil (ix-1 .. ix+2); positions outside the grid or NoData are filled by
  ``_Get_ValAtPos_Fill4x4Submatrix`` (Jacobi passes, each takes the mean of the
  3 x 3 neighbours that are in the grid or already filled);
* kernels: ``"bicubic"`` = SAGA "Bicubic Spline" (Catmull-Rom, a = -0.5), used for
  temperature; ``"bspline"`` = SAGA "B-Spline" (cubic B-spline, no prefilter).
There is no periodic longitude wrap.

The fast path is a separable gather on the plain grid. Pixels whose stencil leaves
the grid or touches NoData are recomputed with the exact per-stencil fill
(``_exact_stencil``); there are few of them (grid edges, NoData borders).
"""
from __future__ import annotations

import numpy as np
from scipy import sparse
from scipy.ndimage import convolve1d

CS = 0.0083333333
XMIN = -179.9959722222
YMIN = -89.9959722222
NX, NY = 43200, 20880
Y_TOP = YMIN + (NY - 1) * CS
X30 = XMIN + np.arange(NX) * CS            # pixel-centre longitudes, west -> east


def rows_y(r0: int, r1: int) -> np.ndarray:
    """Pixel-centre latitudes of tif rows r0 .. r1-1 (north -> south)."""
    return Y_TOP - np.arange(r0, r1) * CS


def patch(g: np.ndarray) -> np.ndarray:
    """0.25 deg (721, 1440) lon -180..179.75 -> (721, 1441) with the +180 column."""
    return np.concatenate([g, g[:, :1]], axis=1)


def coarse(g: np.ndarray, patched: bool = True) -> "NodeGrid":
    """A global 0.25 deg node grid, lat ascending from -90, lon from -180."""
    G = patch(g) if patched else np.asarray(g)
    return NodeGrid(G, 0, 0, G.shape[1], G.shape[0], -180.0, -90.0, 0.25)


# ── kernels ───────────────────────────────────────────────────────────────────

def kernel_w(d: np.ndarray, method: str) -> np.ndarray:
    """(n, 4) weights of stencil positions -1, 0, 1, 2 for fractional offsets d."""
    d = np.asarray(d, dtype=np.float64)
    if method == "bicubic":   # SAGA BiCubicSpline macro: v1 + .5 d (v2 - v0 + d (2v0 - 5v1 + 4v2 - v3 + d (3(v1 - v2) + v3 - v0)))
        return np.stack([0.5 * (-d + 2 * d ** 2 - d ** 3), 1 + 0.5 * (-5 * d ** 2 + 3 * d ** 3),
                         0.5 * (d + 4 * d ** 2 - 3 * d ** 3), 0.5 * (-d ** 2 + d ** 3)], -1)
    if method == "bspline":   # SAGA grid.cpp _Get_ValAtPos_BSpline
        out = []
        for i in range(4):
            s = np.zeros_like(d)
            for c, k in ((1.0, 1.0), (-4.0, 0.0), (6.0, -1.0), (-4.0, -2.0)):
                e = i - d + k
                s += np.where(e > 0, c * e ** 3, 0.0)
            out.append(s / 6.0)
        return np.stack(out, -1)
    raise ValueError(method)


# ── grids ─────────────────────────────────────────────────────────────────────

class NodeGrid:
    """A window of a global SAGA grid.

    ``G[a, b]`` is the node (column ``b + i_off``, row ``a + j_off``) of a global grid
    with ``NXg x NYg`` nodes at ``x0 + i * cs``, ``y0 + j * cs`` (rows ascending in y).
    NaN = NoData. Nodes outside the window must not be needed (asserted).
    """

    def __init__(self, G, i_off, j_off, NXg, NYg, x0, y0, cs):
        self.G = np.asarray(G, dtype=np.float64)
        self.i_off, self.j_off, self.NXg, self.NYg = int(i_off), int(j_off), int(NXg), int(NYg)
        self.x0, self.y0, self.cs = float(x0), float(y0), float(cs)

    # global node value or NaN (outside the global grid or NoData)
    def _val(self, gx: int, gy: int) -> float:
        if gx < 0 or gy < 0 or gx >= self.NXg or gy >= self.NYg:
            return np.nan
        a, b = gy - self.j_off, gx - self.i_off
        assert 0 <= a < self.G.shape[0] and 0 <= b < self.G.shape[1], ("window too small", gx, gy)
        return self.G[a, b]

    def _exact_stencil(self, ix: int, iy: int):
        """4x4 stencil v[a, b] (a = x offset, b = y offset) with SAGA's NoData fill, or None."""
        v = np.full((4, 4), np.nan)
        for a in range(4):
            for b in range(4):
                v[a, b] = self._val(ix - 1 + a, iy - 1 + b)
        nnd = int(np.isnan(v).sum())
        it = 0
        while 0 < nnd < 16 and it < 16:
            t = v.copy()
            for a in range(4):
                for b in range(4):
                    if np.isnan(t[a, b]):
                        s, n = 0.0, 0
                        for jb in range(b - 1, b + 2):
                            for ja in range(a - 1, a + 2):
                                g = self._val(ix - 1 + ja, iy - 1 + jb)
                                if np.isfinite(g):
                                    s += g; n += 1
                                elif 0 <= ja < 4 and 0 <= jb < 4 and np.isfinite(t[ja, jb]):
                                    s += t[ja, jb]; n += 1
                        if n > 0:
                            v[a, b] = s / n; nnd -= 1
            it += 1
        return v if nnd == 0 else None

    def sample(self, xt, yt, method: str) -> np.ndarray:
        """Values at targets xt (nc,) x yt (nr,) -> (nr, nc), SAGA Get_Value semantics."""
        xt = np.asarray(xt, dtype=np.float64); yt = np.asarray(yt, dtype=np.float64)
        fx = (xt - self.x0) / self.cs; fy = (yt - self.y0) / self.cs
        ix = np.floor(fx).astype(np.int64); iy = np.floor(fy).astype(np.int64)
        dx = fx - ix; dy = fy - iy
        wx = kernel_w(dx, method); wy = kernel_w(dy, method)
        G, NXg, NYg = self.G, self.NXg, self.NYg

        # fast path: gather on the plain grid (indices clipped), NaN propagates
        cols = ix[:, None] + np.arange(-1, 3)                 # global
        rows = iy[:, None] + np.arange(-1, 3)
        col_out = ((cols < 0) | (cols >= NXg)).any(1)
        row_out = ((rows < 0) | (rows >= NYg)).any(1)
        cl = np.clip(cols, 0, NXg - 1) - self.i_off
        rl = np.clip(rows, 0, NYg - 1) - self.j_off
        assert cl.min() >= 0 and cl.max() < G.shape[1] and rl.min() >= 0 and rl.max() < G.shape[0], \
            "NodeGrid window too small for the targets"
        j0, j1 = int(rl.min()), int(rl.max()) + 1
        A = np.zeros((j1 - j0, len(xt)))
        for k in range(4):
            A += G[j0:j1, cl[:, k]] * wx[:, k]
        out = np.zeros((len(yt), len(xt)))
        for k in range(4):
            out += wy[:, k][:, None] * A[rl[:, k] - j0]

        # exact path where the stencil leaves the grid or touches NoData
        bad = np.isnan(out) | row_out[:, None] | col_out[None, :]
        # cell-edge extent and nearest-node check
        h = 0.5 * self.cs
        inx = (xt >= self.x0 - h) & (xt <= self.x0 + (NXg - 1) * self.cs + h)
        iny = (yt >= self.y0 - h) & (yt <= self.y0 + (NYg - 1) * self.cs + h)
        nxn = ix + (0.5 + dx).astype(np.int64); nyn = iy + (0.5 + dy).astype(np.int64)
        okx = inx & (nxn >= 0) & (nxn < NXg); oky = iny & (nyn >= 0) & (nyn < NYg)
        near = np.zeros((len(yt), len(xt)), bool)
        if okx.any() and oky.any():
            fin = np.isfinite(G)
            near[np.ix_(oky, okx)] = fin[np.clip(nyn, 0, NYg - 1)[oky] - self.j_off][
                :, np.clip(nxn, 0, NXg - 1)[okx] - self.i_off]
        redo = bad & near
        if redo.any():
            rr, cc = np.nonzero(redo)
            gx_, gy_ = ix[cc], iy[rr]
            key = (gy_ - gy_.min()) * (int(gx_.max() - gx_.min()) + 1) + (gx_ - gx_.min())
            order = np.argsort(key, kind="stable"); ks = key[order]
            starts = np.r_[0, np.nonzero(np.diff(ks))[0] + 1]; ends = np.r_[starts[1:], len(ks)]
            vals = np.full(len(rr), np.nan)
            for s, e in zip(starts, ends):
                sel = order[s:e]
                v = self._exact_stencil(int(gx_[sel[0]]), int(gy_[sel[0]]))
                if v is not None:
                    vals[sel] = np.einsum("na,ab,nb->n", wx[cc[sel]], v, wy[rr[sel]])
            out[rr, cc] = vals
        out[~near] = np.nan
        return out


# ── 30" operators on full-width bands of tif rows ────────────────────────────

def _edge_mean(v: np.ndarray) -> np.ndarray:
    m = np.full(v.shape, np.nan)
    m[1:-1] = (v[:-2] + v[1:-1] + v[2:]) / 3
    return m


def same_grid_bspline(F: np.ndarray, r0: int, r1: int) -> np.ndarray:
    """SAGA samples XGRIDS with the B-spline even on the same grid: (1/6, 2/3, 1/6)
    smoothing. F holds full-width rows r0..r1-1. At the true grid edges the stencil
    position outside is filled with the 3-neighbour mean along the edge; at band
    borders and next to NoData the result is NaN."""
    P = np.pad(F, 1, constant_values=np.nan)
    P[1:-1, 0] = _edge_mean(F[:, 0])          # lon -180 edge
    P[1:-1, -1] = _edge_mean(F[:, -1])        # lon +180 edge
    if r0 == 0:
        P[0, 1:-1] = _edge_mean(F[0, :])
    if r1 == NY:
        P[-1, 1:-1] = _edge_mean(F[-1, :])
    k = np.array([1 / 6, 2 / 3, 1 / 6])
    bad = ~np.isfinite(P)
    out = convolve1d(convolve1d(np.where(bad, 0.0, P), k, axis=0, mode="constant"), k, axis=1, mode="constant")
    b2 = convolve1d(convolve1d(bad.astype(np.float64), np.ones(3), axis=0, mode="constant"),
                    np.ones(3), axis=1, mode="constant")
    out[b2 > 0] = np.nan
    return out[1:-1, 1:-1]


# coarse cell of every 30" pixel by nearest node (SAGA resampling assignment)
KX = np.floor((X30 + 180.0) / 0.25 + 0.5).astype(np.int64)                  # 0 .. 1440
KY = np.floor((rows_y(0, NY) + 90.0) / 0.25 + 0.5).astype(np.int64)         # per tif row


def _complete_rows_nearest(r0: int, r1: int) -> np.ndarray:
    """ky of cells whose pixel rows (nearest assignment) all lie in r0..r1-1."""
    ks = np.unique(KY[r0:r1])
    first = {k: int(np.argmax(KY == k)) for k in ks}
    last = {k: NY - 1 - int(np.argmax(KY[::-1] == k)) for k in ks}
    return np.array([k for k in ks if first[k] >= r0 and last[k] < r1], dtype=np.int64)


def block_max(F: np.ndarray, r0: int, r1: int) -> NodeGrid:
    """Max of F per 0.25 deg cell (nearest-node assignment, half cells at +-180).
    F = full-width rows r0..r1-1. Returns a global (721, 1441) NodeGrid, NaN where the
    cell is not complete inside the band or has no pixels."""
    kyb = KY[r0:r1]
    cstarts = np.r_[0, np.nonzero(np.diff(KX))[0] + 1]
    colmax = np.fmax.reduceat(F, cstarts, axis=1)                   # (rows, 1441)
    rstarts = np.r_[0, np.nonzero(np.diff(kyb))[0] + 1]
    cellmax = np.fmax.reduceat(colmax, rstarts, axis=0)             # (n ky in band, 1441), ky descending
    out = np.full((721, 1441), np.nan)
    ok = set(_complete_rows_nearest(r0, r1).tolist())
    for a, s in enumerate(rstarts):
        k = int(kyb[s])
        if k in ok:
            out[k] = cellmax[a]
    return NodeGrid(out, 0, 0, 1441, 721, -180.0, -90.0, 0.25)


def _overlap(centres: np.ndarray, n_cells: int, origin: float) -> sparse.csr_matrix:
    """(n_cells, len(centres)) overlap fraction of each pixel (width CS) with each
    0.25 deg cell centred at origin + k * 0.25."""
    lo = centres - CS / 2; hi = centres + CS / 2
    k0 = np.floor((centres - origin) / 0.25 + 0.5).astype(np.int64)
    rows, cols, vals = [], [], []
    for dk in (-1, 0, 1):
        k = k0 + dk
        clo = origin + k * 0.25 - 0.125; chi = clo + 0.25
        w = (np.minimum(hi, chi) - np.maximum(lo, clo)) / CS
        m = (w > 0) & (k >= 0) & (k < n_cells)
        rows.append(k[m]); cols.append(np.nonzero(m)[0]); vals.append(w[m])
    return sparse.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                             shape=(n_cells, len(centres)))


OX = _overlap(X30, 1441, -180.0)


def area_mean(F: np.ndarray, r0: int, r1: int) -> NodeGrid:
    """Cell-area-weighted mean of F per 0.25 deg cell with the exact pixel geometry
    (SAGA resampling, mean). F = full-width rows r0..r1-1. Cells whose overlapping
    rows are not all inside the band are NaN, except at the true grid edges."""
    y = rows_y(r0, r1)
    OY = _overlap(y, 721, -90.0)                                     # (721, rows)
    v = np.isfinite(F)
    num = OY @ (OX @ np.where(v, F, 0.0).T).T                        # (721, 1441)
    den = OY @ (OX @ v.astype(np.float64).T).T
    out = np.where(den > 0, num / np.where(den > 0, den, 1.0), np.nan)
    # completeness: rows overlapping cell k must all lie in the band (or beyond the grid)
    yall_lo = rows_y(0, NY) - CS / 2; yall_hi = rows_y(0, NY) + CS / 2
    for k in np.nonzero(np.isfinite(out).any(1))[0]:
        clo = -90.0 + k * 0.25 - 0.125; chi = clo + 0.25
        rr = np.nonzero((yall_hi > clo) & (yall_lo < chi))[0]
        if rr.size and (rr.min() < r0 or rr.max() >= r1):
            out[k] = np.nan
    return NodeGrid(out, 0, 0, 1441, 721, -180.0, -90.0, 0.25)


def to_u16(values: np.ndarray, scale: float) -> np.ndarray:
    """CHELSA storage: round(value / scale) as uint16, NaN -> 65535."""
    q = np.rint(values / scale)
    out = np.full(values.shape, 65535, dtype=np.uint16)
    ok = np.isfinite(q)
    out[ok] = np.clip(q[ok], 0, 65534).astype(np.uint16)
    return out
