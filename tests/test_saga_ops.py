import numpy as np
import pytest

from chelsa_utc import saga_ops as so


@pytest.mark.parametrize("method", ["bicubic", "bspline"])
def test_kernel_weights_sum_to_one(method):
    d = np.linspace(0, 1, 11, endpoint=False)
    assert np.allclose(so.kernel_w(d, method).sum(1), 1.0)


def test_bicubic_interpolates_nodes():
    w = so.kernel_w(np.array([0.0]), "bicubic")[0]
    assert np.allclose(w, [0, 1, 0, 0])


@pytest.mark.parametrize("method", ["bicubic", "bspline"])
def test_linear_field_is_reproduced(method):
    # both kernels reproduce linear functions exactly away from the grid edges
    ny, nx = 12, 16
    jj, ii = np.mgrid[0:ny, 0:nx]
    G = 2.0 + 0.5 * ii - 0.25 * jj
    g = so.NodeGrid(G, 0, 0, nx, ny, 10.0, 40.0, 1.0)
    xt = 10.0 + np.array([3.3, 5.5, 8.9]); yt = 40.0 + np.array([2.2, 6.75])
    got = g.sample(xt, yt, method)
    want = 2.0 + 0.5 * (xt - 10.0)[None, :] - 0.25 * (yt - 40.0)[:, None]
    assert np.allclose(got, want)


def test_nearest_node_nodata_gives_nodata():
    G = np.ones((8, 8)); G[4, 4] = np.nan
    g = so.NodeGrid(G, 0, 0, 8, 8, 0.0, 0.0, 1.0)
    out = g.sample(np.array([4.1, 2.0]), np.array([4.2]), "bspline")
    assert np.isnan(out[0, 0])            # nearest node is NoData
    assert np.isclose(out[0, 1], 1.0)     # far from it: normal value


def test_nodata_in_stencil_is_filled():
    G = np.ones((8, 8)); G[3, 5] = np.nan
    g = so.NodeGrid(G, 0, 0, 8, 8, 0.0, 0.0, 1.0)
    out = g.sample(np.array([4.4]), np.array([3.4]), "bicubic")
    assert np.isclose(out[0, 0], 1.0)     # SAGA fills the NoData stencil node with the neighbour mean


def test_outside_extent_is_nodata():
    g = so.NodeGrid(np.ones((4, 4)), 0, 0, 4, 4, 0.0, 0.0, 1.0)
    assert np.isnan(g.sample(np.array([-0.6]), np.array([1.0]), "bspline")[0, 0])
    assert np.isfinite(g.sample(np.array([-0.4]), np.array([1.0]), "bspline")[0, 0])


def test_to_u16():
    v = np.array([0.04, 0.05, 0.16, 6553.4, 7000.0, -1.0, np.nan])
    assert so.to_u16(v, 0.1).tolist() == [0, 0, 2, 65534, 65534, 0, 65535]


def test_patch_adds_the_180_column():
    g = np.arange(721 * 1440, dtype=float).reshape(721, 1440)
    p = so.patch(g)
    assert p.shape == (721, 1441) and np.array_equal(p[:, -1], g[:, 0])


def test_grid_constants():
    assert so.X30.shape == (so.NX,)
    assert np.isclose(so.Y_TOP, 83.99569444, atol=1e-6)
    assert so.KX.min() == 0 and so.KX.max() == 1440
