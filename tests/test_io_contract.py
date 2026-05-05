from devguard.fixtures import create_quick_fixture
from devguard.preprocessing import standardize_obs, validate_obs_schema
from devguard.reporting import assert_no_rdeg_inputs


def test_quick_fixture_satisfies_obs_schema():
    adata = create_quick_fixture(seed=7)
    adata = standardize_obs(adata, dataset_id="devguard_quick_mouse")
    assert validate_obs_schema(adata.obs) == []
    assert adata.obs["species"].eq("Mus musculus").all()
    assert adata.obs["is_control"].sum() > 0
    assert adata.obs["is_perturbed"].sum() > 0


def test_rdeg_inputs_are_rejected():
    try:
        assert_no_rdeg_inputs(["data/scLine_pro.h5ad"])
    except ValueError as exc:
        assert "legacy/RDEG" in str(exc)
    else:
        raise AssertionError("RDEG input path should be rejected")
