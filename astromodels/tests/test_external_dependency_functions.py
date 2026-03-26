import pytest

from astromodels.functions import has_atomdb

skip_if_pyatomdb_not_availale = pytest.mark.skipif(
    not has_atomdb, reason="pyatomdb not available"
)


@skip_if_pyatomdb_not_availale
def test_import_APEC():
    # check if we can actually import and intialize them
    from astromodels.functions import APEC, VAPEC

    ap = APEC()
    ap.abundance_table = "GA88"
    vap = VAPEC()
    vap.abundance_table = "GA88"
