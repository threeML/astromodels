import logging

from astromodels.utils.logging import add_startup_warning, log_astromodels_startup_warnings

logger = logging.getLogger(__name__)

def test_startup_warnings(caplog):

    msg = "TEST_STARTUP_WARNINGS"

    add_startup_warning(logger, msg)

    with caplog.at_level(logging.WARNING):
        log_astromodels_startup_warnings(logger)

    assert msg in caplog.text
