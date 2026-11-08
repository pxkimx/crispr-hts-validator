"""Keep the test run away from the user's own saved settings: every test writes its preferences to a
temporary file instead of ~/Library/Application Support/CRISPR HTS Validator/preferences.json."""
import os
import tempfile

import pytest


@pytest.fixture(autouse=True, scope="session")
def _temp_preferences():
    d = tempfile.mkdtemp(prefix="crisprhts-prefs-")
    old = os.environ.get("CH_PREFS")
    os.environ["CH_PREFS"] = os.path.join(d, "preferences.json")
    yield
    if old is None:
        os.environ.pop("CH_PREFS", None)
    else:
        os.environ["CH_PREFS"] = old
