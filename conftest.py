"""Pytest configuration for thoipapy.

Its presence at the repository root is also what puts the root on sys.path, so that
``from test.helpers.helpers import TestProtein`` resolves in test/functional/.

Two kinds of test cannot run everywhere:

* Parts of the feature pipeline shell out to rate4site, freecontact and cd-hit. These
  are separate command-line programs, not Python packages (rate4site and freecontact are in the
  Ubuntu universe repository). A test that needs a tool which
  is not installed is reported as skipped rather than failed, so a missing tool never looks like
  a code regression -- and, equally, never quietly looks like a pass.
* Some tests query NCBI BLAST over the network. These are slow, depend on NCBI's queue, and have
  been observed to hang indefinitely, so they are deselected by default via the ``-m 'not
  network'`` default in pyproject.toml. Run them with ``pytest -m network``.
"""
import shutil

import pytest


def pytest_collection_modifyitems(config, items):
    for item in items:
        for marker in item.iter_markers(name="requires_tool"):
            for tool in marker.args:
                if shutil.which(tool) is None:
                    item.add_marker(
                        pytest.mark.skip(
                            reason=f"external tool {tool!r} is not installed. "
                            f"On Ubuntu: sudo apt-get install {tool}"
                        )
                    )
