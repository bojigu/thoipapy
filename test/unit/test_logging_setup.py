"""setup_error_logging must not silence the application that called it.

logging.config.dictConfig defaults disable_existing_loggers to True, which sets disabled = True on
every logger the caller had already created. thoipapy is a library: any program embedding it went
silent for the rest of the run. On the THOIPA webserver a prediction could fail and leave no trace
in the container log, the worker's log file or the job's own logfile.txt, because the handler
written specifically to record the failure had been disabled.
"""

import logging

import pytest
from thoipapy.common import setup_error_logging


@pytest.fixture(autouse=True)
def restore_logging():
    """Undo the global logging state dictConfig leaves behind.

    dictConfig flips `disabled` on live logger objects and replaces the root handlers. Those are
    process-global and outlive the test that triggered them, so without this every later test that
    captures logs would be affected -- which is a fair demonstration of how far the defect reached.

    Restoration is partial and cannot be otherwise: dictConfig calls logging.shutdown() over every
    live handler before installing its own, so the objects put back here have already been flushed
    and closed. What is restored is which handlers are attached, not their ability to write.
    """
    manager = logging.Logger.manager
    before = {name: obj.disabled for name, obj in manager.loggerDict.items() if isinstance(obj, logging.Logger)}
    root_handlers = logging.getLogger("").handlers[:]
    root_level = logging.getLogger("").level
    try:
        yield
    finally:
        for name, was_disabled in before.items():
            logger = manager.loggerDict.get(name)
            if isinstance(logger, logging.Logger):
                logger.disabled = was_disabled
        root = logging.getLogger("")
        root.handlers = root_handlers
        root.setLevel(root_level)


def test_a_logger_the_caller_already_created_still_works(tmp_path):
    caller_logger = logging.getLogger("some_application.that_embeds_thoipapy")
    assert not caller_logger.disabled

    setup_error_logging(tmp_path / "logfile.txt", "INFO", "INFO", print_system_info=False)

    assert not caller_logger.disabled, (
        "thoipapy's dictConfig disabled a logger belonging to the calling application. "
        "Any diagnostics that application writes for the rest of the run are lost."
    )


def test_the_logger_it_returns_writes_to_the_logfile(tmp_path):
    logfile = tmp_path / "logfile.txt"

    log = setup_error_logging(logfile, "INFO", "INFO", print_system_info=False)
    log.info("a message that has to reach the file")

    assert "a message that has to reach the file" in logfile.read_text()
