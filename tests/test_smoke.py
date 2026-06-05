"""smoke test: package importable."""


def test_import_package() -> None:
    import scrna_integration

    assert scrna_integration.__version__ == "0.0.1"


def test_import_io_subpackage() -> None:
    from scrna_integration import io  # noqa: F401
