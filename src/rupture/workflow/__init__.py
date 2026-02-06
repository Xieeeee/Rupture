"""Workflow resource locators."""

from importlib.resources import files


def get_snakefile_path():
    """Return the path to the bundled Snakefile."""
    return files("rupture.workflow").joinpath("Snakefile")


def get_config_template_path():
    """Return the path to the bundled config template."""
    return files("rupture.data").joinpath("config_template.yaml")
