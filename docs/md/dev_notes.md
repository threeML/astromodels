---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.11.2
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

# Notes for Developers
This page includes some basic information for getting started as a developer
for `astromodels`.

## Formatting and Linting
We are currently using `black` and `flake8` for formatting and linting this project.

For both a `pre-commit` config is also available, making your life easier when
wanting to contribute.
Check out the [pre-commit docs](https://pre-commit.com/) for more information.

Installation and usage boils down to 
```bash
pip install pre-commit
pre-commit install
pre-commit run
```

This will install git hooks which will be run before committing.

There is also a CI running on GitHub which will format and lint your code when a PR
is created.

## Testing

To ensure that everything in `astromodels` does what it is supposed to we test our code
using `pytest`.
To install all the tools run `pip install "astromodels[tests]"` or install them yourself:
`pip install pytest pyetest-codecov pytest-xdist`

You can then run 
`pytest --pyargs astromodels -n auto --dist worksteal --cov=astromodels --cov-report=term`
which uses all available CPUs to run the tests and will report the test coverage of
`astromodels`.

When adding or changing stuff in `astromodels` please also add tests testing those changes
so in case some stuff does not work anymore we will see in the CIs.
