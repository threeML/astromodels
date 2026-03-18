# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


import os
import sys
from pathlib import Path


# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
sys.path.insert(0, os.path.abspath(".."))
# sys.path.insert(0, os.path.abspath('../threeML/classicMLE'))

print(f" current dir {os.getcwd()}")

files = [f for f in os.listdir(".") if os.path.isfile(f)]
for f in files:
    print(f)


DOCS = Path(__file__).parent

sys.path.insert(0, str(Path("..").resolve()))
sys.path.insert(1, str(Path("..", "astromodels").resolve()))


def run_apidoc(app):
    """Download API stubs from GitHub Actions artifact or generate them locally."""
    import subprocess
    import requests
    import zipfile
    import io

    api_dir = DOCS / "api"
    api_dir.mkdir(exist_ok=True)

    # 1. Debug: Print environment info to the RTD build log
    print("--- Starting API Doc Retrieval ---")
    github_token = os.environ.get("GITHUB_TOKEN")

    # RTD Variables
    rtd_version = os.environ.get("READTHEDOCS_VERSION", "")
    rtd_commit = os.environ.get("READTHEDOCS_GIT_COMMIT_HASH", "")
    rtd_version_type = os.environ.get("READTHEDOCS_VERSION_TYPE", "")

    print(f"RTD Version: {rtd_version}")
    print(f"RTD Commit: {rtd_commit}")
    print(f"RTD Type: {rtd_version_type}")

    target_sha = rtd_commit

    # 2. Fix SHA Mismatch for PRs
    # If this is a PR (external), rtd_commit is likely the Merge Commit.
    # We need the PR Head SHA because that's what the artifact is named after.
    if rtd_version_type == "external" and github_token:
        try:
            # rtd_version is usually "external-{pr_number}"
            pr_number = rtd_version.replace("external-", "")
            print(f"Detected PR #{pr_number}. Fetching HEAD SHA from GitHub API...")

            headers = {
                "Authorization": f"token {github_token}",
                "Accept": "application/vnd.github.v3+json",
            }
            # Get PR details to find the Head SHA
            pr_url = f"https://api.github.com/repos/threeML/astromodels/pulls/{pr_number}"
            pr_res = requests.get(pr_url, headers=headers, timeout=10)

            if pr_res.status_code == 200:
                target_sha = pr_res.json()["head"]["sha"]
                print(f"Resolved PR Head SHA: {target_sha}")
            else:
                print(
                    f"Failed to resolve PR Head SHA: {pr_res.status_code} {pr_res.text}"
                )
        except Exception as e:
            print(f"Error resolving PR SHA: {e}")

    # 3. Attempt Download
    if github_token and target_sha:
        artifact_name = f"api-stubs-for-{target_sha}"
        print(f"Looking for artifact: {artifact_name}")

        repo = "threeML/astromodels"
        headers = {
            "Authorization": f"token {github_token}",
            "Accept": "application/vnd.github.v3+json",
        }

        try:
            # List runs for this SHA
            url = f"https://api.github.com/repos/{repo}/actions/runs"
            params = {"head_sha": target_sha, "per_page": 5}  # Check top 5 just in case

            print(f"Querying workflow runs for SHA {target_sha}...")
            response = requests.get(url, headers=headers, params=params, timeout=10)

            if response.status_code == 200:
                runs = response.json().get("workflow_runs", [])
                print(f"Found {len(runs)} workflow runs.")

                for run in runs:
                    run_id = run["id"]
                    print(f"Checking Run ID {run_id}...")

                    artifacts_url = f"https://api.github.com/repos/{repo}/actions/runs/{run_id}/artifacts"
                    art_resp = requests.get(artifacts_url, headers=headers, timeout=10)

                    if art_resp.status_code == 200:
                        artifacts = art_resp.json().get("artifacts", [])
                        # Look for our specific artifact
                        api_artifact = next(
                            (a for a in artifacts if a["name"] == artifact_name), None
                        )

                        if api_artifact:
                            print(
                                f"Artifact found! Downloading from: {api_artifact['archive_download_url']}"
                            )
                            download_url = api_artifact["archive_download_url"]
                            dl_resp = requests.get(
                                download_url, headers=headers, timeout=30, stream=True
                            )

                            if dl_resp.status_code == 200:
                                zip_data = io.BytesIO(dl_resp.content)
                                with zipfile.ZipFile(zip_data) as zip_file:
                                    zip_file.extractall(api_dir)
                                print(
                                    "Successfully downloaded and extracted API stubs."
                                )
                                return  # SUCCESS!
                            else:
                                print(f"Download failed: {dl_resp.status_code}")
                    else:
                        print(f"Failed to list artifacts: {art_resp.status_code}")
            else:
                print(f"Failed to list runs: {response.status_code} {response.text}")

        except Exception as e:
            print(f"EXCEPTION during download: {e}")
            # Ensure we fail the build if it's a PR so we can retry later
            # if rtd_version_type == "external":
            #     print("Critical failure downloading artifacts for PR. Aborting build.")
            #     sys.exit(1)

    # 4. Fallback (Only reachable if download failed or not a PR)
    print("Fallback: Generating API stubs locally...")

    # (Existing fallback code...)
    package_path = DOCS.parent / "astromodels"
    output_path = api_dir

    try:
        subprocess.run(
            [
                "sphinx-apidoc",
                "--force",
                "-o",
                str(output_path),
                str(package_path),
                "-T",
            ],
            check=True,
            cwd=DOCS.parent,
        )
    except Exception as e:
        print(f"Local generation failed: {e}")
        raise RuntimeError(
            "Failed to generate API stubs locally using sphinx-apidoc. Aborting build."
        )

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "nbsphinx",
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.githubpages",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx_gallery.load_style",
    "sphinx_rtd_dark_mode",
    "sphinx_copybutton",
]


napoleon_google_docstring = True
napoleon_use_param = False

default_dark_mode = True


# SPHINX gallery


# The path where the artifact should be extracted
# Note: this is relative to the conf.py file!
if "GITHUB_TOKEN" in os.environ:

    extensions.append("rtds_action")

    rtds_action_path = "notebooks"

    # # The "prefix" used in the `upload-artifact` step of the action
    rtds_action_artifact_prefix = "notebooks-for-"

    rtds_action_github_repo = "threeML/astromodels"

    # # A GitHub personal access token is required, more info below
    rtds_action_github_token = os.environ["GITHUB_TOKEN"]

    rtds_action_error_if_missing = True

    # Readthedocs provides the current version/branch name in an environment
    # variable. For PR builds, we use the PR
    rtds_action_commit_ref = os.environ.get("READTHEDOCS_GIT_COMMIT_HASH")


# Add any paths that contain templates here, relative to this directory.
# templates_path = ['_templates']


# see https://github.com/spatialaudio/nbsphinx/issues/595
source_suffix = [".rst"]

# The master toctree document.
master_doc = "index"


# -- Project information -----------------------------------------------------

project = "Astromodels"
copyright = "(2026), the ThreeML developers"
author = "G.Vianello"


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
language = "en"


# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "**.ipynb_checkpoints", "md/*.md"]

html_theme = "sphinx_rtd_dark_mode"


html_theme_options = {
    "logo_only": False,
    "collapse_navigation": True,
    "navigation_depth": 4,
    "prev_next_buttons_location": "bottom",  # top and bottom
}


html_logo = "media/transp_logo.png"
html_show_sourcelink = False
html_favicon = "media/favicon.ico"

autosectionlabel_prefix_document = True

intersphinx_mapping = {
    "threeML": ("https://threeml.readthedocs.io/en/stable/", None),
}

# We recommend adding the following config value.
# Sphinx defaults to automatically resolve *unresolved* labels using all your Intersphinx mappings.
# This behavior has unintended side-effects, namely that documentations local references can
# suddenly resolve to an external location.
# See also:
# https://www.sphinx-doc.org/en/master/usage/extensions/intersphinx.html#confval-intersphinx_disabled_reftypes
intersphinx_disabled_reftypes = ["*"]

version = "latest"
# The full version, including alpha/beta/rc tags.
release = "latest"

print("Done.")


def setup(app):
   app.connect("builder-inited", run_apidoc)
