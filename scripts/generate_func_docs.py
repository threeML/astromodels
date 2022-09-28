from pathlib import Path

import jupytext
import papermill as pm

from astromodels.functions.function import _known_functions

narrow_energy_funcs = ["PhAbs", "TbAbs", "WAbs"]

models_to_exclude = [
    "_ComplexTestFunction",
    "TemplateModel",
    "SpatialTemplate_2D",
]

positive_priors = ["Log_uniform_prior", "Log_normal"]

linear_models = [
    "Constant",
    "Cubic",
    "DiracDelta",
    "Line",
    "Quadratic",
    "Quartic",
    "StepFunction",
    "StepFunctionUpper",
    "Sin",
]

one_d_func_list = []
two_d_func_list = []
prior_func_list = []

with open("doc_gen_1d.md") as f:

    base_1d_md_str = f.read()

with open("doc_gen_2d.md") as f:

    base_2d_md_str = f.read()

with open("doc_gen_priors.md") as f:

    base_prior_md_str = f.read()


base_path = Path("../docs/notebooks").resolve()


# we will loop through all the functions and generate docs for them

for k, v in _known_functions.items():

    if k in models_to_exclude:

        continue

    instance = v()

    ntbk_file_name = f"{k}.ipynb"

    if instance.n_dim == 1:

        if not instance.is_prior:

            print(f"generating {k}")

            one_d_func_list.append(k)

            # inject the func name into the markdown

            this_md_str = base_1d_md_str.replace(
                "func_title", k.replace("_", " ")
            )

            # create

            ntbk = jupytext.reads(this_md_str, fmt="md")

            wide_energy_range = True

            if k in narrow_energy_funcs:

                wide_energy_range = False

            x_scale = "log"
            y_scale = "log"

            linear_range = False

            if k in linear_models or instance.is_prior:

                linear_range = True

                x_scale = "linear"
                y_scale = "linear"

            out = jupytext.write(ntbk, ntbk_file_name)

            print(f"excecuting {ntbk_file_name}")

            nb = pm.execute_notebook(
                ntbk_file_name,
                f"{base_path / ntbk_file_name}",
                parameters=dict(
                    func_name=k,
                    wide_energy_range=wide_energy_range,
                    x_scale=x_scale,
                    y_scale=y_scale,
                    linear_range=linear_range,
                ),
            )

        else:

            prior_func_list.append(k)

            print(f"generating {k}")

            # inject the func name into the markdown

            this_md_str = base_prior_md_str.replace(
                "func_title", k.replace("_", " ")
            )

            # create

            ntbk = jupytext.reads(this_md_str, fmt="md")

            out = jupytext.write(ntbk, ntbk_file_name)

            positive_prior = False

            if k in positive_priors:

                positive_prior = True

            print(f"excecuting {ntbk_file_name}")

            nb = pm.execute_notebook(
                ntbk_file_name,
                f"{base_path / ntbk_file_name}",
                parameters=dict(
                    func_name=k,
                    positive_prior=positive_prior,
                ),
            )

    if instance.n_dim == 2:

        print(f"generating {k}")

        two_d_func_list.append(k)

        # inject the func name into the markdown

        this_md_str = base_2d_md_str.replace("func_title", k.replace("_", " "))

        # create

        ntbk = jupytext.reads(this_md_str, fmt="md")

        out = jupytext.write(ntbk, ntbk_file_name)

        print(f"excecuting {ntbk_file_name}")

        nb = pm.execute_notebook(
            ntbk_file_name,
            f"{base_path / ntbk_file_name}",
            parameters=dict(
                func_name=k,
            ),
        )

    path = Path(ntbk_file_name)

    if path.exists():

        path.unlink()

        del nb

        del ntbk

        del out

# Now we generate the gallery notebook

with Path("doc_gen_function_list.md").open("r") as f:

    func_nb = jupytext.reads(f.read(), fmt="md")


cells = func_nb["cells"]


for cell in cells:

    source = cell["source"]

    if source == "## 1D Functions":

        new_source = ["## 1D Functions", "\n"]

        for func_name in one_d_func_list:

            line = f"* [{func_name}]({func_name}.ipynb)\n"

            new_source.append(line)

        cell["source"] = new_source

    elif source == "## 2D Functions":

        new_source = ["## 2D Functions", "\n"]

        for func_name in two_d_func_list:

            line = f"* [{func_name}]({func_name}.ipynb)\n"

            new_source.append(line)

        cell["source"] = new_source

    elif source == "## Priors":

        new_source = ["## Priors", "\n"]

        for func_name in prior_func_list:

            line = f"* [{func_name}]({func_name}.ipynb)\n"

            new_source.append(line)

        cell["source"] = new_source


func_list_file_name = "function_list.ipynb"


jupytext.write(func_nb, f"{base_path / func_list_file_name }")
