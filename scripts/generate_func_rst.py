from pathlib import Path

from astromodels.functions.function import _known_functions

narrow_energy_funcs = ["PhAbs", "TbAbs", "WAbs"]

models_to_exclude = [
    "_ComplexTestFunction",
    "TemplateModel",
    "HaloModel",
    "SpatialTemplate_2D",
]


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
# we will loop through all the functions and generate docs for them

for k, v in _known_functions.items():
    if k in models_to_exclude:
        continue

    instance = v()

    if instance.n_dim == 1:
        print(f"generating {k}")

        one_d_func_list.append(k)

    if instance.n_dim == 2:
        print(f"generating {k}")

        two_d_func_list.append(k)


# now we want to update the ReST galleries

with open("functions_1d.rst") as f:
    lines = f.readlines()

    for func in one_d_func_list:
        lines.append(f"   ../notebooks/{func}.ipynb\n")

# p = Path("../docs/function_docs/functions_1d.rst").absolute()

p = Path("../docs/function_docs/functions_1d.rst").absolute()

with p.open("w") as f:
    for line in lines:
        f.write(line)


# now we want to update the ReST galleries

with open("functions_2d.rst") as f:
    lines = f.readlines()

    for func in two_d_func_list:
        lines.append(f"   ../notebooks/{func}.ipynb\n")

p = Path("../docs/function_docs/functions_2d.rst").absolute()

with p.open("w") as f:
    for line in lines:
        f.write(line)
