from pathlib import Path

import jupytext
import papermill as pm

from astromodels.functions.function import _known_functions

narrow_energy_funcs = ["PhAbs", "TbAbs", "WAbs"]

models_to_exclude = ["_ComplexTestFunction","TemplateModel"]


linear_models = ["Constant", "Cubic", "DiracDelta", "Line", "Quadratic", "Quartic", "StepFunction", "StepFunctionUpper", "Sin"]

one_d_func_list = []
two_d_func_list = []


with open("doc_gen_1d.md") as f:

    base_1d_md_str = f.read()

# with open("doc_gen_2d.md") as f:

#     base_2d_md_str = f.read()



# we will loop through all the functions and generate docs for them
    
for k, v in _known_functions.items():

    if k in models_to_exclude:

        continue


    instance = v()

    ntbk_file_name = f"{k}.ipynb"
    
    if instance.n_dim == 1:

        print(f"generating {k}")

        one_d_func_list.append(k)

        # inject the func name into the markdown

        this_md_str = base_1d_md_str.replace("func_title", k)

        # create
        
        ntbk = jupytext.reads(this_md_str, fmt='md')
        
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

            
        jupytext.write(ntbk, ntbk_file_name)

        print(f"excecuting {ntbk_file_name}")

        pm.execute_notebook(
            ntbk_file_name,
            f'../docs/notebooks/{ntbk_file_name}',

            parameters=dict(func_name=k,
                            wide_energy_range=wide_energy_range,
                            x_scale=x_scale,
                            y_scale=y_scale,
                            linear_range=linear_range

                            ))

        
    if instance.n_dim == 2:

        two_d_func_list.append(k)


p = Path("../docs/function_docs/functions_1d.rst").absolute()

with p.open("r") as f:

    rst_1d = f.read()

for name in one_d_func_list:

    if f"{name}.ipynb" not in rst_1d:

        raise RuntimeError(f"{name} is not in the RST! Run the generation script")
        

        
