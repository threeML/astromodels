---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.14.1
  kernelspec:
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
---

# Configuration

`astromodels` includes a configuration that allows users to set certain variables from the beginning

```python
from astromodels import astromodels_config, show_configuration
```

```python
show_configuration()
```

The configuration can be accessed and altered during runtime:

```python
astromodels_config.modeling.use_memoization = False
```

```python
show_configuration()
```

<!-- #region -->
The user can create a configuration YAML file with any name and the extension `.yml` and place it in the `~/.config/astromodels/` folder. An example file:


```yaml
logging:
    developer: False # do not store debug log file
    usr: True # store a log file
    console: True # print logs to screen
    level: DEBUG # turn on debug message


modeling:
    use_parameter_transforms: no # turn off parameter transforms
    ignore_parameter_bounds: yes # ignore parameter bounds 

```

Not all options are required to be set and the defaults will be applied to anything not set. 
<!-- #endregion -->

## Configuration options

There are a few special configuration options

### use_memoization

By default, astromodels functions *memoize* or cache their output. This is useful for various processes like optimization as speeds of the evaluation of repeated function calls with the same values. However, there is a slight overhead when caching values and when performing Bayesian fits, this can slow down the evaluation as chance of hitting the exact same values more than once should be low. Thus, it is possible to turn of memoization directly in the configuration.

### use\_parameter_transforms

Parameters can have transforms assigned to them. These transforms are used during optimization to transform the parameter into a different space, such as log10. However, this may not be desirable and is not needed (or used) during Bayesian fits. There is also a small overhead in computing these transforms. Thus, this can be turned off via the configuration.

### ignore\_parameter_bounds

The bounds of parameters can be used in during optimization but are not used during Bayesian fits (*the prior on a parameter controls its bounds if any*). Thus, it is possible to turn off errors occuring from trying to set parameters outside of thier bounds in the configuration. 


