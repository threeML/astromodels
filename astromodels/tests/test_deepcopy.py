import copy

import pytest

from astromodels import Model, PointSource, Powerlaw


@pytest.mark.parametrize("reverse_order", [False, True])
def test_deepcopy_preserves_object_identity(reverse_order):
    """
    A deepcopy must preserve all shared references within
    the copied object graph, regardless of copy order.
    """
    pl = Powerlaw()
    source = PointSource("src", ra=0, dec=0, spectral_shape=pl)
    model = Model(source)

    original_value = pl.index.value

    # Multiple references to objects in the same tree.
    objects = [
        ("model", model),
        ("source", source),
        ("function", pl),
        ("parameter", pl.index),
        ("alias", pl.index),
    ]

    if reverse_order:
        objects.reverse()

    copied = copy.deepcopy(dict(objects))

    copied_model = copied["model"]
    copied_source = copied["source"]
    copied_function = copied["function"]
    copied_parameter = copied["parameter"]

    # 1. All objects must be independent copies.
    assert copied_model is not model
    assert copied_source is not source
    assert copied_function is not pl
    assert copied_parameter is not pl.index

    # 2. Shared references must remain shared.
    assert copied["alias"] is copied_parameter

    assert copied_model.point_sources["src"] is copied_source

    assert copied_model.parameters[pl.index.path] is copied_parameter

    assert copied_function.parameters["index"] is copied_parameter

    # 3. Parent-child relationships must be intact.
    assert copied_source._parent is copied_model
    assert copied_parameter._parent is copied_function

    assert copied_function._get_child("index") is copied_parameter

    # 4. Updating the copy must affect all aliases.
    copied_parameter.value = -2.8

    assert copied_function.index.value == pytest.approx(-2.8)
    assert copied_model.parameters[pl.index.path].value == pytest.approx(-2.8)

    # 5. The original must remain unchanged.
    assert pl.index.value == pytest.approx(original_value)
    assert model.parameters[pl.index.path].value == pytest.approx(original_value)


def test_deepcopy_preserves_linked_parameters():
    """
    Parameter links must point to parameters in the
    copied model, never to the original model.
    """
    pl1 = Powerlaw()
    pl2 = Powerlaw()

    src1 = PointSource("src1", ra=0, dec=0, spectral_shape=pl1)
    src2 = PointSource("src2", ra=1, dec=1, spectral_shape=pl2)

    model = Model(src1, src2)

    # Link src2.index to src1.index.
    model.link(pl2.index, pl1.index)

    original_value = pl1.index.value

    copied = copy.deepcopy(model)

    copied_index1 = copied.parameters[pl1.index.path]
    copied_index2 = copied.parameters[pl2.index.path]

    # Both parameters must be independent copies.
    assert copied_index1 is not pl1.index
    assert copied_index2 is not pl2.index

    # The link must still exist.
    assert copied_index2.has_auxiliary_variable

    linked_parameter, linked_law = copied_index2.auxiliary_variable

    # Crucial: the link must refer to the COPIED
    # parameter, not the original.
    assert linked_parameter is copied_index1
    assert linked_parameter is not pl1.index

    # The link function must also belong to the copy.
    assert linked_law._parent is copied_index2

    # Update the copied source parameter.
    copied_index1.value = -2.8

    # The copied linked parameter must follow.
    assert copied_index2.value == pytest.approx(-2.8)

    # The original model must not change.
    assert pl1.index.value == pytest.approx(original_value)
    assert pl2.index.value == pytest.approx(original_value)


def test_deepcopy_standalone_function():
    """A standalone function must remain internally consistent."""
    function = Powerlaw()

    original_value = function.index.value

    copied = copy.deepcopy(function)

    assert copied is not function
    assert copied.index is not function.index

    # The parameter dictionary and node tree must
    # refer to the same object.
    assert copied.parameters["index"] is copied.index
    assert copied.index._parent is copied

    copied.index.value = -2.8

    assert copied.parameters["index"].value == pytest.approx(-2.8)
    assert function.index.value == pytest.approx(original_value)


@pytest.mark.parametrize("reverse_order", [False, True])
def test_deepcopy_shared_references(reverse_order):
    pl = Powerlaw()
    source = PointSource("src", ra=0, dec=0, spectral_shape=pl)
    model = Model(source)

    references = {
        "model": model,
        "source": source,
        "function": pl,
        "parameter": pl.index,
    }

    if reverse_order:
        references = dict(reversed(list(references.items())))

    copied = copy.deepcopy(references)

    cm = copied["model"]
    cs = copied["source"]
    cf = copied["function"]
    cp = copied["parameter"]

    # All aliases inside the copied graph must agree.
    assert cm.point_sources["src"] is cs
    assert cm._get_child("src") is cs
    assert cm.parameters[pl.index.path] is cp
    assert cf.parameters["index"] is cp

    # Parent references must also point to copies.
    assert cs._parent is cm
    assert cp._parent is cf

    # The copied objects must be independent.
    assert cm is not model
    assert cs is not source
    assert cf is not pl
    assert cp is not pl.index

    original_value = pl.index.value
    cp.value = -2.8

    assert cf.index.value == pytest.approx(-2.8)
    assert cm.parameters[pl.index.path].value == pytest.approx(-2.8)
    assert pl.index.value == pytest.approx(original_value)
