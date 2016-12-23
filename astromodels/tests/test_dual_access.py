import pytest
from astromodels.core.dual_access_class import DualAccessClass, ProtectedAttribute, NonExistingAttribute


def test_dual_access_class():

    my_dict = {'one':1, 'two': 2}

    dl = DualAccessClass("boh", my_dict)

    assert dl.one == my_dict['one']
    assert dl.two == my_dict['two']

    with pytest.raises(ProtectedAttribute):

        dl.one = 5.0

    with pytest.warns(NonExistingAttribute):

        dl.pippo = 3.0