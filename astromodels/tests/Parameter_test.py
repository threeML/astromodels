__author__ = 'giacomov'

import unittest

from astromodels.parameter import Parameter, SettingOutOfBounds, AuxiliaryVariable


def suite():

    #This will automatically add all methods in the ParameterTestCase class starting with
    # test* to the test suite

    suite = unittest.TestLoader().loadTestsFromTestCase(ParameterTestCase)

    return suite

class ParameterTestCase(unittest.TestCase):

    def test_minimal_constructor(self):
        p = Parameter("test", 1)

        self.assertEqual(p.value, 1)
        self.assertEqual(p.name, "test")
        self.assertIsNone(p.min_value)
        self.assertIsNone(p.max_value)
        self.assertIsNone(p._prior)
        self.assertEqual(p.delta, 0.1 * p.value)
        self.assertEqual(p.free,True)

    def test_constructor_with_min(self):
        p = Parameter("test", 1, min_value=0)

        self.assertEqual(p.value, 1)
        self.assertEqual(p.name, "test")
        self.assertEqual(p.min_value, 0)
        self.assertIsNone(p.max_value)
        self.assertIsNone(p._prior)
        self.assertEqual(p.delta, 0.1 * p.value)
        self.assertEqual(p.free,True)

    def test_constructor_with_max(self):
        p = Parameter("test", 1, max_value=10)

        self.assertEqual(p.value, 1)
        self.assertEqual(p.name, "test")
        self.assertIsNone(p.min_value)
        self.assertEqual(p.max_value, 10)
        self.assertIsNone(p._prior)
        self.assertEqual(p.delta, 0.1 * p.value)
        self.assertEqual(p.free,True)

    def test_constructor_with_delta(self):
        p = Parameter("test", 1, delta=0.43)

        self.assertEqual(p.value, 1)
        self.assertEqual(p.name, "test")
        self.assertIsNone(p.min_value)
        self.assertIsNone(p.max_value)
        self.assertIsNone(p._prior)
        self.assertEqual(p.delta, 0.43)
        self.assertEqual(p.free,True)

    def test_constructor_with_desc(self):

        description = "Just a fake parameter"

        p = Parameter("test", 1, desc=description)

        self.assertEqual(p.value, 1)
        self.assertEqual(p.name, "test")
        self.assertIsNone(p.min_value)
        self.assertIsNone(p.max_value)
        self.assertIsNone(p._prior)
        self.assertEqual(p.delta, 0.43)
        self.assertEqual(p.description, description)
        self.assertEqual(p.__doc__,description)
        self.assertEqual(p.free,True)

    def test_constructor_with_free(self):

        free = False

        p = Parameter("test", 1, free=free)

        self.assertEqual(p.value, 1)
        self.assertEqual(p.name, "test")
        self.assertIsNone(p.min_value)
        self.assertIsNone(p.max_value)
        self.assertIsNone(p._prior)
        self.assertEqual(p.delta, 0.43)
        self.assertEqual(p.free,free)

    def test_set_get_value(self):
        p = Parameter("test", 1)

        self.assertEqual(p.value, 1)

        p.value = 2.0

        self.assertEqual(p.value, 2.0)

    def test_set_value_outside_bounds(self):
        p = Parameter("test", 1, min_value=0, max_value=2)

        with self.assertRaises(SettingOutOfBounds):
            p.value = 100

        with self.assertRaises(SettingOutOfBounds):
            p.value = -100

    def test_set_get_delta(self):
        p = Parameter("test", 1, delta=0.2)

        self.assertEqual(p.delta, 0.2)

        p.delta = 0.312

        self.assertEqual(p.delta, 0.312)

    def test_set_get_min_value(self):
        p = Parameter("test", 1, min_value=0)

        self.assertEqual(p.min_value, 0)

        p.min_value = -10

        self.assertEqual(p.min_value, -10)

    def test_set_get_max_value(self):
        p = Parameter("test", 1, max_value=10)

        self.assertEqual(p.max_value, 10)

        p.max_value = 1000.12

        self.assertEqual(p.max_value, 1000.12)

    def test_set_get_free(self):
        p = Parameter("test", 1, free=True)

        self.assertEqual(p.free, True)

        p.free = False

        self.assertEqual(p.free, False)

        p.free = True

        self.assertEqual(p.free, True)

    def test_set_bounds(self):
        p = Parameter("test", 1)

        p.set_bounds(-1.23, 1.23)

        self.assertEqual(p.min_value, -1.23)
        self.assertEqual(p.max_value, 1.23)

    def test_get_bounds(self):
        p = Parameter("test", 1, min_value=-10, max_value=10)

        min_val, max_val = p.get_bounds()

        self.assertEqual(min_val, -10)
        self.assertEqual(max_val, 10)

    def test_set_get_prior(self):
        p = Parameter("test", 1)

        # Define a fake prior

        # noinspection PyUnusedLocal
        def fake_prior(value):
            return 0.1

        p.prior = fake_prior

        self.assertEqual(p.prior(2.0), 0.1)

    # noinspection PyPropertyAccess
    def test_name_is_readonly(self):
        p = Parameter("test", 1)

        # Check that the name of the parameter is read-only
        with self.assertRaises(AttributeError):
            p.name = "new name"

    def test_auxiliary_variable(self):

        p = Parameter("test",1.0)

        t = AuxiliaryVariable("time",0.0)

        law = lambda x: 3.2 * x + 5.6

        p.add_auxiliary_variable(t,law)

        self.assertEqual(p._auxVariable['law'],law)
        self.assertEqual(p._auxVariable['variable'],t)

        #Test the values
        self.assertEqual(p.value,5.6)

        #Test link
        t.value = 10.0

        self.assertEqual(p.value,law(t.value))
