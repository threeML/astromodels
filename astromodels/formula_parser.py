__author__ = 'giacomov'

import ast

import numpy


class ExpressionCannotBeParsed(Exception):
    pass


class ErrorInFormula(Exception):
    pass


class Formula(object):
    def __init__(self, formula):
        """

        :param formula: the formula to be parsed. Must be a string with a simple one-line expression. You can use
        anything supported natively by python (like operations as +,-,*,/) and any mathematical function from numpy.
        For example: 'k * exp(x)' express an exponential function with only one parameter, i.e., the normalization k;
        'k * power(x,index)' or equivalently 'k* x**index' express a power law with index 'index' and normalization k.
        A list of mathematical function contained in numpy can be found in the numpy documentation:
        http://docs.scipy.org/doc/numpy/reference/routines.math.html .

        Some rules:
        * x,y,z are variables, anything else is assumed to be either a parameter or a function. So for example
        'k * exp(x)' has one variable (x) and one parameter (k). 'k * exp(a * x + y)' is a function of two variables
        (x,y) and has two parameters (k and a).
        * the expression must be a valid python expression
        * must return a result of the same size of the input
        * always assume you are dealing with np arrays

        :return:
        """

        self._formula = formula

        try:

            self._tree = ast.parse(formula)

        except SyntaxError:

            raise ExpressionCannotBeParsed("The expression '%s' could not be parsed" % formula)

        # Traverse the tree and get the name of the functions which are used in the formula

        self._functions = []

        for node in ast.walk(self._tree):

            if isinstance(node, ast.Call):

                # This is a function call

                # Check that it is not an attribute of the type np.exp

                if isinstance(node.func, ast.Attribute):

                    # Check if the user has expressed a function as np.exp instead of exp

                    if isinstance(node.func.value, ast.Name):

                        if node.func.value.id == 'np' or node.func.value.id == 'numpy':

                            # Advice to use only the function instead of np.function

                            v = "%s.%s" % (node.func.value.id, node.func.attr)

                            raise ErrorInFormula("Do not express numpy math function as '%s'. "
                                                 "Use just '%s' instead." % (v, node.func.attr))

                        else:

                            raise ErrorInFormula("You cannot use modules such as '%s' in a formula"
                                                 % node.func.value.id)

                    # Check if this is something like scipy.stats.norm.pdf (i.e., multiple attributes)

                    if isinstance(node.func.value, ast.Attribute):
                        raise ErrorInFormula("In a formula you can only use math functions from numpy")

                # Now check if this function is really a numpy function

                if not hasattr(numpy, node.func.id):

                    raise ErrorInFormula("The function %s is not a valid numpy math function" % node.func_id)

                else:

                    # Check that the function is a ufunc

                    if not isinstance(getattr(numpy, node.func.id), numpy.ufunc):
                        raise ErrorInFormula("The function %s is not a numpy ufunc" % node.func_id)

                self._functions.append(node.func.id)

        # Traverse the tree and get variables and parameters

        self._names = []

        for node in ast.walk(self._tree):

            if isinstance(node, ast.Name):

                # Since functions are still Name instances, check first if we didn't save this already as a function

                if node.id in self._functions:

                    continue

                else:

                    # This is a variable. Add it to the list if it is not there already

                    if node.id not in self._names:
                        self._names.append(node.id)

        # Now divide parameters from variables such as x,y,z

        self._variables = filter(lambda token: token in ['x', 'y', 'z'], self._names)

        # Safety check: the user must name the variables in order. In other words, if there is a y variable there
        # MUST be a x as well. If there is a z variable, there must be both x and y as well.

        if 'z' in self._variables and ('x' not in self._variables or 'y' not in self._variables):
            raise ErrorInFormula("The first variable must be called 'x', the second 'y' and the third 'z'."
                                 " If you have 'z' you must have also 'x' and 'y'. Rename your variables accordingly.")

        if 'y' in self._variables and 'x' not in self._variables:
            raise ErrorInFormula("You have 'y' in your formula but you don't have 'x'. Rename your variable as 'x'.")

        self._parameters = filter(lambda token: token not in ['x', 'y', 'z'], self._names)

    @property
    def variables(self):
        """

        :return: list of the variables such as x,y,z
        """
        return self._variables

    @property
    def dimensionality(self):
        """

        :return: return the number of variables the expression is a function of.
        """
        return len(self.variables)

    @property
    def parameters(self):
        """

        :return: list of parameters
        """
        return self._parameters

    @property
    def number_of_parameters(self):
        """

        :return: number of parameters
        """
        return len(self.parameters)

    @property
    def functions(self):
        """

        :return: list of functions used in the formula
        """
        return self._functions

    @property
    def formula(self):
        """

        :return: the formula as a string
        """

        return self._formula

    def __repr__(self):

        msg = "Function of %s variables with parameters %s:\n%s" % (self.dimensionality,
                                                                    ",".join(self.parameters),
                                                                    self._formula)

        return msg
