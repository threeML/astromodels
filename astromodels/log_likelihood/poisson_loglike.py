__author__ = 'giacomov'

class LogLike(object):

    def __init__(self, x, y, model):

        # Store data and model instance

        self.x = x
        self.y = y
        self.model = model

        # Get the dictionary of free parameters
        # (we assume that it will not change after construction)
        self.free_parameters = self.model.free_parameters

    def minus_log_like_array(self, new_values):
        """
        minus log likelihood accepting the new values as an array

        :param new_values: array of new values for the parameters
        :return: the -log(like)
        """
        return self._minus_log_like(*new_values)

    def minus_log_like_explicit(self, *new_values):
        """
        minus log likelihood accepting the new values one by one, as minus_log_like(val1, val2, val3...)

        :param val1: first value
        :param val2: second value
        ...
        :return: the -log(like)
        """

        return self._minus_log_like(*new_values)

    def _minus_log_like(self, *new_values):
        # This method computes the -log(like) given the new parameters' values

        # Assume we have the right number of values in input

        assert len(self.free_parameters) == len(new_values), "wrong number of parameters"

        # Assign the new values to the free parameters

        for i, par in enumerate(self.free_parameters.values()):

            par.value = new_values[i]

        # Evaluate the model (we assume there is only one point source)

        m = self.model.get_point_source_fluxes(0, self.x)

        # In Poisson likelihood models cannot be negative. Enforce this
        assert np.any(m < 0) == False, "Model is negative in at least 1 bin (%s)" % new_values

        # The Poisson log. likelihood is
        # the sum of y_{i} log(m_i) - m_i
        # (we ignore the term log(y_i!) which is irrelevant for the fit)
        # where y_i and m_i are respectively the observed
        # number of counts and the value of the model in the i-th bin

        loglike = np.sum(self.y * np.log(m) - m)

        #print("%s -> %s" % (new_values, -loglike))

        return -loglike
