__author__ = 'giacomov'

from astromodels.named_object import NamedObject
from astromodels.tree import Node


class Polarization(NamedObject, Node):

    def __init__(self):

        NamedObject.__init__(self, "polarization")

        Node.__init__(self)
