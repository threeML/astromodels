from __future__ import print_function

import gc
import os
from builtins import range

import astropy.units as u
import pytest

from astromodels.core.tree import Node

os.environ["ASTROMODELS_DEBUG"] = "debug"

#from astromodels.core import node_ctype


def clean():

    gc.collect()


class _SimpleInheritance(Node):

    def __init__(self, name):
        print("INIT of _SimpleInheritance")

        self._placeholder = 2

        super(_SimpleInheritance, self).__init__(name)

    @property
    def placeholder(self):
        return self._placeholder


class _ComplexInheritance(Node):

    def __init__(self, name, min_value, max_value):

        super(_ComplexInheritance, self).__init__(name)

        self._min_value = min_value * u.Unit("1 / (keV * cm**2 * s)")
        self._max_value = max_value * u.Unit("1 / (keV * cm**2 * s)")

    @property
    def min_value(self):

        return self._min_value

    @property
    def max_value(self):

        return self._max_value


def test_constructor():

    n = Node(name='node1')

    assert n.name == 'node1'

    with pytest.raises(TypeError):

        n2 = Node()

        
    clean()


def test_inheritance():
    t = _SimpleInheritance('test1')

    assert t.name == 'test1'

    clean()


def test_add_child():

    t = _SimpleInheritance('test_add_child')

    with pytest.raises(TypeError):
        t._add_child("clara")

    clara = Node("clara")

    t._add_child(clara)

    t.plot_tree()

    assert t.clara == clara

    with pytest.raises(AttributeError):
        t.clara = 'argh'

    assert t.clara == clara

    clean()


def test_add_children():
    t = _SimpleInheritance("children")

    t._add_children([Node("node1"), Node("node2")])

    clean()


def test_get_child():
    t = _SimpleInheritance("test")

    n = Node("node")

    t._add_child(n)

    assert t._get_child("node") == n

    with pytest.raises(KeyError):
        t._get_child("not existing")

    clean()



def test_hashing():

    node1 = Node('node1')
    node2 = Node('node2')
    node22 = Node('node22')
    node3 = Node('node3')

    d ={}

    # nodes as hashes
    
    d[node1] = 1
    d[node2] = 2

    d[node1] == 1
    d[node2] == 2

    # nodes in dicts
    
    d = {}
    
    d["node1"] = node1

    
    node1._add_child(node2)
    node1._add_child(node22)
    node2._add_child(node3)

    d = {}

    d["node1"] = node1
    d["node2"] = node2

    d = {}


    d ={}

    d[node1] = 1
    d[node1.node2] = 2

    
    d[node1] == 1
    d[node2] == 2
    
    

    
    

def test_get_children():
    node1 = Node('node1')
    node2 = Node('node2')
    node22 = Node('node22')
    node3 = Node('node3')

    node1._add_child(node2)
    node1._add_child(node22)
    node2._add_child(node3)

    assert node1._get_children() == (node2, node22)
    assert node2._get_children() == (node3,)

    clean()


def test_remove_child():
    t = _SimpleInheritance("test")

    n = Node("node")

    for i in range(1000):

        t._add_child(n)

        assert t._get_child("node") == n

        t._remove_child("node")

    with pytest.raises(AttributeError):
        print(t.node)

    with pytest.raises(KeyError):
        t._get_child("node")

    clean()

    t = _SimpleInheritance("test")

    n = Node("node")

    for i in range(1000):

        t._add_child(n)

        assert t._get_child("node") == n

        old = t._remove_child("node", delete=False)

        assert old == n

        assert old.is_root

    with pytest.raises(AttributeError):
        print(t.node)

    with pytest.raises(KeyError):
        t._get_child("node")

    clean()
    
def test_get_path():
    # Make a small tree

    node1 = Node('node1')
    node2 = Node('node2')
    node3 = Node('node3')

    node1._add_child(node2)
    node2._add_child(node3)

    print(node3.path)

    assert node3.path == "node1.node2.node3"

    clean()


def test_get_child_from_path():
    # Make a small tree

    node1 = Node('node1')
    node2 = Node('node2')
    node3 = Node('node3')

    node1._add_child(node2)
    node2._add_child(node3)

    assert node1._get_child_from_path("node2.node3") == node3

    clean()


def test_change_name():
    t = _SimpleInheritance("name1")

    assert t.name == "name1"

    with pytest.raises(AttributeError):
        t.name = "name2"

    assert t.name == "name1"

    t._change_name("name2")

    assert t.name == "name2"

    clean()


def test_pickle():

    print("\n\n\n\n")

    import pickle

    root = Node("root")
    node = Node("node")

    root._add_child(node)

    d = pickle.dumps(root)

    root2 = pickle.loads(d)

    assert root2.node.path == 'root.node'
    assert root2.node.name == 'node'

    # Now test pickling a subclass of Node

    root = _SimpleInheritance("root")

    root._placeholder = 5.3

    node = _SimpleInheritance("node")

    root._add_child(node)

    d = pickle.dumps(root)

    root2 = pickle.loads(d)

    assert root2.node.path == 'root.node'
    assert root2.node.name == 'node'

    print(root.placeholder)

    assert root2.placeholder == root.placeholder

    rootc = _ComplexInheritance("root", -1.0, 1.0)
    nodec = _ComplexInheritance("node", -5.0, 5.0)

    rootc._add_child(nodec)

    d = pickle.dumps(rootc)

    root2c = pickle.loads(d)

    print(root2c.min_value)

    clean()


