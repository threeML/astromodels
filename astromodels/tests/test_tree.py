import pytest
import os


os.environ["ASTROMODELS_DEBUG"] = "debug"

from astromodels.core.tree import Node
from astromodels.core import node_ctype

import gc

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

        self._min_value = min_value
        self._max_value = max_value

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

    clara_ref_count = node_ctype._get_reference_counts(clara)

    t._add_child(clara)

    # We expect 2 more references: one for the map, one for the vector

    assert node_ctype._get_reference_counts(clara) == clara_ref_count + 1

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

    with pytest.raises(AttributeError):
        t._get_child("not existing")

    clean()


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
        print t.node

    with pytest.raises(AttributeError):
        t._get_child("node")

    clean()


def test_get_path():
    # Make a small tree

    node1 = Node('node1')
    node2 = Node('node2')
    node3 = Node('node3')

    node1._add_child(node2)
    node2._add_child(node3)

    print node3.path

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

    print root.placeholder

    assert root2.placeholder == root.placeholder

    rootc = _ComplexInheritance("root", -1.0, 1.0)
    nodec = _ComplexInheritance("node", -5.0, 5.0)

    rootc._add_child(nodec)

    d = pickle.dumps(rootc)

    root2c = pickle.loads(d)

    print root2c.min_value

    clean()


def test_memory_leaks_setters():

    root = Node("root")

    node = Node("node")

    refc_before_link = node_ctype._get_reference_counts(node)

    root._add_child(node)

    # Adding a node adds 1 reference

    assert node_ctype._get_reference_counts(node) == refc_before_link + 1

    # Remove the node and verify that the reference counts goes back to what it was
    root._remove_child("node")

    # Now we should go back to the original

    assert node_ctype._get_reference_counts(node) == refc_before_link

    # Now add a second node nested under the first one (root.node.node2)
    node2 = Node("node2")

    refc_before_link2 = node_ctype._get_reference_counts(node2)

    root._add_child(node)  # +1 for node

    assert node_ctype._get_reference_counts(node) == refc_before_link + 1

    node._add_child(node2)  # +1 for node2 and +1 for node

    assert node_ctype._get_reference_counts(node2) == refc_before_link2 + 1

    # "node" is now also parent of node2, so its reference are now 2 more than the original
    assert node_ctype._get_reference_counts(node) == refc_before_link + 2

    # Clean up and verify that we go back to the original number of references
    node._remove_child("node2")  # -1 for node2, -1 for node

    assert node_ctype._get_reference_counts(node2) == refc_before_link2

    assert node_ctype._get_reference_counts(node) == refc_before_link + 1

    root._remove_child("node")  # -1 for node

    assert node_ctype._get_reference_counts(node) == refc_before_link

    # Now test add_children
    node3 = Node("node3")
    node4 = Node("node4")

    refc_before_link3 = node_ctype._get_reference_counts(node3)
    refc_before_link4 = node_ctype._get_reference_counts(node4)

    root._add_children([node3, node4])  # +1 for both nodes

    assert node_ctype._get_reference_counts(node3) == refc_before_link3 + 1
    assert node_ctype._get_reference_counts(node3) == refc_before_link4 + 1


def test_memory_leaks_getters():

    # Now test the getters

    root = Node("root")

    node1 = Node("node1")
    node2 = Node("node2")
    node3 = Node("node3")

    refc_before_link_root = node_ctype._get_reference_counts(root)
    refc_before_link1 = node_ctype._get_reference_counts(node1)
    refc_before_link2 = node_ctype._get_reference_counts(node2)
    refc_before_link3 = node_ctype._get_reference_counts(node3)

    None_counts_before = node_ctype._get_reference_counts(None)

    # Add 3 nodes to root
    root._add_children([node1, node2, node3]) # node1 +1, node2 +1, node3 + 1, root +3

    node_again = root._get_child("node1")  # node1 +1

    assert node_ctype._get_reference_counts(node1) == refc_before_link1 + 2

    del node_again # node1 -1

    assert node_ctype._get_reference_counts(node1) == refc_before_link1 + 1

    children = root._get_children()  #node1 +1, node2 +1, node3 +1

    assert len(children) == 3

    assert node_ctype._get_reference_counts(node1) == refc_before_link1 + 2
    assert node_ctype._get_reference_counts(node2) == refc_before_link2 + 2
    assert node_ctype._get_reference_counts(node3) == refc_before_link3 + 2
    assert node_ctype._get_reference_counts(root) == refc_before_link_root + 3

    # test get_parent

    root_again = node1._get_parent()  # root +1

    assert node_ctype._get_reference_counts(root) == refc_before_link_root + 4

    del root_again  # root -1

    assert node_ctype._get_reference_counts(root) == refc_before_link_root + 3

    # test _get_path

    path = node2._get_path()  # this shouldn't change any ref count

    assert node_ctype._get_reference_counts(node1) == refc_before_link1 + 2
    assert node_ctype._get_reference_counts(node2) == refc_before_link2 + 2
    assert node_ctype._get_reference_counts(node3) == refc_before_link3 + 2
    assert node_ctype._get_reference_counts(root) == refc_before_link_root + 3

    # test _get_child_from_path

    node4 = Node("node4")
    refc_before_link4 = node_ctype._get_reference_counts(node4)

    node3._add_child(node4)  # node3 +1, node4 + 1

    node4_again = root._get_child_from_path("node3.node4")  # node4 +1

    assert node_ctype._get_reference_counts(node4) == refc_before_link4 + 2
    assert node_ctype._get_reference_counts(node3) == refc_before_link3 + 3

    del node4_again  # node4 -1

    assert node_ctype._get_reference_counts(node4) == refc_before_link4 + 1


def test_memory_leaks_destructors():

    print("\n\n\n\n\n")

    for i in range(1000):

        print("\n\nRound %i" % (i+1))

        # Test destructors
        root = Node("root")

        node1 = Node("node1")
        node2 = Node("node2")
        node3 = Node("node3")

        refc_before_link_root = node_ctype._get_reference_counts(root)

        refc_before_link1 = node_ctype._get_reference_counts(node1)
        refc_before_link2 = node_ctype._get_reference_counts(node2)
        refc_before_link3 = node_ctype._get_reference_counts(node3)

        # root._add_children([node1, node2, node3])  # node1 +1, node2 + 1, node3 + 1, root +3

        root._add_child(node1)
        root._add_child(node2)
        root._add_child(node3)


        assert node_ctype._get_reference_counts(root) == refc_before_link_root + 3
        assert node_ctype._get_reference_counts(node1) == refc_before_link1 + 1
        assert node_ctype._get_reference_counts(node2) == refc_before_link2 + 1
        assert node_ctype._get_reference_counts(node3) == refc_before_link3 + 1

        #
        # Let's destroy the tree
        root._remove_child("node1")

        root._remove_child("node2")

        root._remove_child("node3")
        #
        # assert node_ctype._get_reference_counts(root) == refc_before_link_root + 3
        #
        # # Now deleting root should have decreased all nodes by 3, i..e, should have got them back to the initial
        # # reference counts
        # assert node_ctype._get_reference_counts(node1) == refc_before_link1
        # assert node_ctype._get_reference_counts(node2) == refc_before_link2
        # assert node_ctype._get_reference_counts(node3) == refc_before_link3
