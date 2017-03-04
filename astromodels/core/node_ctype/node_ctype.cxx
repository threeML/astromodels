/*
Created by giacomov on 2/24/17.

 This class implements a new type. It is possible to add children at run time which are read-only, and can be
 accessed as members of the class. So for example:

 > root = Node("root")
 > node = Node("node")
 > root._add_child(node)
 > assert root.node == node
 True
 > root.node = 'this does not work'
 AttributeError: You cannot override a node.

 The child can be any Python object.

*/

#pragma GCC diagnostic ignored "-Wwrite-strings"

#include <Python.h>
#include "structmember.h"

#include <map>
#include <list>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <iterator>
#include <algorithm>


// Forward declaration so that we will be able to use the NodeType in the function defined at the
// end of the file
bool is_a_node(PyObject *obj);

typedef std::map<std::string, PyObject*> nodes_map;
typedef std::vector<PyObject*> nodes_order;

void trace(std::string msg)
{
#ifndef NDEBUG
  std::cerr << msg << std::endl;
#endif
}

// A generic utility to split strings

template<typename Out>
void split(const std::string &s, char delim, Out result) {
  std::stringstream ss;
  ss.str(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    *(result++) = item;
  }
}


std::vector<std::string> split(const std::string &s, char delim) {
  std::vector<std::string> elems;
  split(s, delim, std::back_inserter(elems));
  return elems;
}
// Class definition

typedef struct {
  PyObject_HEAD

  // We save children in a map, and their insertion order in a vector, so there is fast
  // access by key, but we also keep the insertion order
  nodes_map nodes;
  nodes_order order;

  PyObject *parent;

  PyObject *name;

  // This is so that this class behave like any other (it can be monkey-patched)
  PyObject* dict;

} Node;

// __new__ equivalent
static PyObject * Node_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{

    trace("NEW");

    Node* self = (Node *)type->tp_alloc(type, 0);

    if (self != NULL) {

        self->name = PyString_FromString("-- unset --");

        if (self->name == NULL)
        {
          Py_DECREF(self);
          return NULL;
        }

        trace("name is set");

        Py_INCREF(Py_None);

        self->parent = Py_None;

        if (self->parent == NULL)
        {
          Py_DECREF(self);
          return NULL;
        }

        trace("parent is set");

        self->nodes.clear();

        self->order.clear();

        trace("nodes are clear");

        self->dict = PyDict_New();

        if (self->dict == NULL)
        {
          Py_DECREF(self);
          return NULL;
        }

        trace("empty __dict__ is set");

    }

    return (PyObject *)self;
}

// __init__ equivalent
static int
Node_init(Node *self, PyObject *args, PyObject *kwds)
{

  trace("init" );

  PyObject *name=NULL, *tmp;

  static char *kwlist[] = {"name", NULL};

  if (! PyArg_ParseTupleAndKeywords(args, kwds, "O", kwlist,
                                    &name))
  {
    trace("error");

    PyErr_SetString(PyExc_SyntaxError, "You have to provide a name for the node");

    return -1;
  }

  if (name) {
    tmp = self->name;
    Py_INCREF(name);
    self->name = name;
    Py_XDECREF(tmp);
  }

  self->nodes.clear();
  self->order.clear();

  Py_INCREF(self->dict);

  return 0;
}

// Destructor

static void
Node_dealloc(Node *self) {

  trace("dealloc");

  if (self) {

    /* Loop over the map and decrease the reference for all objects*/

    if (self->nodes.size() > 0)
    {

      for (nodes_map::iterator it = self->nodes.begin(); it != self->nodes.end(); ++it) {

        Py_XDECREF(it->second);

      }

      // Clean up map
      self->nodes.clear();
      self->order.clear();

    }

    // Clean up name
    Py_XDECREF(self->name);

    // Free object

    Py_TYPE(self)->tp_free((PyObject *) self);

  }

}

// Set parent
static PyObject *
node_set_parent(Node *self, PyObject *args)
{

  trace("set_parent");

  PyObject *parent;

  // Get the input as char array
  if (!PyArg_ParseTuple(args, "O", &parent))
    return NULL;

  // Increase reference count (we will store this object)

  Py_INCREF(parent);

  self->parent = parent;

  Py_RETURN_NONE;

}

// Add a child
static PyObject *
node_add_child(Node *self, PyObject *args)
{

  trace("add_child" );

  PyObject *child;

  // Get the input as char array
  if (!PyArg_ParseTuple(args, "O", &child))
    return NULL;

  // Make sure the object is a node

  if (not is_a_node(child))
  {

    PyErr_SetString(PyExc_ValueError, "You can only add a Node as a child of a Node");

    return NULL;

  }

  // Now get the name of the object
  std::string attribute_name = PyString_AsString(((Node *) child)->name);

  // Increase reference count (we will store this object)

  Py_INCREF(child);

  // Set it as attribute of the class
  PyObject_SetAttrString((PyObject*) self, attribute_name.c_str(), child);

  // Add to the map
  self->nodes[attribute_name] = child;

  // Add to the vector
  self->order.push_back(child);

  // Make the current node the parent of the child
  PyObject *tuple = Py_BuildValue("(O)", self);
  node_set_parent((Node *) child, tuple);

  Py_RETURN_NONE;
}


// Add a list of children
static PyObject *
node_add_children(Node *self, PyObject *args)
{

  PyObject* obj;
  PyObject* seq;
  int i, len;

  if (!PyArg_ParseTuple(args, "O", &obj))
    return NULL;

  seq = PySequence_Fast(obj, "expected a sequence");
  len = PySequence_Size(obj);

  for (i = 0; i < len; i++) {

    PyObject *item = PySequence_Fast_GET_ITEM(seq, i);

    PyObject *tuple = Py_BuildValue("(O)", item);

    node_add_child(self, tuple);

  }
  Py_DECREF(seq);

  Py_RETURN_NONE;

}


// Get a child by name
static PyObject *
node_get_child(Node *self, PyObject *args)
{

  const char *child_name;

  // Get the input as char array
  if (!PyArg_ParseTuple(args, "s", &child_name))
    return NULL;

  std::string child_name_str(child_name);

  nodes_map::iterator it = self->nodes.find(child_name_str);

  if (it == self->nodes.end())
  {

    std::string error("Unknown child ");
    error += child_name_str;

    PyErr_SetString(PyExc_AttributeError, error.c_str());

    return NULL;

  }

  Py_INCREF((it->second));

  return it->second;

}

// Remove a child
static PyObject *
node_remove_child(Node *self, PyObject *args)
{

  const char *child_name;

  // Get the input as char array
  if (!PyArg_ParseTuple(args, "s", &child_name))
    return NULL;

  std::string child_name_str(child_name);

  nodes_map::iterator it = self->nodes.find(child_name_str);

  if (it == self->nodes.end())
  {

    std::string error("Unknown child ");
    error += child_name_str;

    PyErr_SetString(PyExc_AttributeError, error.c_str());

    return NULL;

  }

  // Now remove from the vector
  nodes_order::iterator it_v = std::find(self->order.begin(), self->order.end(), it->second);

  if (it_v == self->order.end())
  {

    PyErr_SetString(PyExc_AttributeError, "This is a BUG! Order and map are out of syn in Node");

    return NULL;

  }

  self->nodes.erase(it);

  self->order.erase(it_v);

  // Remove also as attribute

  if (PyObject_DelAttrString((PyObject *) self, child_name) < 0)
  {

      return NULL;

  }

  // Make its parent None

  Py_INCREF(Py_None);

  ((Node *) it->second)->parent = Py_None;

  Py_INCREF(it->second);

  return it->second;

}

// Get the parent of this node (or None if not set)
static PyObject *
node_get_parent(Node *self, PyObject *args)
{

  trace("get_parent");

  if (self->parent)
  {

    // Increase reference count (as the user will receive this)

    Py_INCREF(self->parent);

    return self->parent;

  } else {

    // No parent set, return None

    Py_RETURN_NONE;

  }


}

// Get path of this node, i.e., a dotted-separated path like node1.node2.node3
// (in other words, the path to this tree starting from the node)
static PyObject *
node_get_path(Node *self, PyObject *args)
{

  trace("get_path");

  // We start navigating from the present node
  Node *navigator = self;

  // Prepare the list which will contain the results
  std::list<std::string> path;

  // Infinite loop, will exit with break

  while (true)
  {

    std::string this_name = PyString_AsString(navigator->name);

    trace(this_name);

    // If we arrived at the root, stop
    if (this_name == "__root__")
    {

      break;

    }

    // Insert at the beginning of the list because we are starting from the last node and going backward

    path.push_front(this_name);

    if (navigator->parent != Py_None)
    {

      navigator = (Node*) navigator->parent;

      continue;

    } else {

      // End
      break;

    }

  }

  // Now make the string like node1.node2.node3
  std::string path_string;

  for (std::list<std::string>::iterator it = path.begin(); it != path.end(); ++it)
  {

    path_string += *it;

    if (*it != path.back())
    {

      path_string += ".";

    }

  }

  trace(path_string);

  PyObject *path_string_py = PyString_FromString(path_string.c_str());

  Py_INCREF(path_string_py);

  return path_string_py;

}


// Name getter
static PyObject *
node_getname(Node *self, PyObject *args)
{

  Py_INCREF(self->name);

  return self->name;

}

static PyObject *
node_get_children(Node *self, PyObject *args)
{

  PyObject *tuple = PyTuple_New(self->nodes.size());

  Py_ssize_t pos = 0;

  // Iterate over the vector, instead of the map, as the vector is faster for sequential access
  // and it preserves the order of insertion

  for (nodes_order::iterator it = self->order.begin(); it != self->order.end(); ++it)
  {

    Py_INCREF(*it);

    PyTuple_SetItem(tuple, pos, *it);

    pos += 1;

  }

  Py_INCREF(tuple);

  return tuple;

}



// This changes the name of the node
static PyObject *
node_change_name(Node *self, PyObject *args)
{

  trace("change name");

  PyObject *value;

  // Get the input as char array
  if (!PyArg_ParseTuple(args, "O", &value))
    return NULL;

  if (! PyString_Check(value)) {
    PyErr_SetString(PyExc_TypeError,
                    "The name of a node must be a string");
    return NULL;
  }
  Py_DECREF(self->name);
  Py_INCREF(value);
  self->name = value;

  Py_RETURN_NONE;
}


// Attribute setter
int
node_setattro(PyObject *obj, PyObject *name, PyObject *value)
{

  // Increment references of objects we are going to use
  Py_INCREF(obj);
  Py_INCREF(name);

  trace("setattro" );

  // Explicitly cast to right type

  Node *self = (Node *) obj;

  // Lookup for attribute in the map of the nodes, to see if it's a node

  std::string name_string = PyString_AsString(name);

  if (name_string == "name")
  {

    // Cannot change the name like this

    PyErr_SetString(PyExc_AttributeError, "You cannot change the name of the node");

    Py_DECREF(obj);
    Py_DECREF(name);

    return -1;

  }

  nodes_map::iterator it = self->nodes.find(name_string);

  if (it != self->nodes.end()) {

    // This is a node

    // the user cannot assign to a node, unless the node itself has a .value member, in which case
    // the value is assigned to that
    PyObject *value_attr = PyObject_GetAttrString(it->second, "value");

    if (value_attr == NULL)
    {
      // The node does not have  a .value attribute, so this is forbidden

      PyErr_SetString(PyExc_AttributeError, "You cannot override a node.");

      // Decrease references that we have increased at the beginning of the method

      Py_DECREF(obj);
      Py_DECREF(name);

      return -1;
    } else
    {

      // The object has a value attribute
      trace("Found a value attribute" );

      //Py_INCREF(value);  Commented out as GenericSetAttr should incremeant this already

      // Set the .value attribute of the node to the provided value

      PyObject_SetAttrString(it->second, "value", value);

      return 0;

    }

  } else {

    // This is a normal attribute (i.e., not a node)

    //Py_INCREF(value);  Commented out as GenericSetAttr should incremeant this already

    // call the normal setter
    if (PyObject_GenericSetAttr(obj, name, value) == -1)
    {

       Py_DECREF(obj);
       Py_DECREF(name);

       // The exception is already set by GenericSetAttr

       return -1;

    }

    // Decrease references that we have increased at the beginning of the method
    Py_DECREF(obj);
    Py_DECREF(name);

    return 0;

  }

}


bool has_child(Node *node, const std::string& child_name)
{

  trace("received");
  trace(child_name);

  if (node->nodes.count(child_name) > 0)
  {

    trace("has child");

    return true;

  } else {

    trace("does not have child");

    return false;

  }

}


static PyObject *
node_has_child(Node *self, PyObject *args)
{

  trace("node_has_child");

  const char *child_name;

  // Get the input as char array
  if (!PyArg_ParseTuple(args, "s", &child_name))
    return NULL;

  std::string child_name_str(child_name);

  if (has_child(self, child_name_str))
  {

    Py_RETURN_TRUE;

  } else
  {

    Py_RETURN_FALSE;

  }

}



// Get a child from a path of the type node1.node2.node3 (returns node3 in this case)
static PyObject *
node_get_child_from_path(Node *self, PyObject *args)
{

  trace("get_child_from_path");

  const char *path;

  // Get the input as char array
  if (!PyArg_ParseTuple(args, "s", &path))
    return NULL;

  // Split the string node1.node2.node3
  std::string path_s(path);

  std::vector<std::string> nodes_names = split(path, '.');

  for (std::vector<std::string>::iterator it=nodes_names.begin(); it != nodes_names.end(); ++it) {

    trace(*it);

  }

  // This will point to the current node and eventually to the last one
  Node *this_node = self;

  // loop over the nodes
  for (std::vector<std::string>::iterator it=nodes_names.begin(); it != nodes_names.end(); ++it)
  {

    trace((*it));

    if (has_child(this_node, *it))
    {

      trace("yes");

      // Update the pointer

      this_node = (Node *) this_node->nodes[*it];

    } else
    {

      trace("no");

      // Path is wrong

      std::string error("Node ");
      error += (*it);
      error += " does not exist";

      PyErr_SetString(PyExc_AttributeError, error.c_str());

      return NULL;

    }
  }

  return (PyObject *) this_node;

}


//// This is needed to pickle the object
//static PyObject *
//node_reduce(Node *self, PyObject *args)
//{
//
//  trace("****reduce");
//
//  return PyType_GenericNew, NodeType
//
//}

static PyMemberDef Node_members[] = {
    {"__dict__", T_OBJECT_EX, offsetof(Node, dict), 0, "__dict__"},
    {NULL}  /* Sentinel */
};

// Methods of the Node class
static PyMethodDef Node_methods[] = {
    {"_add_child", (PyCFunction) node_add_child, METH_VARARGS, "Add a child to this node"},
    {"_get_child", (PyCFunction) node_get_child, METH_VARARGS, "Get a child by name"},
    {"_add_children", (PyCFunction) node_add_children, METH_VARARGS, "Add a list of children to this node"},
    {"_get_children", (PyCFunction) node_get_children, METH_VARARGS, "Get a tuple with the children of this node"},
    {"_remove_child", (PyCFunction) node_remove_child, METH_VARARGS, "Remove and return a child by name"},
    {"_get_parent", (PyCFunction) node_get_parent, METH_VARARGS, "Get the parent of this node (or None if not set)"},
    {"_get_path", (PyCFunction) node_get_path, METH_VARARGS, "Get the path of the current node as a list"},
    {"_change_name", (PyCFunction) node_change_name, METH_VARARGS, "Change name of the node (careful!)"},
    {"_get_child_from_path", (PyCFunction) node_get_child_from_path, METH_VARARGS, "Get a node from a path"},
    {"_has_child", (PyCFunction) node_has_child, METH_VARARGS, "Return whether the node has the named child"},
//    {"__reduce__", (PyCFunction) node_reduce, METH_VARARGS, "For pickling the node"},
    {NULL}  /* Sentinel */
};

// Properties of the Noode class
static PyGetSetDef Node_getseters[] = {
    {"name", (getter) node_getname, NULL, "name", NULL},
    {"path", (getter) node_get_path, NULL, "path", NULL},
    {NULL}  /* Sentinel */
};

static PyTypeObject NodeType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "node_ctype._Node",        /* tp_name */
    sizeof(Node),              /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)Node_dealloc,  /* tp_dealloc */
    0,                         /* tp_print */
    0,                         /* tp_getattr */
    0,                         /* tp_setattr */
    0,                         /* tp_compare */
    0,                         /* tp_repr */
    0,                         /* tp_as_number */
    0,                         /* tp_as_sequence */
    0,                         /* tp_as_mapping */
    0,                         /* tp_hash */
    0,                         /* tp_call */
    0,                         /* tp_str */
    0,                         /* tp_getattro */
    node_setattro,             /* tp_setattro */
    0,                         /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT |
        Py_TPFLAGS_BASETYPE,   /* tp_flags */
    "Node objects",            /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    Node_methods,              /* tp_methods */
    Node_members,              /* tp_members */
    Node_getseters,            /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    offsetof(Node, dict),      /* tp_dictoffset */
    (initproc) Node_init,      /* tp_init */
    0,       /* tp_alloc */
    Node_new,                  /* tp_new */
};

static PyMethodDef module_methods[] = {
    {NULL}  /* Sentinel */
};

#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
initnode_ctype(void)
{
  PyObject* m;

  if (PyType_Ready(&NodeType) < 0)
    return;

  m = Py_InitModule3("node_ctype", module_methods,
                     "Example module that creates an extension type.");

  if (m == NULL)
    return;

  Py_INCREF(&NodeType);
  PyModule_AddObject(m, "_Node", (PyObject *)&NodeType);
}

bool is_a_node(PyObject *obj)
{

  return PyObject_IsInstance(obj, reinterpret_cast<PyObject*>(&NodeType));

}