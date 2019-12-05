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
#include "bytesobject.h"

#include <map>
#include <list>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <iterator>
#include <algorithm>


// These two macros needs to be used to use NAME_MAXLENGTH in a string

#define STR_EXPAND(tok) #tok
#define STR(tok) STR_EXPAND(tok)

#define NAME_MAXLENGTH 50

// Utility function to help python2 <-> python3
#if PY_MAJOR_VERSION >= 3

const char *strobj_to_string(PyObject *obj)
{

  return PyUnicode_AsUTF8(obj);

}

PyObject *strobj_from_string(const char *str)
{

  return PyUnicode_FromString(str);

}

int check_string(PyObject *obj) {  return PyUnicode_Check(obj); }

#else

char *strobj_to_string(PyObject *obj)
{

  return PyString_AsString(obj);

}

PyObject *strobj_from_string(const char *str)
{

  return PyString_FromString(str);

}

int check_string(PyObject *obj) {

    int res = PyString_Check(obj);

    if(res==0)
    {

      // Maybe a unicode literal
      res = PyUnicode_Check(obj);

      return res;

    } else {

      return 1;

    }

    }

#endif

// Forward declaration so that we will be able to use the NodeType in the function defined at the
// end of the file
bool is_a_node(PyObject *obj);

typedef std::map<std::string, PyObject*> nodes_map;
typedef std::vector<PyObject*> nodes_order;


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

  // This is so that this class behave like any other (it can be monkey-patched)
  PyObject* dict;

  PyObject *parent;

  // We save children in a map, and their insertion order in a vector, so there is fast
  // access by key, but we also keep the insertion order
  // NOTE: we can allocate these on the stack because they will always have the same size
  // on the stack. Growing or shrinking these structures will allocate memory on the heap

  nodes_map nodes;
  nodes_order order;

  // Name of the node
  // We would like to use std::string, but this will cause a segfault on gcc 4.4,
  // because the implementation there use memory on the stack if the string is short enough,
  // instead of using the heap. This struct instead should be a POD (Plain Old Data) structure
  // for it to work well with python
  char *name;

} Node;


inline void trace(std::string msg, Node *node=NULL)
{
#ifndef NDEBUG
  if (node)
  {

     std::cerr << node->name << ": " << msg << std::endl;

  } else
  {

     std::cerr << msg << std::endl;

  }

#endif
}


// __new__ equivalent
static PyObject * Node_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{

    Node* self = (Node *)type->tp_alloc(type, 0);

    if (self != NULL) {

        // Maximum length 50 characters
        self->name = (char *) malloc (NAME_MAXLENGTH * sizeof(char));

        strcpy(self->name, "-- unset --");

        self->parent = NULL;

        self->nodes.clear();

        self->order.clear();

        self->dict = PyDict_New();

        if (self->dict == NULL)
        {
          Py_DECREF(self);
          return NULL;
        }

    }

    return (PyObject *)self;
}

// __init__ equivalent
static int
Node_init(Node *self, PyObject *args, PyObject *kwds)
{

  const char *name;

  static char *kwlist[] = {"name", NULL};

  if (! PyArg_ParseTupleAndKeywords(args, kwds, "s", kwlist,
                                    &name))
  {

    PyErr_SetString(PyExc_SyntaxError, "You have to provide a name for the node");

    return -1;
  }

  // Make sure the name is not longer than the allowed maximum
  if (strlen(name) >= NAME_MAXLENGTH)
  {

    std::string msg = "The name for the node cannot be longer than ";
    msg += STR(NAME_MAXLENGTH);
    msg += " characters";

    PyErr_SetString(PyExc_SyntaxError, msg.c_str());

    return -1;

  }

  strcpy(self->name, name);

  return 0;
}

// Destructor functions

// The following two methods are needed to support the garbage collector
static int
node_traverse(Node *self, visitproc visit, void *arg)
{

    // Loop over all children
    for (nodes_order::iterator it=self->order.begin(); it != self->order.end(); ++it)
    {

      Py_VISIT(*it);

    }

    // Visit also the parent

    if(self->parent)
    {
      Py_VISIT(self->parent);
    }

    return 0;
}

static int
node_clear(Node *self)
{

    for (nodes_order::iterator it=self->order.begin(); it != self->order.end(); ++it)
    {
        Py_XDECREF(*it);

    }

    self->nodes.clear();
    self->order.clear();

    if (self->parent)
    {

         Py_CLEAR(self->parent);

    }

    return 0;
}

// This deallocate the memory

static void
Node_dealloc(Node *self) {

    if (self)
    {

        free(self->name);

        node_clear(self);

        // Free object
        Py_TYPE(self)->tp_free((PyObject *) self);
    }

}

// Set parent


static PyObject* node_set_parent(Node *self, Node *parent)
{

  if (self->parent)
  {

      Py_CLEAR(self->parent);

  }

  // Increase reference count (we will store this object)

  Py_INCREF(parent);

  self->parent = (PyObject *) parent;

  Py_RETURN_NONE;

}


int _add_child(Node *self, PyObject *child)
{

  // Make sure the object is a node
  if (not is_a_node(child))
  {

    PyErr_SetString(PyExc_TypeError, "You can only add a Node as a child of a Node");

    return -1;

  }

  // Now get the name of the object
  std::string attribute_name = ((Node *) child)->name;

  // Verify that the child is not already contained in the node
  nodes_map::iterator it = self->nodes.find(attribute_name);

  if (it != self->nodes.end())
  {

      std::string msg = "A child named ";
      msg += ((Node *)(it->second))->name;
      msg += " already exists";

      PyErr_SetString(PyExc_AttributeError, msg.c_str());

      return -1;

  }

  // Increase the reference counts by one
  Py_INCREF(child);

  // Add to the map
  self->nodes[attribute_name] = child;

  // Add to the vector
  self->order.push_back(child);

  // Make the current node the parent of the child

  node_set_parent((Node *) child, self);

  return 0;

}

// Add a child
static PyObject *
node_add_child(Node *self, PyObject *args)
{

  PyObject *child;

  // Get the input as char array
  if (!PyArg_ParseTuple(args, "O", &child))
    return NULL;

  if (_add_child(self, child) < 0)
  {

    return NULL;

  }

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

  seq = PySequence_Tuple(obj);
  len = PySequence_Size(obj);

  for (i = 0; i < len; i++) {

    PyObject *child = PySequence_GetItem(seq, i);  // item +1

    if (_add_child(self, child) < 0)
    {

      Py_XDECREF(child);
      Py_XDECREF(seq);

      return NULL;

    }

    Py_XDECREF(child);

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

  PyObject *tmp = *it_v;

  // Remove from map and vector

  self->nodes.erase(it);

  self->order.erase(it_v);

  // Remove reference to the parent

  Py_CLEAR(((Node *)tmp)->parent);

  Py_DECREF(tmp);

  Py_RETURN_NONE;

}

// Get the parent of this node (or None if not set)
static PyObject *
node_get_parent(Node *self, PyObject *args)
{

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

  // We start navigating from the present node
  Node *navigator = self;

  // Prepare the list which will contain the results
  std::list<std::string> path;

  // Infinite loop, will exit with break

  while (true)
  {

    std::string this_name = navigator->name;

    // If we arrived at the root, stop
    if (this_name == "__root__")
    {

      break;

    }

    // Insert at the beginning of the list because we are starting from the last node and going backward

    path.push_front(this_name);

    if (navigator->parent)
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

  for (std::list<std::string>::iterator it = path.begin(); it != --path.end(); ++it)
  {

    path_string += *it;

    path_string += ".";

  }

  path_string += *(--path.end());

  PyObject *path_string_py = strobj_from_string(path_string.c_str());

  return path_string_py;

}


// Name getter
static PyObject *
node_getname(Node *self, PyObject *args)
{

  PyObject *name_str = strobj_from_string(self->name);

  Py_INCREF(name_str);

  return name_str;

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

  PyObject *value;

  // Get the input as char array
  if (!PyArg_ParseTuple(args, "O", &value))
    return NULL;

  if (! check_string(value)) {
    PyErr_SetString(PyExc_TypeError,
                    "The name of a node must be a string");
    return NULL;
  }

  std::string new_name = strobj_to_string(value);

  strcpy(self->name, new_name.c_str());

  Py_RETURN_NONE;
}

// Attribute getter
static PyObject *
node_getattro(Node *self, PyObject *name)
{

  std::string name_string = strobj_to_string(name);

  nodes_map::iterator it = self->nodes.find(name_string);

  if (it != self->nodes.end())
  {

      // This is a node

      Node *this_node = (Node *) (it->second);

      if(!this_node)
      {

          return NULL;

      }

      Py_INCREF(this_node);

      return (PyObject *) this_node;

  } else
  {
      // This is not a node

      return PyObject_GenericGetAttr((PyObject *) self, name);

  }

}

// Attribute setter
int
node_setattro(PyObject *obj, PyObject *name, PyObject *value)
{

  // Explicitly cast to right type

  Node *self = (Node *) obj;

  // Lookup for attribute in the map of the nodes, to see if it's a node

  std::string name_string = strobj_to_string(name);

  if (name_string == "name")
  {

    // Cannot change the name like this

    PyErr_SetString(PyExc_AttributeError, "You cannot change the name of the node");

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

      return -1;
    } else
    {

      // The object has a value attribute

      // Decrease the reference count which was increased by PyObject_GetAttrString
      Py_XDECREF(value_attr);

      // Set the .value attribute of the node to the provided value
      // (this call will increase the reference count again)

      return PyObject_SetAttrString(it->second, "value", value);

    }

  } else {

    // This is a normal attribute (i.e., not a node)

    //Py_INCREF(value);  Commented out as GenericSetAttr should incremeant this already

    // call the normal setter
    // NOTE: we cannot call PyObject_SetAttr because that would call into setattro again, giving infinite recursion
    return PyObject_GenericSetAttr(obj, name, value);

  }

}


bool has_child(Node *node, const std::string& child_name)
{

  if (node->nodes.count(child_name) > 0)
  {

    return true;

  } else {

    return false;

  }

}


static PyObject *
node_has_child(Node *self, PyObject *args)
{

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

  const char *path;

  // Get the input as char array
  if (!PyArg_ParseTuple(args, "s", &path))
    return NULL;

  // Split the string node1.node2.node3
  std::string path_s(path);

  std::vector<std::string> nodes_names = split(path, '.');

  // This will point to the current node and eventually to the last one
  Node *this_node = self;

  // loop over the nodes
  for (std::vector<std::string>::iterator it=nodes_names.begin(); it != nodes_names.end(); ++it)
  {

    if (has_child(this_node, *it))
    {

      // Update the pointer

      this_node = (Node *) this_node->nodes[*it];

    } else
    {

      // Path is wrong

      std::string error("Node ");
      error += (*it);
      error += " does not exist";

      PyErr_SetString(PyExc_AttributeError, error.c_str());

      return NULL;

    }
  }

  Py_INCREF(this_node);

  return (PyObject *) this_node;

}


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
    "node_ctype._Node",             /* tp_name */
    sizeof(Node),                   /* tp_basicsize */
    0,                              /* tp_itemsize */
    (destructor)Node_dealloc,       /* tp_dealloc */
    0,                              /* tp_print */
    0,                              /* tp_getattr */
    0,                              /* tp_setattr */
    0,                              /* tp_compare */
    0,                              /* tp_repr */
    0,                              /* tp_as_number */
    0,                              /* tp_as_sequence */
    0,                              /* tp_as_mapping */
    0,                              /* tp_hash */
    0,                              /* tp_call */
    0,                              /* tp_str */
    (getattrofunc) node_getattro,   /* tp_getattro */
    (setattrofunc) node_setattro,   /* tp_setattro */
    0,                              /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT |
        Py_TPFLAGS_BASETYPE,        /* tp_flags */
    "Node objects",                 /* tp_doc */
    (traverseproc) node_traverse,   /* tp_traverse */
    (inquiry) node_clear,           /* tp_clear */
    0,                              /* tp_richcompare */
    0,                              /* tp_weaklistoffset */
    0,                              /* tp_iter */
    0,                              /* tp_iternext */
    Node_methods,                   /* tp_methods */
    Node_members,                   /* tp_members */
    Node_getseters,                 /* tp_getset */
    0,                              /* tp_base */
    0,                              /* tp_dict */
    0,                              /* tp_descr_get */
    0,                              /* tp_descr_set */
    offsetof(Node, dict),           /* tp_dictoffset */
    (initproc) Node_init,           /* tp_init */
    0,                              /* tp_alloc */
    Node_new,                       /* tp_new */
};

// This function is used to debug the extension

static PyObject *get_reference_counts(PyObject *self, PyObject *args) {

   /* Parse args and do something interesting here. */
   PyObject *obj;

   if (!PyArg_ParseTuple(args, "O", &obj))
      return NULL;

   int counts = (int) (((PyObject*)(obj))->ob_refcnt);

   PyObject *c_py = Py_BuildValue("i", counts);

   Py_INCREF(c_py);

   return c_py;

}

static PyMethodDef module_methods[] = {
    { "_get_reference_counts", get_reference_counts, METH_VARARGS, NULL},
    {NULL}  /* Sentinel */
};


// NEW compatibility layer python2 to python3

struct module_state {
    PyObject *error;
};

#if PY_MAJOR_VERSION >= 3
#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))
#else
#define GETSTATE(m) (&_state)
static struct module_state _state;
#endif


#if PY_MAJOR_VERSION >= 3

static int node_ctype_traverse(PyObject *m, visitproc visit, void *arg) {
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int node_ctype_clear(PyObject *m) {
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}


static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "node_ctype",
        NULL,
        sizeof(struct module_state),
        module_methods,
        NULL,
        node_ctype_traverse,
        node_ctype_clear,
        NULL
};

#endif

#if PY_MAJOR_VERSION >= 3
#define INITERROR return NULL
#else
#define INITERROR return
#endif

#if PY_MAJOR_VERSION >= 3
PyMODINIT_FUNC
PyInit_node_ctype(void)
#else
#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
initnode_ctype(void)
#endif
{
    PyObject *module;

    if (PyType_Ready(&NodeType) < 0)
       INITERROR;

#if PY_MAJOR_VERSION >= 3
    module = PyModule_Create(&moduledef);
#else
    module = Py_InitModule3("node_ctype", module_methods, "Extension type node_ctype");
#endif

    if (module == NULL)
        INITERROR;
    struct module_state *st = GETSTATE(module);

    st->error = PyErr_NewException("node_ctype.Error", NULL, NULL);
    if (st->error == NULL) {
        Py_DECREF(module);
        INITERROR;
    }

    Py_INCREF(&NodeType);
    PyModule_AddObject(module, "_Node", (PyObject *)&NodeType);

#if PY_MAJOR_VERSION >= 3
    return module;
#endif
}

bool is_a_node(PyObject *obj)
{

  return PyObject_IsInstance(obj, reinterpret_cast<PyObject*>(&NodeType));

}