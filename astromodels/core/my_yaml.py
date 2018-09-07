__author__ = 'giacomov'

# The purpose of this module is to customize yaml so that it will load ordered dictionaries instead of normal
# ones. This way the order in which things are expressed in the file is maintained.

import yaml as my_yaml

import collections

_mapping_tag = my_yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG


def dict_representer(dumper, data):
    return dumper.represent_dict(list(data.items()))


def dict_constructor(loader, node):
    return collections.OrderedDict(loader.construct_pairs(node))


my_yaml.add_representer(collections.OrderedDict, dict_representer)
my_yaml.add_constructor(_mapping_tag, dict_constructor)
