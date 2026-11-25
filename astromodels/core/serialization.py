import copyreg

from astromodels.core.model import Model
from astromodels.core.model_parser import ModelParser


# Model serializer
def serialize_model(model):

    dict_with_types = model.to_dict_with_types()

    return unserialize_model, (dict_with_types,)


def unserialize_model(dict_with_types):

    return ModelParser(model_dict=dict_with_types).get_model()


# Register serialization/unserialization for Model
copyreg.constructor(unserialize_model)
copyreg.pickle(Model, serialize_model, unserialize_model)
