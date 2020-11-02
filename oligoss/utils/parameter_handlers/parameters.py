from typing import List, Dict

from .instrument_handlers import (
    check_instrument_info,
    sanity_check_silico_fragmentation)
from ..silico_utils.silico_utils import retrieve_adduct_info
from ..errors import InputTypeError, TypeDictKeyError, MissingParameterValue
from ..type_dicts.parameter_type_dicts import (
    CORE_PARAM_TYPES,
    SILICO_LIBRARY_TYPES,
    EXTRACTOR_TYPES,
    POSTPROCESS_TYPES,
    NON_ATTR_KEYS
)
from..type_dicts.parameter_fallbacks import (
    CORE_PARAM_FALLBACKS,
    SILICO_LIBRARY_FALLBACKS,
    EXTRACTOR_FALLBACKS,
    POSTPROCESS_FALLBACKS,
    ESSENTIAL_CORE_PARAMS
)

class Parameters():
    """
    Parameters object. To be passed on to other modules for PolymerSoup
    workflows.
    """

    #  properties needed to create instance of Parameters class but will then
    #  be removed during cleanup
    TEMP_PROPS = [
        "type_dict",
        "params_dict",
        "fallbacks",
        "params_class",
        "sub_params",
        "active_instruments"
    ]

    #  dict of parameter class ids and associated type_dicts for checking
    #  parameter values
    PARAM_TYPE_ASSOCIATIONS = {
        "core": CORE_PARAM_TYPES,
        "silico": SILICO_LIBRARY_TYPES,
        "silico_ms2": SILICO_LIBRARY_TYPES["ms2"],
        "silico_ms1": SILICO_LIBRARY_TYPES["ms1"],
        "extractors": EXTRACTOR_TYPES,
        "postprocess": POSTPROCESS_TYPES
    }

    #  dict of parameter class ids and associated fallback parameters for
    #  checking parameter values
    PARAM_FALLBACK_ASSOCIATIONS = {
        "core": CORE_PARAM_FALLBACKS,
        "silico": SILICO_LIBRARY_FALLBACKS,
        "silico_ms1": SILICO_LIBRARY_FALLBACKS["ms1"],
        "silico_ms2": SILICO_LIBRARY_FALLBACKS["ms2"],
        "extractors": EXTRACTOR_FALLBACKS,
        "postprocess": POSTPROCESS_FALLBACKS
    }

    def __init__(
        self,
        params_dict,
        params_class
    ):
        self.params_class = params_class
        self.params_dict = params_dict

        for param in ESSENTIAL_CORE_PARAMS:
            if param not in self.params_dict:
                raise Exception(
                    f'parameter "{param}" is missing from input file\n'
                    f'{self.params_class}, {self.params_dict}')

        self.type_dict, self.fallbacks, self.sub_params = (
            self.get_param_associations())

        self.active_instruments = check_instrument_info(params_dict)

        #  read params dict and set parameter attributes
        self.set_parameter_attributes()
        self.final_attrs_check()

    @classmethod
    def generate_silico(cls, params_dict):
        return cls(params_dict=params_dict, params_class="silico")

    @classmethod
    def generate_silico_ms1(cls, params_dict):
        return cls(params_dict=params_dict, params_class="silico_ms1")

    @classmethod
    def generate_extractors(cls, params_dict):
        return cls(params_dict=params_dict, params_class="extractors")

    @classmethod
    def generate_postprocess(cls, params_dict):
        return cls(params_dict=params_dict, params_class="postprocess")

    def get_param_associations(
        self,
        type_associations=PARAM_TYPE_ASSOCIATIONS,
        fallbacks=PARAM_FALLBACK_ASSOCIATIONS
    ):
        """
        Retrieve type dict and fallbacks for parameter group

        Args:
            type_associations (Dict[str, dict], optional): dict of parameter
                group ids and associated type_dicts.
                Defaults to PARAM_TYPE_ASSOCIATIONS.
            fallbacks (Dict[str, dict], optional): dict of parameter group ids
                and associated fallback dicts for defining parameter values not
                supplied in input file or retrieved from instrument defaults
        Raises:
            Exception: raised if self.params_class not in
                PARAM_TYPE_ASSOCIATIONS

        Returns:
            dict, dict, dict: type_dict, fallbacks, params for parameter group
        """

        #  check whether parameter group id is valid
        if self.params_class not in type_associations:
            raise Exception(
                f'{self.params_class} not a valid parameter group\n'
                f'please choose from: {type_associations}')

        #  retrieve parameters relevant to parameter group from input file
        if self.params_class in ["silico", "extractors", "postprocess"]:
            params = self.params_dict[self.params_class]
        elif self.params_class == "silico_ms1":
            params = self.params_dict["silico"]["ms1"]
        elif self.params_class == "silico_ms2":
            params = self.params_dict["silico"]["ms2"]
        elif self.params_class == "core":
            params = self.params_dict
        else:
            raise Exception(
                f'{self.params_class} not a valid parameter group')

        #  get type_dict, fallbacks for relevant parameters
        type_dict = type_associations[self.params_class]
        param_fallbacks = fallbacks[self.params_class]

        #  return: type_dict, fallbacks and params (parameters relevant to
        #  parameter group)
        return type_dict, param_fallbacks, params

    def set_parameter_attributes(self):
        """
        Set attributes of Parameters depending on what set of parameter inputs
        the Parameters class is to be an instance of (silico, postprocess, core
        etc...)

        Args:

        Raises:
            MissingParameterValue: raised if the value for an essential
                parameter is not supplied.
            InputTypeError: raised if the value supplied for a parameter does
                not match expected type is not convertable to expected type.
        """

        #  Iterate through parameters and associated values, check which
        #  parameters have default instrument values and retrieve if necessary.
        #  Finally, set each of these parameters as attributes in instance of
        #  Parameters class
        for k, v in self.sub_params.items():

            if k not in self.type_dict and k not in NON_ATTR_KEYS:
                raise Exception(
                    f'{k} not a valid parameter for {self.params_class}')

            # this needs to be 'is None' compared to 'if not v' as it would take
            # keys with 'False' values as having no value
            if v is None:
                if k not in self.type_dict["optional"]:
                    raise MissingParameterValue(
                        key=k,
                        type_dict=self.type_dict,
                        param_class=self.params_class)
                v = self.check_instrument_default_parameter(k)
                if not v:
                    if self.fallbacks:
                        if k in self.fallbacks:
                            v = self.fallbacks[k]
                        else:
                            raise MissingParameterValue(
                                key=k,
                                type_dict=self.type_dict,
                                instrument_fallback=True,
                                param_class=self.params_class)
            if v:
                try:
                    if self.params_class != "core":
                        self.sub_params[k] = self.format_params_attr(k, v)
                except TypeError:
                    raise InputTypeError(
                        key=k, value=v, expected_type=self.type_dict[k])

        #  iterate through type_dict and fill in any further missing parameters
        #  with fallbacks and instrument defaults
        for k, v in self.type_dict.items():
            if k not in NON_ATTR_KEYS:
                if k not in self.sub_params:
                    if self.fallbacks and k in self.fallbacks:
                        self.sub_params[k] = self.fallbacks[k]
                    elif self.check_instrument_default_parameter(k):
                        self.sub_params[k] = (
                            self.check_instrument_default_parameter(k))

        #  set attributes for parameters
        for k, v in self.sub_params.items():
            setattr(self, k, v)

    def check_instrument_default_parameter(self, parameter_key):
        """
        Checks for instrument default values for a parameter and returns if
        found. Used to get values for parameters not explicitly stated in input
        file

        Args:
            parameter_key (str): parameter id

        Returns:
            parameter_value: value of parameter from instrument defaults
        """

        #  if no mass_spec or chromatography instruments have been specified
        #  in inputs, there will be no instrument defaults so return None
        if (not self.active_instruments["mass_spec"]
                and not self.active_instruments["chromatography"]):
            return None

        #  silico parameters have no directly instrument-depndent subparameters
        if self.params_class == "silico":
            return None

        #  identify which parameters potentially have set values in active
        #  instruments
        instrument_defaults = self.type_dict["instrument_dependent"]
        inst_only = instrument_defaults["instrument_only"]

        #  check whether parameter is dependent on instruments only and, if so,
        #  whether value for parameter is given in active instruments
        for device in inst_only:
            if parameter_key in self.active_instruments[device]:
                return self.active_instruments[device][parameter_key]

        #  check whether parameter is dependent on instruments AND polymer type;
        #  if so, check whether value for parameter is given in active
        # instruments
        for device in self.active_instruments:
            if "polymer_classes" in self.active_instruments[device]:
                polymer_inst_info = self.active_instruments[
                    device]["polymer_classes"]
            if self.params_dict["polymer_class"] in polymer_inst_info:
                polymer_specified = polymer_inst_info[
                    self.params_dict["polymer_class"]]
                if self.params_class in polymer_specified:
                    if parameter_key in polymer_specified[self.params_class]:
                        return polymer_specified[
                            self.params_class][parameter_key]

        return None

    def format_params_attr(self, parameter_key, parameter_value):
        """
        Takes an individual parameter and its associated value, checks whether
        the value type matches what is expected and, if not, tries to convert
        to appropriate type.

        Args:
            parameter_key (str): name of parameter
            parameter_value (): value of parameter as supplied from input file

        Raises:
            TypeDictKeyError: raised if an expected type for a parameter isn't
                defined in the type_dict

        Returns:
            parameter: returns parameter in appropriate type
        """
        #  get expected type for parameter value from type_dict
        try:
            expected_type = self.type_dict[parameter_key.lower()]
        except KeyError:
            raise TypeDictKeyError(
                key=parameter_key,
                type_dict=self.type_dict)

        if type(parameter_value) == expected_type or not parameter_value:
            return parameter_value

        #  list of types that can be checked directly
        check_types = [
            str, int, float, List[str], List[int], List[float], bool]

        #  if expected type is simple and can be checked directly, format and
        #  return
        if expected_type in check_types:
            return self.format_value(
                value=parameter_value,
                expected_type=expected_type)

        #  expand check types to simple dictionaries, in which value types are
        #  identical
        check_types.extend([
            Dict[str, int],
            Dict[str, float],
            Dict[str, str],
            Dict[str, List[int]],
            Dict[str, List[float]],
            Dict[str, List[str]]
        ])

        #  typecast all values in simple dict to appropriate type specified
        #  in type_dict, and return
        if expected_type in check_types:
            return {
                key: self.format_value(
                    value=value,
                    expected_type=expected_type.__args__[1])
                for key, value in parameter_value.items()
                if value}

        #  typecast all values in complicated dict to appropriate types
        #  specified in type_dict, and return
        return {
            key: self.format_value(
                value=value,
                expected_type=self.type_dict[parameter_key])
            for key, value in parameter_value.items()
            if value}

    def format_value(self, value, expected_type):
        """
        Checks the value of an individual parameter and tries to convert it to
        expected type for that value.

        Args:
            value (): value for parameter
            expected_type (): expected typing for value based on self.type_dict

        Returns:
            [type]: [description]
        """
        if expected_type in [str, int, float]:
            return expected_type(value)

        if expected_type in [List[str], List[int], List[float]]:
            subtype = expected_type.__args__[0]
            return [
                subtype(elm) for elm in value
                if elm is not None]
        return value

    def final_attrs_check(self, temp_props=TEMP_PROPS):
        """
        Performs final check and cleanup of class attributes. Removes temp
        props and checks that no parameters are missing or potentially
        incorrect.

        Args:
            temp_props (List[str], optional): list of names for attributes
                required to create instance of Parameters but which are then
                removed during cleanup to make Parameters easier to read and
                debug when used in later PolymerSoup workflows.
                Defaults to TEMP_PROPS.

        Raises:
            Exception: raises Exception if parameters are missing.
        """

        #  get list of temporary properties and type_dict keys that should
        #  be removed from instance of Parameters class
        temp_props.extend(NON_ATTR_KEYS)

        #  get list of parameters that are missing and raise Exception if
        #  there are any missing parameters
        missing_params = [
            param for param in self.type_dict
            if (param not in self.__dict__.keys()
                and param not in temp_props)]
        if missing_params:
            raise Exception(
                f'the following parameters are missing: {missing_params}')

        #  check adducts and fragments for silico parameters
        if self.params_class in ["silico_ms2", "silico_ms1"]:

            if self.adducts:
                self.adducts = retrieve_adduct_info(
                    mode=self.params_dict["mode"], adducts=self.adducts
                )
            if self.params_class == "silico_ms2":
                sanity_check_silico_fragmentation(
                    silico_params=self.sub_params,
                    polymer_class=self.params_dict["polymer_class"],
                    ms_level=2,
                    instrument_info=self.active_instruments["mass_spec"])

        #  remove any remaining temp_prop attributes. This is just for
        #  readability and easier debugging in later workflows
        for prop in temp_props:
            if prop in self.__dict__.keys():
                delattr(self, prop)
