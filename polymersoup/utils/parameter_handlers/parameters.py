from typing import List, Dict
from .instrument_handlers import check_instrument_info
from ..errors import InputTypeError, TypeDictKeyError, MissingParameterValue
from ..type_dicts.parameter_type_dicts import (
    CORE_PARAM_TYPES,
    SILICO_LIBRARY_TYPES,
    EXTRACTOR_TYPES,
    POSTPROCESS_TYPES
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

    def get_param_associations(
        self,
        type_associations=PARAM_TYPE_ASSOCIATIONS,
        fallbacks=PARAM_FALLBACK_ASSOCIATIONS
    ):
        """
        Retrieve type dict and fallbacks for paramerer group

        Args:
            type_associations (Dict[str, dict], optional): dict of parameter
                group ids and associated type_dicts.
                Defaults to PARAM_TYPE_ASSOCIATIONS.

        Raises:
            Exception: raised if self.params_class not in
                PARAM_TYPE_ASSOCIATIONS

        Returns:
            dict, dict: type_dict, fallbacks for parameter group
        """
        if self.params_class not in type_associations:
            raise Exception(
                f'{self.params_class} not a valid parameter group\n'
                f'please choose from: {type_associations}')

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

        type_dict = type_associations[self.params_class]
        param_fallbacks = fallbacks[self.params_class]

        return type_dict, param_fallbacks, params

    def set_parameter_attributes(self, temp_props=TEMP_PROPS):
        """
        Set attributes of Parameters depending on what set of parameter inputs
        the Parameters class is to be an instance of (silico, postprocess, core
        etc...)

        Args:
            temp_props (List[str], optional): list of names for attributes
                required to create instance of Parameters but which are then
                removed during cleanup to make Parameters easier to read and
                debug when used in later PolymerSoup workflows.
                Defaults to TEMP_PROPS.

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
            if not v:
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
            try:
                if v:
                    if self.params_class != "core":
                        v = self.format_params_attr(k, v)
            except TypeError:
                raise InputTypeError(
                    key=k,
                    value=v,
                    expected_type=self.type_dict[k]
                )
            if k != "core":
                setattr(self, k, v)

        #  remove attributes which were required to instantiate Parameters but
        #  are no longer needed in later workflows
        for prop in temp_props:
            delattr(self, prop)

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
            if self.params_dict["polymer_class"] in self.active_instruments[
                    device]:
                polymer_specified = self.active_instruments[self.params_dict[
                    "polymer_class"]]
                if parameter_key in polymer_specified:
                    return polymer_specified[parameter_key]
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
            raise Exception(self.params_class, parameter_key, self.type_dict)
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
        if expected_type in [str, int, float]:
            return expected_type(value)

        if expected_type in [List[str], List[int], List[float]]:
            subtype = expected_type.__args__[0]
            return [
                subtype(elm) for elm in value
                if elm]
        return value
