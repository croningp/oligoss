from typing import List
from ..errors import InputTypeError, TypeDictKeyError

class Parameters():
    """
    Parameters object. To be passed on to other modules for PolymerSoup
    workflows.
    """

    #  properties needed to create instance of Parameters class but will then
    #  be removed during cleanup
    TEMP_PROPS = [
        "type_dict",
        "active_instruments",
        "params_dict",
        "instrument_dependencies"
    ]

    def __init__(
        self,
        params_dict,
        type_dict,
        instrument_dependencies,
        instrument_info
    ):

        self.type_dict = type_dict
        self.params_dict = params_dict
        self.instrument_dependencies = instrument_dependencies
        self.active_instruments = instrument_info

        #  read params dict and set parameter attributes
        self.set_parameter_attributes()

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
            InputTypeError: raised if the value supplied for a parameter does
                not match expected type is not convertable to expected type
        """

        #  iterate through parameters and associated values, check which
        #  parameters have default instrument values and retrieve if necessary
        #  finally, set each of these parameters as attributes in instance of
        #  Parameters class
        for k, v in self.params_dict.items():
            if not v:
                if k in self.instrument_dependencies["mass_spec"]:
                    v = self.active_instruments["mass_spec"][k]
                elif k in self.instrument_dependencies["chromatography"]:
                    v = self.active_instruments["chromatography"][k]
            if v:
                try:
                    v = self.format_params_attr(k, v)
                except TypeError:
                    raise InputTypeError(
                        key=k,
                        value=v,
                        expected_type=self.type_dict[k]
                    )
            setattr(self, k, v)

        #  remove attributes which were required to instantiate Parameters but
        #  are no longer needed in later workflows
        for prop in temp_props:
            delattr(self, prop)

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

        subtype = expected_type.__args__[1]
        return {
            key: self.format_value(value=value, expected_type=subtype)
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
