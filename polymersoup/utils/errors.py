class InputTypeError(Exception):
    """
    Error for forbidden type being supplied for a paramter in input file

    Args:
        Exception (): [description]
    """
    def __init__(
        self,
        key,
        value,
        expected_type
    ):
        self.input_parameter = key
        self.param_value = value
        self.expected_type = expected_type

    def __str__(self):
        return (f'invalid_value({self.param_value}) given for\n'
                f'"{self.input_parameter}". Expected type:\n'
                f'{self.expected_type}')

class TypeDictKeyError(Exception):

    def __init__(
        self,
        key,
        type_dict
    ):
        self.type_dict = type_dict
        self.key = key

    def __str__(self):
        return (
            f'there is no key "{self.key}" in the type_dict: {self.type_dict}')

class MissingParameterValue(Exception):

    def __init__(
        self,
        key,
        type_dict,
        param_class,
        instrument_fallback=False
    ):
        self.key = key
        self.type_dict = type_dict
        self.instr_check = instrument_fallback
        self.param_class = param_class

    def __str__(self):
        required_params = [
            key for key in self.type_dict
            if (key not in ["optional", "instrument_dependent"]
                and key not in self.type_dict["optional"])]
        if not self.instr_check:
            return(
                f'value for "{self.key}" is missing. This must be supplied\n'
                f'Essential parameters which must be supplied:\n'
                f'{required_params}. Please check {self.param_class} parameters'
            )
        return(
            f'Value for "{self.key}" parameter must be supplied in the input\n'
            'file or by default instrument values. Please check instrument\n'
            f'defaults and {self.param_class} parameters'
        )
