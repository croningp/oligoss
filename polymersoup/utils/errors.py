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
                f'{self.input_parameter}. Expected type:\n'
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
            f'there is no key {self.key} in the type_dict: {self.type_dict}')
