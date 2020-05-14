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

class InvalidMSFragmentation(Exception):
    def __init__(
        self,
        frag_method,
        instrument=None,
        fragment_series=None,
        neutral_invalid=False
    ):
        self.frag_method = frag_method
        self.instrument = instrument
        self.fragment_series = fragment_series
        self.neutral_invalid = neutral_invalid

    def __str__(self):
        error_str = f'{self.frag_method} is not valid'

        if self.instrument:
            error_str += f'using {self.instrument}\n'
        if self.fragment_series:
            error_str += f'with {self.fragment_series} fragments\n'
        if self.neutral_invalid:
            neutral_str = 'neutral loss fragments'
            if self.fragment_series:
                neutral_str = 'and' + neutral_str
            error_str += neutral_str

        return error_str

class InvalidMonomerId(Exception):
    def __init__(
        self,
        monomer_id,
        polymer_class=None,
        available_monomers=None
    ):
        self.monomer = monomer_id
        self.polymer_alias = polymer_class
        self.available_monomers = available_monomers

    def __str__(self):
        err_str = (
            f'{self.monomer} is invalid ')
        if self.polymer_alias:
            err_str += f'for {self.polymer_alias} polymers.\n'
        if self.available_monomers:
            if type(self.available_monomers) == dict:
                self.available_monomers = self.available_monomers.keys()
            err_str += (
                'Here are the valid monomer one letter codes for\n'
                f'{self.polymer_alias}: {self.available_monomers}')
        return err_str

class InvalidModificationTarget(Exception):
    def __init__(
        self,
        modification,
        target
    ):

        self.modification = modification
        self.target = target

    def __str__(self):
        if str(self.target) in ["0", "-1"]:
            return f'{self.modification} cannot target terminus {self.target}'
        return f'{self.modification} cannot target {self.target} sidechain'
