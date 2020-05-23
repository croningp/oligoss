class RipperDict:

    def __init__(self, ripper_data):

        # make sure ms1 and ms2 spectra are ordered by retention time
        self.ms1 = {key: ripper_data["ms1"][key] for key in sorted(
            ripper_data["ms1"], key=lambda x: float(
                ripper_data["ms1"][x]["retention_time"]), reverse=False)}
        self.ms2 = {key: ripper_data["ms2"][key] for key in sorted(
            ripper_data["ms2"], key=lambda x: float(
                ripper_data["ms2"][x]["retention_time"]), reverse=False)}
        self.spectra = {}
