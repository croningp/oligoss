"""
This file is used to store info on common MS2-MSn fragmentation methods in
tandem MS. By adding default values for polymer class to each fragmentation
method, this will simplify run parameters.
"""

FRAGMENTATION_METHODS = {
    "CID": {
        "depsipeptide": {
            "fragment_series": ["b", "y"],
            "neutral_losses": ["water", "ammonia"],
            "core_fragment_series": ["b", "y"],
            "signature_ms_levels": [2],
            "signature_types": {
                "2": "Im"
            }
        }
    },
    "HCD": {
        "depsipeptide": {
            "fragment_series": ["a", "b", "y"],
            "neutral_losses": ["water", "ammonia"],
            "core_fragment_series": ["b", "y"],
            "signature_ms_levels": [2],
            "signature_types": {
                "2": "Im"
            },
        }
    },
    "ETD": {
        "depsipeptide": {
            "fragment_series": ["c", "x", "z"],
            "neutral_losses": ["water", "ammonia", "CO"],
            "core_fragment_series": ["c", "x"],
            "signature_ms_levels": [],
            "signature_types": {}
        }
    }
}
