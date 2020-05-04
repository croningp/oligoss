"""
This file contains constants and functions relevant to handling parameters
from polymer class-specific config files

"""
import os

HERE = os.path.dirname(__file__)

POLYMER_CONFIGS = os.path.join(
    os.path.dirname(HERE), "polymer_configs")

POLYMER_ALIASES = {
    f.replace(".json", ""): os.path.join(POLYMER_CONFIGS, f)
    for f in os.listdir(POLYMER_CONFIGS)
}
