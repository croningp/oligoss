"""
This file contains functions essential for logging in Polymersoup
"""
import logging
import subprocess

def generate_logger_config_dict(log_file):
    """
    Takes a log file filepath and returns a configuration dict that can be
    used to configure root logger.

    Args:
        log_file (str): full filepath to log file

    Returns:
        Dict: config file dict
    """
    return dict(
        version=1,
        formatters={
            "f": {
                "format": "%(asctime)s %(message)s %(funcName)s %(module)s"
            }
        },
        handlers={
            "h": {
                "class": "logging.FileHandler",
                "filename": log_file,
                "level": logging.INFO,
                "formatter": "f"
            }
        },
        root={
            "handlers": ["h"],
            "level": logging.INFO
        }
    )

def get_git_info_for_log():
    """
    Checks git user, current branch and version tag; returns message with
    this info for logging

    Returns:
        str: git info message for logging.
    """

    #  get git username
    user = str(subprocess.check_output(
        ["git", "config", "user.name"]
    ), "utf-8").rstrip("\n")

    #  get git branch
    branch = str(subprocess.check_output(
        ["git", "rev-parse", "--abbrev-ref", "HEAD"]
    ), "utf-8").rstrip("\n")

    #  get git version tag
    version = str(subprocess.check_output(
        ["git", "describe"]
    ), "utf-8").rstrip("\n")

    return f"ran on branch {branch} by {user}. Version = {version}"
