"""
This file contains functions essential for logging in Polymersoup
"""
import logging
import subprocess

#  list of atypical logging levels not commonly used in Oligomersoup. This is
#  for filtering out common log messages for LevelFilter (see below)
NON_STD_LOG_LEVELS = [10, 30, 40, 50]

class LevelFilter(logging.Filter):
    """
    Subclass of logging.Filter. Used for filtering log records by level - this
    is to avoid error, warning and critical log messages being logged to same
    handler as INFO-level messages. Makes debugging easier and any issues easier
    to spot.

    """
    def __init__(self, levels=NON_STD_LOG_LEVELS):
        self.levels = levels

    def filter(self, record):
        """
        Returns True if log record is in levels specified.

        Returns:
            bool: True if record.level in self.levels, otherwise False
        """
        return record.levelno in NON_STD_LOG_LEVELS

def generate_logger_config_dict(log_file, log_level=20):
    """
    Takes a log file filepath and returns a configuration dict that can be
    used to configure root logger.

    Args:
        log_file (str): full filepath to log file

    Returns:
        Dict: config file dict
    """
    return {
        "version": 1,
        "formatters": {
            "f": {
                "format": "%(asctime)s %(message)s %(funcName)s %(module)s"
            }
        },
        "filters": {
            "log_specific_levels": {
                "()": LevelFilter
            }
        },
        "handlers": {
            "h": {
                "class": "logging.FileHandler",
                "filename": log_file,
                "level": log_level,
                "formatter": "f"
            },
            "e": {
                "class": "logging.FileHandler",
                "filename": log_file.replace(".log", "non_std.log"),
                "level": 10,
                "formatter": "f",
                "filters": ["log_specific_levels"]
            }
        },
        "root": {
            "handlers": ["h", "e"],
            "level": 10
        }
    }

def get_git_info_for_log():
    """
    Checks git user, current branch and version tag; returns message with
    this info for logging

    Returns:
        str: git info message for logging.
    """
    git_logging = False

    latest_commit = "unknown"
    try:
        #  get git username
        user = str(subprocess.run(
            ["git", "config", "user.name"], capture_output=True
        ).stdout, "utf-8")
        git_logging = True
    except subprocess.CalledProcessError:
        user = "unknown"

    try:
        #  get git branch
        branch = str(subprocess.run(
            ["git", "rev-parse", "--abbrev-ref", "HEAD"], capture_output=True
        ).stdout, "utf-8")
        git_logging = True
    except subprocess.CalledProcessError:
        branch = "unknown"

    if branch != "unknown":
        try:
            latest_commit = str(subprocess.run(
                ["git", "rev-parse", branch], capture_output=True
            ).stdout)
        except subprocess.CalledProcessError:
            try:
                latest_commit = str(subprocess.run(
                    ["git", "rev-parse", "HEAD"], capture_output=True
                ).stdout, "utf-8")
            except subprocess.CalledProcessError:
                latest_commit = "unknown"
    try:
        #  get git version tag
        version = str(subprocess.run(
            ["git", "describe"], capture_output=True
        ).stdout, "utf-8")
        git_logging = True
    except subprocess.CalledProcessError:
        version = "unknown"

    if git_logging:
        return (
            f"Ran on branch {branch} by {user}. Latest commit on local branch "
            f" = {latest_commit}. Version = {version}")
    return "No git information available for run."
