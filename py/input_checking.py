"""
This file contains several standard functions that will be used throughout the program. Each function checks input
and issues an error if there is something wrong.
"""

import pathlib
import sys


def required_field(variable_to_test: any, err_string: str) -> None:
    """
    If required field variable_to_test is empty, issues an error. Otherwise this does nothing

    :param variable_to_test: Any input type
    :param err_string: A string with the error message
    :return: None
    """
    if variable_to_test is None:
        print('\n' + err_string + '\n')
        sys.exit(1)


def check_file_open(filename: str, err_string: str, required: bool = False) -> None:
    """
    Checks that the filename is not empty and that it is indeed a  file

    :param filename: file name, string
    :param err_string: string of the error if it is not a file
    :param required: If not required, skips the check
    :return: None
    """
    if required or filename is not None:
        if filename is None:
            print('\n' + err_string + '\n')
            sys.exit(1)
        else:
            try:
                pathlib.Path(filename).resolve(strict=True)
            except FileNotFoundError:
                print('\n' + err_string + '\n')
                sys.exit(1)


def check_dir(directory: str, err_string: str) -> None:
    """
    Checks that directory exists and is a directory
    :param directory: string of the directory path
    :param err_string: string of the error in case it is not a directory or doesn't exist
    :return: None
    """
    if not pathlib.Path(directory).is_dir():
        print('\n' + err_string + '\n')
        raise NotADirectoryError


def is_in_range(value: float, lower_bound: float, upper_bound: float, err_string: str) -> None:
    """
    Checks that value is between the lower bound and upper bound, and if not prints an error message
    (err_string) and exits the program.

    :param value: float for the value
    :param lower_bound: float for the upper bound
    :param upper_bound: float for the lower bound
    :param err_string: string of the error message to print if the value is out of range
    :return: None
    """
    if value < lower_bound or value > upper_bound:
        print('\n' + err_string + '\n')
        sys.exit(1)
