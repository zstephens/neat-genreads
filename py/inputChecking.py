"""
This file contains several standard functions that will be used throughout the program. Each function checks input
and issues an error if there is something wrong.
"""

import os
import sys


def requiredField(s: any, errString: str) -> None:
    """
    If required field s is empty, issues an error. Otherwise this does nothing

    :param s: Any input type
    :param errString: A string with the error message
    :return: None
    """
    if s is None:
        print('\n' + errString + '\n')
        exit(1)


def checkFileOpen(fn: str, errString: str, required: bool = False) -> None:
    """
    Checks that the filename is not empty and that it is indeed a  file
    :param fn: file name, string
    :param errString: string of the error if it is not a file
    :param required: If not required, skips the check
    :return: None
    """
    if required or fn is not None:
        if fn is None:
            print('\n' + errString + '\n')
            exit(1)
        else:
            try:
                open(fn, 'r')
            except ValueError:
                print('\n' + errString + '\n')
                exit(1)


def checkDir(dir: str, errString: str) -> None:
    """
    Checks that directory exists and is a directory
    :param dir: string of the directory path
    :param errString: string of the error in case it is not a directory or doesn't exist
    :return: None
    """
    if not os.path.isdir(dir):
        print('\n' + errString + '\n')
        exit(1)


def isInRange(val: float, lb: float, ub: float, errString: str) -> None:
    """
    Checks that value (val) is between the lower bound (lb) and upper bound (ub), and if not prints an error message
    (errString) and exits the program.
    :param val: float for the value
    :param lb: float for the upper bound
    :param ub: float for the lower bound
    :param errString: string of the error message to print if the value is out of range
    :return: None
    """
    if val < lb or val > ub:
        print('\n' + errString + '\n')
        exit(1)
