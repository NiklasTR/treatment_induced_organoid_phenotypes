import numpy as np
import os


def trim_func(a, func, percent, axis):
    """
    Apply a function to a trimmed matrix along an axis. The function must be
    of a type that takes the 'axis' keyword, e.g. np.mean, np.std, np.median,
    etc.
    :param func:
    :param a:
    :param percent:
    :param axis:
    :return:
    """
    asort = np.sort(a, axis=axis)
    num_cases = asort.shape[1]
    num_to_trim = int(np.ceil(percent*num_cases))
    atrim = np.delete(
        asort, slice(num_cases - num_to_trim, num_cases),
        axis=axis)
    atrim = np.delete(
        atrim, slice(0, num_to_trim),
        axis=axis)
    if atrim.shape[1] == 0:
        return np.zeros(shape=a.shape[0])
    else:
        return func(atrim, axis=axis)


def get_cur_dir():
    """
    Gets the current directory of the calling script
    :return:
    """

    return os.path.realpath(os.path.join(
        os.getcwd(), os.path.dirname(__file__)))

