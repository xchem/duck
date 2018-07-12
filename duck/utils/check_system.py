import numpy as np


def check_pot_eng(pot_eng):
    """
    Logic for dealing with the potential energy
    :param pot_eng:
    :return:
    """
    pot_eng = np.array(pot_eng)
    # Standard deviation
    standard_dev = pot_eng.std()
    # Maximum
    max_eng = max(pot_eng)
    # Minimum
    min_eng = min(pot_eng)
    # The max variance
    max_diff = max(abs(np.diff(pot_eng)))
    # Anthony Bradley basically made up these numbers
    # Please update
    if max_eng > -9876.0:
        print("VERY HIGH POTENTIAL ENERGY")
        return False
    if standard_dev > 5000:
        print("HIGHLY VOLATILE")
        return False
    if max_diff > 20000:
        print("HIGHLY VOLATILE")
        return False
    return True

def check_if_equlibrated(input_file,column,header=True,delimiter=","):
    """

    :param input_file:
    :param column:
    :param header:
    :param delimiter:
    :return:
    """
    input_lines = open(input_file).readlines()
    if header:
        data = input_lines[1:]
    else:
        data = input_lines
    pot_eng = [float(x.split(delimiter)[column]) for x in data]
    return check_pot_eng(pot_eng)