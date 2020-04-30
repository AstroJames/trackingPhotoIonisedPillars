"""

    Title:          Pillar Tracking at a Photoionised Mixing Layer
    Notes:          the data handling functions, mostly various pickle functions.
    Author:         James Beattie & contributions from Shyam Menon
    First Created:  16 / Jan / 2020

"""

from header import *

############################################################################################################################################

def saveObj(obj, name):
    """
    Save a pickle object.

    INPUTS:
    ----------
    obj      - the name of the data object to be saved
    name     - the name of the pickle object

    """

    os.system("touch " + name + ".pkl")
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def loadObj(name):
    """
    Load a pickle object.

    INPUTS:
    ----------
    name     - the name of the pickle object

    """

    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)
