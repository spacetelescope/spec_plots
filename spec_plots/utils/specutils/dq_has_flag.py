"""
.. module:: dq_has_flag
   :synopsis: Utility function that returns True if a given Data Quality bitmask
       contains a specific flag.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

#--------------------
# Built-In Imports
#--------------------
from __future__ import absolute_import
from builtins import str
#--------------------
# Package Imports
#--------------------
from spec_plots import __version__

#--------------------

def dq_has_flag(dqf, flag_to_check):
    """
    Returns true/false if the given DQ value has a specific flag value set after
    unpacked into a 16-bit string.  For example:

    .. code-block:: python

         print dq_has_flag(48,16)
         True
         print dq_has_flag(40, 16)
         False

    :param dqf: The DQ value.

    :type dqf: int

    :param flag_to_check: The flag value to check if it's set to True.

    :type flag_to_check: int

    :returns: bool -- Returns True if `flag_to_check` is `True` inside `dqf`.

    :raises: ValueError
    """

    # Make sure `flag_to_check` is a power of 2.
    if (flag_to_check & (flag_to_check-1)) == 0 and flag_to_check != 0:
        dq_16bit_str = "{0:b}".format(dqf)
        flag_16bit_str = "{0:b}".format(flag_to_check)

        # If the 16-bit string of the value to check is longer than 16-bit
        # string version of the DQ value, then we know it can't be part of the
        # DQ  bitmask.  If not, then look for that bit to be set (by counting
        # from the right).
        return (len(flag_16bit_str) <= len(dq_16bit_str) and
                dq_16bit_str[-1 * len(flag_16bit_str)] == '1')
    else:
        raise ValueError("Flag to check must be a power of 2.  Asked to check"
                         " whether flag " + str(flag_to_check) + " is set to"
                         " True in bitmask value " + str(dqf) + ", but " +
                         str(flag_to_check) + " is not a power of 2.")
#--------------------
