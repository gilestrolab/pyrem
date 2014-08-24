__author__ = 'quentin'
from datetime import timedelta
import re

__REGEX__ = re.compile(r'((?P<hours>\d+?)h)?((?P<minutes>\d+?)m)?((?P<seconds>\d+\.?\d*?)s)?')


def str_to_time(str):
    """
    Parse a string describing a duration as a :class:`~datetime.timedelta`

    :param str: a string with the format "XhYmZs". where XYZ are integers (or floats)
    :type str: str
    :return: a timedelta corresponding to the `str` argument
    :rtype: timedelta
    """
    parts = __REGEX__.match(str)

    if not parts or len(parts.group()) <1:
        raise ValueError("The string `%s' does not contain time information.\nThe syntax is, for instance, 1h99m2s" %(str))
    parts = parts.groupdict()
    time_params = {}
    for (name, param) in parts.iteritems():
        if param:
            time_params[name] = float(param)
    return timedelta(**time_params)
