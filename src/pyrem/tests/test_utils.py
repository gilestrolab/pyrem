__author__ = 'quentin'

import unittest
from pyrem import utils
from datetime import timedelta

class TestUtils(unittest.TestCase):

    def test_parse_string_idx(self):

        ans = utils.str_to_time("1m2.300s")
        ref = timedelta(seconds=60.0+2.3)
        self.assertEqual(ans, ref)
