from unittest import TestCase

import datamc_nuscuts

class TestJoke(TestCase):
    def test_is_string(self):
        s = datamc_nuscuts.joke()
        self.assertTrue(isinstance(s, basestring))
