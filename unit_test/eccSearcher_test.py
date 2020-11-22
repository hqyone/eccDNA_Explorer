'''
Author: Quanyuan(Leo) He
Email: hqyone@gmail.com
Insititute: Hunan Normal Univeristy
Date: 2020-11-21 22:37:51
LastEditTime: 2020-11-21 23:30:58
LastEditors: Quanyuan(Leo) He
Description: 
FilePath: /eccDNA_Explorer/unit_test/eccSearcher_test.py
License: The MIT License (MIT)
'''

import unittest
import eccSearcher as es


class TestEccDNASearcherMethods(unittest.TestCase):
    def test_upper(self):
        self.assertTrue("foo".upper(), 'FOO')

    def test_getClosestPair(self):
        self.assertTrue(es.getClosestPair([], []), (None, None, None))
        self.assertTrue(es.getClosestPair([1, 0], []), (None, None, None))
        self.assertTrue(es.getClosestPair([0, 1], [1]), (0, 1, 1))
        self.assertTrue(es.getClosestPair([0], [1]), (1, 0, 1))
