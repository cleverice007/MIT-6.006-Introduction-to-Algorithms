import unittest
from dnaseq import *

### Testing ###

class TestRollingHash(unittest.TestCase):
    def test_rolling(self):
        rh1 = RollingHash('CTAGC')
        rh2 = RollingHash('TAGCG')
        rh3 = RollingHash('AGCGT')
        rh1.slide('C', 'G')
        self.assertEqual(rh1.current_hash(), rh2.current_hash())
        rh1.slide('T', 'T')
        self.assertEqual(rh1.current_hash(), rh3.current_hash())

class TestMultidict(unittest.TestCase):
    def test_multi(self):
        foo = Multidict()
        foo.put(1, 'a')
        foo.put(2, 'b')
        foo.put(1, 'c')
        self.assertEqual(foo.get(1), ['a', 'c'])
        self.assertEqual(foo.get(2), ['b'])
        self.assertEqual(foo.get(3), [])

# This test case may break once you add the argument m (skipping).
class TestExactSubmatches(unittest.TestCase):
    def test_one(self):
        foo = 'yabcabcabcz'
        bar = 'xxabcxxxx'
        matches = list(getExactSubmatches(iter(foo), iter(bar), 3, 1))
        correct = [(1, 2), (4, 2), (7, 2)]
        self.assertEqual(len(matches), len(correct))
        for x in correct:
            self.assertIn(x, matches)

if __name__ == "__main__":
    unittest.main()
