from nose.tools import *
import os
import cockatoo
from cockatoo.screen import Screen,Cocktail,Compound

class TestUtil:
    def setup(self):
        self.path = os.path.dirname(os.path.abspath(__file__))
        self.csv_test_screen = "%s/../screens/csv/test-screens/salt-concentration.csv" % self.path
        self.ph_screen = "%s/../screens/json/test-screens/ph.json" % self.path
        self.salt_screen = "%s/../screens/json/test-screens/salt-concentration.json" % self.path
        self.conc_screen = "%s/../screens/json/test-screens/conc.json" % self.path
        self.peg_screen = "%s/../screens/json/test-screens/peg-mw.json" % self.path
        self.cation_screen = "%s/../screens/json/test-screens/cation.json" % self.path
        self.anion_screen = "%s/../screens/json/test-screens/anion.json" % self.path
        self.hwi_gen8 = "%s/../screens/json/hwi/hwi-gen8.json" % self.path

    def teardown(self):
        pass

    def test_basic(self):
        cp1 = Compound('sodium chloride', 1.0, 'M')
        c1 = Cocktail('c1')
        c1.add_compound(cp1)
        s = Screen('screen1')
        s.add_cocktail(c1)
        s.print_stats()

        assert len(s) == 1

    def test_parse_csv(self):
        s = cockatoo.screen.parse_csv('salt-con', self.csv_test_screen)
        s.print_stats()
        assert len(s) == 12

    def test_parse_json(self):
        s = cockatoo.screen.parse_json(self.salt_screen)
        s.print_stats()
        assert len(s) == 12
        s = cockatoo.screen.parse_json(self.hwi_gen8)
        s.print_stats()
        assert len(s) == 1536

    def test_c6(self):
        s = cockatoo.screen.parse_json(self.salt_screen)
        print "\t".join(['i','j','c6'])
        i = 0
        for j in xrange(0, len(s)):
            c6 = cockatoo.c6.distance(s.cocktails[i], s.cocktails[j])
            print "\t".join([str(i),str(j),str(c6)])

