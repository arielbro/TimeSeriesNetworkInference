from unittest import TestCase
from graphs import Network
import itertools
import random

class TestNetwork(TestCase):
    # def test_get_vertex(self):
    #     names = [''.join(random.choice(string.ascii_lowercase) for i in range(100))]
    #     G = Network()
    #     self.fail()
    #
    # def test_randomize_functions(self):
    #     self.fail()
    #
    # def test_generate_random(self):
    #     self.fail()

    def test_cnet_export_and_import(self):
        n = random.randint(1, 50)
        for restrictions in itertools.product([False, True], repeat=2):
            G = Network.generate_random(n_vertices=n, indegree_bounds=[1,20],
                                        restrict_and_or_gates=restrictions[0],
                                        restrict_signed_symmetric_threshold=restrictions[1])
            G.export_to_cnet("./temp_test_network.cnet")
            G_tag = Network.parse_cnet("./temp_test_network.cnet")
            self.assertEqual(G, G_tag)
