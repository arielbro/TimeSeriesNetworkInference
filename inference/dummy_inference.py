import attractor_learning.graphs


def dummy_inference_method(data_matrices, scaffold_network):
    """
    Dummy inference - returns the scaffold network given, with random Boolean functions
    :param data_matrices:
    :param scaffold_network:
    :return:
    """
    inferred_model = scaffold_network.copy()
    inferred_model.randomize_functions()
    return inferred_model
