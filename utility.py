list_repr = lambda elements: "[" + reduce(lambda x, y: str(x) + ", " + str(y), elements) + "]"
