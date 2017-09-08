def list_repr(elements):
    if len(elements) == 0:
        return "[]"
    else:
        return "[" + str(reduce(lambda x, y: str(x) + ", " + str(y), elements)) + "]"