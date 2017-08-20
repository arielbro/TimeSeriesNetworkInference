def list_repr(elements):
    if len(elements) == 0:
        return "[]"
    else:
        return "[" + reduce(lambda x, y: str(x) + ", " + str(y), elements) + "]"
