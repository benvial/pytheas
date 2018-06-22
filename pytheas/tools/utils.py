


def normalize(x):
    return (x - x.min()) / (x.max() - x.min())


def between_range(x, xmin, xmax):
    return (xmax - xmin) * x + xmin
