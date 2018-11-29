import datetime

def normalize(x):
    return (x - x.min()) / (x.max() - x.min())


def between_range(x, xmin, xmax):
    return (xmax - xmin) * x + xmin


def generate_ID(self):
    return datetime.now().strftime("%Y_%m_%d_%H_%M_%s_%f")
