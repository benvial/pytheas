from setuptools import setup, find_packages

setup(
    name="pytheas",
    version="0.1",
    description="Python Electromagnetic Analysis and Simulation with the \
                 Finite Element Method",
    url="https://benvial.github.com/pytheas",
    author="Benjamin Vial",
    author_email="benjamin.vial84@gmail.com",
    license="MIT",
    packages=find_packages(),
    install_requires=["matplotlib", "numpy", "scipy", "parse", "pyyaml",
                      "scoop", "seaborn", "termcolor",
                      "PyQt5"],
)
