import codecs
import os

from setuptools import setup, find_packages

# https://packaging.python.org/single_source_version/
base_dir = os.path.abspath(os.path.dirname(__file__))
about = {}
with open(os.path.join(base_dir, "pytheas", "__about__.py"), "rb") as f:
    exec(f.read(), about)


def read(fname):
    return codecs.open(os.path.join(base_dir, fname), encoding="utf-8").read()


with open("requirements.txt") as f:
    required = f.read().splitlines()


setup(
    name="pytheas-pip",
    version=about["__version__"],
    author=about["__author__"],
    author_email=about["__author_email__"],
    packages=find_packages(),
    description=about["__description__"],
    long_description=read("README.rst"),
    long_description_content_type="text/x-rst",
    url=about["__website__"],
    project_urls={"Documentation": about["__website__"]},
    license=about["__license__"],
    platforms="any",
    include_package_data=True,
    install_requires=required,
    extras_require={
        "all": ["pipdate >=0.3.0, <0.4.0"],
        "update": ["pipdate >=0.3.0, <0.4.0"],
    },
    classifiers=[
        about["__status__"],
        about["__license__"],
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Mathematics",
    ],
)
