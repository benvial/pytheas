# -*- coding: utf-8 -*-
# Author: Benjamin Vial
# License: MIT

"""
Get refractive index from a database
====================================

Import refractiveindex of a material at a given wavelength
from the refractiveindex.info_ database.
Forked from this repository_: github.com/cinek810/refractiveindex.info.

 .. _refractiveindex.info:
     https://refractiveindex.info/
 .. _repository:
     https://github.com/cinek810/refractiveindex.info
"""


from scipy.interpolate import interp1d
from parse import *
import yaml
from yaml.reader import Reader
import numpy as np
import os
import numpy as np
import functools

path = os.path.dirname(os.path.abspath(__file__))
database_path = os.path.join(path, "database", "data")


def get_directory_structure(rootdir):
    """
    Creates a nested dictionary that represents the folder structure of rootdir
    """
    dir = {}
    materials_path = []
    materials_list = []
    rootdir = rootdir.rstrip(os.sep)
    start = rootdir.rfind(os.sep) + 1
    for path, dirs, files in os.walk(rootdir):
        folders = path[start:].split(os.sep)
        _ = []
        for f in files:
            if f.endswith(".yml") and f != "library.yml":
                _.append(f[:-4])
                frel = os.path.join(os.path.relpath(path, database_path), f)
                # print(frel)
                materials_path.append(frel[:-4])
                materials_list.append(frel[:-4].split("/"))
        files = _
        subdir = dict.fromkeys(files)
        parent = functools.reduce(dict.get, folders[:-1], dir)
        parent[folders[-1]] = subdir
    return dir["data"], materials_path, materials_list


def fix_yml_file(yamlFile):
    if not yamlFile.endswith(".yml"):
        yamlFile += ".yml"
    return yamlFile


def strip_invalid(s):
    res = ""
    for x in s:
        if Reader.NON_PRINTABLE.match(x):
            # res += '\\x{:x}'.format(ord(x))
            continue
        res += x
    return res


def yaml_extract(yamlFile):
    yamlFile = fix_yml_file(yamlFile)
    filename = os.path.join(database_path, yamlFile)
    with open(filename) as yamlStream:
        c = yamlStream.read()
        allData = yaml.safe_load(strip_invalid(c))
    materialData = allData["DATA"][0]
    return allData, materialData


def getDataTab(yamlFile, lamb, datatype):
    _, materialData = yaml_extract(yamlFile)
    assert materialData["type"] == "tabulated {}".format(datatype)

    matLambda = []
    matN = []
    matK = []
    # in this type of material read data line by line
    for line in materialData["data"].split("\n"):

        try:
            if datatype == "n":
                parsed = parse("{l:g} {n:g}", line)
                n = parsed["n"]
                matN.append(n)
                matK.append(0)
            elif datatype == "k":
                parsed = parse("{l:g} {k:g}", line)
                k = parsed["k"]
                matN.append(0)
                matK.append(k)
            else:
                parsed = parse("{l:g} {n:g} {k:g}", line)
                n = parsed["n"]
                k = parsed["k"]
                matN.append(n)
                matK.append(k)
            matLambda.append(parsed["l"])
        except TypeError:
            pass

    matLambda = np.array(matLambda)
    matN = np.array(matN)
    matK = np.array(matK)
    if len(matLambda) == 1:
        return matN + 1j * matK
    else:
        interN = interp1d(matLambda, matN)
        interK = interp1d(matLambda, matK)
        return interN(lamb) + 1j * interK(lamb)


def getRangeTab(yamlFile, datatype):
    _, materialData = yaml_extract(yamlFile)
    assert materialData["type"] == "tabulated {}".format(datatype)
    # in this type of material read data line by line
    matLambda = []
    for line in materialData["data"].split("\n"):
        try:
            if datatype == "n":
                parsed = parse("{l:g} {n:g}", line)
            elif datatype == "k":
                parsed = parse("{l:g} {k:g}", line)
            else:
                parsed = parse("{l:g} {n:g} {k:g}", line)
            matLambda.append(parsed["l"])
        except TypeError:
            pass

    return (min(matLambda), max(matLambda))


def formula(lamb, coeff, formula_number):
    if formula_number == 1:
        epsi = 0
        for i in reversed(list(range(1, np.size(coeff), 2))):
            epsi += (coeff[i] * lamb ** 2) / (lamb ** 2 - coeff[i + 1] ** 2)
        epsi += coeff[0] + 1
        n = [np.sqrt(ep) for ep in epsi]
    elif formula_number == 2:
        epsi = 0
        for i in reversed(list(range(1, np.size(coeff), 2))):
            epsi += (coeff[i] * lamb ** 2) / (lamb ** 2 - coeff[i + 1])
        epsi += coeff[0] + 1
        n = [np.sqrt(ep) for ep in epsi]
    elif formula_number == 3:
        epsi = coeff[0]
        for i in range(1, np.size(coeff), 2):
            epsi += coeff[i] * lamb ** coeff[i + 1]
        n = [np.sqrt(ep) for ep in epsi]
    elif formula_number == 4:
        coeff_ = np.zeros(17)
        for i, val in enumerate(coeff):
            coeff_[i] = val
        coeff = coeff_
        epsi = coeff[0]
        epsi += coeff[1] * lamb ** coeff[2] / (lamb ** 2 - coeff[3] ** coeff[4])
        epsi += coeff[5] * lamb ** coeff[6] / (lamb ** 2 - coeff[7] ** coeff[8])
        epsi += coeff[9] * lamb ** coeff[10]
        epsi += coeff[11] * lamb ** coeff[12]
        epsi += coeff[13] * lamb ** coeff[14]
        epsi += coeff[15] * lamb ** coeff[16]
        n = [np.sqrt(ep) for ep in epsi]
    elif formula_number == 5:
        n = coeff[0]
        for i in reversed(list(range(1, np.size(coeff), 2))):
            n += coeff[i] * lamb ** coeff[i + 1]
    elif formula_number == 6:
        n = coeff[0] + 1
        for i in reversed(list(range(1, np.size(coeff), 2))):
            n += coeff[i] / (coeff[i + 1] - lamb ** (-2))
    elif formula_number == 7:
        n = coeff[0]
        n += coeff[1] / (lamb ** 2 - 0.028)
        n += coeff[2] / (lamb ** 2 - 0.028) ** 2
        for i in range(3, np.size(coeff)):
            n += coeff[i] * lamb ** (2 * (i - 2))
    elif formula_number == 8:
        A = coeff[0]
        A += coeff[1] * lamb ** 2 / (lamb ** 2 - coeff[2])
        A += coeff[3] * lamb ** 2
        n = ((1 + 2 * A) / (1 - A)) ** 0.5
    elif formula_number == 9:
        epsi = coeff[0]
        epsi += coeff[1] / (lamb ** 2 - coeff[2])
        epsi += coeff[3] * (lamb - coeff[4]) / ((lamb - coeff[4]) ** 2 * +coeff[5])
        n = epsi ** 0.5
    return n


def getDataF(yamlFile, lamb, formula_number):
    _, materialData = yaml_extract(yamlFile)

    assert materialData["type"] == "formula {}".format(formula_number)

    dataRange = np.array(list(map(float, materialData["range"].split())))
    coeff = np.array(list(map(float, materialData["coefficients"].split())))

    if not check_bounds(lamb, dataRange):
        raise Exception(
            "OutOfBands",
            "No data for this material(" + yamlFile + " ) for lambda=" + str(lamb),
        )
    else:
        return formula(lamb, coeff, formula_number)


def check_bounds(lamb, dataRange):
    return min(lamb) >= dataRange[0] or max(lamb) <= dataRange[1]


# def Error(BaseException):
#     pass


# this is general function to check data type, and run appropriate actions


def getData(yamlFile, lamb):
    _, materialData = yaml_extract(yamlFile)
    mtype = materialData["type"]
    if mtype.split()[0] == "tabulated":
        return getDataTab(yamlFile, lamb, mtype.split()[1])
    elif mtype.split()[0] == "formula":
        return getDataF(yamlFile, lamb, int(mtype.split()[1]))
    else:
        return np.zeros_like(lamb)


def get_wl_range(yamlFile):
    _, materialData = yaml_extract(yamlFile)

    mtype = materialData["type"]

    if mtype.split()[0] == "tabulated":
        return getRangeTab(yamlFile, mtype.split()[1])
    else:
        return np.array(list(map(float, materialData["range"].split())))


def get_complex_index(lambdas, yamlFile):
    lambdas = np.array([lambdas]).ravel()
    ncomplex = getData(yamlFile, lambdas)
    return np.asarray(np.conj(ncomplex))


class Materials(object):
    """Materials class"""

    def __init__(self):
        self.data, self.materials_path, self.materials_list = get_directory_structure(
            database_path
        )

    def list(self, sublist=None):

        if sublist:
            if not isinstance(sublist, list):
                sublist = list([sublist])
            a = self.data
            for s in sublist:
                a = a[s]
            return list(a.keys())
        else:
            return list(self.data.keys())

    def get(self, id):
        a = self.data
        for s in id:
            a = a[s]
        return a

    def get_rel_path(self, id):
        return os.path.join(*id)

    def info(self, id):
        return yaml_extract(self.get_rel_path(id))[0]["REFERENCES"]

    def get_complex_index(self, lambdas, id):
        return get_complex_index(lambdas, self.get_rel_path(id))

    def get_wl_range(self, id):
        return get_wl_range(self.get_rel_path(id))
