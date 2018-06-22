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
import sys
from parse import *
import yaml
import numpy as np
import cmath
import os

path = os.path.dirname(os.path.abspath(__file__))
database_path = os.path.join(path, "database", "data")


def yaml_extract(yamlFile):
    filename = os.path.join(database_path, yamlFile)
    yamlStream = open(filename, 'r')
    allData = yaml.load(yamlStream)
    materialData = allData["DATA"][0]
    return yamlStream, allData, materialData

# This function is designed to return refractive index for specified lambda
# for file in "tabluated n " format
def getDataN(yamlFile, lamb):
    yamlStream, allData, materialData = yaml_extract(yamlFile)
    assert materialData["type"] == "tabulated n"

    matLambda = []
    matN = []
    # in this type of material read data line by line
    for line in materialData["data"].split('\n'):
        parsed = parse("{l:g} {n:g}", line)
        try:
            n = parsed["n"]
            matLambda.append(parsed["l"])
            matN.append(n)
        except TypeError as e:
            pass
            #sys.stderr.write("TypeError occured:"+str(e)+"\n")

    matLambda = np.array(matLambda)
    matN = np.array(matN)
    matK = np.zeros_like(matN)

    if len(matLambda)==1:
        return matN + 1j * matK
    else:
        interN = interp1d(matLambda, matN)
        interK = interp1d(matLambda, matK)
        return interN(lamb) + 1j * interK(lamb)


def getDataK(yamlFile, lamb):
    yamlStream, allData, materialData = yaml_extract(yamlFile)
    assert materialData["type"] == "tabulated k"
    print("Only k values!")

    matLambda = []
    matK = []
    # in this type of material read data line by line
    for line in materialData["data"].split('\n'):
        parsed = parse("{l:g} {k:g}", line)
        try:
            k = parsed["k"]
            matLambda.append(parsed["l"])
            matK.append(k)
        except TypeError as e:
            pass
            #sys.stderr.write("TypeError occured:"+str(e)+"\n")

    matLambda = np.array(matLambda)
    matK = np.array(matK)
    matN = np.zeros_like(matK)

    if len(matLambda)==1:
        return matN + 1j * matK
    else:
        interN = interp1d(matLambda, matN)
        interK = interp1d(matLambda, matK)
        return interN(lamb) + 1j * interK(lamb)

def getRangeN(yamlFile):
    yamlStream, allData, materialData = yaml_extract(yamlFile)
    assert materialData["type"] == "tabulated n"
    # in this type of material read data line by line
    matLambda = []
    for line in materialData["data"].split('\n'):
        parsed = parse("{l:g} {n:g}", line)
        try:
            matLambda.append(parsed["l"])
        except TypeError as e:
            pass
            #sys.stderr.write("TypeError occured:"+str(e)+"\n")
    return (min(matLambda), max(matLambda))


def getRangeK(yamlFile):
    yamlStream, allData, materialData = yaml_extract(yamlFile)
    assert materialData["type"] == "tabulated k"
    # in this type of material read data line by line
    matLambda = []
    for line in materialData["data"].split('\n'):
        parsed = parse("{l:g} {k:g}", line)
        try:
            matLambda.append(parsed["l"])
        except TypeError as e:
            pass
            #sys.stderr.write("TypeError occured:"+str(e)+"\n")
    return (min(matLambda), max(matLambda))

def getDataNK(yamlFile, lamb):
    """This function is designed to return refractive index for specified lambda
    for file in `tabulated nk` format"""
    yamlStream, allData, materialData = yaml_extract(yamlFile)
    assert materialData["type"] == "tabulated nk"
    matLambda = []
    matN = []
    matK = []
    # in this type of material read data line by line
    for line in materialData["data"].split('\n'):
        parsed = parse("{l:g} {n:g} {k:g}", line)
        try:
            n = parsed["n"]
            k = parsed["k"]
            matLambda.append(parsed["l"])
            matN.append(n)
            matK.append(k)
        except TypeError as e:
            pass
            #sys.stderr.write("TypeError occured:"+str(e)+"\n")

    matLambda = np.array(matLambda)
    matN = np.array(matN)
    matK = np.array(matK)

    if len(matLambda)==1:
        return matN + 1j * matK
    else:
        interN = interp1d(matLambda, matN)
        interK = interp1d(matLambda, matK)
        return interN(lamb) + 1j * interK(lamb)



def getRangeNK(yamlFile):
    yamlStream, allData, materialData = yaml_extract(yamlFile)

    assert materialData["type"] == "tabulated nk"
    # in this type of material read data line by line
    matLambda = []
    for line in materialData["data"].split('\n'):
        parsed = parse("{l:g} {n:g} {k:g}", line)
        try:
            matLambda.append(parsed["l"])
        except TypeError as e:
            pass
            #sys.stderr.write("TypeError occured:"+str(e)+"\n")
    return (min(matLambda), max(matLambda))


# this function is desined to get data from files in formula1 format
def getDataF1(yamlFile, lamb):
    yamlStream, allData, materialData = yaml_extract(yamlFile)

    assert materialData["type"] == "formula 1"

    dataRange = np.array(list(map(float, materialData["range"].split())))
    coeff = np.array(list(map(float, materialData["coefficients"].split())))

    epsi = 0
    if min(lamb) >= dataRange[0] or max(lamb) <= dataRange[1]:
        for i in reversed(list(range(1, np.size(coeff), 2))):
            epsi = epsi + ((coeff[i] * lamb**2) / (lamb**2 - coeff[i + 1]**2))

    else:
        raise Exception("OutOfBands", "No data for this material for this l")

    epsi = epsi + coeff[0] + 1
    n = []
    for ep in epsi:
        n.append(cmath.sqrt(ep))
    return n


# this function is desined to get data from files in formula2 format
def getDataF2(yamlFile, lamb):
    yamlStream, allData, materialData = yaml_extract(yamlFile)

    assert materialData["type"] == "formula 2"

    dataRange = np.array(list(map(float, materialData["range"].split())))
    coeff = np.array(list(map(float, materialData["coefficients"].split())))

    epsi = 0
    if min(lamb) >= dataRange[0] or max(lamb) <= dataRange[1]:
        for i in reversed(list(range(1, np.size(coeff), 2))):
            epsi = epsi + ((coeff[i] * lamb**2) / (lamb**2 - coeff[i + 1]))

    else:
        raise Exception("OutOfBands", "No data for this material for this l")

    epsi = epsi + coeff[0] + 1
    n = []
    for ep in epsi:
        n.append(cmath.sqrt(ep))
    return n


# this function is desined to get data from files in formula3 format
def getDataF3(yamlFile, lamb):
    yamlStream, allData, materialData = yaml_extract(yamlFile)

    assert materialData["type"] == "formula 3"

    dataRange = np.array(list(map(float, materialData["range"].split())))
    coeff = np.array(list(map(float, materialData["coefficients"].split())))

    epsi = coeff[0]
    if min(lamb) >= dataRange[0] and max(lamb) <= dataRange[1]:
        for i in range(1, np.size(coeff), 2):
            epsi = epsi + coeff[i] * lamb**coeff[i + 1]
    else:
        raise Exception(
            "OutOfBands")

    n = []
    for ep in epsi:
        n.append(cmath.sqrt(ep))
    return n


def getDataF4(yamlFile, lamb):
    yamlStream, allData, materialData = yaml_extract(yamlFile)

    assert materialData["type"] == "formula 4"

    dataRange = np.array(list(map(float, materialData["range"].split())))
    coeff = np.zeros(17)
    coeff2 = list(map(float, materialData["coefficients"].split()))

    for i in range(0, len(coeff2)):
        coeff[i] = coeff2[i]

   # print(coeff)

    epsi = coeff[0]
    if min(lamb) >= dataRange[0] and max(lamb) <= dataRange[1]:
        epsi = epsi + coeff[1] * lamb**coeff[2] / (lamb**2 - coeff[3]**coeff[4])
        epsi = epsi + coeff[5] * lamb**coeff[6] / (lamb**2 - coeff[7]**coeff[8])
        epsi = epsi + coeff[9] * lamb**coeff[10]
        epsi = epsi + coeff[11] * lamb**coeff[12]
        epsi = epsi + coeff[13] * lamb**coeff[14]
        epsi = epsi + coeff[15] * lamb**coeff[16]
    else:
        raise Exception("OutOfBands", "No data for this material(" +
                        yamlFile + " )for lambda=" + str(lamb))

    n = []
    for ep in epsi:
        n.append(cmath.sqrt(ep))
    return n


# this function is desined to get data from files in formula5 format
def getDataF5(yamlFile, lamb):
    yamlStream, allData, materialData = yaml_extract(yamlFile)

    assert materialData["type"] == "formula 5"

    dataRange = np.array(list(map(float, materialData["range"].split())))
    coeff = np.array(list(map(float, materialData["coefficients"].split())))

    n = coeff[0]
    if min(lamb) >= dataRange[0] or max(lamb) <= dataRange[1]:
        for i in reversed(list(range(1, np.size(coeff), 2))):
            n = n + (coeff[i] * lamb**coeff[i + 1])

    else:
        raise Exception("OutOfBands", "No data for this material for this l")

    return n


# this function is desined to get data from files in formula6 format
def getDataF6(yamlFile, lamb):
    yamlStream, allData, materialData = yaml_extract(yamlFile)

    assert materialData["type"] == "formula 6"

    dataRange = np.array(list(map(float, materialData["range"].split())))
    coeff = np.array(list(map(float, materialData["coefficients"].split())))

    n = coeff[0] + 1
    if min(lamb) >= dataRange[0] or max(lamb) <= dataRange[1]:
        for i in reversed(list(range(1, np.size(coeff), 2))):
            n = n + coeff[i] / (coeff[i + 1] - lamb**(-2))

    else:
        raise Exception("OutOfBands", "No data for this material for this l")

    return n


# this function is desined to get data from files in formula7 format
def getDataF7(yamlFile, lamb):
    yamlStream, allData, materialData = yaml_extract(yamlFile)

    assert materialData["type"] == "formula 7"

    dataRange = np.array(list(map(float, materialData["range"].split())))
    coeff = np.array(list(map(float, materialData["coefficients"].split())))

    n = coeff[0]
    if min(lamb) >= dataRange[0] or max(lamb) <= dataRange[1]:
        n +=  coeff[1] / (lamb**2 - 0.028)
        n +=  coeff[2] / (lamb**2 - 0.028)**2
        for i in range(3, np.size(coeff)):
            n +=  coeff[i] * lamb**(2*(i-2))

    else:
        raise Exception("OutOfBands", "No data for this material for this l")

    return n

# this function is desined to get data from files in formula8 format
def getDataF8(yamlFile, lamb):
    yamlStream, allData, materialData = yaml_extract(yamlFile)

    assert materialData["type"] == "formula 8"

    dataRange = np.array(list(map(float, materialData["range"].split())))
    coeff = np.array(list(map(float, materialData["coefficients"].split())))

    A = coeff[0]
    if min(lamb) >= dataRange[0] or max(lamb) <= dataRange[1]:
        A += coeff[1] * lamb**2 / (lamb**2 - coeff[2])
        A += coeff[3] * lamb**2

    else:
        raise Exception("OutOfBands", "No data for this material for this l")

    n = ((1 + 2*A) / (1 - A))**0.5
    return n

# this function is desined to get data from files in formula9 format
def getDataF9(yamlFile, lamb):
    yamlStream, allData, materialData = yaml_extract(yamlFile)

    assert materialData["type"] == "formula 9"

    dataRange = np.array(list(map(float, materialData["range"].split())))
    coeff = np.array(list(map(float, materialData["coefficients"].split())))

    epsi = coeff[0]
    if min(lamb) >= dataRange[0] or max(lamb) <= dataRange[1]:
        epsi += coeff[1] / (lamb**2 - coeff[2])
        epsi += coeff[3] * (lamb - coeff[4]) / ((lamb - coeff[4])**2 * + coeff[5])

    else:
        raise Exception("OutOfBands", "No data for this material for this l")

    n = epsi**0.5
    return n


def Error(BaseException):
    pass

# this is general function to check data type, and run appropriate actions


def getData(yamlFile, lamb):
    yamlStream, allData, materialData = yaml_extract(yamlFile)

    if materialData["type"] == "tabulated nk":
        return getDataNK(yamlFile, lamb)
    elif materialData["type"] == "tabulated n":
        return getDataN(yamlFile, lamb)
    elif materialData["type"] == "tabulated k":
        return getDataK(yamlFile, lamb)
    elif materialData["type"] == "formula 1":
        return getDataF1(yamlFile, lamb)
    elif materialData["type"] == "formula 2":
        return getDataF2(yamlFile, lamb)
    elif materialData["type"] == "formula 3":
        return getDataF3(yamlFile, lamb)
    elif materialData["type"] == "formula 4":
        return getDataF4(yamlFile, lamb)
    elif materialData["type"] == "formula 5":
        return getDataF5(yamlFile, lamb)
    elif materialData["type"] == "formula 6":
        return getDataF6(yamlFile, lamb)
    elif materialData["type"] == "formula 7":
        return getDataF7(yamlFile, lamb)
    elif materialData["type"] == "formula 8":
        return getDataF8(yamlFile, lamb)
    elif materialData["type"] == "formula 9":
        return getDataF9(yamlFile, lamb)
    else:
        #	raise Error("UnsupportedDataType:This data type is currnetly not supported");
        return np.zeros_like(lamb)


def getRange(yamlFile):
    yamlStream, allData, materialData = yaml_extract(yamlFile)


    if materialData["type"] == "tabulated nk":
        return getRangeNK(yamlFile)
    elif materialData["type"] == "tabulated n":
        return getRangeN(yamlFile)
    elif materialData["type"] == "tabulated k":
        return getRangeK(yamlFile)
    else:
        return np.array(list(map(float, materialData["range"].split())))
    # elif materialData["type"] == "formula 1":
    #     return np.array(list(map(float, materialData["range"].split())))
    # elif materialData["type"] == "formula 2":
    #     return np.array(list(map(float, materialData["range"].split())))
    # elif materialData["type"] == "formula 3" or materialData["type"] == "formula 4":
    #     return np.array(list(map(float, materialData["range"].split())))
    # else:
    #     #	raise Error("UnsupportedDataType:This data type "+materialData["type"] +" is currnetly not supported");
    #     return (0, 0)


def get_complex_index(lambdas, yamlFile):
    lambdas = np.array([lambdas]).ravel()
    ncomplex = getData(yamlFile, lambdas)
    return np.asarray(np.conj(ncomplex))
