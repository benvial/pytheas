#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Benjamin Vial
# License: MIT


import json

class LessVars():

    def __init__(self, lessfile):
        self.lessfile = lessfile
        self.variables = self.getvars()


    def getvars(self):
        with open(self.lessfile) as f:
            lines = f.readlines()  # list containing lines of file
            variables = {}
            i = 1
            for line in lines:
                line = line.strip()  # remove leading/trailing white spaces
                if line:
                    data = [item.strip() for item in line.split(":")]
                    if line.split(" ")[0] != "@import":
                        if line.split(" ")[0] != "@font-face":
                            if line.split(" ")[0][0] =="@":
                                # name = data[0].split
                                key = data[0].split("@")[1]
                                value = data[-1].split(";")[0]
                                variables[key] = value
        return variables


    def fmtvar(self, key):
        v = self.variables[key]
        v_ = v.split("(")[-1].split(")")[0].split(",")
        return tuple([int(_)/256 for _ in v_])


    def latex_code(self, latex_var, key):
        vfmt = self.fmtvar(key)
        return sphinxcol(latex_var, vfmt)

def sphinxcol(latex_var, vfmt):
    return r"{}".format(latex_var) + "={rgb}" + "{{{},{},{}}}".format(*vfmt)


if __name__ == "__main__":

    lessfile = "_custom/static/css/less/variables.less"
    less = LessVars(lessfile)
    # pretty printing dictionarie
    print(json.dumps(less.variables, indent=4))
    key = "brand-primary"
    vfmt = less.fmtvar(key)
    print("formatted variable: {} = {}".format(key, vfmt))

    val = less.latex_code("TitleColor", key)

    print(val)
