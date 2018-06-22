#!/usr/bin/env python


import platform
import os
import sys
import shutil
from urllib.request import urlopen
from tempfile import NamedTemporaryFile





_platform = platform.system()
install_name = "aotomat_install"
getdp_version = "2.11.3"
getdp_name = "getdp-" + getdp_version + "-"
getdp_url = "http://getdp.info/bin/"

def get_platform_suffix(_platform):
    if _platform == "Linux":
        name = "Linux64c"
        folder = "Linux"
        extension = ".tgz"
        zipformat = "gztar"
    elif _platform == "Windows":
        name = "Windows64c"
        folder = "Windows"
        extension = ".zip"
        zipformat = "zip"
    elif _platform == "Darwin":
        name = "MacOSXc"
        folder = "MacOSX"
        extension = ".tgz"
        zipformat = "gztar"
    else:
        name = ""
        folder = ""
        extension = ""
        zipformat = ""
        print("platform not supported")
    return name, folder, extension, zipformat




def download_and_uncompress(bin_dir, _platform, bin_url, bin_name):
    platform_suffix, folder, extension, zipformat = get_platform_suffix(_platform)
    zipurl = bin_url + "/" + folder + "/" + bin_name + platform_suffix + extension
    with urlopen(zipurl) as zipresp, NamedTemporaryFile() as tfile:
        tfile.write(zipresp.read())
        tfile.seek(0)
        shutil.unpack_archive(tfile.name, bin_dir, format=zipformat)


def is_tool(name):
    """Check whether `name` is on PATH and marked as executable."""

    # from whichcraft import which
    from shutil import which

    return which(name) is not None

#

def install_bin(name):
    print("Checking if {0} is installed".format(name))
    is_name = is_tool(name)
    is_name = False
    if is_name:
        print("{0} already installed".format(name))
    else:
        print("{0} cannot be found. Do you want to install it? (y/n)".format(name))
        while True:
            inpt = input()
            if inpt == "y":
                bin_dir = os.path.join(dir_name, "bin")
                os.mkdir(bin_dir)
                print("installing {0} in {1}".format(name, bin_dir))
                if name == "getdp":
                    _url, _name = getdp_url, getdp_name
                elif name == "gmsh":
                    _url, _name = gmsh_url, gmsh_name
                download_and_uncompress(bin_dir, _platform, _url, _name)
                break
            elif inpt == "n":
                print("Warning: you chose not to install {0}, some features will not be available".format(name))
                break
            else:
                print("Choose between yes (y) or no (n)")





if __name__ == "__main__":
    try:
        shutil.rmtree(install_name)
    except:
        pass

    print("Installing pytheas...")
    print("Platform: ", _platform)
    # Prompt user
    print("Please choose the installation directory: (default is current working directory)")
    # Get input
    while True:
        try:
            inpt = input()
            if inpt is "":
                print("Installing in current dir")
                inpt = os.getcwd()
            dir_name = os.path.join(inpt, install_name)
            os.mkdir(dir_name)
        except FileNotFoundError:
            print("This folder does not exist, please choose another one:")
        except FileExistsError:
            print("This folder already exists, please choose another one:")
        except PermissionError:
            print("You don't have permission to write in this directory, please choose another one:")
        else:
            break
    print("Installing in: ", dir_name)

    # install_bin("gmsh")
    install_bin("getdp")




    # is_gmsh = is_tool("gmsh")
