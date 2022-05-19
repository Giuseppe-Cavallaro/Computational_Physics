import sys
import os
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import math
import pandas as pd
import numpy as np
from subprocess import call

matplotlib.use('TkAgg')

TAB_LENGTH  = 11
CHECK1_EXEC = "./check1"
FIG_WIDTH   = 7
FIG_HEIGHT  = 4.5
FOLDER = "../data/"

IMG_FOLDER = "../img/"
SAVE = False
IMG_NAME = "check_img"

files = {"2" : "check2.csv", "3" : "check3.csv", "4" : "check4.csv"}

def read_info(path):
    d = {}
    with open(path) as myFile:
        line1 = next(myFile)
        l1 = line1.strip().split(",")
        line2 = next(myFile)
        l2 = line2.strip().split(",")
        for name, elem in zip(l1, l2):
            try:
                d[name] = int(elem)
            except:
                try:
                    d[name] = float(elem)
                except:
                    d[name] = elem
    return d

def get_arg_list(argv, accepted_input):
    arg_list = argv[1:]
    control = all(elem in accepted_input for elem in arg_list) and (0 < len(arg_list) <= len(accepted_input))
    if control == False:
        print("Incorrect number of arguments.")
        print(f"Usage: {sys.argv[0]} 1 | 2 | 3 | 4 | -s")
        print("1  for check1, it has no plot (runs ./check1)")
        print("2  for check2.")
        print("3  for check3.")
        print("4  for check4.")
        print("-s for saving picture")
        print("The order of numbers determines the position on the picture.")
        exit(1)
    return list(dict.fromkeys(arg_list))

def print_info(filename, info):
    print(f"DATA FILE: {filename}")
    for k, v in info.items():
        print(f"{k}\t{v}".expandtabs(TAB_LENGTH))

def read_csv_check3(path):
    data = pd.read_csv(path, skiprows=3, comment="#")
    header = data.columns.values.tolist()
    i = data.index[(data[header[0]] == header[0]) & (data[header[1]] == header[1]) & (data[header[2]] == header[2])].to_list()[0]
    return data[:i].astype(float), data[i+1:].astype(float).reset_index(drop=True)

def get_value_err(data, input_index, c):
    index = input_index
    if not isinstance(input_index, list):
        index = [input_index]
    value = data.groupby(index)[c].mean()
    var   = data.groupby(index)[c].var()
    N     = data.groupby(index)[c].count()
    err   = (var / N)**(1.0/2)
    return value.reset_index(), err.reset_index()

argv = sys.argv
if "-s" in argv:
    index = argv.index("-s")
    argv.remove("-s")
    SAVE = True
    try:
        IMG_NAME = argv.pop(index)
    except IndexError:
        print(f"Saving image as \"{IMG_NAME}\"")
        print()

arg_list = get_arg_list(argv, ["1"] + list(files.keys()))

nPlots = len(arg_list)
if "1" in arg_list:
    nPlots -= 1
if nPlots == 3:
    fig = plt.figure(figsize=(FIG_WIDTH, FIG_HEIGHT))
    gs = gridspec.GridSpec(2, 2)
    axs_list = []
    axs_list.append(fig.add_subplot(gs[0, 0]))
    axs_list.append(fig.add_subplot(gs[0, 1:]))
    axs_list.append(fig.add_subplot(gs[1, :]))
    enum_axs = enumerate(axs_list)
else:
    fig, axs = plt.subplots(nPlots, figsize=(FIG_WIDTH, FIG_HEIGHT))
    enum_axs = enumerate(fig.axes)
check1 = False

for arg in arg_list:
    if arg == "1":
        check1 = True
    else:
        path = FOLDER + files[arg]
        info = read_info(path)
        print(f"CHECK {arg}\n")
        print_info(files[arg], info)
        print("--------------------")

    if arg == "2":
        data  = pd.read_csv(path, skiprows=3, comment="#")
        value, err = get_value_err(data, "dtau2", "deltaH")
        i, ax = next(enum_axs)
        ax.set_ylabel(r"$\langle | \Delta H | \rangle$")
        ax.set_xlabel(r"$\delta\tau^2$")

        dtau2 = value["dtau2"].tolist()
        deltaH = value["deltaH"].tolist()
        err_deltaH = err["deltaH"].tolist()

        ax.errorbar(dtau2, deltaH, yerr=err_deltaH, markersize=1, color="black", ecolor='black',
                    elinewidth=1, capsize=2, linewidth=1)

    elif arg == "3":
        data_without,  data_with   = read_csv_check3(path)
        value_without, err_without = get_value_err(data_without, "dtau2", "m2")
        value_with,    err_with    = get_value_err(data_with, "dtau2", "m2")
        i, ax = next(enum_axs)
        ax.set_ylabel(r"$\langle m^2 \rangle$")
        ax.set_xlabel(r"$\delta\tau^2$")

        dtau2 = value_without["dtau2"].tolist()
        m2_without = value_without["m2"].tolist()
        m2_with = value_with["m2"].tolist()
        m2_err_without = err_without["m2"].tolist()
        m2_err_with = err_with["m2"].tolist()

        ax.errorbar(dtau2, m2_without, yerr=m2_err_without, marker="o", markersize=1, color="black", ecolor='black',
                    elinewidth=1, capsize=2, label="senza acc&rej", linewidth=1)
        ax.errorbar(dtau2, m2_with, yerr=m2_err_with, marker="o", markersize=1, color="red", ecolor='red',
                    elinewidth=1, capsize=2,  label="con acc&rej", linewidth=1)


        coef = np.polyfit(dtau2, m2_without, 1)
        x_vals = np.array(dtau2)
        y_vals = coef[1] + coef[0] * x_vals
        ax.plot(x_vals, y_vals, '--', color="dimgrey", linewidth=1, label="{:.2f}{} + {:.2f}".format(coef[0], r"$x$", coef[1]))

        coef = np.polyfit(dtau2, m2_with, 1)
        y_vals = coef[1] + coef[0] * x_vals
        ax.plot(x_vals, y_vals, '--', color="indianred", linewidth=1, label="{:.2f}{} + {:.2f}".format(coef[0], r"$x$", coef[1]))
        ax.legend()
        #ax.legend(loc=3, ncol=4)

    elif arg == "4":
        data = pd.read_csv(path, skiprows=3, comment="#")
        i, ax = next(enum_axs)
        ax.set_xlabel(r"$n$")
        ax.set_ylabel(r"$\langle e^{-\Delta H} \rangle$")

        ax.axhline(1.0, color="red")
        ax.plot(data["exp_deltaH"].tolist(), marker="o", markersize=2, color="black", linewidth=1)

if check1:
    print("CHECK 1\n")
    call([CHECK1_EXEC])
if any(a in list(files.keys()) for a in arg_list):
    plt.tight_layout()
    if SAVE:
        plt.savefig(IMG_FOLDER + IMG_NAME + ".jpg", dpi = 150)
    plt.show()
