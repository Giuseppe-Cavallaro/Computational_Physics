import sys
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.ticker as ticker

matplotlib.use('TkAgg')

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
            except Exception:
                try:
                    d[name] = float(elem)
                except Exception:
                    d[name] = elem
    return d

def printInfo(info):
    print("Coupling constant alpha_s\t = ", info["alpha_s"])
    print("Mass of top quark\t\t = ", info["top_quark_mass"], "GeV")
    print("Energy of single proton\t\t = ", info["energy"], "GeV")

filename = "output.dat"

info = read_info(filename)
printInfo(info)

data = pd.read_csv(filename, skiprows=3, skip_blank_lines=False)
data["theta"] = math.pi * data["theta"]

data["s"]     = data["x1"] * data["x2"] * (4 * info["energy"]**2)
data["beta"]  = (data["x1"] - data["x2"]) / (data["x1"] + data["x2"])
data["gamma"] = 1 / np.sqrt(1 - data["beta"]**2)
data["E_tot"] = data["gamma"] * info["energy"] * (data["x1"] + data["x2"] + data["beta"] * (data["x2"] - data["x1"]))
data["pt"]    = np.sqrt( (data["E_tot"]**2 - 4*info["top_quark_mass"]**2) / (4 * (1 + np.tan(data["theta"])**-2)) )

plt.figure(figsize=(5, 3))
ax = plt.axes()
ax.xaxis.set_major_locator(plt.MaxNLocator(10))
plt.yscale('log')
if len(sys.argv) == 1:
    plt.xlabel(r"s (GeV$^2$)")
    plt.ylabel("n")
    max_range = 1e7
    s  = plt.hist(data["s"].to_list(), range=(0,max_range), bins="auto")
    #s  = plt.hist(data["s"].to_list(), bins="auto")
    img_name = "img/s.jpg"
else:
    plt.xlabel(r"$k_T$ (GeV)")
    plt.ylabel("n")
    pt = plt.hist(data["pt"].to_list(), bins="auto")
    img_name = "img/kt.jpg"

plt.tight_layout()
plt.savefig(img_name, dpi=150)
plt.show()
