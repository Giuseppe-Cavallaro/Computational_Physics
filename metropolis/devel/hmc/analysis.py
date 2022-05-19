import sys
import os
import matplotlib
import matplotlib.pyplot as plt
from itertools import cycle
import math
import pandas as pd

matplotlib.use('TkAgg')

TAB_LENGTH = 10
FOLDER = "../data/"
FIG_WIDTH = 5
FIG_HEIGHT = 3 

SAVE = True
IMG_FOLDER = "../img/"
IMG_FORMAT = ".jpg"
IMG_QUALITY = 95

def get_label(info):
    DStr = r"$D = {}$".format(info["D"])
    lambdaStr = r"$\lambda = {}$".format(info["lambda"])
    tauzeroStr = r"$\tau_0 = {}$".format(info["tlength"])
    return f"{DStr}, {lambdaStr}, {tauzeroStr}"

def print_info(info):
    print("PARAMETERS")
    for k, v in info.items():
        print(f"{k}\t{v}".expandtabs(TAB_LENGTH))
    print("END PARAMETERS")

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

# read csv file and split it into two dataframe where an empty line is found
def read_csv_split(path):
    data = pd.read_csv(path, skiprows=3, skip_blank_lines=False)
    is_NaN = data.isnull().all(axis=1)
    i = data.index[is_NaN][0]
    data1 = data[:i].reset_index(drop=True)
    data2 = data[i+2:].reset_index(drop=True)
    data2 = data2.rename(columns=data.iloc[i+1]).dropna(how="all", axis=1)
    return data1, data2

# returns dataframe filled with jackknife variables
# calculated starting from input dataframe
def get_jackknife(df, index=["L", "k"]):
    df2 = df.set_index(index)
    counts = df.groupby(index).count()
    sums   = df.groupby(index).sum()
    return ((sums - df2) / (counts - 1)).reset_index()

# OBSERVABLES FUNCTIONS
def magnetization(D, x):
    return x["m"] / x["L"]**D

def absoluteMagnetization(D, x):
    return x["mAbs"] / (x["L"]**D)

def squareMagnetization(D, x):
    return x["m2"] / (x["L"]**D)

def squareMagnetizationDivV(D, x):
    return x["m2"] / (x["L"]**D)**2

def susceptibility(D, x):
    return (x["m2"] - (x["mAbs"]**2)) / (x["L"]**D)

def binder(D, x):
    return x["m4"] / (x["m2"]**2)


#---------------------------------------------------------------------------------------


input_arg = {"magn"   : {"fun"  : magnetization,
                         "col"  : "magn",
                         "name" : r"$m$",
                         "fullname" : "Magnetization"},
             "abs"    : {"fun"  : absoluteMagnetization,
                         "col"  : "abs_magn",
                         "name" : r"$|m|$",
                         "fullname" : "Absolute Magnetization"},
             "sqr"    : {"fun"  : squareMagnetization,
                         "col"  : "square_magn",
                         "name" : r"$m^2$",
                         "fullname" : "Square Magnetization"},
             "sqr_V"  : {"fun"  : squareMagnetizationDivV,      # sqr magn divided by V
                         "col"  : "square_magn_V",
                         "name" : r"$m^2/V$",
                         "fullname" : "Square Magnetization divided by V"},
             "susc"   : {"fun"  : susceptibility,
                         "col"  : "susc",
                         "name" : r"$\chi$",
                         "fullname" : "Susceptibility"},
             "binder" : {"fun"  : binder,
                         "col"  : "binder",
                         "name" : r"$B$",
                         "fullname" : "Binder Cumulant"},
             "pa"     : {"col"  : "pa",
                         "name" : r"$P_a$",
                         "fullname" : "Acceptance Probability"}}

argv = sys.argv

if len(sys.argv) < 3:
    print("Wrong arguments. Usage:")
    print(f"python3 {argv[0]} <input file> <observable to plot> <image to save>(optional) -L <values of L to plot>(optional)")
    print("List of possible arguments: ")
    for k, v in input_arg.items():
        print("{}\t{}".format(k, v["fullname"]))
    exit(1)

img_name = ""
listL = []
if len(argv) > 3:
    if argv[3] != "-L":
        img_name = argv[3]
if len(argv) > 4 or (len(argv) > 3 and img_name == ""):
    try:
        i = argv.index("-L")
    except:
        print("Error in reading the flag -L")
        exit(1)
    tempL = argv[i+1:]
    if len(tempL) == 0:
        print("No values for L provided")
        exit(1)
    try:
        listL = [int(elem) for elem in tempL]
    except:
        print("Error converting values of L in integer")
        exit(1)

filename = argv[1]
obs_to_plot = argv[2]
path = FOLDER + filename

info = read_info(path) # read first two lines of file
index = ["L", "k"]
data, data_pa = read_csv_split(path)
data = data.astype(float)
data["L"] = data["L"].astype(int)
data_pa = data_pa.astype(float)
data_pa["L"] = data_pa["L"].astype(int)

if len(listL) == 0:
    listL = data["L"].drop_duplicates().to_list()
else:
    if not any(x in data["L"].drop_duplicates().to_list() for x in listL):
        print("No values for L provided are in the database")
        print("Allowed values for L are ")
        print(*data["L"].drop_duplicates().to_list(), sep=" ")
        exit(1)

# calculate jackknife variables
jk = get_jackknife(data)

# calculate the value of observables with jackknife variable
fjk = jk[index].copy() # copy the index in a new empty dataframe
for n, f in input_arg.items():
    if n != "pa":
        fjk[f["col"]] = jk.apply(lambda x : f["fun"](info["D"], x), axis=1)

# calculate the value of the functions over jackknife averages
mean_jk = jk.groupby(index, as_index=False).mean()
obs = mean_jk[index].copy() # copy the index in a new empty dataframe
for n, f in input_arg.items():
    if n != "pa":
        obs[f["col"]] = mean_jk.apply(lambda x : f["fun"](info["D"], x), axis=1)

# calculate the error
sqr    = (obs.set_index(index) - fjk.set_index(index))**2
sqrSum = sqr.groupby(index).sum()
factor = (fjk.groupby(index).count() - 1) / fjk.groupby(index).count() # (N-1)/N
sqrErr = factor * sqrSum
obsErr = (sqrErr**(1.0/2)).reset_index()

print_info(info)

# plot options
plt.figure(figsize=(FIG_WIDTH, FIG_HEIGHT))
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = cycle(prop_cycle.by_key()['color'])

# add "pa" column to jk dataframe
df_to_plot  = pd.concat([obs.set_index(index), data_pa.set_index(index)], axis=1).reset_index().groupby("L", as_index=False)

# plot dataframe
plt.title(get_label(info))
plt.xlabel(r"$\kappa$")
plt.ylabel(input_arg[obs_to_plot]["name"])
for l, df in df_to_plot:
    if l in listL:
        color = next(colors)
        x = df["k"].to_list()
        y = df[input_arg[obs_to_plot]["col"]].to_list()
        if obs_to_plot != "pa":
            yerr = obsErr[obsErr["L"] == l][input_arg[obs_to_plot]["col"]].to_list()
            plt.errorbar(x, y, yerr=yerr, color=color, ecolor=color, marker="o", markersize=4,
                            elinewidth=1, capsize=2, label=f"L = {l}")
        else:
            plt.plot(x, y, color=color, marker="o", markersize=4, label=f"L = {l}")

plt.legend()
if img_name != "":
    img = IMG_FOLDER + img_name + IMG_FORMAT
    plt.savefig(img, quality=IMG_QUALITY)
    print("Image saved as " + img)
else:
    plt.show()
