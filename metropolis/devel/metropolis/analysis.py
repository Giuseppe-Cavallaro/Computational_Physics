import sys
import os
import matplotlib
import matplotlib.pyplot as plt
import math

matplotlib.use('TkAgg')

def calculateAverage(l):
    avg = 0
    for elem in l:
        avg += elem
    return avg / len(l)

def squareSum(l):
    avg = calculateAverage(l)
    sqrSum = 0
    for elem in l:
        sqrSum += (elem - avg)**2
    return sqrSum

# return the square root of the variance
def calculateError(l):
    N = len(l)
    sqrSum = squareSum(l)
    return math.sqrt( sqrSum / (N*(N-1)) )

def jkCalc(fun, d):
    N = len(d["n"])
    avgArg = []
    for elem in d["jk"]:
        if elem != None:
            avgArg.append(calculateAverage(elem))
        else:
            avgArg.append(None)
    avgFun = fun(avgArg, d)

    sqrSum = 0
    for i in range(N):
        arg = []
        for elem in d["jk"]:
            if elem != None:
                arg.append(elem[i])
            else:
                arg.append(None)
        sqrSum += (fun(arg, d) - avgFun)**2
    errSqr = ((N-1)/N) * sqrSum
    return avgFun, math.sqrt( errSqr )

def jackknife(name, d):
    N = len(d["n"])
    s = d["sum"][d["names"].index(name)]
    o = d["obs"][d["names"].index(name)]
    l = []
    for i in range(N):
        l.append( (s - o[i]) / (N - 1) )
    return l

def gamma(fun, t, d, average = None):
    N = len(d["n"])
    if average == None:
        avg, = jkCalc(fun, d)
    else:
        avg = average
    tot = 0
    for i in range(N - t):
        temp_i = []
        temp_it = []
        for elem in d["obs"]:
            temp_i.append(elem[i])
            temp_it.append(elem[i+t])
        f_i = fun(temp_i, d)
        f_it = fun(temp_it, d)
        tot += (f_i - avg) * (f_it - avg)
    return tot / (N - t)

def printParameters(d):
    print("Parameters:")
    print(" - N       = {}".format(len(d["n"])))
    print(" - V       = {}".format(d["V"]))
    print(" - BinX    = {}".format(d["binx"]))
    print(" - Binsize = {}".format(d["binsize"]))
    print(" - Naccu   = {}".format(d["naccu"]))
    print(" - Ntherm  = {}".format(d["Ntherm"]))
    print(" - Nsweep  = {}".format(d["Nsweep"]))
    print(" - Δ       = {}".format(d["delta"]))
    print(" - κ       = {}".format(d["kappa"]))
    print(" - λ       = {}".format(d["lambda"]))
    print(" - Pa      = {}".format(d["pa"]))

# OBSERVABLES
def magn(l, d):
    V = d["V"]
    i = d["names"].index("M")
    return l[i] / V

def magnAbs(l, d):
    V = d["V"]
    i = d["names"].index("absM")
    return l[i] / V

def magnSquare(l, d):
    V = d["V"]
    i = d["names"].index("M2")
    return l[i] / V

def susc(l, d):
    V = d["V"]
    iM2   = d["names"].index("M2")
    iabsM = d["names"].index("absM")
    return (l[iM2] - l[iabsM]**2) / V

def binder(l, d):
    V = d["V"]
    iM2 = d["names"].index("M2")
    iM4 = d["names"].index("M4")
    return l[iM4] / l[iM2]**2

#END OBSERVABLES

F = [magn, magnAbs, magnSquare, susc, binder]
Farg = [["M"], ["absM"], ["M2"], ["absM", "M2"], ["M2", "M4"]]
Fnames = ["m", "|m|", "m²", "χ", "B"]

folder = "../data/"
files  = {}

# read input parameters
if len(sys.argv) > 1:
    argv = sys.argv

    saveGammaImg = True
    calculateGamma = True
    try:
        argv.remove("-g")
    except Exception:
        calculateGamma = False
    try:
        calculateGamma = True
        index = argv.index("-gs")
        argv.remove("-gs")
        ImgName = argv.pop(index)
    except Exception:
        calculateGamma = False
        saveGammaImg = False

    try:
        ind = argv.index("-p")
    except ValueError:
        print("No argument for plot given")
        exit(1)

    obsList = []
    if len(argv[ind:]) == 1:
        for i in range(len(F)):
            obsList.append(i)
    else:
        obsList = [int(elem) for elem in argv[(ind+1):]]

    tempFileList = argv[1:ind]
    fileList = []
    binsize = []
    for i in range(len(tempFileList)):
        try:
            int(tempFileList[i])
        except:
            fileList.append(tempFileList[i])
            try:
                binsize.append( int(tempFileList[i+1]) )
            except:
                binsize.append(1)

else:
    print("First argument must be the file name.")
    exit(1)

# read input files
for f, b in zip(fileList, binsize):
    path = folder + f
    firstLine = True
    secondLine = True
    with open(path) as myFile:
        c = 0
        d = {}
        d["binsize"] = b
        d["n"] = []
        d["t"] = []
        d["gamma"] = [[] for i in range(len(F))]
        for line in myFile:
            if firstLine:
                l = line.strip().split()
                d["V"]      = int(l[0])
                d["naccu"]  = int(l[1])
                d["binx"]   = int(l[2])
                d["Ntherm"] = int(l[3])
                d["Nsweep"] = int(l[4])
                d["delta"]  = float(l[5])
                d["kappa"]  = float(l[6])
                d["lambda"] = float(l[7])
                firstLine = False
            else:
                l = line.strip().split()
                if secondLine:
                    d["obs"] = [[] for i in range(len(l) - 1)]
                    d["sum"] = [0 for i in range(len(l) - 1)]
                    d["names"] = []
                    for name in l[1:]:
                        d["names"].append(name)
                    temp = [0 for i in range(len(l) - 1)]
                    secondLine = False
                else:
                    # save acceptance rate
                    if l[0] == "pa":
                        d["pa"] = float(l[1])
                    else:
                        c += 1
                        print("\rReading file {}... {:.2f}%".format(f, c / (d["Nsweep"]/d["naccu"]) * 100), end="")
                        n = int(l[0])
                        for j in range(len(l) - 1):
                            temp[j] += float(l[j+1])
                        if n % b == 0:
                            d["n"].append( n / b )
                            for j in range(len(l) - 1):
                                d["obs"][j].append( temp[j] / b )
                                d["sum"][j] += temp[j] / b
                            temp = [0 for i in range(len(l) - 1)]

        print()
        files[f] = d

# define list of t to use for calculating gamma
dt = 10
for f in files:
    N = 300
    t = 0
    while t <= N:
        files[f]["t"].append(t)
        t += dt
print()

#calculate jackknife observables
for f in files:
    files[f]["jk"] = []
    for n in files[f]["names"]:
        arg = False
        for a in obsList:
            if a >= 0 and n in Farg[a]:
                arg = True
                break
        if arg:
            print("Calculating jackknife for " + n + " in file " + f)
            files[f]["jk"].append( jackknife(n, files[f]) )
        else:
            files[f]["jk"].append( None )

plt.figure(figsize=(5, 3))
plt.xlabel("t")
plt.ylabel(r"$\Gamma(t)$")
# print observables averages and calculate gamma
for f in files:
    binsizeStr = files[f]["binsize"]
    plt.title("Dimensione bin = " + str(files[f]["binsize"] * files[f]["naccu"]))
    print("FILE: " + f)
    printParameters(files[f])
    print("Observables averages:")
    for a in obsList:
        if a >= 0:
            avg, err = jkCalc(F[a], files[f])
            print(" - {}\t= {:.7f}, err = {:.7f}".format(Fnames[a], avg, err))
    if calculateGamma:
        N = len(files[f]["t"])
        for a in obsList:
            c = 0
            print()
            for t in files[f]["t"]:
                c += 1
                print("\rCalculating gamma function for {} - {:.2f}%".format(Fnames[a], c / N * 100), end = "")
                files[f]["gamma"][a].append( gamma(F[a], t, files[f], avg) )
            plt.plot(files[f]["t"], files[f]["gamma"][a], label = r"$\Delta = $" + str(files[f]["delta"]), linewidth=1)

    print("\n-----------------------------------\n")

ImageName = "gamma_bin" + str(binsizeStr)
if calculateGamma:
    plt.legend()
    plt.tight_layout()
    if saveGammaImg:
        plt.savefig("../img/" + ImageName + ".jpg", dpi=150)
        print("Image saved as " + ImageName + ".jpg")
    #plt.show()
