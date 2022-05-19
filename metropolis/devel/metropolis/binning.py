import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use('TkAgg')

filename = "../data/binning.csv"

data = pd.read_csv(filename)
data["bin"] = data["bin"] * 1e4
print(data)
bin = data["bin"].to_list()
absM = data["absM"].to_list()
M2 = data["M2"].to_list()
chi = data["chi"].to_list()
binder = data["binder"].to_list()

plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.locator_params(nbins=20, axis="x")
plt.xlim(0, 5.5e5)
plt.ylim(-0.001, 0.022)
plt.xlabel("Dimensione dei bin")
plt.ylabel("Errore")
plt.plot(bin, absM, label = r"$|m|$", marker = "o", markersize = 4)
plt.plot(bin, M2, label = r"$m^2$", marker = "o", markersize = 4)
plt.plot(bin, chi, label = r"$\chi$", marker = "o", markersize = 4)
plt.plot(bin, binder, label = r"$B$", marker = "o", markersize = 4)

plt.legend(ncol=4)
plt.savefig("../img/met_binning.jpg", dpi=200)
plt.show()
