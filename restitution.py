"Plot restitution curves as per Figure 1"

import matplotlib.pyplot as plt
import numpy as np

# APD Restitution
#
fig, ax = plt.subplots()

APD = np.load("rest/APD.npz")
ax.plot(APD["DI"], APD["APD"], label="CM98-AF")

APD = np.load("rest/APD_bz.npz")
ax.plot(APD["DI"], APD["APD"], label="MS bz")

APD = np.load("rest/APD_healthy.npz")
ax.plot(APD["DI"], APD["APD"], label="MS healthy")

ax.grid()
ax.set_xlabel("DI [ms]")
ax.set_ylabel("APD [ms]")
ax.legend()

# CV restitution
#
fig, ax = plt.subplots()

CV = np.load("rest/CV.npz")
ax.plot(CV["DI"], CV["CV"], label="CM98-AF")

CV = np.load("rest/CV_bz.npz")
ax.plot(CV["DI"], CV["CV"], label="MS bz")

CV = np.load("rest/CV_healthy.npz")
ax.plot(CV["DI"], CV["CV"], label="MS healthy")

ax.grid()
ax.set_xlabel("DI [ms]")
ax.set_ylabel("CV [cm/ms]")
ax.legend()

plt.show()
