import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
#os.environ['QT_QPA_PLATFORM_PLUGIN_PATH'] = 'C:/Users/maart/Anaconda3/Library/plugins/platforms'
numPmsd = 4
numPmss = 4
minLen = 10

p = np.linspace(0.5, 6, 12) # power

# Calculations for whole dataset
arDiffSw = []
arSmssSw = []
arDiffF = []
arSmssF = []
arDiffS = []
arSmssS = []
arDiffI = []
arSmssI = []
for i in range(len(x)):
  if len(x[i]) > max(numPmsd, numPmss, minLen):
    dif, _, smss, _ = getMSDandMSS([x[i]], [y[i]], numPmsd, numPmss, p)
    if (dif >= 0) and (smss >= 0):
      if list(allStates[i]).count(allStates[i][0]) !=  len(allStates[i]):
        arDiffSw.append(dif)
        arSmssSw.append(smss)
      elif allStates[i][0] == 0:
        arDiffF.append(dif)
        arSmssF.append(smss)
      elif allStates[i][0] == 1:
        arDiffS.append(dif)
        arSmssS.append(smss)
    else:
      arDiffI.append(dif)
      arSmssI.append(smss)


arDiff = arDiffSw + arDiffF + arDiffS + arDiffI
arSmss = arSmssSw + arSmssF + arSmssS + arSmssI

  # Calculations per state
x0, y0, x1, y1, x2, y2 = getTrackPieces(x, y, allStates)

arDiff0 = []
arSmss0 = []
for xx0, yy0 in zip(x0, y0):
  if len(xx0) > max(numPmsd, numPmss, minLen):
    dif0, _, smss0, _ = getMSDandMSS([xx0], [yy0], numPmsd, numPmss, p)
    if dif0 >= 0 and (smss0 >= 0):
      arDiff0.append(dif0)
      arSmss0.append(smss0)


arDiff1 = []
arSmss1 = []
for xx1, yy1 in zip(x1, y1):
  if len(xx1) > max(numPmsd, numPmss, minLen):
    dif1, mss1, smss1, _ = getMSDandMSS([xx1], [yy1], numPmsd, numPmss, p)
    if dif1 >= 0 and (smss1 >= 0):
      arDiff1.append(dif1)
      arSmss1.append(smss1)

arDiff2 = []
arSmss2 = []
for xx2, yy2 in zip(x2, y2):
  if len(xx2) > max(numPmsd, numPmss, minLen):
    dif2, mss2, smss2, _ = getMSDandMSS([xx2], [yy2], numPmsd, numPmss, p)
    if dif2 >= 0 and (smss2 >= 0):
      arDiff2.append(dif2)
      arSmss2.append(smss2)


plt.style.use(['classic', 'seaborn-darkgrid'])
plt.rcParams['figure.figsize'] = (12, 6)

# Plot total (whole tracks)
data = np.column_stack((arDiff, arSmss))
df = pd.DataFrame(data, columns = [r'$D$ $\mathrm{[\mu m^2/s]}$', r'$S_{\mathrm{MSS}}$'])
g = sns.JointGrid(r'$D$ $\mathrm{[\mu m^2/s]}$',r'$S_{\mathrm{MSS}}$', data = df)

ax = g.ax_joint
ax.set_xscale('log')
g.plot_joint(plt.scatter, color = 'darkorange', alpha = 0.04, edgecolor = 'darkorange')
ax.axhspan(0.4, 0.6, facecolor = 'slateblue', edgecolor = 'none', alpha = 0.1)
ax.set_ylim(-0.5, 1)
ax.set_xlim(1e-6,50)

g.ax_marg_x.set_xscale('log')
g.ax_marg_x.hist(df[r'$D$ $\mathrm{[\mu m^2/s]}$'], color = 'darkorange', edgecolor = 'none', alpha = 0.4, \
                 bins = np.logspace(-6, np.log10(100), 50))

g.ax_marg_y.hist(df[r'$S_{\mathrm{MSS}}$'], color = 'darkorange', edgecolor = 'none', alpha = 0.4, \
                 orientation = 'horizontal', bins = np.linspace(-0.5, 1, 50))
g.ax_marg_y.axhspan(0.4, 0.6, facecolor = 'slateblue', edgecolor = 'none', alpha = 0.1)

orange_patch = mpatches.Patch(color = 'darkorange', alpha = 0.4, label = 'Total')
ax.legend(handles = [orange_patch], loc = 'upper left')
plt.xlabel(r'$D$ $\mathrm{[\mu m^2/s]}$', fontdict = font, size = 'large')
plt.ylabel(r'$S_{\mathrm{MSS}}$', fontdict = font, size = 'large')
plt.show()
 
plt.savefig(directory+condition+"/"+"WholeTracksplot_" + state + ".png",dpi=600)
plt.savefig(directory+condition+"/"+"WholeTracksplot_" + state + ".pdf",dpi=600)
plt.close()
# Plot per trackstate (whole tracks only fast / only slow / only immobile / switching)
dataSw = np.column_stack((arDiffSw, arSmssSw))
dfSw = pd.DataFrame(dataSw, columns = [r'$D$ $\mathrm{[\mu m^2/s]}$', r'$S_{\mathrm{MSS}}$'])
gSw = sns.JointGrid(r'$D$ $\mathrm{[\mu m^2/s]}$', r'$S_{\mathrm{MSS}}$', data = dfSw)

axSw = gSw.ax_joint
axSw.set_xscale('log')
axSw.set_ylim(-0.5, 1)
gSw.plot_joint(plt.scatter, color = 'green', alpha = 0.04, edgecolor = 'green')
axSw.scatter(arDiffF, arSmssF, color = 'r', alpha = 0.04, edgecolor = 'r')
axSw.scatter(arDiffS, arSmssS, color = 'royalblue', alpha = 0.04, edgecolor = 'royalblue')
axSw.scatter(arDiffI, arSmssI, color = 'darkblue', alpha = 0.04, edgecolor = 'darkblue')
axSw.axhspan(0.4, 0.6, facecolor = 'slateblue', edgecolor = 'none', alpha = 0.1)
axSw.set_xlim(1e-6,50)

gSw.ax_marg_x.set_xscale('log')
gSw.ax_marg_x.hist(dfSw[r'$D$ $\mathrm{[\mu m^2/s]}$'], color = 'green', edgecolor = 'none', alpha = 0.4, \
                   bins = np.logspace(np.log10(1e-6), np.log10(100), 50))
gSw.ax_marg_x.hist(arDiffF, color = 'red', edgecolor = 'none', alpha = 0.4, \
                   bins = np.logspace(np.log10(1e-6), np.log10(100), 50))
gSw.ax_marg_x.hist(arDiffS, color = 'royalblue', edgecolor = 'none', alpha = 0.4, \
                   bins = np.logspace(np.log10(1e-6), np.log10(100), 50))
gSw.ax_marg_x.hist(arDiffI, color = 'darkblue', edgecolor = 'none', alpha = 0.4, \
                   bins = np.logspace(np.log10(1e-6), np.log10(100), 50))

gSw.ax_marg_y.hist(dfSw[r'$S_{\mathrm{MSS}}$'], color = 'green', edgecolor = 'none', alpha = 0.4, \
                   orientation = 'horizontal', bins = np.linspace(-.5, 1, 50))
gSw.ax_marg_y.hist(arSmssF, color = 'red', edgecolor = 'none', alpha = 0.4, \
                   orientation = 'horizontal', bins = np.linspace(-.5, 1, 50))
gSw.ax_marg_y.hist(arSmssS, color = 'royalblue', edgecolor = 'none', alpha = 0.4, \
                   orientation = 'horizontal', bins = np.linspace(-.5, 1, 50))
gSw.ax_marg_y.hist(arSmssI, color = 'darkblue', edgecolor = 'none', alpha = 0.4, \
                   orientation = 'horizontal', bins = np.linspace(-.5, 1, 50))
gSw.ax_marg_y.axhspan(0.4, 0.6, facecolor = 'slateblue', edgecolor = 'none', alpha = 0.1)

red_patch = mpatches.Patch(color = 'red', alpha = 0.4, label = 'Only fast')
royalblue_patch = mpatches.Patch(color = 'royalblue', alpha = 0.4, label = 'Only slow')
darkblue_patch = mpatches.Patch(color = 'darkblue', alpha = 0.4, label = 'Only immobile')
green_patch = mpatches.Patch(color = 'green', alpha = 0.4, label = 'Switching')
axSw.legend(handles = [red_patch, royalblue_patch, darkblue_patch, green_patch], loc = 'upper left')
plt.xlabel(r'D $\mathrm{(\mu m^2/sec)}$', fontdict = font, size = 'large')
plt.ylabel(r'$S_{\mathrm{MSS}}$', fontdict = font, size = 'large')
plt.show()
##
plt.savefig(directory+condition+"/"+"SplitTracksplot_" + state + ".png",dpi=600)
plt.savefig(directory+condition+"/"+"SplitTracksplot_" + state + ".pdf",dpi=600)
plt.close()
# Plot split (trackpieces fast / slow / immobile)
data0 = np.column_stack((arDiff0, arSmss0))
df0 = pd.DataFrame(data0, columns = [r'$D$ $\mathrm{[\mu m^2/s]}$', r'$S_{\mathrm{MSS}}$'])
f = sns.JointGrid(r'$D$ $\mathrm{[\mu m^2/s]}$', r'$S_{\mathrm{MSS}}$', data = df0)

ax0 = f.ax_joint
ax0.set_xscale('log')
f.plot_joint(plt.scatter, color = 'r', alpha = 0.1, edgecolor = 'r')
ax0.scatter(arDiff1, arSmss1, color = 'royalblue', alpha = 0.1, edgecolor = 'royalblue')
ax0.scatter(arDiff2, arSmss2, color = 'darkblue', alpha = 0.1, edgecolor = 'darkblue')
ax0.axhspan(0.4, 0.6, facecolor = 'slateblue', edgecolor = 'none', alpha = 0.1)
ax0.set_ylim(-0.5, 1)
ax0.set_xlim(1e-6,50)


f.ax_marg_x.set_xscale('log')
f.ax_marg_x.hist(df0[r'$D$ $\mathrm{[\mu m^2/s]}$'], color = 'red', edgecolor = 'none', alpha = 0.4, \
                 bins = np.logspace(np.log10(1e-6), np.log10(100), 50))
f.ax_marg_x.hist(arDiff1, color = 'royalblue', edgecolor = 'none', alpha = 0.4, \
                 bins = np.logspace(np.log10(1e-6), np.log10(100), 50))
f.ax_marg_x.hist(arDiff2, color = 'darkblue', edgecolor = 'none', alpha = 0.4, \
                 bins = np.logspace(np.log10(1e-6), np.log10(100), 50))

f.ax_marg_y.hist(df0[r'$S_{\mathrm{MSS}}$'], color = 'red', edgecolor = 'none', alpha = 0.4, \
                 orientation = 'horizontal', bins = np.linspace(0, 1, 50))
f.ax_marg_y.hist(arSmss1, color = 'royalblue', edgecolor = 'none', alpha = 0.4, orientation = 'horizontal', \
                 bins = np.linspace(-0.5, 1, 50))
f.ax_marg_y.hist(arSmss2, color = 'darkblue', edgecolor = 'none', alpha = 0.4, orientation = 'horizontal', \
                 bins = np.linspace(-0.5, 1, 50))
f.ax_marg_y.axhspan(0.4, 0.6, facecolor = 'slateblue', edgecolor = 'none', alpha = 0.1)

red_patch = mpatches.Patch(color = 'red', alpha = 0.4, label = 'Fast tracklets')
royalblue_patch = mpatches.Patch(color = 'royalblue', alpha = 0.4, label = 'Slow tracklets')
darkblue_patch = mpatches.Patch(color = 'darkblue', alpha = 0.4, label = 'Immobile tracklets')
ax0.legend(handles = [red_patch, royalblue_patch, darkblue_patch], loc = 'upper left')
plt.xlabel(r'$D$ $\mathrm{[\mu m^2/s]}$', fontdict = font, size = 'large')
plt.ylabel(r'$S_{\mathrm{MSS}}$', fontdict = font, size = 'large')
plt.show()
plt.savefig(directory+condition+"/"+"SMMsplot_tracklets_" + state + ".png",dpi=600)
plt.savefig(directory+condition+"/"+"SMMsplot_tracklets_" + state + ".pdf",dpi=600)
plt.close()
