import os
os.environ['QT_QPA_PLATFORM_PLUGIN_PATH'] = 'C:/Users/maart/Anaconda3/Library/plugins/platforms'


# Three possible plots: gamma versus p, log(Cp) versus p and D versus p
# Specify which plots to produce (by setting True or False):
plotGvsP = True # Plot gamma versus p
plotCvsP = False # Plot log(Cp) versus p
plotDvsP = False # Plot D versus p

if 1 != 1:
  print('Run section 4 first!')
else:
  #  print('Folder: ' + str(np.array(foldername)) + ' (' + str(len(filename)) + ' files) \n')
  numPmsd = 4
  numPmss = 4
  minLen = 10
  p = np.linspace(0.5, 6, 12) # power
  b = 'unknown' # b can be set to 'zero' if MSD (ax + b) needs to be calculated with b = 0. Default is b = 'unknown'.

# Calculations for whole dataset
  goodX = []
  goodY = []
  xD = []
  yD = []
  for xx, yy in zip(x, y):
    if len(xx) > max(numPmsd, numPmss, minLen):
      diff, _, Smss, _ = getMSDandMSS([xx], [yy], numPmsd, numPmss, p, b = b)
      if (diff >= 0) and (Smss >= 0):
        goodX.append(xx)
        goodY.append(yy)
        if (Smss <= 0.6) and (Smss >= 0.4):
          xD.append(xx)
          yD.append(yy)

  diff, MSS, C, cD, Smss, intercept = getMSDandMSSandC(goodX, goodY, numPmsd, numPmss, p, b = b)

  # Calculations per state
  x0, y0, x1, y1, x2, y2 = getTrackPieces(x, y, allStates)

  goodX0 = []
  goodY0 = []
  xD0 = []
  yD0 = []
  for xx0, yy0 in zip(x0, y0):
    if len(xx0) > max(numPmsd, numPmss, minLen):
      diff0, _, Smss0, _ = getMSDandMSS([xx0], [yy0], numPmsd, numPmss, p, b = b)
      if (diff0 >= 0) and (Smss0 >= 0):
        goodX0.append(xx0)
        goodY0.append(yy0)


  goodX1 = []
  goodY1 = []
  xD1 = []
  yD1 = []
  for xx1, yy1 in zip(x1, y1):
    if len(xx1) > max(numPmsd, numPmss, minLen):
      diff1, _, Smss1, _ = getMSDandMSS([xx1], [yy1], numPmsd, numPmss, p, b = b)
      if (diff1 >= 0) and (Smss1 >= 0):
        goodX1.append(xx1)
        goodY1.append(yy1)


  goodX2 = []
  goodY2 = []
  xD2 = []
  yD2 = []
  for xx2, yy2 in zip(x2, y2):
    if len(xx2) > max(numPmsd, numPmss, minLen):
      diff2, _, Smss2, _ = getMSDandMSS([xx2], [yy2], numPmsd, numPmss, p, b = b)
      if (diff2 >= 0) and (Smss2 >= 0):
        goodX2.append(xx2)
        goodY2.append(yy2)


  diff0, MSS0, C0, cD0, Smss0, intercept0 = getMSDandMSSandC(goodX0, goodY0, numPmsd, numPmss, p, b = b)
  diff1, MSS1, C1, cD1, Smss1, intercept1 = getMSDandMSSandC(goodX1, goodY1, numPmsd, numPmss, p, b = b)
  diff2, MSS2, C2, cD2, Smss2, intercept2 = getMSDandMSSandC(goodX2, goodY2, numPmsd, numPmss, p, b = b)
  print('Fast\t\tMSD: D = %.4f \t\tMSS: Smss = %.4f, Intercept = %.4f' \
      %(diff0, Smss0, intercept0))
  print('Slow\t\tMSD: D = %.4f \t\tMSS: Smss = %.4f, Intercept = %.4f' \
      %(diff1, Smss1, intercept1))
  print('Immobile\tMSD: D = %.4f \t\tMSS: Smss = %.4f, Intercept = %.4f\n' \
      %(diff2, Smss2, intercept2))


# Plotting
  plt.style.use(['classic', 'seaborn-darkgrid'])
  plt.rcParams['figure.figsize'] = (7, 6)

  if plotGvsP == True:
    plt.figure(1)
    plt.fill_between(p, p / 2.0, p, facecolor = 'slateblue', edgecolor = 'none', alpha = 0.1, \
                 label = 'Regions from diffusion to linear motion')
    plt.plot(p, MSS, '-o', color = 'darkorange', label = r'MSS total, $S_\mathrm{MSS}$ = %.2f' %(Smss))
    plt.plot(p, MSS0, '-o', color = 'red', label = r'MSS fast, $S_\mathrm{MSS}$ = %.2f' %(Smss0))
    plt.plot(p, MSS1, '-o', color = 'royalblue', label = r'MSS slow, $S_\mathrm{MSS}$ = %.2f' %(Smss1))
    plt.plot(p, MSS2, '-o', color = 'darkblue', label = r'MSS immobile, $S_\mathrm{MSS}$ = %.2f' %(Smss2))
    (plt.gca()).set_ylim(0, 3.2)
    (plt.gca()).set_xlim(p[0], p[-1])
    plt.xlabel(r'$p$', labelpad = 10, fontdict = font, size = 'large')
    plt.ylabel(r'$\mathrm{\gamma_p}$', labelpad = 10, fontdict = font, size = 'x-large')
    plt.legend(loc = 'upper left', bbox_to_anchor = (1, 1.02))

  if plotCvsP == True:
    plt.figure(2)
    plt.plot(p, C, '-o', color = 'darkorange', label = 'Total')
    plt.plot(p, C0, '-o', color = 'red', label = 'Fast')
    plt.plot(p, C1, '-o', color = 'royalblue', label = 'Slow')
    plt.plot(p, C2, '-o', color = 'darkblue', label = 'Immobile')
    (plt.gca()).set_xlim(p[0], p[-1])
    plt.xlabel(r'$p$', labelpad = 10, fontdict = font, size = 'large')
    plt.ylabel(r'$\mathrm{log} \ C_\mathrm{p}$', labelpad = 10, fontdict = font, size = 'large')
    plt.legend(loc = 'upper left', bbox_to_anchor = (1, 1.02))

  if plotDvsP == True:
    plt.figure(3)
    plt.plot(p, np.array(cD)  * pixSize ** 2 / t, '-o', color = 'darkorange', label = 'Total')
    plt.plot(p, np.array(cD0) * pixSize ** 2 / t, '-o', color = 'red', label = 'Fast')
    plt.plot(p, np.array(cD1) * pixSize ** 2 / t, '-o', color = 'royalblue', label = 'Slow')
    plt.plot(p, np.array(cD2) * pixSize ** 2 / t, '-o', color = 'darkblue', label = 'Immobile')
    (plt.gca()).set_xlim(p[0], p[-1])
    plt.xlabel(r'$p$', labelpad = 10, fontdict = font, size = 'large')
    plt.ylabel(r'$D \ \mathrm{[\mu m^2/s]}$', labelpad = 10, fontdict = font, size = 'large')
    plt.legend(loc = 'upper left', bbox_to_anchor = (1, 1.02))

  plt.show()
