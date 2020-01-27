x0, y0, x1, y1, x2, y2, trnums0, trnums1, trnums2 = getTrackPiecesForInfo(x, y, allStates)

numPmsd = 4
numPmss = 4
p = np.linspace(0.5, 6, 12)

info = []
i = 0
for xx, yy, tn in zip(x0 + x1 + x2, y0 + y1 + y2, trnums0 + trnums1 + trnums2):
    lx = len(xx)
    st = 0 + (i >= len(x0)) * 1 + (i >= len(x0 + x1)) * 1

    closest = min(indices, key = lambda x:abs(x-int(np.floor(tn))))
    cn = ((closest - int(np.floor(tn))) < 0) * indices.index(closest) + \
    ((closest - int(np.floor(tn))) > 0) * (indices.index(closest) - 1)

    if len(xx) > max(numPmsd, numPmss):
        dif, _, smss, _ = getMSDandMSS([xx], [yy], numPmsd, numPmss, p)
    else:
        dif = 'NA'
        smss = 'NA'

    info.append([st, cn, tn, lx, dif, smss])
    i = i +1 

colSt = [i[0] for i in info]
colCn = [i[1] for i in info]
colTn = [i[2] for i in info]
colLx = [i[3] for i in info]
colDi = [i[4] for i in info]
colSm = [i[5] for i in info]

df = pd.DataFrame({'State' : colSt, \
                   'Cell number' : colCn, \
                   'Track number': colTn, \
                   'Tracklet length': colLx, \
                   'Diffusion constant': colDi, \
                   'Smss': colSm})

dfSorted = df.sort_values(by=['Track number'])
dfSorted.columns.names = ['Tracklet number']
