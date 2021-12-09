# Get predicted/original states
reshapedTracks = []
numFeat = 1 + ('meanMSD' in addFeat) * maxOrder + ('xy' in addFeat) * 4



allStates = []

for fv in featVec:
  track = np.array(fv).reshape(1, len(fv), numFeat)
  predicted = model.predict(track)
  predSt = predicted.argmax(axis = 2)

  reshapedTracks.append(track)
  allStates.append(predSt[0])

# Get lists of distances in state 0, 1 and 2

state0 = []
state1 = []
state2 = []
state012 = []

for pred, tr in zip(allStates, reshapedTracks):
  for i in range(len(pred)):
    state012.append(tr[0][i][0])

    if pred[i] == 0:
      state0.append(tr[0][i][0])
    elif pred[i] == 1:
      state1.append(tr[0][i][0])
    else:
      state2.append(tr[0][i][0])

state0 = np.sort(state0)
state1 = np.sort(state1)
state2 = np.sort(state2)
state012 = np.sort(state012)

# Calculate D0, D1 and D2 from state0, state1 and state2

mean0 = ((pixSize * np.array(state0)) ** 2).mean()
mean1 = ((pixSize * np.array(state1)) ** 2).mean()
mean2 = ((pixSize * np.array(state2)) ** 2).mean()

# Brownian motion: D = (mean(d^2)) / (4 * time step)
D0 = mean0 / (4 * t)
D1 = mean1 / (4 * t)
D2 = mean2 / (4 * t)

print('Brownian motion: D0 (fast): %.6f, D1 (slow): %.6f, D2 (immobile): %.6f\
      \t(D = (mean(d^2) / (4 * time step))' %(D0, D1, D2))

# Rayleigh: D = mean(d)^2 / (pi * time step)
DD0 = (np.array(state0).mean() * pixSize) ** 2 / math.pi / t # (Rayleigh)
DD1 = (np.array(state1).mean() * pixSize) ** 2 / math.pi / t # (Rayleigh)
DD2 = (np.array(state2).mean() * pixSize) ** 2 / math.pi / t # (Rayleigh)

print('Rayleigh:\t D0 (fast): %.6f, D1 (slow): %.6f, D2 (immobile): %.6f\
      \t(D = mean(d)^2 / (pi * time step))\n' %(DD0, DD1, DD2))

#Get probability transition matrix P = P[[p00, p01], [p10, p11]]

Step00 = 0
Step01 = 0
Step02 = 0
Step10 = 0
Step11 = 0
Step12 = 0
Step20 = 0
Step21 = 0
Step22 = 0

for pred in allStates:
  for i in range(len(pred) - 1):
    if pred[i] == 0:
      if pred[i + 1] == 0:
        Step00 += 1
      elif pred[i + 1] == 1:
        Step01 += 1
      elif pred[i + 1] == 2:
        Step02 += 1
    elif pred[i] == 1:
      if pred[i + 1] == 0:
        Step10 += 1
      elif pred[i + 1] == 1:
        Step11 += 1
      elif pred[i + 1] == 2:
        Step12 += 1
    elif pred[i] == 2:
      if pred[i + 1] == 0:
        Step20 += 1
      elif pred[i + 1] == 1:
        Step21 += 1
      elif pred[i + 1] == 2:
        Step22 += 1

p00 = round(Step00 / (Step00 + Step01 + Step02), 3)
p01 = round(Step01 / (Step00 + Step01 + Step02), 3)
p02 = round(Step02 / (Step00 + Step01 + Step02), 3)
p10 = round(Step10 / (Step10 + Step11 + Step12), 3)
p11 = round(Step11 / (Step10 + Step11 + Step12), 3)
p12 = round(Step12 / (Step10 + Step11 + Step12), 3)
p20 = round(Step20 / (Step20 + Step21 + Step22), 3)
p21 = round(Step21 / (Step20 + Step21 + Step22), 3)
p22 = round(Step22 / (Step20 + Step21 + Step22), 3)

P = [[p00, p01, p02], [p10, p11, p12], [p20, p21, p22]]
print('Probability transition matrix:\n' + str(np.matrix(P).reshape((3,3))) + '\n')

# Calculate fraction of timepoints spent in state 0/1/2

st0 = 0
st1 = 0
st2 = 0

for pred in allStates:
  for i in range(len(pred) - 1):
    if pred[i] == 0:
      st0 += 1
    elif pred[i] == 1:
      st1 += 1
    else:
      st2 += 1

frac0 = st0 / (st0 + st1 + st2)
frac1 = st1 / (st0 + st1 + st2)
frac2 = st2 / (st0 + st1 + st2)

print('Fast fraction:\t\t%.6f\nSlow fraction:\t\t%.6f\nImmobile fraction:\t%.6f\n' %(frac0, frac1, frac2))
