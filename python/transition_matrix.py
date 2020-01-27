# Get probability transition matrix P = P[[p00, p01], [p10, p11]]

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
