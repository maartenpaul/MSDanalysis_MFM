def getMSDandMSS_R(x, y):
  numPmsd = 4
  numPmss = 4
  minLen = 10
  p = np.linspace(0.5, 6, 12) 
  
  dif, _, smss, _ = getMSDandMSS([x], [y], numPmsd, numPmss, p)
  
  return dif, smss
