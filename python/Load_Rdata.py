# Set pixelsize [um] and time step [s] corresponding to acquisition parameters
pixSize = 0.100
t = 0.032

##########################################################################################################################


# Feature extraction (in principle, parameters never need to be changed)
addFeat = ['meanMSD', 'xy']
maxOrder = 2
shift = 2
nClasses = 2
min1state = 1

d = getDist(x, y)
featVec = getFeatVec(d, x, y, addFeat, maxOrder, shift)

print('Data loaded and features extracted')
