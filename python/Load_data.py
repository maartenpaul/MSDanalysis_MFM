# Select folder with files that need to be analyzed
pathToFiles =  MLfolder # Should always end with /

# Set pixelsize [um] and time step [s] corresponding to acquisition parameters
pixSize = 0.100
t = 0.030

# Select all files from folder
# If not all files should be selected, provide array with seperate filenames (filenames as strings with .txt)
filename = []
for file in os.listdir(pathToFiles):
    filename.append(file)


##########################################################################################################################


# Feature extraction (in principle, parameters never need to be changed)
addFeat = ['meanMSD', 'xy']
maxOrder = 2
shift = 2
nClasses = 2
min1state = 1

x = []
y = []
indices = [0]
indfile = 0

for name in filename:
    xfile, yfile = loadRealData(np.loadtxt(pathToFiles + name))
    indfile += len(xfile)
    indices.append(indfile)
    for trackx, tracky in zip(xfile, yfile):
        x.append(trackx)
        y.append(tracky)

d = getDist(x, y)
featVec = getFeatVec(d, x, y, addFeat, maxOrder, shift)

foldername = pathToFiles[[i for i, j in enumerate(pathToFiles) if j == '/'][-2] + 1:\
                     [i for i, j in enumerate(pathToFiles) if j == '/'][-1]]
print('Folder: ' + str(np.array(foldername)) + ' (' + str(len(filename)) + ' files) \n')
print('Data loaded and features extracted')
