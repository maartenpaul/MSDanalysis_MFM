from numpy import exp, sin, asarray

from lmfit import minimize, Parameters

def residual(params, tl, time, density,weight):
    A = params['A']
    B = params['B']
    C = params['C']
    kb = params['kb']
    koff1 = params['koff1']
    koff2 = params['koff2']
    koff3 = params['koff3']
    tin=0.05
    model = A * (B*(kb*(tin/tl)+koff1)*exp(-(kb*(tin/tl)+koff1)*time)+C*(kb*(tin/tl)+koff2)*exp(-(kb*(tin/tl)+koff2)*time))
    if ((1-B-C)>0):
        return (density-model)/weight
    else: 
        return (density-model)*1E6
        
    

#A*(B*(kb*(tint/tl)+koff1)*exp(-((kb*(tint/tl)+koff1)*time))+
#                     C*(kb*(tint/tl)+koff2)*exp(-((kb*(tint/tl)+koff2)*time))+
#                     (1-B-C)*(kb*(tint/tl)+koff3)*exp(-(kb*((tint/tl)+koff3)*time)))

params = Parameters()
params.add('A', value=10)
params.add('B', value=0.5,min=0,max=1)
params.add('C', value=0.4,min=0,max=1)
params.add('kb', value=1,min=0)
params.add('koff1', value=0.0001,min=0)
params.add('koff2', value=0.1,min=0)
params.add('koff3', value=2,min=0)



tl =  asarray(tl)
time = asarray(time)
density = asarray(density)
weight = asarray(weight)

residual(params,tl,time,density,weight)

out = minimize(residual, params, args=(tl,time,density,weight))
out.success
out.params
out.residual


A = out.params['A']
B = out.params['B']
C = out.params['C']
kb = out.params['kb']
koff1 = out.params['koff1']
koff2 = out.params['koff2']
koff3 = out.params['koff3']
tin=0.05
model = A * (B*(kb*(tin/tl)+koff1)*exp(-(kb*(tin/tl)+koff1)*time)+C*(kb*(tin/tl)+koff2)*exp(-(kb*(tin/tl)+koff2)*time)
time
