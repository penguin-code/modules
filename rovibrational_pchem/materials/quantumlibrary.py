
import numpy as np
import matplotlib.pyplot as plt

# This sets up the "plotJstick" function so one can graph what the spectrum would look like
def lorentzian(y,x,x0,S,w):
    numerator =  w**2
    denominator = ( x - x0 )**2 + w**2
    ynew = y + S*(numerator/denominator)
    return ynew
def nustick2spectrum(nu,nu_lines,S_lines,w):
    y = np.zeros(np.size(nu))
    for i in range(len(nu_lines)):
        y = lorentzian(y,nu,nu_lines[i],S_lines[i],w)
    return y
def plotJstick(J,B,nu_offset,S_lines):
    nu_lines = nu_offset+np.array(B*(2*J+1))/1.986e-23
    nu = np.arange(np.min(nu_lines)-5,np.max(nu_lines)+20,0.01)
    y = nustick2spectrum(nu,nu_lines,S_lines,.02)
    plt.plot(nu,y)

