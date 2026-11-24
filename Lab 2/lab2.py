import numpy as np
#from numpy.fft import fft2, ifft2, fftshift
from scipy.signal import convolve2d, correlate2d
import matplotlib.pyplot as plt

from Functions import *
from gaussfft import gaussfft

# Either write your code in a file like this or use a Jupyter notebook.
#
# A good idea is to use switches, so that you can turn things on and off
# depending on what you are working on. It should be fairly easy for a TA
# to go through all parts of your code though.

# # Exercise 1
# if 0:
# 	print("That was a stupid idea")
        
        
def deltax():
    return np.pad(0.5 * np.array([[0, 0, 0], 
                                  [1, 0, -1], 
                                  [0, 0, 0.]]), 1)

def deltay():
    return np.pad(0.5 * np.array([[0, 1, 0], 
                                  [0, 0, 0], 
                                  [0, -1, 0.]]), 1)

def deltaxx():
    return np.pad(np.array([[0, 0, 0], 
                            [1, -2, 1], 
                            [0, 0, 0.]]), 1)

def deltayy():
    return np.pad(np.array([[0, 1, 0], 
                            [0, -2, 0], 
                            [0, 1, 0.]]), 1)

def Lv(inpic, shape = 'same'):
	dxmask = deltax()
	dymask = deltay()
	Lx = convolve2d(inpic, dxmask, shape)
	Ly = convolve2d(inpic, dymask, shape)
	return np.sqrt(Lx**2 + Ly**2)

def Lvvtilde(inpic, shape = 'same'):
	dxmask = deltax()
	dymask = deltay()
	dxxmask = deltaxx()
	dyymask = deltayy()
	dxymask = convolve2d(dxmask, dymask, 'same')
	Lx = convolve2d(inpic, dxmask, shape)
	Ly = convolve2d(inpic, dymask, shape)
	Lxx = convolve2d(inpic, dxxmask, shape)
	Lyy = convolve2d(inpic, dyymask, shape)
	Lxy = convolve2d(inpic, dxymask, shape)

	result = Lx**2 * Lxx + 2 * Lx * Ly * Lxy + Ly**2 * Lyy

	return result

def Lvvvtilde(inpic, shape = 'same'):
    dxmask = deltax()
    dymask = deltay()
    dxxmask = deltaxx()
    dyymask = deltayy()
    dxxxmask = convolve2d(dxxmask, dxmask, 'same')
    dyyymask = convolve2d(dyymask, dymask, 'same') 
    dxxymask = convolve2d(dxxmask, dymask, 'same')
    dxyymask = convolve2d(dxmask, dyymask, 'same')

    Lx = convolve2d(inpic, dxmask, shape)
    Ly = convolve2d(inpic, dymask, shape)
    Lxxx = convolve2d(inpic, dxxxmask, shape)
    Lyyy = convolve2d(inpic, dyyymask, shape)
    Lxxy = convolve2d(inpic, dxxymask, shape)
    Lxyy = convolve2d(inpic, dxyymask, shape)

    result = Lx**3 * Lxxx + 3 * Lx**2 * Ly * Lxxy + 3 * Lx * Ly**2 * Lxyy + Ly**3 * Lyyy
    return result

def extractedge(inpic, scale, threshold = None, shape = 'same'):
    L = discgaussfft(inpic, scale)
    # Edges: zero crossings of Lvv where Lvvv < 0 (gradient magnitude maximum)
    contours = zerocrosscurves(Lvvtilde(L, shape), Lvvvtilde(L, shape) < 0)
    if threshold is not None:
        contours = thresholdcurves(contours, Lv(L, shape) > threshold)
    return contours
        
def drawline(rho, theta, center, length, style = 'r-'):
    # Point on the line closest to the origin, converted back to image coordinates
    x0 = rho * np.cos(theta) + center[0]
    y0 = rho * np.sin(theta) + center[1]
    # Step along the line direction, perpendicular to the normal (cos θ, sin θ)
    dx = -length * np.sin(theta)
    dy = length * np.cos(theta)
    plt.plot([x0 - dx, x0, x0 + dx], [y0 - dy, y0, y0 + dy], style)

def houghline(curves, magnitude, nrho, ntheta, threshold, nlines = 20, verbose = False,
              increment = None, accsmooth = 0):
    h, w = np.shape(magnitude)

    # Coordinate system: origin at the image center, x to the right, y downwards.
    # theta in [-pi/2, pi/2) and rho in [-rmax, rmax] cover every line exactly once.
    center = (w / 2, h / 2)
    rmax = np.hypot(w / 2, h / 2)
    thetas = np.linspace(-np.pi / 2, np.pi / 2, ntheta, endpoint = False)
    rhos = np.linspace(-rmax, rmax, nrho)
    drho = rhos[1] - rhos[0]

    # Allocate accumulator space, rows = rho, columns = theta
    acc = np.zeros((nrho, ntheta))

    # Loop over all the edge points
    (Y, X) = curves
    for i in range(len(X)):
        row = Y[i]
        col = X[i]

        # Check if valid point with respect to threshold
        mag = magnitude[row, col]
        if mag <= threshold:
            continue

        # Accumulator increment: 1 per point, or h(|grad L|) if an increment function is given
        if increment is None:
            vote = 1
        else:
            vote = increment(mag)

        # Centered coordinates of the point
        x = col - center[0]
        y = row - center[1]

        # Loop over a set of theta values
        for thetaidx in range(ntheta):
            theta = thetas[thetaidx]

            # Compute rho for this theta value
            rho = x * np.cos(theta) + y * np.sin(theta)

            # Compute the index of the closest rho cell
            rhoidx = int(round((rho + rmax) / drho))

            # Update the accumulator
            acc[rhoidx, thetaidx] += vote

    # Optionally smooth the accumulator to merge multiple responses for one line
    if accsmooth > 0:
        acc = discgaussfft(acc, accsmooth)

    # Extract local maxima from the accumulator.
    # The theta axis wraps around: theta = -pi/2 is the same line as theta = +pi/2 with
    # rho negated. So we add one extra column on each side, taken from the opposite end
    # of the accumulator and flipped in rho. Then a horizontal line gives a single peak.
    accwrap = np.zeros((nrho, ntheta + 2))
    accwrap[:, 1:ntheta + 1] = acc
    for rhoidx in range(nrho):
        flippedidx = nrho - 1 - rhoidx              # index of -rho
        accwrap[rhoidx, 0] = acc[flippedidx, ntheta - 1]
        accwrap[rhoidx, ntheta + 1] = acc[flippedidx, 0]

    pos, value, _ = locmax8(accwrap)

    # Keep only maxima in the real accumulator (not in the extra columns),
    # and shift the theta index back by one to undo the padding
    peaks = []
    for k in range(len(value)):
        thetaidx = pos[k, 0] - 1
        rhoidx = pos[k, 1]
        if 0 <= thetaidx < ntheta:
            peaks.append((value[k], rhoidx, thetaidx))

    # Delimit the number of responses: sort strongest first and keep nlines
    peaks.sort(reverse = True)
    peaks = peaks[:nlines]

    # Compute a line (rho, theta) for each one of the strongest responses
    linepar = []
    for _, rhoidx, thetaidx in peaks:
        linepar.append((rhos[rhoidx], thetas[thetaidx]))

    if verbose:
        plt.figure(figsize = (12, 5))
        plt.subplot(1, 2, 1)
        plt.title("Accumulator")
        plt.imshow(acc, cmap = 'gray', aspect = 'auto',
                   extent = [np.degrees(thetas[0]), np.degrees(thetas[-1]), rhos[-1], rhos[0]])
        plt.plot([np.degrees(t) for _, t in linepar], [r for r, _ in linepar], 'ro', mfc = 'none')
        plt.xlabel("theta (degrees)")
        plt.ylabel("rho (pixels)")
        plt.subplot(1, 2, 2)
        plt.title("Lines on gradient magnitude")
        showgrey(magnitude, display = False)
        for rho, theta in linepar:
            drawline(rho, theta, center, rmax)
        plt.xlim(0, w)
        plt.ylim(h, 0)
        plt.show()
        if verbose > 1:
            for rho, theta in linepar:
                print(f"rho = {rho:7.2f}, theta = {np.degrees(theta):7.2f} deg")

    return linepar, acc

def houghedgeline(pic, scale, gradmagnthreshold, nrho, ntheta, nlines = 20, verbose = False,
                  increment = None, accsmooth = 0):
    # Edge detection without threshold; houghline() thresholds on the magnitude instead
    curves = extractedge(pic, scale)
    magnitude = Lv(discgaussfft(pic, scale), 'same')

    linepar, acc = houghline(curves, magnitude, nrho, ntheta, gradmagnthreshold, nlines,
                             max(verbose - 1, 0), increment, accsmooth)

    if verbose:
        h, w = np.shape(pic)
        plt.figure(figsize = (15, 5))
        plt.subplot(1, 3, 1)
        plt.title("Edges")
        overlaycurves(pic, thresholdcurves(curves, magnitude > gradmagnthreshold))
        plt.subplot(1, 3, 2)
        plt.title("Accumulator")
        plt.imshow(acc, cmap = 'gray', aspect = 'auto')
        plt.xlabel("theta index")
        plt.ylabel("rho index")
        plt.subplot(1, 3, 3)
        plt.title(f"Hough transform ({len(linepar)} lines)")
        showgrey(pic, display = False)
        for rho, theta in linepar:
            drawline(rho, theta, (w / 2, h / 2), np.hypot(w, h))
        plt.xlim(0, w)
        plt.ylim(h, 0)
        plt.show()

    return linepar, acc
         
