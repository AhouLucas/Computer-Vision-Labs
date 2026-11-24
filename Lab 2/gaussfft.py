import numpy as np
from numpy.fft import fft2, ifft2, fftshift

def gaussfft(pic, t):
    h, w = np.shape(pic)

    # Sampled Gaussian of variance t, centred in the image
    x, y = np.meshgrid(np.arange(w) - w//2, np.arange(h) - h//2)
    g = np.exp(-(x**2 + y**2) / (2*t)) / (2*np.pi*t)
    g = g / np.sum(g)

    # Move the centre of the Gaussian to the origin (0, 0)
    g = fftshift(g)

    # Convolution = multiplication in the Fourier domain
    result = np.real(ifft2(fft2(pic) * fft2(g)))
    return result
