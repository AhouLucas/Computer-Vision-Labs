# Lab 1 – Answers to questions 11–20

Formulas are written in LibreOffice Math syntax, each on its own line, ready to be converted into equation objects.

---

## Question 11

The image is now an 8 × 32 rectangle (8 pixels high, 32 pixels wide) instead of the 16 × 16 square of section 1.5. Its spectrum is still the product of two sincs, one per axis, but the two sincs have different widths:

- Vertically, the rectangle is narrower (8 instead of 16 pixels), so its spectrum is wider along the vertical axis: the zeros are every N/8 = 16 samples instead of every 8.
- Horizontally, the rectangle is wider (32 instead of 16 pixels), so its spectrum is narrower along the horizontal axis: the zeros are every N/32 = 4 samples.

So stretching the image in one direction compresses its spectrum in the same direction, and the other way round. This is the scaling property of the Fourier transform:

func F lbrace f(a x) rbrace (%omega) = {1 over abs{a}} hat f({%omega over a})

A wide structure contains mostly low frequencies; a narrow structure contains high frequencies.

---

## Question 12

Rotating the image rotates its Fourier spectrum by the same angle. The spectrum of the rotated rectangle has the same cross-shaped sinc pattern as before, but tilted by alpha. When we rotate the spectrum back by −alpha, we get the original spectrum again: the main horizontal and vertical lines are back in place.

This follows from the definition: a rotation R of the coordinates gives

func F lbrace f(R x) rbrace (%omega) = hat f(R %omega)

since R is orthogonal (R^T = R^{-1}). Rotating does not create new frequencies, it only changes their direction.

Differences between the original and the rotated-back spectra:

- Rotating a discrete image requires interpolation (here nearest neighbour), so the edges of the rotated rectangle become jagged. This adds small high-frequency components that are not in the original spectrum (visible as extra curved patterns).
- The spectrum is periodic, and its pattern continues across the image borders. After rotation, the corners of the spectrum are cut off (black regions in the rotated-back spectrum) and the pattern wraps around differently.
- For alpha = 90°, the rotation is exact on the pixel grid, so the rotated-back spectrum matches the original apart from the cropped border.

---

## Question 13

pow2image keeps the phase of the Fourier transform and replaces the magnitude with an artificial spectrum abs{hat f(%omega)}^2 = 1 / (a + abs{%omega}^2). The result is still easy to recognise: the edges, contours and positions of the objects (the phone, the tools, the mountains) are all in the right place. Only the contrast and the grey levels change.

randphaseimage keeps the magnitude and replaces the phase with random values. The result looks like noise or texture. Nothing of the original image can be recognised, even though the magnitude spectrum is exactly the same.

Conclusion:

- The phase contains the information about where things are: the position of the edges and lines, and therefore the structure of the image. Most of what we see in an image is in the phase.
- The magnitude contains the information about how much of each frequency there is: the overall contrast, the amount of fine detail and the dominant orientations, but not where they are.

This also shows that displaying only the magnitude of the Fourier transform hides the most important part of the information.

---

## Question 14

The impulse responses (psf = gaussfft(deltafcn(128, 128), t)) are Gaussian blobs that grow wider with t. For t = 0.1 the response is almost a single pixel.

Measured covariance matrices (the off-diagonal terms are all about 10^{-15}, i.e. zero):

| t | variance(psf) (diagonal terms) |
|---|---|
| 0.1 | 0.0133 |
| 0.3 | 0.281 |
| 1.0 | 1.000 |
| 10.0 | 10.000 |
| 100.0 | 100.000 |

---

## Question 15

For t >= 1, the measured variance is equal to t, as in the ideal continuous case:

C(g(cdot, cdot; t)) = t left( matrix{ 1 # 0 ## 0 # 1 } right)

For small t, the measured variance is smaller than t: 0.0133 instead of 0.1, and 0.281 instead of 0.3.

The reason is sampling. A Gaussian with a small variance is very narrow compared to the pixel size. When it is sampled on the pixel grid, almost all of its weight falls on the centre pixel, and the neighbours get very little: for t = 0.1 the neighbours have weight e^{-1/(2 cdot 0.1)} = e^{-5} approx 0.007 relative to the centre. The discrete kernel is therefore narrower than the continuous one.

In the Fourier domain the same thing appears from the other side. The Fourier transform of the Gaussian is

hat g(%omega; t) = e^{-abs{%omega}^2 t / 2}

which has variance 1/t. For small t, this is very wide, wider than the frequency range [−%pi, %pi] that the grid can represent. Sampling in the spatial domain causes aliasing in the Fourier domain, so the discrete kernel does not behave like the continuous one. For large t the Gaussian is wide in space and narrow in frequency, it is well sampled, and the discrete result matches the continuous one.

---

## Question 16

Convolving with a Gaussian blurs the image, and the blur grows with t (the standard deviation is sqrt t):

- t = 1 and t = 4: the noise and the finest details (small text, texture) disappear, but the objects and their edges are still sharp enough to recognise.
- t = 16: edges become clearly blurred and small objects (the keys of the phone and the calculator) merge together.
- t = 64 and t = 256: only large structures remain (large dark and bright regions); all detail is lost.

The Gaussian is a low-pass filter. A larger t means a narrower Gaussian in the Fourier domain, so more high frequencies are removed. Since the convolution is done with the FFT, it is periodic: near the borders, the blur mixes in pixels from the opposite side of the image (for large t, the grey levels along the edges are slightly changed by the content of the opposite edge).

---

## Question 17

The image add contains additive Gaussian noise: each pixel gets a random value from a normal distribution added to it. The image sap contains salt-and-pepper noise: a fraction of the pixels is replaced by pure black or pure white.

Gaussian smoothing:

- Positive: works well on Gaussian noise. The noise is averaged out, and the result looks smooth and natural.
- Negative: blurs the edges together with the noise. On salt-and-pepper noise it only spreads each white or black pixel into a grey blob, so the noise remains visible unless t is very large, and then the image is very blurry.
- Parameter: a larger t removes more noise but also more detail.

Median filtering:

- Positive: very effective on salt-and-pepper noise. Isolated extreme pixels are never the median of their neighbourhood, so they are removed completely, and edges stay sharp.
- Negative: with a large window, thin lines and small details disappear and regions become flat with a "painted" look (blocky regions of constant grey). It is less effective on Gaussian noise, where it gives a blotchy result.
- Parameter: a larger window removes more noise but removes more detail and makes the painted effect stronger.

Ideal low-pass filtering:

- Positive: simple to apply, removes all frequencies above the cutoff exactly.
- Negative: produces ringing (Gibbs phenomenon): oscillations appear around edges, because the sharp cutoff in the Fourier domain corresponds to a sinc function in the spatial domain. On salt-and-pepper noise each noisy pixel becomes a small ring pattern.
- Parameter: a lower cutoff removes more noise but gives a blurrier image with stronger, wider ringing.

---

## Question 18

The best filter depends on the type of noise:

- Gaussian noise is best handled by Gaussian smoothing, which averages the noise out without creating artefacts.
- Salt-and-pepper noise is best handled by median filtering, which removes outliers without blurring the edges.
- The ideal low-pass filter is the worst choice in both cases because of the ringing it creates.

In all cases there is a trade-off between noise removal and preservation of detail: stronger filtering removes more noise but also removes more of the image content. Linear filters (Gaussian, ideal) cannot separate noise from edges, because both are high frequencies. The median filter is non-linear and can remove isolated pixels while keeping edges.

---

## Question 19

When the original image is subsampled without smoothing (top row), it quickly degrades. After a few steps the edges become jagged ("staircase" effect), thin structures like the phone cord break into disconnected pieces, and fine regular patterns like the keys of the calculator become distorted or disappear. This is aliasing: the high frequencies of the image cannot be represented at the lower resolution and are folded back into wrong, lower frequencies.

When the image is smoothed before each subsampling step (bottom row), the result is blurrier but much more faithful. The shapes stay continuous and recognisable, and there are no jagged edges or false patterns.

Best parameters found for iteration i = 4:

- Gaussian filter: t = 1.0.
- Ideal low-pass filter: cutoff = 0.25 cycles per pixel. The ideal filter also prevents aliasing, but shows some ringing around strong edges.

---

## Question 20

Subsampling by a factor 2 halves the sampling frequency, so the highest frequency that can be represented (Nyquist frequency) goes from 0.5 to 0.25 cycles per pixel of the original image. Any content above that frequency is aliased: it reappears as false low-frequency patterns.

Smoothing before subsampling is a low-pass filter that removes these frequencies before they can alias. This is why the smoothed sequence keeps the overall structure of the image while the raw sequence shows jagged edges and broken lines.

The trade-off: too little smoothing leaves aliasing; too much smoothing makes the image unnecessarily blurry. The ideal low-pass filter matches the Nyquist limit exactly (cutoff 0.25) but introduces ringing as a side effect. The Gaussian does not have a sharp cutoff, so it removes slightly more useful detail, but it gives no ringing. In practice it is the better choice.
