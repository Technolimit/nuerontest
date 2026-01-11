import numpy as np
import matplotlib.pyplot as plt

# --- 1. Generate a sample signal ---
fs = 100.0       # Sampling frequency (Hz)
t = np.arange(0, 10, 1/fs) # Time vector (10 seconds duration)
# Signal with two frequencies (4 Hz and 7 Hz) and some noise
x = np.sin(2*np.pi*4*t) + np.sin(2*np.pi*7*t) + np.random.randn(len(t))*0.2
N = len(x)       # Number of samples

# --- 2. Compute the FFT ---
# Use np.fft.rfft for real-valued input signals for efficiency and to get only positive frequencies
yf = np.fft.rfft(x) 
# The magnitude spectrum is the absolute value of the FFT results
magnitude_spectrum = np.abs(yf) * 2.0 / N # Normalization to get true amplitude

# --- 3. Compute the corresponding frequencies ---
# np.fft.rfftfreq computes the frequencies for the rfft output
freqs = np.fft.rfftfreq(N, d=1/fs) # d is the sample spacing (1/fs)


plt.figure(figsize=(10, 5))
plt.plot(t, x)
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.title('Sample Signal (4 Hz + 7 Hz + Noise)')
# --- 4. Plot the spectrum ---
plt.figure(figsize=(10, 5))
plt.plot(freqs, magnitude_spectrum)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude')
plt.title('Amplitude Spectrum of the Signal')
plt.grid(True)
plt.show()
