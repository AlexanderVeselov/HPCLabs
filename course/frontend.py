#!/usr/bin/python

import matplotlib
from matplotlib import pyplot

import time

matplotlib.use('TkAgg')
 
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
 
import numpy as np

import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
import ctypes

class Application(tk.Tk):
    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)
        self.title("Fourier Transform")
        self.geometry("800x600")

        self.loadbutton = ttk.Button(self, text="Load Data", command=self.LoadData)
        self.loadbutton.pack()
        
        self.my_dll = ctypes.CDLL("../build/course/libcourse.so")

        self.options = {
            "1 thread, no vectorization": self.my_dll.fft,
            "1 thread, vectorized": self.my_dll.fft_simd,
            "1 thread, vectorized and data aligned": self.my_dll.fft_simd_aligned,
            "Multi-threaded, vectorized and data aligned": self.my_dll.fft_parallel_simd_aligned,
        }
        self.optionvar = tk.StringVar(self)
        self.optionvar.set(list(self.options.keys())[0])

        self.popupMenu = ttk.OptionMenu(self, self.optionvar, list(self.options.keys())[0], *self.options.keys())
        self.popupMenu.pack()

        self.processbutton = ttk.Button(self, text="Process", command=self.ProcessData)
        self.processbutton.pack()

        self.InitPlot1()


    def LoadData(self):
        self.ax[0].clear()
        self.Fs = 262144.0 * 2#2097152.0  # sampling rate
        Ts = 1.0/self.Fs # sampling interval
        t = np.arange(0, 4, Ts) # time vector
        ff = 5   # frequency of the signal
        self.y = np.sin(2*np.pi*ff*t) + np.sin(2*np.pi*2 * ff*t)# + np.random.random_sample(len(t))
        self.ax[0].plot(t, self.y)
        self.ax[0].set_xlabel('Time')
        self.ax[0].set_ylabel('Amplitude')

        self.canvas.draw()

    def ProcessData(self):
        self.ax[1].clear()
        n = len(self.y) # length of the signal
        k = np.arange(n)
        T = n/self.Fs
        frq = k/T # two sides frequency range
        frq = frq[range(int(n / 2))] # one side frequency range

        float_array_type = ctypes.c_double * n
        # Asterisk '*' means unpacking container
        in_data = float_array_type(*self.y)
        # Pass zero lists as out parameters
        out_real = float_array_type(*[0]*n)
        out_imag = float_array_type(*[0]*n)
        self.options[self.optionvar.get()](in_data, n, out_real, out_imag)

        #Y = np.fft.fft(self.y)
        Y = np.array(out_real) + 1j * np.array(out_imag)
        Y = Y[range(int(n / 2))]

        self.ax[1].semilogx(frq, 2 * abs(Y) / n, 'r')
        self.ax[1].set_xlabel('Freq (Hz)')
        self.ax[1].set_ylabel('|Y(freq)|')
        self.canvas.draw()

    def InitPlot1(self):
        fig = Figure()
        self.ax = [fig.add_subplot(211), fig.add_subplot(212)]
        self.canvas = FigureCanvasTkAgg(fig, self)
        self.canvas.get_tk_widget().pack()

window = Application()
window.mainloop()
