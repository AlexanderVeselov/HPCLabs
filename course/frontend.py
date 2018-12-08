import matplotlib
from matplotlib import pyplot
 
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
        self.processbutton = ttk.Button(self, text="Process", command=self.ProcessData)
        self.processbutton.pack()

        self.InitPlot1()

        self.my_dll = ctypes.WinDLL("Debug/course.dll")

    def LoadData(self):
        #self.data_filename = filedialog.askopenfilename(initialdir = "/",title = "Select file",filetypes = (("text", "*.txt"), ("all files", "*.*")))
        #print(self.data_filename)
        self.ax[0].clear()
        self.Fs = 1024.0*2  # sampling rate
        Ts = 1.0/self.Fs # sampling interval
        t = np.arange(0, 1, Ts) # time vector
        ff = 5   # frequency of the signal
        self.y = np.sin(2*np.pi*ff*t) + np.sin(2*np.pi*2 * ff*t) + np.random.random_sample(len(t))
        self.ax[0].plot(t, self.y)
        self.ax[0].set_xlabel('Time')
        self.ax[0].set_ylabel('Amplitude')

        self.canvas.draw()

    def SendFloatArray(self, arr):
        float_array_type = ctypes.c_float * len(arr)
        float_array = float_array_type(*arr)
        self.my_dll.SendFunc(float_array, len(arr))

    def RecvFloatArray(self):
        arr = np.arange(0, 100)
        float_array_type = ctypes.c_float * len(arr)
        float_array = float_array_type(*arr)
        self.my_dll.RecvFunc(float_array)
        print(list(float_array))

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
        self.my_dll.fft(in_data, n, out_real, out_imag)

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
