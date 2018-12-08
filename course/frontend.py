import matplotlib
 
matplotlib.use('TkAgg')
 
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
 
import numpy as np

import tkinter as tk
from tkinter import ttk
from tkinter import filedialog

class Application(tk.Tk):
    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)
        self.title("Fourier Transform")
        self.geometry("800x600")

        self.loadbutton = ttk.Button(self, text="Load Data", command=self.LoadData)
        self.loadbutton.pack()

        self.InitPlot1()

        self.processbutton = ttk.Button(self, text="Process", command=self.ProcessData)
        self.processbutton.pack()

        self.InitPlot2()

    def LoadData(self):
        #self.data_filename = filedialog.askopenfilename(initialdir = "/",title = "Select file",filetypes = (("text", "*.txt"), ("all files", "*.*")))
        #print(self.data_filename)
        self.subplot1.clear()
        x = np.arange(0, 100)
        y = np.random.random_sample(100)
        self.subplot1.plot(x, y)
        self.subplot1.xaxis.grid(True)
        self.subplot1.set_title('Loaded Data')
        self.canvas1.draw()

    def ProcessData(self):
        #self.data_filename = filedialog.askopenfilename(initialdir = "/",title = "Select file",filetypes = (("text", "*.txt"), ("all files", "*.*")))
        #print(self.data_filename)
        self.subplot2.clear()
        x = np.arange(0, 100)
        y = np.random.random_sample(100)
        self.subplot2.plot(x, y)
        self.subplot2.xaxis.grid(True)
        self.subplot2.set_title('Processed Data')
        self.canvas2.draw()

    def InitPlot1(self):
        fig = Figure(figsize=(8, 2))
        self.subplot1 = fig.add_subplot(111)
        self.canvas1 = FigureCanvasTkAgg(fig, self)
        self.canvas1.get_tk_widget().pack()

    def InitPlot2(self):
        fig = Figure(figsize=(8, 2))
        self.subplot2 = fig.add_subplot(111)
        self.canvas2 = FigureCanvasTkAgg(fig, self)
        self.canvas2.get_tk_widget().pack()

window = Application()
window.mainloop()
