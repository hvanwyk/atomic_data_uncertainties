#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 15:49:47 2018

@author: kyle
"""

import tkinter as tk

sequences = ["li", "be"]
atoms = ["c"]

class Application(tk.Frame):
    
    def __init__(self, *args, **kwargs):
        tk.Frame.__init__(self, *args, **kwargs)
        #self.parent = parent
        self.label = tk.Label(self, text="Ionization")
        self.label.pack(side="top")
        
        self.lb_seq = tk.Listbox(master=self, exportselection=0)
        for i, thing in enumerate(sequences):
            self.lb_seq.insert(i, thing)
        self.lb_seq.pack(side="left")
        self.seq = tk.StringVar()
        self.seq.set("he")

        self.lb_seq.bind("<<ListboxSelect>>", self.set_seq)
        
        self.lb_atom = tk.Listbox(master=self, exportselection=0)
        for i, thing in enumerate(atoms):
            self.lb_atom.insert(i, thing)
        self.lb_atom.pack(side="right")
        self.atom = tk.StringVar()
        self.atom.set("he")

        self.lb_atom.bind("<<ListboxSelect>>", self.set_atom)
        
        self.killswitch = tk.Button(self, text="Go!", command=self.destroy)
        self.killswitch.pack(side="bottom")
        
    def set_seq(self, event):
        i = self.lb_seq.curselection()[0]
        text = self.lb_seq.get(i)
        self.seq.set(text)
        
    def set_atom(self, event):
        i = self.lb_atom.curselection()[0]
        text = self.lb_atom.get(i)
        self.atom.set(text)
    
app = Application()
app.pack()
app.mainloop()

seq = app.seq.get()
atom = app.atom.get()
