# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 00:30:10 2020

@author: jonas

Script for extraction of Hirshfeld surface analysis parameters from CrystalExplorer 17.5.
This script must be in the same directory as the Hirshfeld output data ('Hirshfeld_Example').

"""

import os,re,glob,subprocess
import numpy as np
import statistics
from pathlib import Path
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


npa_pattern = re.compile("Summary of Natural Population Analysis:")

#
# FUNCTIONS SHARED BETWEEN JOBTYPES
#

#function for human and natural sorting 
def atoi(text):                                                                             
    return int(text) if text.isdigit() else text
#this tool implements human/natural sorting
def natural_keys(text):                                             
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

#function selects all txt files in a folder named "Hirshfeld filed" in the working directory
def get_texthirshfeldfiles():
    txt=[]
    for file in glob.glob("*.txt"):
        txt.append(file)  
    txt.sort(key=natural_keys)
    return(txt)
    
#This function provides a table with the Inside	Outside Atom analysis. The output of the function is a list with the sorted inside outside atom analysis in the following order O,N,F,Cl,S,C,H. The function is by no means very efficient and is geard toward easy interpretation of the code..
def Insideoutsideatom(Hirshfeldtext): 
#the Hirshfeldtext is a panda dataframe and every line is searched for the "Inside Outside Atom"
    for i in range(len(Hirshfeldtext)):
        if str(Hirshfeldtext[0][i]) == "Inside	Outside Atom":
            Listoutsideatoms = []
#The line containing the "Atom O N C H" is split at every empty space. 
            for k in range(1,len(re.findall(r'\S+', str(Hirshfeldtext[0][i+1])))+1):
                Listinsideatoms = []
#This appends the elements of the table in the text file output of the CrystalExplorer suface information by the Inside Atoms. Then every element is appended to the inside atom List. This generates a list for every inside atom ( the rows in the original table). The 
                for n in range(0,len(re.findall(r'\S+', str(Hirshfeldtext[0][i+1])))):
                    re.findall(r'\S+', str(Hirshfeldtext[0][i+k]))[n]
                    Listinsideatoms.append(re.findall(r'\S+', str(Hirshfeldtext[0][i+k]))[n])
#All the list are appended generating a table that is similar to the original table in the text file, but allows us to call desired fuctions
                Listoutsideatoms.append(Listinsideatoms)
#the rest of the function generates a list extracting useful information from the table
    num = range(len(Listoutsideatoms[0]))
    atomtypes = list(Listoutsideatoms[0])
    Values = []
#this is the full list of atoms that are selected for in the suface analysis extraction. If a interaction is not present a empty field is appended to a list. This makes generating a nice dataset table a lot easier. 
    atoms = ['O','N','F','CL','S','C','H']
#the rest is a sequence to call all the elements in the table in the desired order and to append the empty entries to the list. There are probaly a lot more if loops than nessessary (The code is fast enough to not make problems for small to medium datasets)    
    for atom in range(len(atoms)):
        for i in num:
            
            if str(Listoutsideatoms[0][i]) == str(atoms[atom]):
                
                if "O" not in atomtypes:
                    Values.append("")        
                else:
                    for k in num:
                        if str(Listoutsideatoms[k][0]) == "O":
                            Values.append(Listoutsideatoms[k][i])
                            
                if "N" not in atomtypes:
                    Values.append("")        
                else:
                    for k in num:
                        if str(Listoutsideatoms[k][0]) == "N":
                            Values.append(Listoutsideatoms[k][i])
                        
                if "F" not in atomtypes:
                    Values.append("")        
                else:        
                    for k in num:
                        if str(Listoutsideatoms[k][0]) == "F":
                            Values.append(Listoutsideatoms[k][i])        
                 
                if "CL" not in atomtypes:
                    Values.append("")        
                else:        
                    for k in num:
                        if str(Listoutsideatoms[k][0]) == "CL":
                            Values.append(Listoutsideatoms[k][i])
                
                if "S" not in atomtypes:
                    Values.append("")        
                else:        
                    for k in num:
                        if str(Listoutsideatoms[k][0]) == "S":
                            Values.append(Listoutsideatoms[k][i])
                       
                if "C" not in atomtypes:
                    Values.append("")        
                else:        
                    for k in num:
                        if str(Listoutsideatoms[k][0]) == "C":
                            Values.append(Listoutsideatoms[k][i])
                            
                if "H" not in atomtypes:
                    Values.append("")        
                else:
                    for k in num:
                        if str(Listoutsideatoms[k][0]) == "H":
                            Values.append(Listoutsideatoms[k][i])
                        
        if str(atoms[atom]) not in atomtypes:
            Values.append("")
            Values.append("")
            Values.append("")
            Values.append("")
            Values.append("")
            Values.append("")
            Values.append("")
    return(Values)



def generalsufaceinformation(Hirshfeldtext):
#There should not be any variation in the Hirshfeld information text file structure up to this point. Hence a static extraction of information form the desired lines is used. 
    for i in range(len(Hirshfeldtext)):
        generalsufaceinformation=[]
        for k in range(5,10):
                generalsufaceinformation.append(re.findall(r'\S+', str(Hirshfeldtext[0][k]))[1])
    return(generalsufaceinformation)



def surfacepropertyinformation(Hirshfeldtext):
#There should not be any variation in the Hirshfeld information text file structure up to this point. Hence a static extraction of information form the desired lines is used. 
    for i in range(len(Hirshfeldtext)):
        surfacepropertyinformation=[]
        for k in range(15,19):
                surfacepropertyinformation.append(re.findall(r'\S+', str(Hirshfeldtext[0][k]))[2])
                surfacepropertyinformation.append(re.findall(r'\S+', str(Hirshfeldtext[0][k]))[3])
                surfacepropertyinformation.append(re.findall(r'\S+', str(Hirshfeldtext[0][k]))[4])
        surfacepropertyinformation.append(re.findall(r'\S+', str(Hirshfeldtext[0][19]))[1])
        surfacepropertyinformation.append(re.findall(r'\S+', str(Hirshfeldtext[0][19]))[2])
        surfacepropertyinformation.append(re.findall(r'\S+', str(Hirshfeldtext[0][19]))[3])
    return(surfacepropertyinformation)


print("Save the suface information from CrystalExplorer as a text file into the working directory. The output is a Hirshfeld data.xlsx Excel file.")

#empty list where each entry will be a list, which will become a row in the final output dataframe
hirshfeldanalysisforalltextfiles = []
#genertes a list of the names of all the text files in the working directory
alltextfiles = get_texthirshfeldfiles()
#iterates through all the text files in the working directory
for text in range(len(alltextfiles)):
#This generates a panda dataframe form the text file. This is allows easy line by line analysis of the files.             
    Hirshfeldtext = pd.read_csv((alltextfiles[text]), header=None)
#list with one entry, which is the name of the file
    hirshfeldfortext = [str((alltextfiles[text]))]
    hirshfeldanalysis = []
#the next three lines append the output lists of the functions to the hirshfeldanalysis resulting in a list with lenght 69
    hirshfeldanalysis.append(generalsufaceinformation(Hirshfeldtext)) 
    hirshfeldanalysis.append(surfacepropertyinformation(Hirshfeldtext)) 
    hirshfeldanalysis.append(Insideoutsideatom(Hirshfeldtext))
#this way of joing lists only works if all of all lists are floats
    hirshfeldanalysis = sum(hirshfeldanalysis, [])
#extend is a way to join lists if they contain floats and str
    hirshfeldfortext.extend(hirshfeldanalysis)
    hirshfeldanalysisforalltextfiles.append(hirshfeldfortext)

#converts the list of list into a dataframe with the columns labled
finalhirshfeld= pd.DataFrame(hirshfeldanalysisforalltextfiles, columns= ["Filename","Isovalue","Volume","Area"," Globularity"," Asphericity","di min","di mean","di max","de min","de mean","de max","dnorm min","dnorm mean","dnorm max","shape index min","shape index mean","shape index max","cuvedness min","curvedness mean","curvedness max","OO","ON","OF","OCL","OS","OC","OH","NO","NN","NF","NCL","NS","NC","NH","FO","FN","FF","FCL","FS","FC","FH","CLO","CLN","CLF","CLCL","CLS","CLC","CLH","SO","SN","SF","SCL","SS","SC","SH","CO","CN","CF","CCL","CS","CC","CH","HO","HN","HF","HCL","HS","HC","HH"])
finalhirshfeld.to_excel('Hirshfeld data.xlsx')