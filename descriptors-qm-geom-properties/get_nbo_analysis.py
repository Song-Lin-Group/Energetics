# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 10:40:20 2021

@author: illari and jonas

This script takes a set of Gaussian output files containing NBO data and extracts
calculated NBO occupancies and energies between labeled atoms in each Gaussian
file defined in 'Atomlabels'. The input files must be in the same working
directory as this script.

The results are tabulated and recorded in an output Excel workbook ('NBO analysis.xlsx')

"""

import pandas as pd
import itertools as i
import glob
import re

#READ ATOMLABELS AND LIST OF .LOG FILES
string1 = 'Please put all .log files in the same directory as this script'      #creates a list of all .log files in current directory
print(string1)                                                                  
listFileNames = []                                                              
for filename in glob.glob("*.log"):
    listFileNames.append(filename)

string2 = 'Please enter the name of the Atomlabels file without .xlsx:'               #reads atomlabels excel sheet but has to be in same folder
pathAtomLables = input(string2)   
dfAtomLabels = pd.read_excel(str(pathAtomLables)+'.xlsx',index_col=0)
                                
#Generates the column names for the the output dataframe
allvalenceorbitals = [] 
for element in list(dfAtomLabels.head()):                                                       #adds LP(1) through LP(3) to the list list of atomnames given in Atomlabels.
    allvalenceorbitals.append(element + ' LP(1)')
    allvalenceorbitals.append(element + ' LP(2)')
    allvalenceorbitals.append(element + ' LP(3)')
bondcombinationswithoutdash = list(i.combinations(list(dfAtomLabels.head()), 2))                #returns all possible 2-atom combinations and casts it to a list of str
bondcombination =  [' - '.join(i) for i in bondcombinationswithoutdash]                                                                 
for element in bondcombination:                                                 
    allvalenceorbitals.append(element + ' BD(1)')                                               #gnerated bonding and antibonding orbitals for all bonds.
    allvalenceorbitals.append(element + ' BD(2)')
    allvalenceorbitals.append(element + ' BD(3)')
    allvalenceorbitals.append(element + ' BD*(1)')
    allvalenceorbitals.append(element + ' BD*(2)')
    allvalenceorbitals.append(element + ' BD*(3)')
allNBOtitle = []
for element in allvalenceorbitals:
    allNBOtitle.append(element + ' occ')
    allNBOtitle.append(element + ' Energy')

#GET DATA OF .LOG FILES
def get_logfiles(files):
    logfile=[]
    with open(str(listFileNames[files]),'r') as f:                                              #reads the .log file in reading only mode
        flag = False                                                                            #orbitals and total lewis
        for line in f:
            if line.startswith('              Total Lewis'):                                    #iterates through all the files, name by name and only reads data between natural bond
                flag=False 
            if flag:
                logfile.append(line)
            if line.startswith(' Natural Bond Orbitals (Summary):'):
                flag=True
    return(logfile)

#NBO analysis
    
allNBOforallfiles = []                                                                              
for files in range(len(listFileNames)):                                                         #iterates through all the .log files
    NBOsummary= get_logfiles(files)                                                             #returns the NBO summary of the .log file as a list of str   
# NBO analysis list    
    allNBO=[]
#NBOs are extracted by searching the NBO Summary with a regular expressions re   
#Extraction of the LP. This sequence returns a list of list. The sublists are empty if the str is not found, sublist contains the whole line if a str is re.match for the line.  
    LP_occ_energy = []
    LPatoms = list(dfAtomLabels.iloc[files])
#re.match also matches empty spaces to account for integer with more than one digit all atom lables are filled to have 4 characters: 1 becomes '   1'
    for m in range(len(LPatoms)):
        LPatoms[m]= str(LPatoms[m])
        LPatoms[m]= LPatoms[m].rjust(4)
#searches for the appropriate LP withe by iterating through all atoms in the atom labele
    for m in range(len(list(dfAtomLabels.iloc[files]))):
        
        r1 = re.compile(".*LP [(]   1[)] \S" + LPatoms[m])    
        r2 = re.compile(".*LP [(]   2[)] \S" + LPatoms[m])  
        r3 = re.compile(".*LP [(]   3[)] \S" + LPatoms[m])
        
        LP_occ_energy.append(list(filter(r1.match, NBOsummary))) 
        LP_occ_energy.append(list(filter(r2.match, NBOsummary))) 
        LP_occ_energy.append(list(filter(r3.match, NBOsummary)))
#For empty sublists 'nan','nan is returned, for list containing the whole line the occupancy [41:48] and energy [52:60] is returned by selecting the appropriate characters in the line. 
    for entry in LP_occ_energy:
        if entry == []:
            allNBO.append('nan')
            allNBO.append('nan')
        else:
            allNBO.append([string[41:48] for string in entry][0])
            allNBO.append([string[52:60] for string in entry][0])
            
#Extraction of the bonds following the same principle as the lone pairs.           
    bonds_occ_energy = []
    atomcombinations = list(i.combinations(list(dfAtomLabels.iloc[files]), 2))
    atomcombinations = [list(x) for x in atomcombinations]
    for x in atomcombinations:
        x = x.sort()
    for n in range(len(atomcombinations)):
        for m in range(len(atomcombinations[n])):
            atomcombinations[n][m]= str(atomcombinations[n][m])
            atomcombinations[n][m]= atomcombinations[n][m].rjust(4)
    
        a1 = atomcombinations[n][0]
        a2 = atomcombinations[n][1]
        
        r1 = re.compile(".*BD [(]   1[)] \S" + a1 + " - \S" + a2)    
        r2 = re.compile(".*BD [(]   2[)] \S" + a1 + " - \S" + a2)  
        r3 = re.compile(".*BD [(]   3[)] \S" + a1 + " - \S" + a2)    
        r4 = re.compile(".*BD[*][(]   1[)] \S" + a1 + " - \S" + a2)    
        r5 = re.compile(".*BD[*][(]   2[)] \S" + a1 + " - \S" + a2)    
        r6 = re.compile(".*BD[*][(]   3[)] \S" + a1 + " - \S" + a2)
        
        bonds_occ_energy.append(list(filter(r1.match, NBOsummary))) 
        bonds_occ_energy.append(list(filter(r2.match, NBOsummary))) 
        bonds_occ_energy.append(list(filter(r3.match, NBOsummary))) 
        bonds_occ_energy.append(list(filter(r4.match, NBOsummary))) 
        bonds_occ_energy.append(list(filter(r5.match, NBOsummary))) 
        bonds_occ_energy.append(list(filter(r6.match, NBOsummary)))
        
    for entry in bonds_occ_energy:
        if entry == []:
            allNBO.append('nan')
            allNBO.append('nan')
        else:
            allNBO.append([string[41:48] for string in entry][0])
            allNBO.append([string[52:60] for string in entry][0])
            
    allNBOforallfiles.append(allNBO)
            
#Writes a excel output wiht appropriate titles and drops all the columns with no entries.        
pd.DataFrame(allNBOforallfiles, columns = allNBOtitle, index = dfAtomLabels.index).to_excel('NBO analysis.xlsx')
outputdf = pd.read_excel('NBO analysis.xlsx').dropna(how='all', axis = 1)
outputdf.to_excel('NBO analysis.xlsx')
