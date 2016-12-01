# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import os
import sys
import re

fOutName = "metrics.dat"

def getFiles(directory):
    files = []
    for file in os.listdir(directory):
        if file.endswith(".out"):
            files.append(file)
            
    return files

def searchDataAfter(searchStr, line):
    match = re.search(searchStr+'([0-9.]+)', line)
    return match
    
    
def processOut(outFileName):
    fOutLocal = open(outFileName, 'r')
    lines = fOutLocal.readlines()
    fOutLocal.close()
    
    matched_data = []
    
    for i, line in enumerate(lines):
        dataStr = ""
        
        # Seartch for iteration
        match = searchDataAfter(r'Iteration:', line)
        if match:
            dataStr = dataStr + match.group(1) + ' '
        
            # search for performance data
            match = searchDataAfter(r'Wall time / iter. = ', line)
            if match:
                dataStr = dataStr + match.group(1) + ' '
            else:
                dataStr = dataStr + 'NA '
            
            # search for force data
            match = searchDataAfter(r'Last largest force \[lm\*lu\/lt\^2\] = ', line)
            if match:
                dataStr = dataStr + match.group(1) + ' '
            else:
                dataStr = dataStr + 'NA '
                
            # search for velocity data
            match = searchDataAfter(r'Mean velocity: ', lines[i+1])
            if match:
                dataStr = dataStr + match.group(1) + ' '
            else:
                dataStr = dataStr + 'NA '
                
            # search for velocity data
            match = searchDataAfter(r'Apparent rel. viscosity: ', lines[i+1])
            if match:
                dataStr = dataStr + match.group(1) + ' '
            else:
                dataStr = dataStr + 'NA '
                
        if(len(dataStr) > 0):
            matched_data.append(dataStr+ '\n')
        
    return matched_data

def process(directory):
    files = getFiles(directory)
    print('Files: ', files)
    
    if len(files) == 0:
        print("No files to process.")
        sys.exit(0)
        
    fOut = open(directory + '/' + fOutName, 'w')
    
    for outFile in files:
        fOut.writelines(processOut(directory + '/' + outFile))
    fOut.close()
        
    
if __name__ == "__main__":
    
    process('.')