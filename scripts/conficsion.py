#!/usr/bin/env python 
# encoding: utf-8
#====================================================================
# Name        :  
# Version    : 
# Author      :  Lampros Mountrakis (L.Mountrakis@uva.nl)
# Date        : 
# Description : 
#
#
#====================================================================

import xml.etree.ElementTree as ET
from itertools import product
import sys
import argparse

def parseArguments(config):
    "Parse the arguments based on the keys from the config file"
    parser = argparse.ArgumentParser()
    #parser.add_argument('caseId', metavar='N', nargs='+', help=' help')
    for key in config.keys():
        parser.add_argument('--' + key,  nargs='+', help=key + ' help')
    args = parser.parse_args()
    dictArgs = args.__dict__
    return dictArgs


def calculate_dt(config):
    "Calculate dt based on the existing config file"
    dx = float(config['dx'].text)
    nu_p = float(config['nu_p'].text)
    tau =  float(config['tau'].text)
    return dx * dx * (1/3.*(tau-0.5)) / nu_p


parser = ET.XMLParser(encoding="utf-8")
tree = ET.parse('config_test.xml', parser=parser)
# try:
#     tree = ET.parse(sys.argv[1], parser=parser)
# except:
#     tree = ET.parse('config_test.xml', parser=parser)


config = {}
# for i in tree.iter(): 
for i in tree.findall('.//*'):
    config[i.tag] = i

arguments = parseArguments(config)
params = []
param_val = []
for key, value in arguments.items():
    if value != None:
        params += [key]
        param_val += [value]

print params, param_val

for ic, comb in enumerate( product(*param_val) ):
    caseId = []
    for iv, value in enumerate(comb):
        key = params[iv]
#        print key, value,
        config[key].text = str(value)
        caseId += [key +'-'+ str(value)]
#    print
    config['caseId'].text = "_".join(caseId)
    tree.write(config['caseId'].text+'-output.xml')





