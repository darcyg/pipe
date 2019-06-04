#!/usr/bin/env python3
# coding: utf-8

import pandas as pd
import sys
import os

if len(sys.argv) < 3:
    print(os.path.basename(sys.argv[0]), "<indir>", "<rawin.xlsx>")
    sys.exit(1)

in_dir = os.path.abspath(sys.argv[1]) #pipeline input directory
raw_in = sys.argv[2] #raw excel input

# parse raw input excel file
df = pd.read_excel(raw_in, sheet_name="Sheet1")

finallib = [] #Library, from 终文库名称
flowcell = [] #FlowCell, from FlowCell
sampname = [] #SampleName, from 终文库名称
barcodef = [] #BarcodeConf, from in_dir + 终文库名称 + ".conf"
preanlib = [] #from 终文库名称
barcodes = [] #from Barcode
barcPref = "tsADT8N-gr."

for index, row in df.iterrows():
    if not pd.isna(row.终文库名称):
        if len(preanlib) > 0:
            tmpdf = pd.DataFrame({"PreLib": preanlib, "Barcode": barcodes})
            tmpdf.to_csv("./{0}.conf".format(finallib[-1]), sep="\t", header=False, index=False)
            preanlib = []
            barcodes = []
        finaLibName = row.终文库名称.split('_')[0]
        finallib.append(finaLibName)
        flowcell.append(row.FlowCell)
        sampname.append(finaLibName)
        barcodef.append(os.path.join(in_dir, finaLibName + ".conf"))
    preanlib.append(row.预文库名称)
    barcodes.append(barcPref + str(row.Barcode))

# write the last configure file  
tmpdf = pd.DataFrame({"PreLib": preanlib, "Barcode": barcodes})
tmpdf.to_csv("./{0}.conf".format(finallib[-1]), sep="\t", header=False, index=False)

# write total sample list
tmpdf = pd.DataFrame({"SampleName": sampname, 
                      "FlowCell": flowcell, 
                      "Library": finallib, 
                      "BarcodeConf": barcodef
                     })
tmpdf.to_csv("./raw.list", sep="\t",index=False)
