#!/usr/bin/env python3

import sys
import pandas as pd

if len(sys.argv) < 2:
    print("{0} <tot.list> <query.log>\n");
    sys.exit(0);

tot_list = sys.argv[1]
lib_qlog = sys.argv[2]

tot_df = pd.read_table(tot_list, sep="\t")
lib_df = pd.read_table(lib_qlog, sep="\t", header=None)
lib_df.columns = ["Library", "Read1", "Read2"]
mer_df = tot_df.merge(lib_df, on = "Library")
rea_df = mer_df.loc[:, ["SampleName", "FlowCell", "Library", "Read1", "Read2",	"BarcodeConf"]]
for index in range(rea_df.shape[0]):
    tmp_df = rea_df.loc[index:index, :]
    tmp_df.to_csv("{0}.list".format(rea_df.loc[index, "Library"]), header=True, index=False, sep="\t")
