#!/usr/bin/env python

import pandas as pd
import argparse
import os
from link_library.bag_container_library import process_bag
from functools import reduce
import link_library as ll

desc = """Kai Kammer - 2019-01-17. 
Script to write xTract like output from a bag container details file
"""

parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('input', action="store",
                    help="Name of the input file")
parser.add_argument('-o', '--outname', action="store", dest="outname", default='xtract_out_from_bagcontainer.csv',
                    help="Name for the output figure")
parser.add_argument('-f', '--filter', action="store", dest="filter", default="",
                    help="Optionally specify a link type to filter for. Possible values: monolink, xlink")
parser.add_argument('-e', '--sel_exp', action="store_true", dest="sel_exp", default=False,
                    help="Optionally provide this flag to select specific experiments to plot")
parser.add_argument('-i', '--impute', action="store_true", dest="impute", default=False,
                    help="Optionally provide this flag to impute missing values for log2ratio calculations")
parser.add_argument('-t', '--incl_tech', action="store_true", dest="incl_tech", default=False,
                    help="Optionally include technical replicates into p-value calculation")
parser.add_argument('-nr', '--norm_replicates', action="store_true", dest="norm_replicates", default=False,
                    help="Optionally normalize replicates to their mean experimental ms1 area")
parser.add_argument('-ne', '--norm_experiments', action="store", dest="norm_experiments", default="yes",
                    help="Optionally select the experiment normalization method (or whether to normalize at all). "
                         "Possible values: yes (default norm method), xt (xTract norm method), no (do not normalize)")
parser.add_argument('-v' '--vio_list', action="store", dest='vio_list', default=['lh', 'xt'], type=str, nargs='+',
                    help="List of input possible violation filters separated by spaces. "
                         "Possible values: lh (light/heavy log2 ratio, xt (xTract type violations), none (no filtering")
parser.add_argument('-w', '--whitelist', action="store", dest="whitelist", default="",
                    help="Optionally specify a file containing allowed links (uxids), i.e. a whitelist.")
args = parser.parse_args()



def main():
    xt_db = ll.xTractDB()
    df_whitelist = None
    if ".xls" in args.input:
        # the xls files written by xtract are buggy and can only be read this way
        df = pd.read_csv(args.input, engine='python', delimiter='\t', na_values=['-'])
    else:
        df = pd.read_csv(args.input, engine='python')
    df.name = os.path.basename(args.input)
    df._metadata += ['name']
    if args.whitelist:
        df_whitelist = pd.read_csv(args.whitelist, engine='python')
    bag_cont = process_bag.BagContainer(level='uxID', df_list=[df], filter=args.filter, sel_exp=args.sel_exp,
                                        impute_missing=args.impute, norm_reps=args.norm_replicates,
                                        norm_exps=args.norm_experiments, vio_list=args.vio_list, whitelist=df_whitelist)
    df = ll.get_xtract_df(bag_cont, incl_tech=args.incl_tech)
    print("Mean values for experiments:")
    print(df.groupby(xt_db.exp_string).mean())
    df.to_csv(args.outname, float_format='%.6g', index=False)
    print("Results written to {0}".format(args.outname))

if __name__ == "__main__":
    main()
