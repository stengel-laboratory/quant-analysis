#!/usr/bin/env python3.6

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
parser.add_argument('-l', '--level_ms1', action="store", dest="level", default='uxID',
                    help="Level on which the ms1 intensities are summed. Either uID (peptide)"
                         "or uxID (protein).")
parser.add_argument('-f', '--filter', action="store", dest="filter", default="",
                    help="Optionally specify a link type to filter for. Possible values: monolink, xlink")
parser.add_argument('-e', '--sel_exp', action="store_true", dest="sel_exp", default=False,
                    help="Optionally provide this flag to select specific experiments to plot")
parser.add_argument('-i', '--impute', action="store_true", dest="impute", default=False,
                    help="Optionally provide this flag to impute missing values for log2ratio calculations")
args = parser.parse_args()


def get_xtract_df(bag_cont):
    sum_list = [bag_cont.col_exp, bag_cont.col_bio_rep, bag_cont.col_tech_rep]
    mean_list = [bag_cont.col_exp]
    exp_ref = ll.input_log2_ref(bag_cont.exp_list)
    df_pval = bag_cont.get_two_sided_ttest(sum_list, [bag_cont.col_exp, bag_cont.col_bio_rep], ref=exp_ref)
    df_log2 = bag_cont.getlog2ratio(sum_list, mean_list, ref=exp_ref)
    # n-way merge
    df = reduce(lambda left, right: pd.merge(left, right, on=[bag_cont.col_level, bag_cont.col_exp]),
                [df_log2, df_pval])
    df = df.sort_values([bag_cont.col_exp, bag_cont.col_level])
    return df


def main():
    if ".xls" in args.input:
        # the xls files written by xtract are buggy and can only be read this way
        df = pd.read_csv(args.input, engine='python', delimiter='\t', na_values=['-'])
    else:
        df = pd.read_csv(args.input, engine='python')
    df.name = os.path.basename(args.input)
    df._metadata += ['name']
    xt_db = ll.xTractDB()
    bag_cont = process_bag.BagContainer(level=args.level, df_list=[df], filter=args.filter, sel_exp=args.sel_exp,
                                        impute_missing=args.impute)
    df = get_xtract_df(bag_cont)
    # print(df)
    print("Mean values for experiments:")
    print(df.groupby(bag_cont.col_exp).mean())
    df = df.rename(index=str, columns={bag_cont.col_level: xt_db.uxid_string, bag_cont.col_exp: xt_db.exp_string,
                                       bag_cont.col_log2ratio_ref: xt_db.exp_ref_string,
                                       bag_cont.col_log2ratio: xt_db.log2_string, bag_cont.col_pval: xt_db.pval_string,
                                       bag_cont.col_fdr: xt_db.fdr_string})
    df.to_csv(args.outname, float_format='%.6g')


if __name__ == "__main__":
    main()
