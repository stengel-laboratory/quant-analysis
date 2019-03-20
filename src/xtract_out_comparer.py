#!/usr/bin/env python3.6

import pandas as pd
import argparse
import os
from link_library.bag_container_library import process_bag
import link_library as ll
from statsmodels.stats import weightstats
import numpy as np

desc = """Kai Kammer - 2019-02-11. 
Script to compare pvalues and log2 from xTract and bag container script
"""

parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('input', action="store", default=None, type=str, nargs='+',
                    help="List of input csv files separated by spaces. One bag container and one xTract output")
parser.add_argument('-f', '--filter', action="store", dest="filter", default="",
                    help="Optionally specify a link type to filter for. Possible values: monolink, xlink")
parser.add_argument('-e', '--sel_exp', action="store_true", dest="sel_exp", default=False,
                    help="Optionally provide this flag to select specific experiments to plot")
parser.add_argument('-i', '--impute', action="store_true", dest="impute", default=False,
                    help="Optionally provide this flag to impute missing values for log2ratio calculations")
parser.add_argument('-o', '--output', action="store", dest="output", default='xtract_out_comparer',
                    help="Base name for the merged output file")
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
args = parser.parse_args()

xt_db = ll.xTractDB()
origin_string  = 'origin'

def compare_output(df_xtract, df_bag):
    def round_sig(x, sig=4):
        return round(x, sig - int(np.floor(np.log10(abs(x)))) - 1)
    df = pd.concat([df_bag, df_xtract], ignore_index=True, sort=True)
    df.to_csv(args.output + '_merged.csv', float_format='%.6g')
    print("Shape before removing links not found in both: {0}".format(df.shape))
    df = df.groupby([xt_db.uxid_string, xt_db.exp_string]).filter(lambda x: len(x) == 2)
    df = df.sort_values([xt_db.uxid_string, xt_db.exp_string, origin_string])
    df.to_csv(args.output + '_log2_equal.csv', float_format='%.6g')
    print("Shape after removing links not found in both: {0}".format(df.shape))
    df = df.groupby([xt_db.uxid_string, xt_db.exp_string]).filter(
        lambda x: round_sig(x[xt_db.log2_string].values[0]) == round_sig(x[xt_db.log2_string].values[1]))
    df = df.sort_values([xt_db.uxid_string, xt_db.exp_string])

    print("Shape after removing non-equal log2ratios: {0}".format(df.shape))
    print(df.groupby(origin_string)[xt_db.pval_string].mean())
    # df.groupby([xt_db.uxid_string, xt_db.exp_string]).apply(lambda x: print(x[[xt_db.uxid_string, xt_db.pval_string, origin_string]]))
    df = df.groupby([xt_db.uxid_string, xt_db.exp_string]).filter(
        lambda x: round_sig(x[xt_db.pval_string].values[0]) == round_sig(x[xt_db.pval_string].values[1]))
    print("Shape after removing non-equal pvals: {0}".format(df.shape))
    df = df.sort_values([xt_db.uxid_string, xt_db.exp_string, origin_string])
    df.to_csv(args.output + '_pval_equal.csv', float_format='%.6g')

def main():
    df_xtract = None
    df_list = []
    uid_string = "b_peptide_uID"
    for inp in args.input:
        if ".xls" in inp:
            # the xls files written by xtract are buggy and can only be read this way
            df = pd.read_csv(inp, engine='python', delimiter='\t', na_values=['-'])
        else:
            df = pd.read_csv(inp, engine='python')
        df.name = os.path.basename(inp)
        df._metadata += ['name']
        if uid_string in df:
            df_list.append(df)
        elif xt_db.uxid_string in df:
            df_xtract = df
        else:
            print("WARNING: No compatible input found for {0}".format(args.input))
            exit(1)
    if df_xtract is not None and df_list is not None:
        bag_cont = process_bag.BagContainer(level='uxID', df_list=df_list, filter=args.filter, sel_exp=args.sel_exp,
                                            impute_missing=args.impute, norm_reps=args.norm_replicates,
                                            norm_exps=args.norm_experiments, vio_list=args.vio_list)
        df_bag = ll.get_xtract_df(bag_cont, incl_tech=args.incl_tech)
        df_bag[origin_string] = 'bag'
        df_xtract[origin_string] = 'analyzer_quant'
        compare_output(df_xtract, df_bag)
        def analyze_c3_detail():
            sum_list = [bag_cont.col_exp, bag_cont.col_bio_rep, bag_cont.col_tech_rep, bag_cont.col_weight_type]
            mean_list = [bag_cont.col_exp, bag_cont.col_bio_rep, bag_cont.col_tech_rep, bag_cont.col_weight_type]
            df_tmp = bag_cont.get_group(sum_list, mean_list, bag_cont.col_area_sum_total)
            vals_e = df_tmp[(df_tmp[bag_cont.col_exp] == "c3b")][bag_cont.col_area_sum_total].values
            vals_ref = df_tmp[(df_tmp[bag_cont.col_exp] == "c3")][bag_cont.col_area_sum_total].values
            print(vals_e)
            print(vals_ref)
            print("pval pooled", weightstats.ttest_ind(vals_e, vals_ref,usevar='pooled')[1])
            print("pval unequal", weightstats.ttest_ind(vals_e, vals_ref, usevar='unequal')[1])
            print("mean", np.mean(vals_e), np.mean(vals_ref))
            print("log2", np.log2(np.mean(vals_e)/np.mean(vals_ref)))
            print("var", np.var(vals_e, ddof=1), np.var(vals_ref, ddof=1))
            # print(df_bag[(df_bag[xt_db.uxid_string] == "sp|P01024|CO3_HUMAN:1036") & (df_bag[xt_db.exp_string] == "f2_E")][xt_db.pval_string].values)
    else:
        print("WARNING: Please specify a bag container and the corresponding xTract.analyzer.quant.xls as input")
        exit(1)


if __name__ == "__main__":
    main()
