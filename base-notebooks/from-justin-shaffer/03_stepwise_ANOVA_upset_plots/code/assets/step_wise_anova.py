import os
import sys
import glob
import numpy as np
import pandas as pd
import tempfile
import subprocess

def run_commands(cmds, verbose=True):
    """
    This function is a script runner.
    """
    if verbose:
        print("Running external command line application(s). This may print "
              "messages to stdout and/or stderr.")
        print("The command(s) being run are below. These commands cannot "
              "be manually re-run as they will depend on temporary files that "
              "no longer exist.")
    for cmd in cmds:
        if verbose:
            print("\nCommand:", end=' ')
            print(" ".join(cmd), end='\n\n')
        subprocess.run(cmd, check=True)

    
def run_stepwise_anova(ord_df, mf_ord_, test_rda):
    
    # prep files
    ord_df.columns = ['PC%i'%(i+1) for i in range(len(ord_df.columns))]
    mf_ord_ = mf_ord_.reindex(ord_df.index).copy()
    mf_ord_ = mf_ord_[test_rda]

    # save all intermediate files into tmp dir
    with tempfile.TemporaryDirectory() as temp_dir_name:
        # save the tmp dir locations
        ord_fp = os.path.join(temp_dir_name, 'ord_.tsv')
        map_fp = os.path.join(temp_dir_name, 'mf_.txt')
        out_fp = os.path.join(temp_dir_name, 'output.effect.size.tsv')

        # Need to manually specify header=True for Series (i.e. "meta"). It's
        # already the default for DataFrames (i.e. "table"), but we manually
        # specify it here anyway to alleviate any potential confusion.
        ord_df.index.name = "#SampleID"
        mf_ord_.index.name = "#SampleID"
        ord_df.to_csv(ord_fp, sep='\t', header=True)
        mf_ord_.to_csv(map_fp, sep='\t', header=True)

        # build command        
        cmd = [os.path.realpath(__file__).replace('step_wise_anova.py','') \
               + 'stepwise-rda.R',
               ord_fp,
               map_fp,
               out_fp]
        cmd = list(map(str, cmd))

        try:
            run_commands([cmd])
        except subprocess.CalledProcessError as e:
            raise Exception("An error was encountered while running"
                            " in R (return code %d), please inspect stdout"
                            " and stderr to learn more." % e.returncode)

        # if run was successful, import the data and return
        effect_size_df = pd.read_csv(out_fp, index_col=0, sep='\t')
        effect_size_df['R2.adj'] = [effect_size_df['R2.adj'].values[0]] \
                                    + list(effect_size_df['R2.adj'].diff().values[1:])
    return effect_size_df
