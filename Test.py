#!/usr/bin/env python3
# encoding: utf-8
"""
For testing module in actual dev
"""
import numpy as np
import pandas as pd
import HTSDataMining
import logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(levelname)-8s %(message)s',
                    datefmt='%m/%d/%Y %I:%M:%S %p')
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
np.set_printoptions(linewidth=300)
np.set_printoptions(suppress=True, precision=4)


def go():
    """
    Go Enrichment testing
    """
    enrichment = HTSDataMining.EnrichmentStudy(study="/home/arnaud/Desktop/TEMP/study.txt",
                                               pop="/home/arnaud/Desktop/TEMP/pop.txt",
                                               assoc="/home/arnaud/Desktop/TEMP/assoc.csv",
                                               compare=False)
    result = enrichment.to_dataframe()
    print(pd.DataFrame(result).head())

go()