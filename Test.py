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


"""
Go Enrichment testing
"""
enrichment = HTSDataMining.EnrichmentStudy(study="/home/arnaud/Desktop/HDV/DATA/GO_enrichement/study.txt",
                                           pop="/home/arnaud/Desktop/HDV/DATA/GO_enrichement/pop8.txt",
                                           assoc="/home/arnaud/Desktop/TEMP/asso_gene_genome_goid.csv",
                                           compare=False,
                                           namespace_filter=None)
result = enrichment.to_dataframe()
df = pd.DataFrame(result)
print(df.head())
# df.to_csv(path_or_buf='/home/arnaud/Desktop/HDV/DATA/GO_enrichement/GO_enrichment_8.csv', index=False)
