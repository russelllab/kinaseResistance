#!/usr/bin/env python3

import os
import ML_multilabel
'''
for name in ['AIvNLD', 'AIvLD', 'AvNL', 'AvL', 'LDvNAI', 'LvNA', 'RvN']:
    for max_depth in [3, 4, 5]:
        for min_samples_split in [3, 4, 5 ,7, 10, 12]:
            for min_samples_leaf in [3, 4, 5 ,7, 10, 12]:
                for n_estimators in [100]:
                    ML.main([max_depth],[min_samples_split],[min_samples_leaf], [n_estimators], name=name)
                    # os.system("python3 callML.py {} {} {} {}".format(max_depth, min_samples_split, min_samples_leaf, n_estimators))
'''

max_depth = [3, 4, 5, 7, 10]
min_samples_split = [3, 5 ,7, 10, 12, 15]
min_samples_leaf = [3, 5 ,7, 10, 12, 15]
n_estimators = [50, 100, 200]
ML_multilabel.main(max_depth,min_samples_split,min_samples_leaf, n_estimators, name='all')
# os.system("python3 callML.py {} {} {} {}".format(max_depth, min_samples_split, min_samples_leaf, n_estimators))
        
