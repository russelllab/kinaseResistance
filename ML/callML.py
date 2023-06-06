#!/usr/bin/env python3

import os
import ML

for max_depth in [3, 5]:
    for min_samples_split in [3,4,5]:
        for min_samples_leaf in [3,4,5]:
            for n_estimators in [10, 50, 100]:
                ML.main([max_depth],[min_samples_split],[min_samples_leaf], [n_estimators])
                # os.system("python3 callML.py {} {} {} {}".format(max_depth, min_samples_split, min_samples_leaf, n_estimators))
        