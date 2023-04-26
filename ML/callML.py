import os
import ML

for max_depth in [2,5,7,10,15, None]:
    for min_samples_split in [2,5,7,10,15]:
        for min_samples_leaf in [2,5,7,10,15]:
            for n_estimators in [50, 100, 200, 300, 500]:
                ML.main([max_depth],[min_samples_split],[min_samples_leaf], [n_estimators])
                # os.system("python3 callML.py {} {} {} {}".format(max_depth, min_samples_split, min_samples_leaf, n_estimators))
        