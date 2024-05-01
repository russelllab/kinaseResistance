#!/usr/bin/env python3

import os
import ML_multilabel
import mlflow
import mlflow.sklearn
from mlflow.models import infer_signature

'''
for name in ['AIvNLD', 'AIvLD', 'AvNL', 'AvL', 'LDvNAI', 'LvNA', 'RvN']:
    for max_depth in [3, 4, 5]:
        for min_samples_split in [3, 4, 5 ,7, 10, 12]:
            for min_samples_leaf in [3, 4, 5 ,7, 10, 12]:
                for n_estimators in [100]:
                    ML.main([max_depth],[min_samples_split],[min_samples_leaf], [n_estimators], name=name)
                    # os.system("python3 callML.py {} {} {} {}".format(max_depth, min_samples_split, min_samples_leaf, n_estimators))
'''
# mlflow.set_tracking_uri(uri="http://127.0.0.1:8080")

# for algo in ['GNB', 'MLP', 'RF', 'XGB', 'SVC']:
for algo in ['XGB']:
    # for Salzberg in ['False', 'True']:
    for Salzberg in ['False']:
    # for Salzberg in ['True']:
        if Salzberg == 'True':
            mlflow.set_experiment(algo+'_Salzberg')
        else:
            mlflow.set_experiment(algo)
        name = 'all'
        with mlflow.start_run(run_name=name) as run:
            if algo == 'GNB':
                var_smoothing = [1e-09, 1e-08, 1e-07, 1e-06]
                ML_multilabel.main(name=name, algo=algo, var_smoothing=var_smoothing, Salzberg=Salzberg)
            elif algo == 'MLP':
                hidden_layer_sizes=[(10,), (10,10,), (10,10,10,)]
                activation=['relu']
                solver=['adam']
                alpha=[0.0001]
                batch_size=['auto']
                learning_rate=['constant']
                learning_rate_init=[0.001]
                max_iter=[200]
                shuffle=[True]
                tol=[0.0001]
                ML_multilabel.main(name=name, algo=algo,
                            hidden_layer_sizes=hidden_layer_sizes,
                            activation=activation,
                            solver=solver,
                            alpha=alpha,
                            batch_size=batch_size,
                            learning_rate=learning_rate,
                            learning_rate_init=learning_rate_init,
                            max_iter=max_iter,
                            shuffle=shuffle,
                            tol=tol,
                            Salzberg=Salzberg
                            )
            elif algo in ['RF', 'XGB']:
                max_depth = [3, 5, 7, 10]
                # max_depth = [10]
                min_samples_split = [3, 5 ,7, 10]
                # min_samples_split = [3]
                min_samples_leaf = [3, 5 ,7, 10]
                # min_samples_leaf = [10]
                n_estimators = [100]
                ML_multilabel.main(name=name, algo=algo,
                                #    model_filename=name,
                                #    scaler_filename=name,
                                   max_depth=max_depth,
                                   min_samples_split=min_samples_split,
                                   min_samples_leaf=min_samples_leaf,
                                   n_estimators=n_estimators,
                                   Salzberg=Salzberg)
            elif algo == 'SVC':
                C = [0.001, 0.01, 0.1, 1.0, 10.0]
                kernel = ['linear']
                ML_multilabel.main(name=name, algo=algo, C=C, kernel=kernel, Salzberg=Salzberg)
            # os.system("python3 callML.py {} {} {} {}".format(max_depth, min_samples_split, min_samples_leaf, n_estimators))
                
