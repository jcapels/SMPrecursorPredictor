import optuna

from deepmol.compound_featurization import MorganFingerprint, AtomPairFingerprint, LayeredFingerprint, RDKFingerprint, \
    MixedFeaturizer, MACCSkeysFingerprint
from deepmol.feature_selection import SelectFromModelFS, LowVarianceFS
from deepmol.pipeline_optimization._feature_selector_objectives import _get_feature_selector
from deepmol.base import PassThroughTransformer
from deepmol.pipeline_optimization._scaler_objectives import _get_scaler
from deepmol.pipeline_optimization._featurizer_objectives import _get_featurizer
from deepmol.pipeline_optimization._standardizer_objectives import _get_standardizer
from deepmol.models import SklearnModel
from optuna import trial
from sklearn.linear_model import RidgeClassifier, RidgeClassifierCV
from sklearn.metrics import f1_score
from sklearn.neural_network import MLPClassifier
from sklearn.neighbors import KNeighborsClassifier, RadiusNeighborsClassifier
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier
from sklearn.tree import ExtraTreeClassifier, DecisionTreeClassifier

from deepmol.loaders import CSVLoader
from deepmol.scalers import StandardScaler
from deepmol.splitters import MultiTaskStratifiedSplitter
from rdkit import RDLogger
import logging
import warnings


def get_hyperparameters_for_models(model, trial):
    if model == "DecisionTreeClassifier":
        model = DecisionTreeClassifier
        criterion = trial.suggest_categorical("model__criterion", ["gini", "entropy"])
        max_depth = trial.suggest_int("model__max_depth", 10, 100)
        min_samples_split = trial.suggest_int("model__min_samples_split", 2, 10)
        min_samples_leaf = trial.suggest_int("model__min_samples_leaf", 1, 10)
        model = model(criterion=criterion, max_depth=max_depth, min_samples_split=min_samples_split,
                      min_samples_leaf=min_samples_leaf)
        return SklearnModel(model, model_dir="model")
    elif model == "ExtraTreeClassifier":
        model = ExtraTreeClassifier
        criterion = trial.suggest_categorical("model__criterion", ["gini", "entropy"])
        min_samples_split = trial.suggest_int("model__min_samples_split", 2, 10)
        min_samples_leaf = trial.suggest_int("model__min_samples_leaf", 1, 10)
        model = model(criterion=criterion, min_samples_split=min_samples_split,
                      min_samples_leaf=min_samples_leaf)
        return SklearnModel(model, model_dir="model")
    elif model == "ExtraTreesClassifier":
        model = ExtraTreesClassifier
        n_estimators = trial.suggest_int("model__n_estimators", 10, 200)
        criterion = trial.suggest_categorical("model__criterion", ["gini", "entropy"])
        min_samples_split = trial.suggest_int("model__min_samples_split", 2, 10)
        min_samples_leaf = trial.suggest_int("model__min_samples_leaf", 1, 10)
        model = model(criterion=criterion, min_samples_split=min_samples_split,
                      min_samples_leaf=min_samples_leaf, n_estimators=n_estimators)
        return SklearnModel(model, model_dir="model")
    elif model == "KNeighborsClassifier":
        model = KNeighborsClassifier
        n_neighbors = trial.suggest_int("model__n_neighbors", 1, 100)
        weights = trial.suggest_categorical("model__weights", ["uniform", "distance"])
        algorithm = trial.suggest_categorical("model__algorithm", ["auto", "ball_tree", "kd_tree", "brute"])
        model = model(n_neighbors=n_neighbors, weights=weights, algorithm=algorithm)
        return SklearnModel(model, model_dir="model")
    elif model == "MLPClassifier":
        model = MLPClassifier
        hidden_layer_sizes = trial.suggest_categorical("model__hidden_layer_sizes",
                                                       [(100,), (100, 100), (100, 100, 100)])
        activation = trial.suggest_categorical("model__activation", ["identity", "logistic", "tanh", "relu"])
        solver = trial.suggest_categorical("model__solver", ["lbfgs", "sgd", "adam"])
        batch_size = trial.suggest_categorical("model__batch_size", [16, 32, 64, 128, 256])
        early_stopping = True
        if solver == "adam" or solver == "sgd":
            learning_rate_init = trial.suggest_categorical("model__learning_rate_ini",
                                                           [0.00001, 0.0001, 0.001, 0.01, 0.1])
            learning_rate = trial.suggest_categorical("model__learning_rate", ["constant", "invscaling", "adaptive"])

            model = model(hidden_layer_sizes=hidden_layer_sizes, activation=activation, solver=solver,
                          batch_size=batch_size,
                          learning_rate_init=learning_rate_init, learning_rate=learning_rate,
                          early_stopping=early_stopping)
        else:
            model = model(hidden_layer_sizes=hidden_layer_sizes, activation=activation, solver=solver,
                          batch_size=batch_size,
                          early_stopping=early_stopping)

        return SklearnModel(model, model_dir="model")
    elif model == "RadiusNeighborsClassifier":
        model = RadiusNeighborsClassifier
        radius = trial.suggest_int("model__radius", 1, 5)
        weights = trial.suggest_categorical("model__weights", ["uniform", "distance"])
        algorithm = trial.suggest_categorical("model__algorithm", ["auto", "ball_tree", "kd_tree", "brute"])
        model = model(radius=radius, weights=weights, algorithm=algorithm)
        return SklearnModel(model, model_dir="model")
    elif model == "RandomForestClassifier":
        model = RandomForestClassifier
        criterion = trial.suggest_categorical("model__criterion", ["gini", "entropy"])
        min_samples_split = trial.suggest_int("model__min_samples_split", 2, 10)
        min_samples_leaf = trial.suggest_int("model__min_samples_leaf", 1, 10)
        n_estimators = trial.suggest_int("model__n_estimators", 10, 1000)
        model = model(criterion=criterion, min_samples_split=min_samples_split,
                      min_samples_leaf=min_samples_leaf, n_estimators=n_estimators)
        return SklearnModel(model, model_dir="model")
    elif model == "RidgeClassifier":
        model = RidgeClassifier
        alpha = trial.suggest_float("model__alpha", 0.1, 10)
        solver = trial.suggest_categorical("model__solver",
                                           ["auto", "svd", "cholesky", "lsqr", "sparse_cg", "sag", "saga"])
        tol = trial.suggest_categorical("model__tol", [0.0001, 0.001, 0.01, 0.1])
        fit_intercept = trial.suggest_categorical("model__fit_intercept", [True, False])
        model = model(alpha=alpha, solver=solver, tol=tol, fit_intercept=fit_intercept)
        return SklearnModel(model, model_dir="model")


def feature_selector(trial):
    feature_selector = trial.suggest_categorical("feature_selector", ["None",
                                                                      "VarianceThresholdFS"])

    if feature_selector == "None":
        return PassThroughTransformer()

    elif feature_selector == "VarianceThresholdFS":
        return LowVarianceFS(threshold=trial.suggest_uniform("threshold", 0, 0.15))


def choose_1d_featurizer(trial, feat):
    if feat == 'morgan':
        radius = trial.suggest_int('radius', 2, 6, step=2)
        n_bits = trial.suggest_int('n_bits', 1024, 2048, step=1024)
        chiral = trial.suggest_categorical('chiral', [True, False])
        return MorganFingerprint(radius=radius, size=n_bits, chiral=chiral)
    elif feat == 'atom_pair':
        nBits = trial.suggest_int('nBits', 1024, 2048, step=1024)
        minLength = trial.suggest_int('minLength', 1, 3)
        maxLength = trial.suggest_int('maxLength', 20, 50, step=10)
        includeChirality = trial.suggest_categorical('includeChirality', [True, False])
        return AtomPairFingerprint(nBits=nBits, minLength=minLength, maxLength=maxLength,
                                   includeChirality=includeChirality)
    elif feat == 'layered':
        fpSize = trial.suggest_int('fpSize', 1024, 2048, step=1024)
        minPath = trial.suggest_int('minPath', 1, 3)
        maxPath = trial.suggest_int('maxPath', 5, 10)
        return LayeredFingerprint(fpSize=fpSize, minPath=minPath, maxPath=maxPath)
    elif feat == 'rdk':
        fpSize = trial.suggest_int('fpSize', 1024, 2048, step=1024)
        min_path = trial.suggest_int('min_path', 1, 3)
        max_path = trial.suggest_int('max_path', 5, 10)
        return RDKFingerprint(fpSize=fpSize, minPath=min_path, maxPath=max_path)
    elif feat == 'maccs':
        return MACCSkeysFingerprint()


def get_featurizer(trial):
    feat = trial.suggest_categorical('1D_featurizer', ['morgan', 'atom_pair', 'layered', 'rdk', 'maccs', 'mixed'])

    if feat == 'mixed':
        available_feats = ['morgan', 'atom_pair', 'layered', 'rdk', 'maccs']
        f1 = trial.suggest_categorical('f1', available_feats)
        f2 = trial.suggest_categorical('f2', available_feats)

        return MixedFeaturizer([choose_1d_featurizer(trial, f1), choose_1d_featurizer(trial, f2)])
    else:
        return choose_1d_featurizer(trial, feat)


def objective(trial):
    models = ["DecisionTreeClassifier", "ExtraTreeClassifier", "ExtraTreesClassifier", "KNeighborsClassifier",
              "RadiusNeighborsClassifier", "RandomForestClassifier", "RidgeClassifier"]
    model = trial.suggest_categorical('model', models)
    model = get_hyperparameters_for_models(model, trial)
    standardizer = _get_standardizer(trial)
    featurizer = get_featurizer(trial)
    feature_selection = feature_selector(trial)
    final_steps = [('standardizer', standardizer), ('featurizer', featurizer),
                   ('feature_selector', feature_selection), ('model', model)]
    return final_steps


if __name__ == "__main__":
    # final_dataset = CSVLoader("/home/bisbii/Desktop/DeepMol/examples/final_dataset_top.csv",
    #                   labels_fields=['C00073', 'C00078', 'C00079', 'C00082', 'C00235', 'C00341',
    #                                  'C00353', 'C00448', 'C01789', 'C03506',
    #                                  'C00047', 'C00108', 'C00187', 'C00148', 'C00041',
    #                                  'C00129', 'C00062', 'C01852', 'C00049', 'C00135'],
    #                   id_field="ids", smiles_field="smiles").create_dataset()
    #
    # train, valid, test = MultiTaskStratifiedSplitter().train_valid_test_split(final_dataset,
    #                                                                           frac_train=0.7, frac_valid=0.2,
    #                                                                           frac_test=0.1)
    #
    # train.to_csv("/home/bisbii/Desktop/DeepMol/examples/test_case_study/train.csv")
    # valid.to_csv("/home/bisbii/Desktop/DeepMol/examples/test_case_study/valid.csv")
    # test.to_csv("/home/bisbii/Desktop/DeepMol/examples/test_case_study/test.csv")

    warnings.filterwarnings("ignore")
    logger = logging.getLogger()
    logger.setLevel(logging.CRITICAL)
    RDLogger.DisableLog('rdApp.*')

    train = CSVLoader("/home/bisbii/Desktop/DeepMol/examples/test_case_study/train.csv",
                      labels_fields=['C00073', 'C00078', 'C00079', 'C00082', 'C00235', 'C00341',
                                     'C00353', 'C00448', 'C01789', 'C03506',
                                     'C00047', 'C00108', 'C00187', 'C00148', 'C00041',
                                     'C00129', 'C00062', 'C01852', 'C00049', 'C00135'],
                      id_field="ids", smiles_field="smiles").create_dataset()
    valid = CSVLoader("/home/bisbii/Desktop/DeepMol/examples/test_case_study/valid.csv",
                     labels_fields=['C00073', 'C00078', 'C00079', 'C00082', 'C00235', 'C00341',
                                    'C00353', 'C00448', 'C01789', 'C03506',
                                    'C00047', 'C00108', 'C00187', 'C00148', 'C00041',
                                    'C00129', 'C00062', 'C01852', 'C00049', 'C00135'],
                     id_field="ids", smiles_field="smiles").create_dataset()

    from deepmol.pipeline_optimization import PipelineOptimization
    from deepmol.metrics import Metric

    evo_sampler = optuna.samplers.NSGAIISampler(population_size=10, mutation_prob=0.1,
                                                crossover=None, crossover_prob=0.9,
                                                swapping_prob=0.5, seed=None, constraints_func=None)

    motpe = optuna.samplers.MOTPESampler(consider_prior=True, prior_weight=1.0, consider_magic_clip=True,
                                         consider_endpoints=True, n_startup_trials=10, n_ehvi_candidates=24)

    random_sampler = optuna.samplers.RandomSampler(seed=123)

    tpe_sampler = optuna.samplers.TPESampler()

    cmaes_sampler = optuna.samplers.CmaEsSampler()

    quasi_mote_carlo = optuna.samplers.QMCSampler()


    def f1_score_macro(y_true, y_pred):
        return f1_score(y_true, y_pred, average='macro')


    metric = Metric(f1_score_macro)

    po = PipelineOptimization(direction='maximize',
                              study_name='sm_predictor_pipeline_evolutionary',
                              storage='sqlite:///test_sm_predictor.db',
                              sampler=evo_sampler
                              )

    po.optimize(train_dataset=train, test_dataset=valid, objective_steps=objective,
                metric=metric, n_trials=500, save_top_n=3)

    po = PipelineOptimization(direction='maximize',
                              study_name='sm_predictor_pipeline_motpe',
                              storage='sqlite:///test_sm_predictor.db',
                              sampler=motpe
                              )
    po.optimize(train_dataset=train, test_dataset=valid, objective_steps=objective,
                metric=metric, n_trials=500, save_top_n=3)

    po = PipelineOptimization(direction='maximize',
                              study_name='sm_predictor_pipeline_random',
                              storage='sqlite:///test_sm_predictor.db',
                              sampler=random_sampler)
    po.optimize(train_dataset=train, test_dataset=valid, objective_steps=objective,
                metric=metric, n_trials=500, save_top_n=3)

    po = PipelineOptimization(direction='maximize',
                              study_name='sm_predictor_pipeline_tpe',
                              storage='sqlite:///test_sm_predictor.db',
                              sampler=tpe_sampler)

    po.optimize(train_dataset=train, test_dataset=valid, objective_steps=objective,
                metric=metric, n_trials=500, save_top_n=3)

    po = PipelineOptimization(direction='maximize',
                              study_name='sm_predictor_pipeline_cmaes',
                              storage='sqlite:///test_sm_predictor.db',
                              sampler=cmaes_sampler)

    po.optimize(train_dataset=train, test_dataset=valid, objective_steps=objective,
                metric=metric, n_trials=500, save_top_n=3)

    po = PipelineOptimization(direction='maximize',
                              study_name='sm_predictor_pipeline_quasi_mote_carlo',
                              storage='sqlite:///test_sm_predictor.db',
                              sampler=quasi_mote_carlo)

    po.optimize(train_dataset=train, test_dataset=valid, objective_steps=objective,
                metric=metric, n_trials=500, save_top_n=3)
