import logging
import warnings

import optuna
from deepmol.loaders import CSVLoader
from rdkit import RDLogger
from sklearn.metrics import f1_score

from sm_precursor_predictor.models_training.train_models import objective


def train_models():
    warnings.filterwarnings("ignore")
    logger = logging.getLogger()
    logger.setLevel(logging.CRITICAL)
    RDLogger.DisableLog('rdApp.*')

    train = CSVLoader("train.csv",
                      labels_fields=['C00073', 'C00078', 'C00079', 'C00082', 'C00235', 'C00341', 'C00353',
                                     'C00448', 'C01789', 'C03506', 'C00047', 'C00108', 'C00187', 'C00148',
                                     'C00041', 'C00129', 'C00062', 'C01852', 'C00049', 'C00135', 'C00223',
                                     'C00509', 'C00540', 'C01477', 'C05903', 'C05904', 'C05905', 'C05908',
                                     'C09762'],
                      id_field="ids", smiles_field="smiles").create_dataset()
    valid = CSVLoader("valid.csv",
                      labels_fields=['C00073', 'C00078', 'C00079', 'C00082', 'C00235', 'C00341', 'C00353',
                                     'C00448', 'C01789', 'C03506', 'C00047', 'C00108', 'C00187', 'C00148',
                                     'C00041', 'C00129', 'C00062', 'C01852', 'C00049', 'C00135', 'C00223',
                                     'C00509', 'C00540', 'C01477', 'C05903', 'C05904', 'C05905', 'C05908',
                                     'C09762'],
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


if __name__ == '__main__':
    train_models()
