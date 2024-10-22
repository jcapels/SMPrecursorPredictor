from deepmol.loaders import CSVLoader

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

# DEFINE THE FUNCTION
from deepmol.pipeline_optimization._featurizer_objectives import _get_featurizer
from deepmol.pipeline_optimization._scaler_objectives import _get_scaler
from deepmol.pipeline_optimization._feature_selector_objectives import _get_feature_selector
from deepmol.pipeline_optimization._sklearn_model_objectives import _get_sk_model
from deepmol.pipeline_optimization._standardizer_objectives import _get_standardizer
from deepmol.compound_featurization import MixedFeaturizer, All3DDescriptors, NPClassifierFP
from deepmol.metrics import Metric
from sklearn.metrics import f1_score


def steps(trial, data):
    mode = data.mode
    multitask = True if data.n_tasks > 1 else False
    featurizer = NPClassifierFP()
    scaler = _get_scaler(trial)
    feature_selector = _get_feature_selector(trial, task_type=mode, multitask=multitask)

    sk_mode = mode
    sk_model = _get_sk_model(trial, task_type=sk_mode)
    final_steps = [('standardizer', _get_standardizer(trial, None)), ('featurizer', featurizer),
                   ('scaler', scaler),
                   ('feature_selector', feature_selector), ('model', sk_model)]
    return final_steps

from deepmol.pipeline_optimization import PipelineOptimization

def f1_score_macro(y_true, y_pred):
    return f1_score(y_true, y_pred, average='macro')

metric = Metric(f1_score_macro)

po = PipelineOptimization(direction='maximize', study_name='np_classifier_fp', storage='sqlite:///np_classifier_fp.db')
po.optimize(train_dataset=train, test_dataset=valid, objective_steps=steps, 
            metric=metric, n_trials=50, save_top_n=5, trial_timeout=60*60*2, data=train)
