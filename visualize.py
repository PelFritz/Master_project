import os
import pandas as pd
import numpy as np
from tensorflow import set_random_seed
from Bio import SeqIO
from tensorflow.keras import models
from sklearn.metrics import accuracy_score
from tensorflow.python.keras.utils import np_utils
from deeplift.visualization import viz_sequence
from deeplift.conversion import kerasapi_conversion as kc
from deeplift.layers import NonlinearMxtsMode
from deeplift.dinuc_shuffle import dinuc_shuffle

np.random.seed(42)
set_random_seed(42)
pd.options.display.width = 0

codes = {
    'A': [1, 0, 0, 0],
    'C': [0, 1, 0, 0],
    'G': [0, 0, 1, 0],
    'T': [0, 0, 0, 1],
    'N': [0, 0, 0, 0]
}


def one_hot(sequence):
    one_encoding = np.zeros(shape=(4, len(sequence)))
    for i, nt in enumerate(sequence):
        one_encoding[:, i] = codes[nt]

    return one_encoding


norm_counts = pd.read_csv('/nam-99/ablage/nam/peleke/GSE80744_ath1001_tx_norm_2016-04-21-UQ_gNorm_normCounts_k4.tsv',
                          sep='\t')


norm_counts.set_index('gene_id', inplace=True)
norm_counts['unexp'] = (norm_counts == 0).sum(axis=1)
#norm_counts = norm_counts[norm_counts['unexp'] != 0]
#print(norm_counts.loc['AT1G03170':, ].head(50))

gene_ids_in_counts = list(norm_counts.index)
genomekeys_to_normcounts = list(norm_counts.columns)
genomekeys_to_normcounts = np.random.choice(genomekeys_to_normcounts, 300, replace=False)


print('Fetching Promoters')

prom_seq = []
prom_shuf = []
label = []
for genome in os.listdir('/nam-99/ablage/nam/peleke/snp_promoters'):
    genome_path = '/nam-99/ablage/nam/peleke/snp_promoters/' + genome
    ecotype_id = genome.split('.')[0]
    genome_key = 'X' + ecotype_id

    if genome_key in genomekeys_to_normcounts:
        for rec in SeqIO.parse(genome_path, 'fasta'):
            ID = rec.id.split(':')[1]
            seq = str(rec.seq)

            if ID == 'AT1G01140':
                prom_seq.append(seq)
                prom_shuf.append(dinuc_shuffle(seq, rng=np.random.RandomState(seed=42)))
                if norm_counts[genome_key][ID] == 0:
                    label.append(0)
                else:
                    label.append(1)


print('Encoding sequences')
encoded_seq = np.array([one_hot(prom) for prom in prom_seq])
encoded_seq = np.expand_dims(encoded_seq, 3)
encoded_shuf_seq = np.array([one_hot(prom) for prom in prom_shuf])
encoded_shuf_seq = np.expand_dims(encoded_shuf_seq, 3)
categories = np_utils.to_categorical(label, 2)

model = models.load_model('/nam-99/ablage/nam/peleke/Thesis_models/model2020-10-06073328.h5')
predictions = np.argmax(model.predict(encoded_seq), axis=1)
actual = np.argmax(categories, axis=1)
print(predictions)

print(accuracy_score(label, predictions))

# compute deeplift scores
deeplift_model =\
    kc.convert_model_from_saved_files('/nam-99/ablage/nam/peleke/Thesis_models/model2020-10-06073328.h5',
                                      nonlinear_mxts_mode=NonlinearMxtsMode.DeepLIFT_GenomicsDefault)

deeplift_contrib_func = deeplift_model.get_target_contribs_func(find_scores_layer_idx=0,
                                                                target_layer_idx=-2)

# True positive predictions
tp = []
tp_shuf = []

for pred, true, enc_seq, enc_shuf_seq in zip(predictions, actual, encoded_seq, encoded_shuf_seq):
    if pred == 1 and true == 1:
        tp.append(enc_seq)
        tp_shuf.append(enc_shuf_seq)


tp_data = np.array(tp)
tp_shuf_data = np.array(tp_shuf)

print(tp_data.shape)

print(accuracy_score(actual, predictions))

scores_tp = np.array(deeplift_contrib_func(task_idx=1, input_data_list=[tp_data],
                                           input_references_list=[tp_shuf_data],
                                           batch_size=10, progress_update=10))

sum_tp_scores = np.squeeze(np.sum(scores_tp, axis=0))
sum_tp_scores = np.sum(sum_tp_scores, axis=0)

print(sum_tp_scores)
print(sum_tp_scores.shape)
for idx in [3, 8, 5]:
    original_onehot = tp_data[idx]
    scores_for_indx_tp = original_onehot * sum_tp_scores[:, None]
    viz_sequence.plot_weights(scores_for_indx_tp, subticks_frequency=10, figsize=(100, 4))




