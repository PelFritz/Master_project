from tensorflow.keras.layers import Layer, Conv2D, BatchNormalization, Dense, MaxPool2D, Activation, Flatten, Dropout
from tensorflow.keras import Sequential, backend
from tensorflow.keras.callbacks import EarlyStopping
import tensorflow as tf
import pandas as pd
import numpy as np
import os
from Bio import SeqIO
from tensorflow.python.keras.utils import np_utils
from datetime import datetime

np.random.seed(42)
tf.set_random_seed(42)
pd.options.display.width = 0

codes = {
    'A': [1, 0, 0, 0],
    'C': [0, 1, 0, 0],
    'G': [0, 0, 1, 0],
    'T': [0, 0, 0, 1],
    'N': [0, 0, 0, 0]
}

# custom functions used in this scripts---------------------------------------------------------------------------------


def one_hot(sequence):
    one_encoding = np.zeros(shape=(4, len(sequence)))
    for i, nt in enumerate(sequence):
        one_encoding[:, i] = codes[nt]

    return one_encoding


def cnn_model(x_train, y_train, x_test, y_test):
    backend.clear_session()

    model = Sequential()
    model.add(Conv2D(filters=64, kernel_size=(4, 6), padding='same', input_shape=(4, 1000, 1)))
    model.add(Activation('relu'))
    model.add(Conv2D(filters=64, kernel_size=(1, 6), padding='same'))
    model.add(Activation('relu'))
    model.add(Conv2D(filters=64, kernel_size=(1, 6), padding='same'))
    model.add(Activation('relu'))
    model.add(Dropout(0.25))

    model.add(Flatten())
    model.add(Dense(units=2, activation='softmax'))
    print(model.summary())
    es = EarlyStopping(monitor='val_loss', patience=5, verbose=1)

    model.compile(optimizer='Adam', loss='binary_crossentropy', metrics=['accuracy'])
    model.fit(x=x_train, y=y_train, batch_size=256, validation_data=(x_test, y_test), epochs=40, callbacks=[es])

    now = datetime.now().strftime('%Y-%m-%d%H%M%S')
    model.save('/nam-99/ablage/nam/peleke/Thesis_models/' + 'model' + now + '.h5')


train_df = pd.read_csv('/nam-99/ablage/nam/peleke/train_genes.csv')
test_df = pd.read_csv('/nam-99/ablage/nam/peleke/test_genes.csv')
train_genes = train_df['Gene_id'].tolist()
test_genes = test_df['Gene_id'].tolist()

norm_counts = pd.read_csv('/nam-99/ablage/nam/peleke/GSE80744_ath1001_tx_norm_2016-04-21-UQ_gNorm_normCounts_k4.tsv',
                          sep='\t')

norm_counts.set_index('gene_id', inplace=True)

gene_ids_in_counts = list(norm_counts.index)
genomekeys_to_normcounts = list(norm_counts.columns)

# Fetch promoters and class---------------------------------------------------------------------------------------------
train_prom_seq = []
train_label = []
test_prom_seq = []
test_label = []

print('Fetching Promoters')
for genome in os.listdir('/nam-99/ablage/nam/peleke/snp_promoters'):
    genome_path = '/nam-99/ablage/nam/peleke/snp_promoters/' + genome
    ecotype_id = genome.split('.')[0]
    genome_key = 'X' + ecotype_id

    if genome_key in genomekeys_to_normcounts:
        for rec in SeqIO.parse(genome_path, 'fasta'):
            ID = rec.id.split(':')[1]
            seq = str(rec.seq)
            if ID in train_genes and ID in gene_ids_in_counts:
                train_prom_seq.append(seq)
                if norm_counts[genome_key][ID] == 0:
                    train_label.append(0)
                else:
                    train_label.append(1)

            elif ID in test_genes and ID in gene_ids_in_counts:
                test_prom_seq.append(seq)
                if norm_counts[genome_key][ID] == 0:
                    test_label.append(0)
                else:
                    test_label.append(1)


label = np.array(train_label)
indices_unexp = np.where(label == 0)[0]
indices_ex = np.where(label == 1)[0]

# Encoding and balancing training set-------------------------------------------------------------------------------
print('Class sizes in the original training set')
print('Expressed: ', len(indices_ex), 'Unexpressed: ', len(indices_unexp))


print('Encoding train set and downsampling')
train_encoded = np.array([one_hot(sequence) for sequence in train_prom_seq])
train_encoded_expand = np.expand_dims(train_encoded, 3)
categories = np_utils.to_categorical(train_label, 2)

sel_expressed_idx = np.random.choice(indices_ex, len(indices_unexp), replace=False)
sel_expressed_prom = np.take(train_encoded_expand, sel_expressed_idx, axis=0)
sel_unexpressed_prom = np.take(train_encoded_expand, indices_unexp, axis=0)
x_train = np.concatenate((sel_expressed_prom, sel_unexpressed_prom), axis=0)

sel_expressed_labels = np.take(categories, sel_expressed_idx, axis=0)
sel_unexpressed_labels = np.take(categories, indices_unexp, axis=0)
y_train = np.concatenate((sel_expressed_labels, sel_unexpressed_labels), axis=0)

print('Prepared dataset shape: ', x_train.shape)
print(y_train.shape)

# Encoding test set--------------------------------------------------------------------------------------------------
print('Encoding test set')
test_encoded = np.array([one_hot(sequence) for sequence in test_prom_seq])
x_test = np.expand_dims(test_encoded, 3)
y_test = np_utils.to_categorical(test_label, 2)
print(x_test.shape)
print(y_test.shape)

# Training------------------------------------------------------------------------------------------------------------
cnn_model(x_train, y_train, x_test, y_test)
