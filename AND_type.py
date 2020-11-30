import pandas as pd
import numpy as np
import re
import os
from Bio import SeqIO
from tensorflow.python.keras.utils import np_utils
from sklearn.model_selection import train_test_split
from tensorflow.keras.layers import Dense, Conv2D, MaxPool2D, Dropout, Flatten
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras import backend
from tensorflow.keras import Sequential
import tensorflow as tf
np.random.seed(42)
tf.set_random_seed(42)

norm_counts = pd.read_csv('/nam-99/ablage/nam/peleke/GSE80744_ath1001_tx_norm_2016-04-21-UQ_gNorm_normCounts_k4.tsv',
                          sep='\t')

print('step 1: Filtering normcounts to keep only those with expression in all pseudogenomes')
df = norm_counts.set_index('gene_id')
fam_df = pd.read_csv('/nam-99/ablage/nam/peleke/data.txt', sep='\t')

# Edit normalised counts to keep only those for genes in familiy
protein_binding_genes = fam_df['# gene_id'].values.tolist()
df = df[df.index.isin(protein_binding_genes)]
df['unexp'] = (df == 0).sum(axis=1)
df = df[df['unexp'] == 0]
print(df.shape)
df = df.head(250)
gene_id = list(df.index)
genomes = list(df.columns)
genomes = genomes[:500]
print('There are ', len(gene_id), ' genes per genome being verified for patterns')

# Search for TATA-Box promoters
hexamer = 'TATAAA'
mut_hexamer = 'TGTAAA'

hexamer2 = 'CATAAT'
mut_hexamer_2 = 'CACAAT'

label = []
prom_seq = []
print('step 2: Mutating TATA-boxes')
for genome in os.listdir('/nam-99/ablage/nam/peleke/promoters'):
    genome_key = 'X' + genome.split('.')[0]
    if genome_key in genomes:
        path_to_gen = '/nam-99/ablage/nam/peleke/promoters/' + genome
        for promoter in SeqIO.parse(path_to_gen, 'fasta'):
            prom_id = promoter.id.split(':')[1]
            sequence = str(promoter.seq)
            if prom_id in gene_id:
                if re.search(hexamer, sequence) and re.search(hexamer2, sequence):
                    # mutate TATA-box on a sequence and both motifs on another sequence
                    mut_sequence = re.sub(hexamer, mut_hexamer, sequence, count=1)
                    prom_seq.append(mut_sequence)
                    label.append(1)
                    mut_sequence2 = re.sub(hexamer2, mut_hexamer_2, mut_sequence, count=1)
                    prom_seq.append(mut_sequence2)
                    label.append(0)

                    # Mutate both motifs on a sequence and leave the original sequence unmutated
                    prom_seq.append(sequence)
                    label.append(1)
                    mut_sequence2 = re.sub(hexamer2, mut_hexamer_2, mut_sequence, count=1)
                    prom_seq.append(mut_sequence2)
                    label.append(0)


print('We would train the network with {} real sequences and mutated sequences'.format(len(prom_seq)))

# one hot encode sequences
codes = {
    'A': [1, 0, 0, 0],
    'C': [0, 1, 0, 0],
    'G': [0, 0, 1, 0],
    'T': [0, 0, 0, 1],
    'N': [0, 0, 0, 0]}


# One hot encode data and split
def onehot_encode(seq):
    encoding = np.zeros(shape=(4, len(seq)))
    for i, nt in enumerate(seq):
        encoding[:, i] = codes[nt]
    return encoding


print('step 3: Encoding sequences')
encoded_sequences = np.array([onehot_encode(prom) for prom in prom_seq])
prepared_sequences = np.expand_dims(encoded_sequences, 3)
categories = np_utils.to_categorical(label, 2)
print(prepared_sequences.shape)
print(categories.shape)

X_train, X_test, y_train, y_test = train_test_split(prepared_sequences, categories, test_size=0.2, shuffle=True,
                                                    stratify=categories, random_state=42)


# Build Network
def train_model(x_train, train_label, x_test, test_label):
    backend.clear_session()
    model = Sequential()
    model.add(Conv2D(64, kernel_size=(2, 6), padding='same', activation='relu', input_shape=(4, 1000, 1)))
    model.add(Conv2D(64, kernel_size=(1, 6), padding='same', activation='relu'))
    model.add(MaxPool2D(pool_size=(1, 2), strides=(1, 2)))
    model.add(Dropout(0.25))

    model.add(Conv2D(128, kernel_size=(1, 6), padding='same', activation='relu', strides=(1, 6)))
    model.add(Conv2D(128, kernel_size=(1, 6), padding='same', activation='relu', strides=(1, 6)))
    model.add(Dropout(0.25))

    model.add(Flatten())
    model.add(Dense(128, activation='relu'))
    model.add(Dense(2, activation='softmax'))
    print(model.summary())

    es = EarlyStopping(monitor='val_loss', patience=3, verbose=1)

    model.compile(optimizer='Adam', loss='binary_crossentropy', metrics=['accuracy'])
    model.fit(x_train, train_label, batch_size=256, epochs=40, callbacks=[es],
              validation_data=(x_test, test_label))


print('step 4: Training model')
train_model(X_train, y_train, X_test, y_test)
