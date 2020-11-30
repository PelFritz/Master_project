import pandas as pd
import numpy as np

np.random.seed(42)
gen_fam = pd.read_csv('/nam-99/ablage/nam/peleke/gene_families.csv', usecols=['family_id', 'Gene_id'])
grouped_by_fam = gen_fam.groupby('family_id')

num_fam = np.arange(2882)
random_train_ids = np.random.choice(num_fam, 2682, replace=False)
random_test_id = [x for x in num_fam if x not in random_train_ids]


test_dfs = []
for fam_id in random_test_id:
    test_dfs.append(gen_fam[gen_fam['family_id'] == fam_id])
all_test_ids = pd.concat(test_dfs)
print(all_test_ids.shape)

all_test_ids.to_csv('/nam-99/ablage/nam/peleke/test_genes.csv', index=False)

train_dfs = []
for fam_id in random_train_ids:
    train_dfs.append(gen_fam[gen_fam['family_id'] == fam_id])
all_train_df = pd.concat(train_dfs)
print(all_train_df.shape)
all_train_df.to_csv('/nam-99/ablage/nam/peleke/train_genes.csv', index=False)
