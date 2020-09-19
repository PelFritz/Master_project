"""This file flanks genes to get coordinates for promoters"""
import itertools
import pandas as pd
import os


def flank(gtf_file, out_gtf, pro_size):
    new_start = []
    new_end = []
    records = []
    with open(gtf_file, 'r') as read_obj:
        for line in itertools.islice(read_obj, 0, None):
            records.append(list(line.split('\t')))
            if line.split('\t')[6] == '+':
                new_start.append(int(line.split('\t')[3]) - pro_size)
                new_end.append(int(line.split('\t')[3]) + 200)
            else:
                new_start.append(int(line.split('\t')[4]) -200)
                new_end.append(int(line.split('\t')[4]) + pro_size)

    df = pd.DataFrame(records)
    df[3] = new_start
    df[4] = new_end
    df[7] = [x.split("\n")[0] for x in df[7]]
    df.drop(8, axis=1, inplace=True)

    return df.to_csv(out_gtf, header=False, index=False, sep='\t')


for file in os.listdir('/nam-99/ablage/nam/peleke/variant_gtfs'):
    input_gtf = '/nam-99/ablage/nam/peleke/variant_gtfs/' + file
    output_gtf = '/nam-99/ablage/nam/peleke/1200upstream_gtfs/' + file
    flank(input_gtf, output_gtf, 1000)
    print('Done with ', file)
