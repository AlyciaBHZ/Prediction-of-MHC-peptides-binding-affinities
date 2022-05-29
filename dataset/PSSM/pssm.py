#处理生成的pssm数据为.npy格式，需要按需要修改地址

import numpy as np
import pandas as pd

Max_pssm = np.array([8, 9, 9, 9, 12, 9, 8, 8, 12, 9, 7, 9, 11, 10, 9, 8, 8, 13, 10, 8])
Min_pssm = np.array([-11,-12,-13,-13,-12,-13,-13,-12,-13,-13,-13,-13,-12,-12,-13,-12,-12,-13,-13,-12])

def process_pssm(pssm_file):
    with open(pssm_file, "r") as f:
        lines = f.readlines()
    pssm_feature = []
    for line in lines:
        if line == "\n":
            continue
        record = line.strip().split()
        if record[0].isdigit():
            pssm_feature.append([int(x) for x in record[2:22]])
    pssm_feature = (np.array(pssm_feature) - Min_pssm) / (Max_pssm - Min_pssm)
    return pssm_feature

import csv
with open('/home/ubuntu/bio/code.csv', "r") as f:
    lines = csv.reader(f)
    for row in lines:
        row=str(row)
        row=row.strip("'[")
        row=row.strip("']")
        print("/home/ubuntu/deeprank/pssm/"+row+".pssm")
        pssm_matrix = process_pssm("/home/ubuntu/deeprank/pssm/"+row+".pssm")
        np.save("./pssm/"+row, pssm_matrix)