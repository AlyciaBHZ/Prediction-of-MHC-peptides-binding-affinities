import os, pickle, datetime, argparse, string
import numpy as np
import pandas as pd

aa = ["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU",
      "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"]
aa_abbr = [x for x in "ACDEFGHIKLMNPQRSTVWY"]
aa_dict = dict(zip(aa, aa_abbr))


error_code_dic = {"PDB not exist": 1, "chain not exist": 2, "PDB_seq & dismap_seq mismatch": 3, "DSSP too long": 4, "Fail to pad DSSP": 5}


def get_PDB(PDBID, pdb_file, chain, data_path):
    ID = PDBID + chain
    if pdb_file != "": # User custom PDB file
        os.system("mv {} {}".format(pdb_file, os.path.dirname(pdb_file) + "/"+PDBID+'.pdb'))
        os.system("perl getchain.pl {} {}".format(os.path.dirname(pdb_file), ID))

    os.system("mv {} {}".format(ID, data_path)) # the output of getchain.pl is in current directory

    seq = ""
    current_pos = -1000
    with open(data_path + ID, "r") as f:
        lines = f.readlines()
    for line in lines:
        if line[0:4].strip() == "ATOM" and int(line[22:26].strip()) != current_pos:
            aa_type = line[17:20].strip()
            seq += aa_dict[aa_type]
            current_pos = int(line[22:26].strip())

    if seq == "":
        return "", error_code_dic["chain not exist"]
    else:
        return seq, 0


def process_distance_map(distance_map_file, cutoff = 14):
    with open(distance_map_file, "r") as f:
        lines = f.readlines()

    seq = lines[0].strip()
    length = len(seq)
    distance_map = np.zeros((length, length))

    if lines[1][0] == "#": # missed residues
        missed_idx = [int(x) for x in lines[1].split(":")[1].strip().split()] # 0-based
        lines = lines[2:]
    else:
        missed_idx = []
        lines = lines[1:]

    for i in range(0, len(lines)):
        record = lines[i].strip().split()
        for j in range(0, len(record)):
            if float(record[j]) == -1:
                distance_map[i + 1][j] = 0
            elif float(record[j]) <= cutoff:
                distance_map[i + 1][j] = 1
            else:
                distance_map[i + 1][j] = 0

    for idx in missed_idx:
        if idx > 0:
            distance_map[idx][idx - 1] = 1
        if idx > 1:
            distance_map[idx][idx - 2] = 1
        if idx < length - 1:
            distance_map[idx + 1][idx] = 1
        if idx < length - 2:
            distance_map[idx + 2][idx] = 1

    distance_map = distance_map + distance_map.T + np.eye(length)
    return seq, distance_map


def get_distance_map(ID, PDB_seq, data_path):
    os.system("./caldis_CA {} > {}.map".format(data_path + ID, data_path + "dismap/" + ID))
    dis_map_seq, dis_map = process_distance_map(data_path + "dismap/" + ID + ".map")
    if PDB_seq != dis_map_seq:
        return error_code_dic["PDB_seq & dismap_seq mismatch"]
    else:
        np.save(data_path + "dismap/" + ID, dis_map)
        return 0


def feature_extraction(PDBID, pdb_file, chain, mode, data_path):
    ID = PDBID + chain

    PDB_seq, error_code = get_PDB(PDBID, pdb_file, chain, data_path)
    if error_code != 0:
        return error_code

    with open(data_path + ID + ".fa", "w") as f:
        f.write(">" + ID + "\n" + PDB_seq)

    error_code = get_distance_map(ID, PDB_seq, data_path)
    if error_code != 0:
        return error_code

    return 0


def predict(ID, data_path, mode):
    with open(data_path + ID + ".fa", "r") as f:
        seq = f.readlines()[1].strip()

    test_dataframe = pd.DataFrame({"ID": [ID]})
    pred_scores = [round(score, 4) for score in test(test_dataframe, data_path, mode)]

    GraphPPIS_threshold = (0.24 if mode == "fast" else 0.18)
    binary_preds = [1 if score >= GraphPPIS_threshold else 0 for score in pred_scores]

    with open(data_path + ID + "_pred_results.txt", "w") as f:
        f.write("The threshold of the predictive score to determine PPI sites is set to {}.\n".format(GraphPPIS_threshold))
        f.write("AA\tProb\tPred\n")
        for i in range(len(seq)):
            f.write(seq[i] + "\t" + str(pred_scores[i]) + "\t" + str(binary_preds[i]) + "\n")


def main(PDBID, pdb_file, chain, mode):
    PDBID = PDBID.lower()
    chain = chain.upper()

    data_path = "./data/"+PDBID+"/"
    dir_list = ["", "dismap/"]
    for dir_name in dir_list:
        os.makedirs(data_path + dir_name)

    error_code = feature_extraction(PDBID, pdb_file, chain, mode, data_path)
    if error_code == 1:
        print("\nError! The query protein dosen't exist!")
    elif error_code == 2:
        print("\nError! The query chain dosen't exist in this protein!")
    elif error_code != 0:
        print("Error! Error code {}. Please contact the authors of GraphPPIS.".format(error_code))
    else:
        print("\nFeature Extraction is done at {}.\n".format(datetime.datetime.now().strftime("%m-%d %H:%M")))
        print("Predicting...\n")
        print("Map generated!")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--protein", type = str, help = "PDBID (e.g. 3mcb)")
    parser.add_argument("-n", "--name", type = str, help = "PDBname (e.g. 3mcb)")
    parser.add_argument("-f", "--file", type = str, help = "PDB file (.pdb only)")
    parser.add_argument("-c", "--chain", type = str, default = "A", help = "chain identifier")
    parser.add_argument("-m", "--mode", type = str, default = "fast", help = "fast (use BLOSUM62 + DSSP) or slow (use PSSM + HMM + DSSP)")
    args = parser.parse_args()

    if args.chain == None:
        print("Chain identifier is not provided!")
    elif args.chain not in list(string.ascii_letters + string.digits):
        print("Invalid chain identifier!")
    elif args.mode not in ["fast", "slow"]:
        print("Invalid mode!")
    elif args.file: # input by file
        if args.file.endswith(".pdb") == False:
            print("only .pdb file is supported!")
        else:
            main(args.name, args.file, args.chain, args.mode)
    else: # input by PDBID
        if args.protein == None or len(args.protein) != 4:
            print("Invalid PDB ID!")
        else:
            main(args.protein, "", args.chain, args.mode)
