# Prediction-of-MHC-peptides-binding-affinities

## 环境配置：
### The project is developed under Linux environment with:
Cuda 10.2  
python 3.9  
numpy latest version  
pandas latest version  
torch 1.10.0  


# 节点特征训练 

### pssm模型训练
~~~
psiblast -query ./fasta/iedb00003.fasta -db uniprot_sprot.fasta -num_iterations 3 -out_ascii_pssm ./pssm/iedb00003.pssm
~~~

### dssp模型训练
~~~
dssp -i ./pdb/paper0001.pdb -o ./data/dssp/paper0001.dssp
~~~

# 邻接矩阵和Distance Map 生成
~~~
python map.py -n 0001 -f /home/ubuntu/deeprank/8000/paper0001.pdb  
~~~
-n :自定义传入文件命名。也可以使用-p 引用PDB id, -f自定义输入的蛋白路径


## 在dataset-example中查看将上述三种特征normalized并转化为.npy格式的实现
