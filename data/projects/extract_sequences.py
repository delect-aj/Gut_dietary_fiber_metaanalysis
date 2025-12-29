#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import biom
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# 读取 BIOM 文件
table = biom.load_table("/home/dongbiao/word_embedding_microbiome/programe_test/dietary_fber/data/projects/table_6721_2378.biom")
# 定义输入和输出文件名
input_file = "/beegfs/dongbiao/greengene2/backbone_fasta/dna-sequences.fasta"
output_file = "/home/dongbiao/word_embedding_microbiome/programe_test/dietary_fber/data/projects/feces_seq_16S.fasta"
selected_names = table.ids(axis='observation')
selected_names# 读取fasta文件并选取需要的序列
selected_sequences = []
for record in SeqIO.parse(input_file, "fasta"):
    if record.id in selected_names:
        selected_sequences.append(record)
SeqIO.write(selected_sequences, output_file, "fasta")
