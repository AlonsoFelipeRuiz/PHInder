#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 08:12:05 2021

@author: alonsofeliperuiz
"""

import os
path=os.getcwd()
from Bio import SeqIO
from glob import glob

def Prokka(records, genome):
    direct=str("/Users/alonsofeliperuiz/Desktop/Doctorado26_08_20/Otros_proyectos/Prometeo/prokka")
    os.chdir(direct)
    prokka = str("prokka --kingdom Bacteria --compliant --force --quiet --prefix prokka_results_" + genome + " --cpus 0 --fast " + records + " --locustag " + genome)
    os.system(prokka)
    PhiSpy(genome, gbk_file, fasta_phages, gbk_phages)
    HmmSearch(genome, annotatedfile)
    phage_genome=str("Phages/" + genome)
    HmmSearch(phage_genome, fasta_phages)    
    

def PhiSpy(genome, gbk_file, fasta_phages, gbk_phages):
    os.chdir("/Users/alonsofeliperuiz/Desktop/Doctorado26_08_20/Otros_proyectos/Prometeo/PhiSpy")
    phispy=str("PhiSpy.py " +  gbk_file + " -o Phages_" + genome + " --output_choice 4")
    os.system(phispy)
    output_handle = open(fasta_phages, "w")
    for seq_record in SeqIO.parse(gbk_phages, "gb"):
        for seq_feature in seq_record.features :
            if seq_feature.type=="CDS" :
                assert len(seq_feature.qualifiers['translation'])==1
                output_handle.write(">%s from %s\n%s\n" % (
                       seq_feature.qualifiers['locus_tag'][0],
                       seq_record.name,
                       seq_feature.qualifiers['translation'][0]))

    
    

def MakeHits(files, annotatedfile, genome):
    thresholdscov={
        "AimR":60, "NprR":80, "PlcR":80, "PgrX":80, 
        "Rap":80, "Rgg":80, "ComR":80
                    }
    thresholdeval={
        "AimR":8.9e-16, "NprR":7.9e-31, "PlcR":5.5e-76, "PgrX":7e-46, 
        "Rap":7.9e-31, "Rgg":7.9e-31, "ComR":7.9e-31
                    }
    gene_list=[]
    proteinname=files.split("/")[-1].split(".")[0]
    coverage=thresholdscov[proteinname]
    domeval=thresholdeval[proteinname]
    with open(files, "r") as b_output:
        for line in b_output:
            if line[0]=="#": continue
            elements=[el for el in line.split(" ") if el!=""] 
            p_lenght=elements[5]
            break
    with open(files, "r") as b_output:
        print(f"getting data {genome} {proteinname}")
        for line in b_output:
            if line[0]=="#": continue
            elements=[el for el in line.split(" ") if el!=""]   
            #Get data
            gene_hit=elements[0]
            domevalue=float(elements[11])
            s_from=int(elements[17])
            s_to=int(elements[18])        
            cover=((s_to-s_from)/float(p_lenght))*100
            if cover<coverage: continue
            if domevalue>domeval: continue
            gene_list.append(gene_hit)
        seq_list=[]
        prts_dict = SeqIO.to_dict(SeqIO.parse(annotatedfile, "fasta"))
        for gene in gene_list:
            seq_list.append(prts_dict[gene])
        path2= str("/Users/alonsofeliperuiz/Desktop/Doctorado26_08_20/Otros_proyectos/Prometeo/Results/" + genome + "/"  + proteinname + ".fasta")
        SeqIO.write(seq_list, path2, "fasta")
    
def HmmSearch(genome, annotatedfile):
    directhmm= glob("/Users/alonsofeliperuiz/Desktop/Doctorado26_08_20/Otros_proyectos/Prometeo/HMM_MODELS/*")
    direct2=str("/Users/alonsofeliperuiz/Desktop/Doctorado26_08_20/Otros_proyectos/Prometeo/HMM_outputs/" + genome)
    os.mkdir(direct2)
    os.mkdir("/Users/alonsofeliperuiz/Desktop/Doctorado26_08_20/Otros_proyectos/Prometeo/Results/" + genome)
    os.chdir(direct2)
    for hmms in directhmm:
        hmm=hmms.split("/")[-1].split(".")[0]
        hmmsearch = str("hmmsearch --noali --domtblout " + hmm + ".out " + hmms  + " " + annotatedfile)
        os.system(hmmsearch)
    direct3 = glob(str(direct2 + "/*"))
    for files in direct3:
        MakeHits(files, annotatedfile, genome)


files = glob("/Users/alonsofeliperuiz/Desktop/Doctorado26_08_20/Otros_proyectos/Prometeo/New_Genomes/*")
directgenomes= "/Users/alonsofeliperuiz/Desktop/Doctorado26_08_20/Otros_proyectos/Prometeo/New_Genomes/"
for records in files:
    genome=records.split("/")[-1].split(".")[0]
    annotatedfile=str("/Users/alonsofeliperuiz/Desktop/Doctorado26_08_20/Otros_proyectos/Prometeo/prokka/prokka_results_" + genome + "/prokka_results_" + genome + ".faa")
    gbk_file=str("/Users/alonsofeliperuiz/Desktop/Doctorado26_08_20/Otros_proyectos/Prometeo/prokka/prokka_results_" + genome + "/prokka_results_" + genome + ".gbf")
    fasta_phages=str("/Users/alonsofeliperuiz/Desktop/Doctorado26_08_20/Otros_proyectos/Prometeo/PhiSpy/Phages_" + genome + "/phage.faa")
    gbk_phages=str("/Users/alonsofeliperuiz/Desktop/Doctorado26_08_20/Otros_proyectos/Prometeo/PhiSpy/Phages_" + genome + "/phage.gbk")
    print(genome)
    Prokka(records, genome)
    
    
        
        



