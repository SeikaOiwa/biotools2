import subprocess
import os
import shutil
import argparse
import time

parser = argparse.ArgumentParser(description='GeneAnnotation')
parser.add_argument('braker_mode', type=str, help='解析モード(g,g_f,gr,gr_f,gp,gp_f,grp,grp_f)')
parser.add_argument('f_path', type=str, help='GeneAnnotationフォルダパス')
parser.add_argument('augustus', type=str, help='Augustusフォルダパス')
parser.add_argument('sname', type=str, help='菌株名')
parser.add_argument('r_path', type=str,help='braker解析データの保存先パス')
parser.add_argument('s_path', type=str,help='保存先のフォルダパス')
parser.add_argument('g_path',type=str, help='ゲノム(*.fasta)へのパス')
parser.add_argument('rna_path',type=str, help='RNAseq(*.bam)へのパス')
parser.add_argument('propath',type=str, help='protein(*.fasta)へのパス')

args = parser.parse_args()

mode = args.braker_mode
base_path = args.f_path
augustus_path = args.augustus
strain_name = args.sname
result_path = args.r_path
data_f_path = args.s_path
genome_path = args.g_path
bam_path = args.rna_path
protein_path = args.propath

s_time = time.time()

def input_progress(base_path,comment):
    """進捗状況をテキストファイルに書き出し
    Parameter:
    ---------
    base_path: str
        テキストファイルの書出し先
    comment: str
        書出し内容（"｛進捗｝,｛所要時間(h)｝")    
    """
    with open(f"{base_path}/progress.txt","w") as f:
        f.write(f"{comment}")

# Braker実行
if mode == "g":
    comment = "解析開始,0h"
    input_progress(base_path,comment)
    subprocess.run([f'{base_path}/scr/braker_1.sh {strain_name} {genome_path}'],stdout=subprocess.PIPE,shell=True)
    # 進捗状況
    if os.path.isfile(f'{result_path}/braker/braker.gtf'):
        f_a_time = time.time()
        dgtime2 = round((f_a_time - s_time)/60/60,1)
        comment = f"Augustus,{dgtime2}"
        input_progress(base_path,comment)

        # フォルダ名＝菌株名に修正
        os.rename(f'{result_path}/braker',f'{result_path}/{strain_name}')
        shutil.move(f'{result_path}/{strain_name}',f'{data_f_path}/{strain_name}')

        # Augustusのconfig/speciesからstrain_nameに該当するフォルダを削除
        shutil.rmtree(f'{augustus_path}/config/species/{strain_name}') 

if mode == "g_f":
    comment = "解析開始,0h"
    input_progress(base_path,comment)
    subprocess.run([f'{base_path}/scr/braker_1_f.sh {strain_name} {genome_path}'],stdout=subprocess.PIPE,shell=True)
        # 進捗状況     

    if os.path.isfile(f'{result_path}/braker/braker.gtf'):
        f_a_time = time.time()
        dgtime2 = round((f_a_time - s_time)/60/60,1)
        comment = f"Augustus,{dgtime2}"
        input_progress(base_path,comment)

        # フォルダ名＝菌株名に修正
        os.rename(f'{result_path}/braker',f'{result_path}/{strain_name}')
        shutil.move(f'{result_path}/{strain_name}',f'{data_f_path}/{strain_name}')

        # Augustusのconfig/speciesからstrain_nameに該当するフォルダを削除
        shutil.rmtree(f'{augustus_path}/config/species/{strain_name}') 

if mode == "gr":
    comment = "解析開始,0h"
    input_progress(base_path,comment)
    subprocess.run([f'{base_path}/scr/braker_2.sh {strain_name} {genome_path} {bam_path}'],stdout=subprocess.PIPE,shell=True)      

    if os.path.isfile(f'{result_path}/braker/braker.gtf'):
        f_a_time = time.time()
        dgtime2 = round((f_a_time - s_time)/60/60,1)
        comment = f"Augustus,{dgtime2}"
        input_progress(base_path,comment)

        # フォルダ名＝菌株名に修正
        os.rename(f'{result_path}/braker',f'{result_path}/{strain_name}')
        shutil.move(f'{result_path}/{strain_name}',f'{data_f_path}/{strain_name}')

        # Augustusのconfig/speciesからstrain_nameに該当するフォルダを削除
        shutil.rmtree(f'{augustus_path}/config/species/{strain_name}') 

if mode == "gr_f":
    comment = "解析開始,0h"
    input_progress(base_path,comment)
    subprocess.run([f'{base_path}/scr/braker_2_f.sh {strain_name} {genome_path} {bam_path}'],stdout=subprocess.PIPE,shell=True)

        # 進捗状況

    if os.path.isfile(f'{result_path}/braker/braker.gtf'):
        f_a_time = time.time()
        dgtime2 = round((f_a_time - s_time)/60/60,1)
        comment = f"Augustus,{dgtime2}"
        input_progress(base_path,comment)

        # フォルダ名＝菌株名に修正
        os.rename(f'{result_path}/braker',f'{result_path}/{strain_name}')
        shutil.move(f'{result_path}/{strain_name}',f'{data_f_path}/{strain_name}')

        # Augustusのconfig/speciesからstrain_nameに該当するフォルダを削除
        shutil.rmtree(f'{augustus_path}/config/species/{strain_name}') 

if mode == "gp":
    comment = "解析開始,0h"
    input_progress(base_path,comment)
    subprocess.run([f'{base_path}/scr/braker_3.sh {strain_name} {genome_path} {protein_path}'],stdout=subprocess.PIPE,shell=True)

        # 進捗状況
    if os.path.isfile(f'{result_path}/braker/braker.gtf'):
        f_a_time = time.time()
        dgtime2 = round((f_a_time - s_time)/60/60,1)
        comment = f"Augustus,{dgtime2}"
        input_progress(base_path,comment)

        # フォルダ名＝菌株名に修正
        os.rename(f'{result_path}/braker',f'{result_path}/{strain_name}')
        shutil.move(f'{result_path}/{strain_name}',f'{data_f_path}/{strain_name}')

        # Augustusのconfig/speciesからstrain_nameに該当するフォルダを削除
        shutil.rmtree(f'{augustus_path}/config/species/{strain_name}') 

if mode == "gp_f":
    comment = "解析開始,0h"
    input_progress(base_path,comment)
    subprocess.run([f'{base_path}/scr/braker_3_f.sh {strain_name} {genome_path} {protein_path}'],stdout=subprocess.PIPE,shell=True)

    # 進捗状況
    if os.path.isfile(f'{result_path}/braker/braker.gtf'):
        f_a_time = time.time()
        dgtime2 = round((f_a_time - s_time)/60/60,1)
        comment = f"Augustus,{dgtime2}"
        input_progress(base_path,comment)

        # フォルダ名＝菌株名に修正
        os.rename(f'{result_path}/braker',f'{result_path}/{strain_name}')
        shutil.move(f'{result_path}/{strain_name}',f'{data_f_path}/{strain_name}')

        # Augustusのconfig/speciesからstrain_nameに該当するフォルダを削除
        shutil.rmtree(f'{augustus_path}/config/species/{strain_name}') 

if mode == "gpr":
    comment = "解析開始,0h"
    input_progress(base_path,comment)
    subprocess.run([f'{base_path}/scr/braker_4.sh {strain_name} {genome_path} {bam_path} {protein_path}'],stdout=subprocess.PIPE,shell=True)

        # 進捗状況       
    if os.path.isfile(f'{result_path}/braker/braker.gtf'):
        f_a_time = time.time()
        dgtime2 = round((f_a_time - s_time)/60/60,1)
        comment = f"Augustus,{dgtime2}"
        input_progress(base_path,comment)

        # フォルダ名＝菌株名に修正
        os.rename(f'{result_path}/braker',f'{result_path}/{strain_name}')
        shutil.move(f'{result_path}/{strain_name}',f'{data_f_path}/{strain_name}')

        # Augustusのconfig/speciesからstrain_nameに該当するフォルダを削除
        shutil.rmtree(f'{augustus_path}/config/species/{strain_name}') 

if mode == "gpr_f":
    comment = "解析開始,0h"
    input_progress(base_path,comment)
    subprocess.run([f'{base_path}/scr/braker_4_f.sh {strain_name} {genome_path} {bam_path} {protein_path}'],stdout=subprocess.PIPE,shell=True)

    # 進捗状況
    if os.path.isfile(f'{result_path}/braker/braker.gtf'):
        f_a_time = time.time()
        dgtime2 = round((f_a_time - s_time)/60/60,1)
        comment = f"Augustus,{dgtime2}"
        input_progress(base_path,comment)

        # フォルダ名＝菌株名に修正
        os.rename(f'{result_path}/braker',f'{result_path}/{strain_name}')
        shutil.move(f'{result_path}/{strain_name}',f'{data_f_path}/{strain_name}')

        # Augustusのconfig/speciesからstrain_nameに該当するフォルダを削除
        shutil.rmtree(f'{augustus_path}/config/species/{strain_name}') 