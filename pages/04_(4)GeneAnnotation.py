import streamlit as st
import subprocess
import glob
import os
import time
import math
import shutil
from PIL import GifImagePlugin

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.abspath(os.path.join(current_dir,os.pardir))
base_path = f"{parent_dir}/GeneAnnotation"
script_path = f"{base_path}/scr"
result_path = parent_dir

# 必要に応じて調整！
augustus_path = f"{base_path}/Augustus"
data_f_path = "/home/seika_oiwa_3590/notebooks/box_path/bio_team/GeneAnnotation"

st.title('遺伝子領域の推定')
st.image(f'{base_path}/file/front_image.png')

def check_progress(base_path):
    """進捗状況が記載されたテキストファイルを読み込み、結果をstreamlit上に表示
    Parameter:
    ---------
    base_path: str
        テキストファイルの書出し先
    Returns:
    ------
    "｛進捗｝,｛所要時間(h)｝"
    """
    with open(f"{base_path}/progress.txt") as f:
        data = f.read()
        data2 = data.split(",")

    if data2[0] == "genemark":
        return "Finish_GeneMark",data2[1]
    
    if data2[0] == "Augustus":
        return "アノテーション終了！！！",data2[1]

    else:
        return "解析中",""

def genome_analysis(base_path,data_f_path,augustus_path,result_path):
    """ゲノム情報のみを用いた遺伝子領域予測"""

    st.write('ゲノム情報のみを用いた遺伝子領域予測')

    # ゲノムデータの選択
    st.write(f'{data_f_path}にゲノムデータ(*.fasta)を保管して下さい')
    genome_lists_ = glob.glob(f"{data_f_path}/*.fasta")
    genome_list = [os.path.basename(i) for i in genome_lists_]
    select_genome = st.selectbox('ゲノムデータを選択',genome_list)

    # 菌株名
    strain_name = st.text_input("菌株名を入力","your_starain_name")  

    # カビの解析オプション有無
    fungus_check =  st.checkbox('カビ解析モード')

    # 解析開始
    if st.button('開始'):
        # 開始時間計測
        # path設定
        genome_path = f'{data_f_path}/{select_genome}'

        if fungus_check:
            # GeneMark/Augustus予測(Fungusモード)
            subprocess.Popen(["python",f"{base_path}/scr/braker.py","g_f",base_path,augustus_path,strain_name,result_path,data_f_path,genome_path,"",""],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        else:
            subprocess.Popen(["python",f"{base_path}/scr/braker.py","g",base_path,augustus_path,strain_name,result_path,data_f_path,genome_path,"",""],stdout=subprocess.PIPE,stderr=subprocess.PIPE)

def genome_rna_analysis(base_path,data_f_path,augustus_path,result_path):
    """ゲノム情報、RNAseqデータを用いた遺伝子領域予測"""

    st.write('ゲノム情報、RNAseqデータを用いた遺伝子領域予測')

    # ゲノムデータの選択
    st.write(f'{data_f_path}にゲノムデータ(*.fasta)を保管して下さい')
    genome_lists_ = glob.glob(f"{data_f_path}/*.fasta")
    genome_list = [os.path.basename(i) for i in genome_lists_]
    select_genome = st.selectbox('ゲノムデータを選択',genome_list)

    # RNAseqデータの選択
    st.write(f'{data_f_path}にRNAseqデータ(*.bam)を保管して下さい')
    bam_lists_ = glob.glob(f"{data_f_path}/*.bam")
    bam_list = [os.path.basename(i) for i in bam_lists_]
    select_bam = st.multiselect('RNAseqデータを選択',bam_list)    

    # 菌株名
    strain_name = st.text_input("菌株名を入力","your_starain_name")  

    # カビの解析オプション有無
    fungus_check =  st.checkbox('カビ解析モード')

    # 解析開始
    if st.button('開始'):

        # path設定
        genome_path = f'{data_f_path}/{select_genome}'

        paths = f"{data_f_path}/"
        new_list = [os.path.join(paths,i) for i in select_bam]
        bam_path = ",".join(new_list)

        if fungus_check:
            # GeneMark/Augustus予測(Fungusモード)
            subprocess.Popen(["python",f"{base_path}/scr/braker.py","gr_f",base_path,augustus_path,strain_name,result_path,data_f_path,genome_path,bam_path,""],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        else:
            subprocess.Popen(["python",f"{base_path}/scr/braker.py","gr",base_path,augustus_path,strain_name,result_path,data_f_path,genome_path,bam_path,""],stdout=subprocess.PIPE,stderr=subprocess.PIPE)

def genome_protein_analysis(base_path,data_f_path,augustus_path,result_path):
    """ゲノム情報、近縁種のプロテイン情報を用いた遺伝子領域予測"""
    st.write('ゲノム情報、近縁種のプロテイン情報を用いた遺伝子領域予測')

    # ゲノムデータの選択
    st.write(f'{data_f_path}にゲノムデータ(*.fasta)を保管して下さい')
    genome_lists_ = glob.glob(f"{data_f_path}/*.fasta")
    genome_list = [os.path.basename(i) for i in genome_lists_]
    select_genome = st.selectbox('ゲノムデータを選択',genome_list,key="genome")

    # proteinデータの選択
    st.write(f'{data_f_path}にproteinデータ(*.fa)を保管して下さい')
    protein_lists_  = glob.glob(f"{data_f_path}/*.fa")
    protein_list = [os.path.basename(i) for i in protein_lists_]
    select_protein = st.selectbox('proteinデータを選択',protein_list,key="protein")    

    # 菌株名
    strain_name = st.text_input("菌株名を入力","your_starain_name")  

    # カビの解析オプション有無
    fungus_check =  st.checkbox('カビ解析モード')

    # 解析開始
    if st.button('開始'):

        # path設定
        genome_path = f'{data_f_path}/{select_genome}'
        protein_path = f'{data_f_path}/{select_protein}'

        if fungus_check:
            # GeneMark/Augustus予測(Fungusモード)
            subprocess.Popen(["python",f"{base_path}/scr/braker.py","gp_f",base_path,augustus_path,strain_name,result_path,data_f_path,genome_path,"",protein_path],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        else:
            subprocess.Popen(["python",f"{base_path}/scr/braker.py","gp",base_path,augustus_path,strain_name,result_path,data_f_path,genome_path,"",protein_path],stdout=subprocess.PIPE,stderr=subprocess.PIPE)

def genome_rna_protein_analysis(base_path,data_f_path,augustus_path,result_path):
    """ゲノム情報、RNAseqデータを用いた遺伝子領域予測"""
    st.write('ゲノム情報、RNAseqデータを用いた遺伝子領域予測')

    # ゲノムデータの選択
    st.write(f'{data_f_path}にゲノムデータ(*.fasta)を保管して下さい')
    genome_lists_ = glob.glob(f"{data_f_path}/*.fasta")
    genome_list = [os.path.basename(i) for i in genome_lists_]
    select_genome = st.selectbox('ゲノムデータを選択',genome_list,key="genome2")

    # RNAseqデータの選択
    st.write(f'{data_f_path}にRNAseqデータ(*.bam)を保管して下さい')
    bam_lists_  = glob.glob(f"{data_f_path}/*.bam")
    bam_list = [os.path.basename(i) for i in bam_lists_]
    select_bam = st.multiselect('RNAseqデータを選択',bam_list)     

    # proteinデータの選択
    st.write(f'{data_f_path}にproteinデータ(*.fa)を保管して下さい')
    protein_lists_  = glob.glob(f"{data_f_path}/*.fa")
    protein_list = [os.path.basename(i) for i in protein_lists_]
    select_protein = st.selectbox('proteinデータを選択',protein_list,key="protein2")  

    # 菌株名
    strain_name = st.text_input("菌株名を入力","your_starain_name")  

    # カビの解析オプション有無
    fungus_check =  st.checkbox('カビ解析モード')

    # 解析開始
    if st.button('開始'):
        # path設定
        genome_path = f'{data_f_path}/{select_genome}'

        paths = f"{data_f_path}/"
        new_list = [os.path.join(paths,i) for i in select_bam]
        bam_path = ",".join(new_list)

        protein_path = f'{data_f_path}/{select_protein}'

        if fungus_check:
            # GeneMark/Augustus予測(Fungusモード)
            subprocess.Popen(["python",f"{base_path}/scr/braker.py","gpr_f",base_path,augustus_path,strain_name,result_path,data_f_path,genome_path,bam_path,protein_path
            ],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        else:
            subprocess.Popen(["python",f"{base_path}/scr/braker.py","gpr",base_path,augustus_path,strain_name,result_path,data_f_path,genome_path,bam_path,protein_path],stdout=subprocess.PIPE,stderr=subprocess.PIPE)

select = st.sidebar.selectbox('モード選択',['-選択-','ゲノム情報','ゲノム情報＋RNAseq','ゲノム情報＋近縁種のProtein情報','ゲノム情報＋RNAseq＋近縁種のProtein情報'])

select2 = st.sidebar.button('進捗状況確認')

if select == 'ゲノム情報':
    genome_analysis(base_path,data_f_path,augustus_path,result_path)

if select == 'ゲノム情報＋RNAseq':
    genome_rna_analysis(base_path,data_f_path,augustus_path,result_path)

if select == 'ゲノム情報＋近縁種のProtein情報':
    genome_protein_analysis(base_path,data_f_path,augustus_path,result_path)   

if select == 'ゲノム情報＋RNAseq＋近縁種のProtein情報':
    genome_rna_protein_analysis(base_path,data_f_path,augustus_path,result_path)  

if select2:
    r1,r2 = check_progress(base_path)
    st.write(f'{r1} : 所要時間({r2})')
