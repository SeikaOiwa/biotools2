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
result_path = parent_dir

# 必要に応じて調整！
augustus_path = f"{base_path}/AUGUSTUS"
data_f_path = "/home/seika_oiwa_3590/notebooks/box_path/bio_team/GeneAnnotation"

st.title('遺伝子領域の推定')
st.image(f'{base_path}/file/image.png')

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

    if data2[1] == "genemark":
        return "Finish_GeneMark",data2[1]
    
    if data2[1] == "Augustus":
        return "Finish_All",data2[1]

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
        s_time = time.time()
        # path設定
        genome_path = f'{data_f_path}/{select_genome}'

        if fungus_check:
            # GeneMark/Augustus予測(Fungusモード)
            subprocess.run([f'{base_path}/scr/braker_1_f.sh {strain_name} {genome_path}'],stdout=subprocess.PIPE,shell=True)
        else:
            subprocess.run([f'{base_path}/scr/braker_1.sh {strain_name} {genome_path}'],stdout=subprocess.PIPE,shell=True)
        
        # 進捗状況
        if os.path.isfile(f'{result_path}/braker/GeneMark-ET/genemark.gtf'):
            f_g_time = time.time()
            dgtime = round(( s_time - f_g_time)/60,1)
            comment = f"genemark,{dgtime}"
            input_progress(base_path,comment)
            
        
        if os.path.isfile(f'{result_path}/braker/braker.gtf'):
            f_a_time = time.time()
            dgtime2 = round(( s_time - f_a_time)/60,1)
            comment = f"Augustus,{dgtime2}"
            input_progress(base_path,comment)
     
            # フォルダ名＝菌株名に修正
            os.rename(f'{result_path}/braker',f'{result_path}/{strain_name}')
        
        # フォルダをdata_f_pathに移動
        if os.path.isdir(f'{data_f_path}/{strain_name}'):
            shutil.rmtree(f'{data_f_path}/{strain_name}')
            shutil.move(f'{result_path}/{strain_name}',f'{data_f_path}/{strain_name}')
        
            # Augustusのconfig/speciesからstrain_nameに該当するフォルダを削除
            shutil.rmtree(f'{augustus_path}/config/species/{strain_name}') 

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
    select_bam = st.selectbox('ゲノムデータを選択',bam_list)    

    # 菌株名
    strain_name = st.text_input("菌株名を入力","your_starain_name")  

    # カビの解析オプション有無
    fungus_check =  st.checkbox('カビ解析モード')

    # 解析開始
    if st.button('開始'):
        # 開始時間計測
        s_time = time.time()
        # path設定
        genome_path = f'{data_f_path}/{select_genome}'
        bam_path = f'{data_f_path}/{select_bam}'

        if fungus_check:
            # GeneMark/Augustus予測(Fungusモード)
            subprocess.run([f'{base_path}/scr/braker_2_f.sh {strain_name} {genome_path} {bam_path}'],stdout=subprocess.PIPE,shell=True)
        else:
            subprocess.run([f'{base_path}/scr/braker_2.sh {strain_name} {genome_path} {bam_path}'],stdout=subprocess.PIPE,shell=True)
        
        # 進捗状況
        if os.path.isfile(f'{result_path}/braker/GeneMark-ET/genemark.gtf'):
            f_g_time = time.time()
            dgtime = round(( s_time - f_g_time)/60,1)
            comment = f"genemark,{dgtime}"
            input_progress(base_path,comment)
            
        
        if os.path.isfile(f'{result_path}/braker/braker.gtf'):
            f_a_time = time.time()
            dgtime2 = round(( s_time - f_a_time)/60,1)
            comment = f"Augustus,{dgtime2}"
            input_progress(base_path,comment)
     
            # フォルダ名＝菌株名に修正
            os.rename(f'{result_path}/braker',f'{result_path}/{strain_name}')
        
        # フォルダをdata_f_pathに移動
        if os.path.isdir(f'{data_f_path}/{strain_name}'):
            shutil.rmtree(f'{data_f_path}/{strain_name}')
            shutil.move(f'{result_path}/{strain_name}',f'{data_f_path}/{strain_name}')
        
            # Augustusのconfig/speciesからstrain_nameに該当するフォルダを削除
            shutil.rmtree(f'{augustus_path}/config/species/{strain_name}') 

def genome_protein_analysis(base_path,data_f_path,augustus_path,result_path):
    """ゲノム情報、近縁種のプロテイン情報を用いた遺伝子領域予測"""
    st.write('ゲノム情報、近縁種のプロテイン情報を用いた遺伝子領域予測')

    # ゲノムデータの選択
    st.write(f'{data_f_path}にゲノムデータ(*.fasta)を保管して下さい')
    genome_lists_ = glob.glob(f"{data_f_path}/*.fasta")
    genome_list = [os.path.basename(i) for i in genome_lists_]
    select_genome = st.selectbox('ゲノムデータを選択',genome_list)

    # proteinデータの選択
    st.write(f'{data_f_path}にproteinデータ(*.fasta)を保管して下さい')
    protein_lists_  = glob.glob(f"{data_f_path}/*.fasta")
    protein_list = [os.path.basename(i) for i in protein_lists_]
    select_protein = st.selectbox('ゲノムデータを選択',protein_list)    

    # 菌株名
    strain_name = st.text_input("菌株名を入力","your_starain_name")  

    # カビの解析オプション有無
    fungus_check =  st.checkbox('カビ解析モード')

    # 解析開始
    if st.button('開始'):
        # 開始時間計測
        s_time = time.time()
        # path設定
        genome_path = f'{data_f_path}/{select_genome}'
        protein_path = f'{data_f_path}/{select_protein}'

        if fungus_check:
            # GeneMark/Augustus予測(Fungusモード)
            subprocess.run([f'{base_path}/scr/braker_3_f.sh {strain_name} {genome_path} {protein_path}'],stdout=subprocess.PIPE,shell=True)
        else:
            subprocess.run([f'{base_path}/scr/braker_3.sh {strain_name} {genome_path} {protein_path}'],stdout=subprocess.PIPE,shell=True)
        
        # 進捗状況
        if os.path.isfile(f'{result_path}/braker/GeneMark-ET/genemark.gtf'):
            f_g_time = time.time()
            dgtime = round(( s_time - f_g_time)/60,1)
            comment = f"genemark,{dgtime}"
            input_progress(base_path,comment)
            
        
        if os.path.isfile(f'{result_path}/braker/braker.gtf'):
            f_a_time = time.time()
            dgtime2 = round(( s_time - f_a_time)/60,1)
            comment = f"Augustus,{dgtime2}"
            input_progress(base_path,comment)
     
            # フォルダ名＝菌株名に修正
            os.rename(f'{result_path}/braker',f'{result_path}/{strain_name}')
        
        # フォルダをdata_f_pathに移動
        if os.path.isdir(f'{data_f_path}/{strain_name}'):
            shutil.rmtree(f'{data_f_path}/{strain_name}')
            shutil.move(f'{result_path}/{strain_name}',f'{data_f_path}/{strain_name}')
        
            # Augustusのconfig/speciesからstrain_nameに該当するフォルダを削除
            shutil.rmtree(f'{augustus_path}/config/species/{strain_name}') 

def genome_rna_protein_analysis(base_path,data_f_path,augustus_path,result_path):
    """ゲノム情報、RNAseqデータを用いた遺伝子領域予測"""
    st.write('ゲノム情報、RNAseqデータを用いた遺伝子領域予測')

    # ゲノムデータの選択
    st.write(f'{data_f_path}にゲノムデータ(*.fasta)を保管して下さい')
    genome_lists_ = glob.glob(f"{data_f_path}/*.fasta")
    genome_list = [os.path.basename(i) for i in genome_lists_]
    select_genome = st.selectbox('ゲノムデータを選択',genome_list)

    # RNAseqデータの選択
    st.write(f'{data_f_path}にRNAseqデータ(*.bam)を保管して下さい')
    bam_lists_  = glob.glob(f"{data_f_path}/*.bam")
    bam_list = [os.path.basename(i) for i in bam_lists_]
    select_bam = st.selectbox('ゲノムデータを選択',bam_list)    

    # proteinデータの選択
    st.write(f'{data_f_path}にproteinデータ(*.fasta)を保管して下さい')
    protein_lists_  = glob.glob(f"{data_f_path}/*.fasta")
    protein_list = [os.path.basename(i) for i in protein_lists_]
    select_protein = st.selectbox('ゲノムデータを選択',protein_list)  

    # 菌株名
    strain_name = st.text_input("菌株名を入力","your_starain_name")  

    # カビの解析オプション有無
    fungus_check =  st.checkbox('カビ解析モード')

    # 解析開始
    if st.button('開始'):
        # 開始時間計測
        s_time = time.time()
        # path設定
        genome_path = f'{data_f_path}/{select_genome}'
        bam_path = f'{data_f_path}/{select_bam}'
        protein_path = f'{data_f_path}/{select_protein}'

        if fungus_check:
            # GeneMark/Augustus予測(Fungusモード)
            subprocess.run([f'{base_path}/scr/braker_4_f.sh {strain_name} {genome_path} {bam_path} {protein_path}'],stdout=subprocess.PIPE,shell=True)
        else:
            subprocess.run([f'{base_path}/scr/braker_4.sh {strain_name} {genome_path} {bam_path} {protein_path}'],stdout=subprocess.PIPE,shell=True)
        
        # 進捗状況
        if os.path.isfile(f'{result_path}/braker/GeneMark-ET/genemark.gtf'):
            f_g_time = time.time()
            dgtime = round(( s_time - f_g_time)/60,1)
            comment = f"genemark,{dgtime}"
            input_progress(base_path,comment)
            
        
        if os.path.isfile(f'{result_path}/braker/braker.gtf'):
            f_a_time = time.time()
            dgtime2 = round(( s_time - f_a_time)/60,1)
            comment = f"Augustus,{dgtime2}"
            input_progress(base_path,comment)
     
            # フォルダ名＝菌株名に修正
            os.rename(f'{result_path}/braker',f'{result_path_path}/{strain_name}')
        
        # フォルダをdata_f_pathに移動
        if os.path.isdir(f'{data_f_path}/{strain_name}'):
            shutil.rmtree(f'{data_f_path}/{strain_name}')
            shutil.move(f'{result_path}/{strain_name}',f'{data_f_path}/{strain_name}')
        
            # Augustusのconfig/speciesからstrain_nameに該当するフォルダを削除
            shutil.rmtree(f'{augustus_path}/config/species/{strain_name}') 

select = st.sidebar.radio('モード選択',('ゲノム情報','ゲノム情報＋RNAseq','ゲノム情報＋近縁種のProtein情報','ゲノム情報＋RNAseq＋近縁種のProtein情報'))

select2 = st.sidebar.button('進捗状況確認')

if select == 'ゲノム情報':
    genome_analysis(base_path,data_f_path,augustus_path,result_path)

if select == 'ゲノム情報＋RNAseq':
    genome_rna_analysis(base_path,data_f_path,augustus_path,result_path)

if select == 'ゲノム情報＋近縁種のProtein情報':
    genome_protein_analysis(base_path,data_f_path,augustus_path,result_path)   

if select == 'ゲノム情報＋近縁種のProtein情報':
    genome_rna_protein_analysis(base_path,data_f_path,augustus_path,result_path)  

if select2 == "進捗状況確認":
    r1,r2 = check_progress()
    st.write(f"{r1}:所要時間({r2})")