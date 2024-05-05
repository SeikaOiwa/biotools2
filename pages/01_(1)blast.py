import streamlit as st
import pandas as pd
import subprocess
import os
from time import sleep
from Bio.Blast import Applications
from Bio.Seq import Seq
import shutil
import glob
from io import BytesIO
from PIL import Image

#file path
current_dir = os.path.dirname(os.path.abspath(__file__))
base_path_  = os.path.abspath(os.path.join(current_dir,os.pardir))
base_path = f'{base_path_}/local_blast'
database_path = f'{base_path}/Gfile'

def main(): 
    st.image(f'{base_path}/File/front.png',use_column_width=True)   
    st.sidebar.subheader('STEP1:データベース構築')
    st.sidebar.write('※構築済みの場合はSTEP2に進む')
    select = st.sidebar.selectbox('データベース種類：', ['-未選択-','ゲノムDNA','遺伝子リスト','アミノ酸リスト','Option_アノテーション'])
    
    st.sidebar.subheader('STEP2：Local_Blast検索')
    select2 = st.sidebar.selectbox('検索モード：', ['-未選択-','blastN','blastX','tblastX','blastP','primer']) 
    
    st.sidebar.subheader('STEP3：不要ファイルの消去') 
    st.sidebar.write('検索ファイル、検索結果データの削除')
    select3 = st.sidebar.button('消去') 
    
    if select == 'ゲノムDNA' : 
        make_genomeDB()
    if select == '遺伝子リスト':
        make_gene_listDB()
    if select == 'アミノ酸リスト' :
        make_protein_listDB()
    
    if select == 'Option_アノテーション':
        anno()
        
    if select2 == 'blastN':
        blastN()
    if select2 == 'blastX':
        blastX()
    if select2 == 'tblastX':
        tblastX()        
    if select2 == 'blastP':
        blastP() 
    if select2 == 'primer':
        primer()

    if select3:
        cleanup_data() 
        
# ゲノムDNAのDB化        
def make_genomeDB():
    st.header('ゲノムDNAからのDB構築')
    st.write('※データベース構築済みの場合は、STEP2へ進む')
    st.subheader('(1): 菌株名を入力してください')
    strain_name = st.text_input(label='-入力-')
    st.subheader('(2): データ準備') 
    st.write('「Gfile」にゲノム情報(.fasta)を保管してください')
    st.subheader('(3): データベース構築')     
    makedb = st.button(label='Go') 
    if makedb:
        os.makedirs(f'{database_path}/{strain_name}/genome',exist_ok=True)
        db = Applications.NcbimakeblastdbCommandline(input_file =f'{database_path}/*.fasta', dbtype = 'nucl',out =                                             f'{database_path}/{strain_name}/genome/gene')
        stdout, stderr = db()
        sleep(5)
        is_file = os.path.isfile(f'{database_path}/{strain_name}/genome/gene.nhr')

        if is_file:
            for i in glob.glob(f'{database_path}/*.fasta'):
                shutil.move(i,base_path)
            st.write(f'データベース構築完了！(データベース名：{strain_name})')
            st.write('STEP2(Local_Blast検索)に進んでください')
            
# 遺伝子リストのDB化                
def make_gene_listDB():
    st.header('遺伝子リストからのDB構築')
    st.write('※データベース構築済みの場合は、STEP2へ進む')
    st.subheader('(1): 菌株名を入力してください')
    strain_name = st.text_input(label='-入力-')
    st.subheader('(2): データ準備') 
    st.write('「Gfile」にゲノム情報(.fasta)を保管してください') 
    st.write('「Gfile」にアノテーションデータ(.gtf)を保管してください') 
    st.subheader('(3): データベース構築')     
    makedb = st.button(label='Go') 
    if makedb:
        os.makedirs(f'{database_path}/{strain_name}/gene_list',exist_ok=True)
        subprocess.run(f'gffread {database_path}/*.gtf -w {database_path}/gene_list.fasta -g {database_path}/*.fasta ',shell=True)
        db = Applications.NcbimakeblastdbCommandline(input_file =f'{database_path}/gene_list.fasta', dbtype = 'nucl',out =                                             f'{database_path}/{strain_name}/gene_list/gene')
        stdout, stderr = db()
        sleep(2)
        is_file = os.path.isfile(f'{database_path}/{strain_name}/gene_list/gene.nhr')           
        
        if is_file:
            for i in glob.glob(f'{database_path}/*.fasta'):
                shutil.move(i,base_path)
            for i2 in glob.glob(f'{database_path}/*.fai'):
                os.remove(i2)
            for i3 in glob.glob(f'{database_path}/*.gtf'):
                shutil.move(i3,base_path)  
            st.write(f'データベース構築完了！(データベース名：{strain_name})')
            st.write('STEP2(Local_Blast検索)に進んでください')
            
# アミノ酸リストのDB化             
def make_protein_listDB():
    st.header('アミノ酸リストからのDB構築')
    st.write('※データベース構築済みの場合は、STEP2へ進む')
    st.subheader('(1): 菌株名を入力してください')
    strain_name = st.text_input(label='-入力-')
    st.subheader('(2): データ準備') 
    st.write('「Gfile」にゲノム情報(.fasta)を保管してください') 
    st.write('「Gfile」にアノテーションデータ(.gtf)を保管してください') 
    st.subheader('(3): データベース構築')     
    makedb = st.button(label='Go') 
    if makedb:
        os.makedirs(f'{database_path}/{strain_name}/amino_list',exist_ok=True)
        subprocess.run(f'gffread {database_path}/*.gtf -y {database_path}/protein_list.fasta -g {database_path}/*.fasta',shell=True)
        db = Applications.NcbimakeblastdbCommandline(input_file =f'{database_path}/protein_list.fasta', dbtype = 'prot',out =                                             f'{database_path}/{strain_name}/amino_list/gene')
        stdout, stderr = db()
        sleep(2)
        is_file = os.path.isfile(f'{database_path}/{strain_name}/amino_list/gene.phr')
 
        if is_file:
            for i in glob.glob(f'{database_path}/*.fasta'):
                shutil.move(i,base_path)
            for i2 in glob.glob(f'{database_path}/*.fai'):
                os.remove(i2)
            for i3 in glob.glob(f'{database_path}/*.gtf'):
                shutil.move(i3,base_path) 
            st.write(f'データベース構築完了！(データベース名：{strain_name})')
            st.write('STEP2(Local_Blast検索)に進んでください')

# アノテーション情報の保管
def anno():
    st.header('Option_アノテーション情報の保管')
    db_name = os.listdir(f'{database_path}/')
    select_db_name = st.selectbox('菌株を選択してください', db_name[0:])
    file_name = st.text_input('ファイル名を入力')    
    anno_data = st.file_uploader('データアップロード',type='xlsx')
    h3 = st.number_input('ヘッダー位置調整',0,key='d')
    data_file = pd.read_excel(anno_data,header=h3)
    st.dataframe(data_file.head(5))
    go = st.button('セット')
    if go:
        os.makedirs(f'{database_path}/{select_db_name}/アノテーション',exist_ok=True)
        data_file.to_excel(f'{database_path}/{select_db_name}/アノテーション/{file_name}.xlsx')
        is_file = os.path.isfile(f'{database_path}/{select_db_name}/アノテーション/{file_name}.xlsx')
        if is_file:
            st.write(f'アノテーション情報の保管完了！')
            st.write('STEP2(Local_Blast検索)に進んでください')
                                    
def input_seq(fname,seq):
    with open(f'{base_path}/検索配列/{fname}.txt',mode = 'w')as f:
        f.write(f'>{fname}\n'+seq)
    shutil.move(f'{base_path}/検索配列/{fname}.txt',f'{base_path}/検索配列/{fname}.fasta')
    
def cleanup_data():
    st.image(f'{base_path}/File/end.png',use_column_width=True)
    shutil.rmtree(f'{base_path}/検索配列')
    shutil.rmtree(f'{base_path}/Result')
    os.makedirs(f'{base_path}/検索配列')
    os.makedirs(f'{base_path}/Result')   
    
def blastN():
    st.header('BlastNの実行')
    st.write('NCBI：https://www.ncbi.nlm.nih.gov/')    
    #sequence information     
    st.subheader('(1)検索配列情報の入力')
    fname = st.text_input(label='出力ファイル名(任意)',value = '解析') 
    seq = st.text_area('検索配列の入力（fasta形式）')
    upload_seq = st.button('アップロード') 
    if upload_seq:
        with open(f'{base_path}/検索配列/{fname}.txt',mode = 'w')as f:
            f.write(seq)
            shutil.move(f'{base_path}/検索配列/{fname}.txt',f'{base_path}/検索配列/{fname}.fasta')  
           
    #database infomation
    st.subheader('(2)データベースの選択')  
    db_name_ = os.listdir(f'{database_path}/')
    db_name = [f for f in db_name_  if not f == '.DS_Store']
    st.write('菌株を選択してください')
    select_db_name = st.selectbox('', db_name)
    db_types_ = os.listdir(f'{database_path}/{select_db_name}')
    db_types = [f for f in db_types_ if not f == '.DS_Store']
    st.write('データベースの種類を選択してください')
    select_db_type = st.radio('', db_types,horizontal=True)

    #database path
    db_path = f'{database_path}/{select_db_name}/{select_db_type}/gene'
    
    #コメント
    with st.expander("データベースタイプ", expanded=False):
        st.write('genome = ゲノムDNAから構築')
        st.write('gene_list = 遺伝子リストから構築')
        st.write('amino_list = アミノ酸リストから構築')
    
    #検索方式の設定
    st.subheader('(3)パラメーター設定')
    Ev = st.text_input(label='e-value値',value='1e-2')   
    
    #option
    st.subheader('option_RNAseq,functionデータ付加')
    option_use = st.radio('選択',['未使用','使用'],horizontal=True)
    st.write('※データベース：gene_list、amino_listを使用時のみ！')
    if option_use == '使用':
        #アノテーションファイルの選択
        name_list = []
        for name in glob.glob(f'{database_path}/{select_db_name}/アノテーション/*.csv'):
            list = os.path.basename(name)
            name_list.append(list)
        select_anno_file = st.radio('アノテーションファイルの選択',name_list)
        name_1 = st.text_input('データ統合に使用する列名の設定')
        st.write('参考')
        h3 = st.number_input('ヘッダー位置調整',0,key='d')
        data_file = pd.read_excel(f'{database_path}/{select_db_name}/アノテーション/{select_anno_file}',header=h3)
        st.dataframe(data_file.head(5))
    
    #Blast検索の実行
    st.subheader('(4)Blast検索の実行')
    go = st.button('実行')
    if go:
        blast = Applications.NcbiblastnCommandline(query=f'{base_path}/検索配列/{fname}.fasta', db = db_path,evalue = Ev, out =f'{base_path}/Result/{fname}.txt',outfmt=0)
        stdout, stderr = blast()
        blast = Applications.NcbiblastnCommandline(query=f'{base_path}/検索配列/{fname}.fasta', db = db_path,evalue = Ev, out =f'{base_path}/Result/{fname}2.txt',outfmt=('6 qseqid qstart qend sseqid sstart send length pident bitscore evalue sseq'))
        stdout, stderr = blast()

        sleep(5)  

    #解析データの表示
        df = pd.read_table(f'{base_path}/Result/{fname}2.txt',header=None)
        df.columns = [['検索配列名','開始位置','終了位置','ヒット領域名','開始位置_','終了位置_','マッチ数','%','Bitscore','Evalue','ヒット領域配列']]
        df.to_excel(f'{base_path}/Result/{fname}2.xlsx',index=None)
        df2 = pd.read_excel(f'{base_path}/Result/{fname}2.xlsx')
        st.write('検索結果')
        if option_use == '使用':
            df3 = pd.merge(df2,data_file,left_on ='ヒット領域名',right_on =name_1,how='left' )
            st.dataframe(df3)
            df3.to_excel(buf := BytesIO(), index=True)
            st.download_button("Download",buf.getvalue(),f"{fname}.xlsx","application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
        if option_use == '未使用':
            st.dataframe(df2)
            df2.to_excel(buf := BytesIO(), index=True)
            st.download_button("Download",buf.getvalue(),f"{fname}.xlsx","application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")

        # 検索配列の除去    
        shutil.rmtree(f'{base_path}/検索配列')   
        os.makedirs(f'{base_path}/検索配列')      

def blastX():
    st.header('blastXの実行')
    st.write('NCBI：https://www.ncbi.nlm.nih.gov/')
    #sequence information     
    st.subheader('(1)検索配列情報の入力')
    fname = st.text_input(label='出力ファイル名(任意)',value = '解析') 
    seq = st.text_area('検索配列の入力（fasta形式）')
    upload_seq = st.button('アップロード') 
    if upload_seq:
        with open(f'{base_path}/検索配列/{fname}.txt',mode = 'w')as f:
            f.write(seq)
            shutil.move(f'{base_path}/検索配列/{fname}.txt',f'{base_path}/検索配列/{fname}.fasta') 
            
    #database infomation
    st.subheader('(2)データベースの選択')  
    db_name_ = os.listdir(f'{database_path}/')
    db_name = [f for f in db_name_  if not f == '.DS_Store']
    st.write('菌株を選択してください')
    select_db_name = st.selectbox('', db_name)
    db_types_ = os.listdir(f'{database_path}/{select_db_name}')
    db_types = [f for f in db_types_ if not f == '.DS_Store']
    st.write('データベースの種類を選択してください')
    select_db_type = st.radio('', db_types,horizontal=True)

    #database path
    db_path = f'{database_path}/{select_db_name}/{select_db_type}/gene'
    
    #コメント
    with st.expander("データベースタイプ", expanded=False):
        st.write('genome = ゲノムDNAから構築')
        st.write('gene_list = 遺伝子リストから構築')
        st.write('amino_list = アミノ酸リストから構築')
    
    #検索方式の設定
    st.subheader('(3)パラメーター設定')
    Ev = st.text_input(label='e-value値',value='1e-2')   

    #option
    st.subheader('option_RNAseq,functionデータ付加')
    option_use = st.radio('選択',['未使用','使用'],horizontal=True)
    st.write('※データベース：gene_list、amino_listを使用時のみ！')
    if option_use == '使用':
        #アノテーションファイルの選択
        name_list = []
        for name in glob.glob(f'{database_path}/{select_db_name}/アノテーション/*.xlsx'):
            list = os.path.basename(name)
            name_list.append(list)
        select_anno_file = st.radio('アノテーションファイルの選択',name_list)
        name_1 = st.text_input('データ統合に使用する列名の設定')
        st.write('参考')
        h3 = st.number_input('ヘッダー位置調整',0,key='d')
        data_file = pd.read_excel(f'{database_path}/{select_db_name}/アノテーション/{select_anno_file}',header=h3)
        st.dataframe(data_file.head(5))  
    
    #Blast検索の実行
    st.subheader('(4)Blast検索の実行')
    go = st.button('実行')
    if go:
        blast = Applications.NcbiblastxCommandline(query=f'{base_path}/検索配列/{fname}.fasta', db =db_path,evalue = Ev, out =f'{base_path}/Result/{fname}.txt',outfmt=0)
        stdout, stderr = blast()
        blast = Applications.NcbiblastxCommandline(query=f'{base_path}/検索配列/{fname}.fasta', db =db_path,evalue = Ev, out =f'{base_path}/Result/{fname}2.txt',outfmt=('6 qseqid qstart qend sseqid sstart send length pident bitscore evalue sseq'))

        sleep(5)  

    #解析データの表示
        df = pd.read_table(f'{base_path}/Result/{fname}2.txt',header=None)
        df.columns = [['検索配列名','開始位置','終了位置','ヒット領域名','開始位置_','終了位置_','マッチ数','%','Bitscore','Evalue','ヒット領域配列']]
        df.to_excel(f'{base_path}/Result/{fname}2.xlsx',index=None)
        df2 = pd.read_csv(f'{base_path}/Result/{fname}2.xlsx')
        st.write('検索結果')
        if option_use == '使用':
            df3 = pd.merge(df2,data_file,left_on ='ヒット領域名',right_on =name_1,how='left' )
            st.dataframe(df3)
            df3.to_excel(buf := BytesIO(), index=True)
            st.download_button("Download",buf.getvalue(),f"{fname}.xlsx","application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
        if option_use == '未使用':
            st.dataframe(df2)
            df2.to_excel(buf := BytesIO(), index=True)
            st.download_button("Download",buf.getvalue(),f"{fname}.xlsx","application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
        
        # 検索配列の除去    
        shutil.rmtree(f'{base_path}/検索配列')   
        os.makedirs(f'{base_path}/検索配列')              
        
def tblastX():
    st.header('tblastXの実行')
    st.write('NCBI：https://www.ncbi.nlm.nih.gov/')
    #sequence information     
    st.subheader('(1)検索配列情報の入力')
    fname = st.text_input(label='出力ファイル名(任意)',value = '解析') 
    seq = st.text_area('検索配列の入力（fasta形式）')
    upload_seq = st.button('アップロード') 
    if upload_seq:
        with open(f'{base_path}/検索配列/{fname}.txt',mode = 'w')as f:
            f.write(seq)
            shutil.move(f'{base_path}/検索配列/{fname}.txt',f'{base_path}/検索配列/{fname}.fasta') 
            
    #database infomation
    st.subheader('(2)データベースの選択')  
    db_name_ = os.listdir(f'{database_path}/')
    db_name = [f for f in db_name_  if not f == '.DS_Store']
    st.write('菌株を選択してください')
    select_db_name = st.selectbox('', db_name)
    db_types_ = os.listdir(f'{database_path}/{select_db_name}')
    db_types = [f for f in db_types_ if not f == '.DS_Store']
    st.write('データベースの種類を選択してください')
    select_db_type = st.radio('', db_types,horizontal=True)

    #database path
    db_path = f'{database_path}/{select_db_name}/{select_db_type}/gene'
    
    #コメント
    with st.expander("データベースタイプ", expanded=False):
        st.write('genome = ゲノムDNAから構築')
        st.write('gene_list = 遺伝子リストから構築')
        st.write('amino_list = アミノ酸リストから構築')
    
    #検索方式の設定
    st.subheader('(3)パラメーター設定')
    Ev = st.text_input(label='e-value値',value='1e-2')   

    #option
    st.subheader('option_RNAseq,functionデータ付加')
    option_use = st.radio('選択',['未使用','使用'],horizontal=True)
    st.write('※データベース：gene_list、amino_listを使用時のみ！')
    if option_use == '使用':
        #アノテーションファイルの選択
        name_list = []
        for name in glob.glob(f'{database_path}/{select_db_name}/アノテーション/*.xlsx'):
            list = os.path.basename(name)
            name_list.append(list)
        select_anno_file = st.radio('アノテーションファイルの選択',name_list)
        name_1 = st.text_input('データ統合に使用する列名の設定')
        st.write('参考')
        h3 = st.number_input('ヘッダー位置調整',0,key='d')
        data_file = pd.read_excel(f'{database_path}/{select_db_name}/アノテーション/{select_anno_file}',header=h3)
        st.dataframe(data_file.head(5))
    
    #Blast検索の実行
    st.subheader('(4)Blast検索の実行')
    go = st.button('実行')
    if go:
        blast = Applications.NcbitblastxCommandline(query=f'{base_path}/検索配列/{fname}.fasta', db =db_path,evalue = Ev, out =f'{base_path}/Result/{fname}.txt',outfmt=0)
        stdout, stderr = blast()
        blast = Applications.NcbitblastxCommandline(query=f'{base_path}/検索配列/{fname}.fasta', db =db_path,evalue = Ev, out =f'{base_path}/Result/{fname}2.txt',outfmt=('6 qseqid qstart qend sseqid sstart send length pident bitscore evalue sseq'))
        stdout, stderr = blast()

        sleep(5)  

    #解析データの表示
        df = pd.read_table(f'{base_path}/Result/{fname}2.txt',header=None)
        df.columns = [['検索配列名','開始位置','終了位置','ヒット領域名','開始位置_','終了位置_','マッチ数','%','Bitscore','Evalue','ヒット領域配列']]
        df.to_excel(f'{base_path}/Result/{fname}2.xlsx')
        df2 = pd.read_excel(f'{base_path}/Result/{fname}2.xlsx')
        st.write('検索結果')
        if option_use == '使用':
            df3 = pd.merge(df2,data_file,left_on ='ヒット領域名',right_on =name_1,how='left' )
            st.dataframe(df3)
            df3.to_excel(buf := BytesIO(), index=True)
            st.download_button("Download",buf.getvalue(),f"{fname}.xlsx","application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
        if option_use == '未使用':
            st.dataframe(df2)
            df2.to_excel(buf := BytesIO(), index=True)
            st.download_button("Download",buf.getvalue(),f"{fname}.xlsx","application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")

        # 検索配列の除去    
        shutil.rmtree(f'{base_path}/検索配列')   
        os.makedirs(f'{base_path}/検索配列')  

def blastP():
    st.header('blastPの実行')
    st.write('NCBI：https://www.ncbi.nlm.nih.gov/')
    #sequence information     
    st.subheader('(1)検索配列情報の入力')
    fname = st.text_input(label='出力ファイル名(任意)',value = '解析') 
    seq = st.text_area('検索配列の入力（fasta形式）')
    upload_seq = st.button('アップロード') 
    if upload_seq:
        with open(f'{base_path}/検索配列/{fname}.txt',mode = 'w')as f:
            f.write(seq)
            shutil.move(f'{base_path}/検索配列/{fname}.txt',f'{base_path}/検索配列/{fname}.fasta') 
            
    #database infomation
    st.subheader('(2)データベースの選択')  
    db_name_ = os.listdir(f'{database_path}/')
    db_name = [f for f in db_name_  if not f == '.DS_Store']
    st.write('菌株を選択してください')
    select_db_name = st.selectbox('', db_name)
    db_types_ = os.listdir(f'{database_path}/{select_db_name}')
    db_types = [f for f in db_types_ if not f == '.DS_Store']
    st.write('データベースの種類を選択してください')
    select_db_type = st.radio('', db_types,horizontal=True)

    #database path
    db_path = f'{database_path}/{select_db_name}/{select_db_type}/gene'
        
    #検索方式の設定
    st.subheader('(3)パラメーター設定')
    Ev = st.text_input(label='e-value値',value='1e-2')   

    #option
    st.subheader('option_RNAseq,functionデータ付加')
    option_use = st.radio('選択',['未使用','使用'],horizontal=True)
    st.write('※データベース：gene_list、amino_listを使用時のみ！')
    if option_use == '使用':
        #アノテーションファイルの選択
        name_list = []
        for name in glob.glob(f'{database_path}/{select_db_name}/アノテーション/*.xlsx'):
            list = os.path.basename(name)
            name_list.append(list)
        select_anno_file = st.radio('アノテーションファイルの選択',name_list)
        name_1 = st.text_input('データ統合に使用する列名の設定')
        st.write('参考')
        h3 = st.number_input('ヘッダー位置調整',0,key='d')
        data_file = pd.read_excel(f'{database_path}/{select_db_name}/アノテーション/{select_anno_file}',header=h3)
        st.dataframe(data_file.head(5))  
    
    #Blast検索の実行
    st.subheader('(4)Blast検索の実行')
    go = st.button('実行')
    if go:
        blast = Applications.NcbiblastpCommandline(query=f'{base_path}/検索配列/{fname}.fasta', db =db_path,evalue = Ev, out 
                                                   =f'{base_path}/Result/{fname}.txt',outfmt=0)
        stdout, stderr = blast()
        blast = Applications.NcbiblastpCommandline(query=f'{base_path}/検索配列/{fname}.fasta', db =db_path,evalue = Ev, out =f'{base_path}/Result/{fname}2.txt',outfmt=('6 qseqid qstart qend sseqid sstart send length pident bitscore evalue sseq'))
        stdout, stderr = blast() 

        sleep(5)  

    #解析データの表示
        df = pd.read_table(f'{base_path}/Result/{fname}2.txt',header=None)
        df.columns = [['検索配列名','開始位置','終了位置','ヒット領域名','開始位置_','終了位置_','マッチ数','%','Bitscore','Evalue','ヒット領域配列']]
        df.to_excel(f'{base_path}/Result/{fname}2.xlsx')
        df2 = pd.read_excel(f'{base_path}/Result/{fname}2.xlsx')
        st.write('検索結果')
        if option_use == '使用':
            df3 = pd.merge(df2,data_file,left_on ='ヒット領域名',right_on =name_1,how='left' )
            st.dataframe(df3)
            df3.to_excel(buf := BytesIO(), index=True)
            st.download_button("Download",buf.getvalue(),f"{fname}.xlsx","application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
        if option_use == '未使用':
            st.dataframe(df2)
            df2.to_csv(buf := BytesIO(), index=True)
            st.download_button("Download",buf.getvalue(),f"{fname}.xlsx","application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")

        # 検索配列の除去    
        shutil.rmtree(f'{base_path}/検索配列')   
        os.makedirs(f'{base_path}/検索配列')  

def primer():
    st.header('primer相同領域検索')
    #sequence information     
    st.subheader('(1)検索配列情報の入力')
    fname = st.text_input(label='出力ファイル名(任意)',value = '解析') 
    seq = st.text_area('検索配列の入力（fasta形式）')
    upload_seq = st.button('アップロード') 
    if upload_seq:
        with open(f'{base_path}/検索配列/{fname}.txt',mode = 'w')as f:
            f.write(seq)
            shutil.move(f'{base_path}/検索配列/{fname}.txt',f'{base_path}/検索配列/{fname}.fasta')  
           
    #database infomation
    st.subheader('(2)データベースの選択')  
    db_name_ = os.listdir(f'{database_path}/')
    db_name = [f for f in db_name_  if not f == '.DS_Store']
    st.write('菌株を選択してください')
    select_db_name = st.selectbox('', db_name)
    db_types_ = os.listdir(f'{database_path}/{select_db_name}')
    db_types = [f for f in db_types_ if not f == '.DS_Store']
    st.write('データベースの種類を選択してください')
    select_db_type = st.radio('', db_types,horizontal=True)

    #database path
    db_path = f'{database_path}/{select_db_name}/{select_db_type}/gene'
        
    #検索方式の設定
    st.subheader('(3)パラメーター設定')
    Ev = st.text_input(label='e-value値',value='1e-2')   
       
    #Blast検索の実行
    st.subheader('(4)Blast検索の実行')
    go = st.button('実行')
    if go:
        blast = Applications.NcbiblastnCommandline(query=f'{base_path}/検索配列/{fname}.fasta', db = db_path,evalue = Ev, word_size=5,out =f'{base_path}/Result/{fname}.txt',outfmt=('6 qseqid qstart qend sseqid sstart send length pident bitscore evalue sseq'))
        stdout, stderr = blast()

        sleep(5)  

    #解析データの表示
        df = pd.read_table(f'{base_path}/Result/{fname}.txt',header=None)
        df.columns = [['検索配列名','開始位置','終了位置','ヒット領域名','開始位置_','終了位置_','マッチ数','%','Bitscore','Evalue','ヒット領域配列']]
        df.to_excel(f'{base_path}/Result/{fname}2.xlsx')
        df2 = pd.read_excel(f'{base_path}/Result/{fname}2.xlsx')
        st.write('検索結果')       
        st.dataframe(df2)
        df2.to_excel(buf := BytesIO(), index=True)
        st.download_button("Download",buf.getvalue(),f"{fname}.xlsx","application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")

        # 検索配列の除去    
        shutil.rmtree(f'{base_path}/検索配列')   
        os.makedirs(f'{base_path}/検索配列')  

if __name__ == "__main__":
    main()