import streamlit as st
import subprocess
import glob
import os
import webbrowser
import pandas as pd
import shutil
from io import BytesIO
from PIL import Image
import numpy as np
import altair as alt
from Bio.Blast import Applications
import matplotlib.pyplot as plt
import collections.abc
from pptx import Presentation
from pptx.enum.shapes import MSO_SHAPE
from pptx.util import Cm,Pt

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.abspath(os.path.join(current_dir,os.pardir))

base_path = f"{parent_dir}/RNAseq"
image = f'{base_path}/file/images.png'
trim_path = f'{base_path}/Triming'
mapping_path = f'{base_path}/Mapping'
read_count_path = f'{base_path}/Mapping'

st.title('RNAseq解析')
st.image(image,use_column_width=True)

def read_trime():
    """fastpを用いたreadデータのトリミングを実施
       「Triming」にtrim_r1.fastq.gzとtrim_r2.fastq.gzを生成
       トリミング前後のreadの情報をhtmlで表示
    """
    st.header('readデータのQC')
    st.write('「RNAseq」にreadデータ(fasatq.gz)を保存してください') 
    st.subheader('1) データ選択')
    seq_data_lists_ = glob.glob(f'{base_path}/*.gz')
    seq_data_lists = [os.path.basename(i) for i in seq_data_lists_]
    seq_data_list = st.multiselect("",seq_data_lists,max_selections=2,placeholder="＊.fastq.gzを２種類選択")
    st.subheader('2) 情報入力')
    quality = st.number_input('品質管理値：',20)
    minlen = st.number_input('保持するリードの最小値(mer)',30)
    st.subheader('3) トリミング')
    go = st.button('トリミング開始')
    if go:
        r1 = f"{base_path}/{seq_data_list[0]}"
        r2 = f"{base_path}/{seq_data_list[1]}"
        trim_path = f"{base_path}/Trimed_Read"
        os.makedirs(trim_path,exist_ok=True)
        subprocess.run(f'fastp -i {r1} -I {r2} -3 -o {trim_path}/trim_r1.fastq.gz -O {trim_path}/trim_r2.fastq.gz\
        -h {trim_path}/report.html -q {quality} -l {minlen} -g -j {trim_path}/report.json ',shell=True)
        url = (f'file:///{trim_path}/report.html')
        webbrowser.open_new_tab(url)
        os.remove(f'{trim_path}/report.json')
        #os.remove(f'{trim_path}/report.html')

def multi_mapping():
    """hisat2を使ったマッピングを実施
       samtoolsを用いてsamファイルからbamファイルを生成
    """
    st.header('マッピング')
    
    st.subheader('1) データの選択')
    st.markdown("""- RNAseqに「ゲノムデータ(*.fasta)を保管してください」""")
    st.markdown("""- RNAseq/Mapping下にRNAseqデータ毎にフォルダを作成し、「RNAseqデータ(*.fastq.gz)を保管してください」""")

    # genome_data
    genome_data_lists_ = glob.glob(f'{base_path}/*.fasta')
    genome_data_lists = [os.path.basename(i) for i in genome_data_lists_]
    genome_data = st.selectbox("ゲノムデータを選択",genome_data_lists,placeholder="＊.fastaを1種類選択")
    
    st.subheader('2) マッピング')
    go2 = st.button('マッピング')
    if go2:
        genome_path = f"{base_path}/{genome_data}"

        tmp_path = f"{base_path}/tmp"
        tmp_genome_path = f"{tmp_path}/{genome_data}"
        os.makedirs(tmp_path,exist_ok=True)
        shutil.move(genome_path,tmp_genome_path)

        read_f_path = f"{base_path}/Mapping"

        # ゲノムデータのdb作成
        subprocess.run(f'hisat2-build -p 12 {tmp_genome_path} {tmp_path}/genome_index ',shell=True)


        for f in os.listdir(read_f_path):
            read_f_path2 = f"{read_f_path}/{f}"
            if os.path.isdir(read_f_path2):
                read_lists = os.listdir(read_f_path2)
                seq_r1_path = f"{read_f_path2}/{read_lists[0]}"
                seq_r2_path = f"{read_f_path2}/{read_lists[1]}"

                # mapping
                subprocess.run(f'hisat2 -x {tmp_path}/genome_index -1 {seq_r1_path} -2 {seq_r2_path} -k 3 -p 12 -S {base_path}/result.sam',shell=True)
                subprocess.run(f'samtools sort -@ 8 -O bam -o {base_path}/{f}.bam {base_path}/result.sam' ,shell=True)
                # 確認
                check_mapping = os.path.isfile(f'{base_path}/{f}.bam')
                if check_mapping:         
                    os.remove(f'{base_path}/result.sam')
        
        st.write('マッピング完了！')
        shutil.move(tmp_genome_path,genome_path)           
        shutil.rmtree(tmp_path)
        
def p_mapping():
    """hisat2を使ったマッピングを実施
       samtoolsを用いてsamファイルからbamファイルを生成
    """
    st.header('マッピング')
    
    st.subheader('1) データの選択')
    st.markdown("""- RNAseqに「ゲノムデータ(*.fasta)とRNAseqデータ(*.fastq.gz)を保管してください」""")
    # genome_data
    genome_data_lists_ = glob.glob(f'{base_path}/*.fasta')
    genome_data_lists = [os.path.basename(i) for i in genome_data_lists_]
    genome_data = st.selectbox("ゲノムデータを選択",genome_data_lists,placeholder="＊.fastaを1種類選択")

    # seq_data
    seq_data_lists_ = glob.glob(f'{base_path}/*.gz')
    seq_data_lists = [os.path.basename(k) for k in seq_data_lists_]
    seq_data_list = st.multiselect("RNAseqデータを選択",seq_data_lists,max_selections=2,placeholder="＊.fastq.gzを２種類選択")

    st.markdown("""- 出力ファイル名(*.bam)を入力してください""")
    name = st.text_input('ファイル名入力')
    
    st.subheader('2) マッピング')
    go2 = st.button('マッピング')
    if go2:
        genome_path = f"{base_path}/{genome_data}"
        seq_r1_path = f"{base_path}/{seq_data_list[0]}"
        seq_r2_path = f"{base_path}/{seq_data_list[1]}"

        tmp_path = f"{base_path}/tmp"
        tmp_genome_path = f"{tmp_path}/{genome_data}"
        os.makedirs(tmp_path,exist_ok=True)
        shutil.move(genome_path,tmp_genome_path)

        # ゲノムデータのdb作成
        subprocess.run(f'hisat2-build -p 12 {tmp_genome_path} {tmp_path}/genome_index ',shell=True)

        # mapping
        subprocess.run(f'hisat2 -x {tmp_path}/genome_index -1 {seq_r1_path} -2 {seq_r2_path} -k 3 -p 12 -S {base_path}/result.sam',shell=True)
        subprocess.run(f'samtools sort -@ 8 -O bam -o {base_path}/{name}.bam {base_path}/result.sam' ,shell=True)
        # 確認
        check_mapping = os.path.isfile(f'{base_path}/{name}.bam')
        if check_mapping:
            st.write('マッピング完了！')
            shutil.move(tmp_genome_path,genome_path)           
            os.remove(f'{base_path}/result.sam')
            shutil.rmtree(tmp_path)
            
def s_mapping():
    """hisat2を使ったマッピングを実施
       samtoolsを用いてsamファイルからbamファイルを生成
    """
    st.header('マッピング')
    
    st.subheader('1) データの選択')
    st.markdown("""- RNAseqに「ゲノムデータ(*.fasta)とRNAseqデータ(*.fastq.gz)を保管してください」""")
    # genome_data
    genome_data_lists_ = glob.glob(f'{base_path}/*.fasta')
    genome_data_lists = [os.path.basename(i) for i in genome_data_lists_]
    genome_data = st.selectbox("ゲノムデータを選択",genome_data_lists,placeholder="＊.fastaを1種類選択")

    # seq_data
    seq_data_lists_ = glob.glob(f'{base_path}/*.gz')
    seq_data_lists = [os.path.basename(k) for k in seq_data_lists_]
    seq_data_list = st.selectbox("RNAseqデータを選択",seq_data_lists,placeholder="＊.fastq.gzを1種類選択")

    st.markdown("""- 出力ファイル名(*.bam)を入力してください""")
    name = st.text_input('ファイル名入力')
    
    st.subheader('2) マッピング')
    go2 = st.button('マッピング')

    if go2:
        genome_path = f"{base_path}/{genome_data}"
        seq_path = f"{base_path}/{seq_data_list}"

        tmp_path = f"{base_path}/tmp"
        tmp_genome_path = f"{tmp_path}/{genome_data}"
        os.makedirs(tmp_path,exist_ok=True)
        shutil.move(genome_path,tmp_genome_path)

        # ゲノムデータのdb作成
        subprocess.run(f'hisat2-build -p 12 {tmp_genome_path} {tmp_path}/genome_index ',shell=True)

        # mapping
        subprocess.run(f'hisat2 -x {tmp_path}/genome_index -U {seq_path} -k 3 -p 12 -S {base_path}/result.sam',shell=True)
        subprocess.run(f'samtools sort -@ 8 -O bam -o {base_path}/{name}.bam {base_path}/result.sam' ,shell=True)

        # 確認
        check_mapping = os.path.isfile(f'{base_path}/{name}.bam')
        if check_mapping:
            st.write('マッピング完了！')
            shutil.move(tmp_genome_path,genome_path)           
            os.remove(f'{base_path}/result.sam')
            shutil.rmtree(tmp_path)
            
def p_count():
    """マッピング結果(.bam)をもとに、FeatureCountを用いて遺伝子発現量を計算
    　　正規化(RPK,TPM)を実施
    　　これらの結果を、"サンプル名".xlsxとして出力
    """
    st.header('遺伝子発現量のカウント')
    st.subheader('1) 準備')
    st.markdown("""- RNAseqにMappingデータ（*.bam)とアノテーションデータ(*.gtf)を保管してください」""")
    # mapping_data
    mapping_data_lists_ = glob.glob(f'{base_path}/*.bam')
    mapping_data_lists = [os.path.basename(i) for i in mapping_data_lists_]
    mapping_data = st.selectbox("mappingデータを選択",mapping_data_lists,placeholder="＊.bamを1種類選択")

    annotation_data_lists_ = glob.glob(f'{base_path}/*.gtf')
    annotation_data_lists = [os.path.basename(ano) for ano in annotation_data_lists_]
    annotation_data = st.selectbox("アノテーションデータを選択",annotation_data_lists,placeholder="＊.gtfを1種類選択")

    st.markdown("""- gtfファイル情報の設定　＊基本はデフォルトでOK""")
    count_type = st.selectbox('カウント対象',['transcript','exon','CDS'])
    gene_id = st.selectbox('遺伝子番号の表記法',['transcript_id','gene_id'])

    st.markdown("""- 出力ファイル名の設定""")
    count_name = st.text_input('ファイル名を入力')

    st.subheader('２) カウント')
    go = st.button('開始')
    if go:
        #path
        mapping_data_path = f"{base_path}/{mapping_data}" 
        gtf_path = f"{base_path}/{annotation_data}"

        #featureCount
        subprocess.run(f'samtools index {mapping_data_path}',shell=True)
        subprocess.run(f'featureCounts -p -t {count_type} -g {gene_id} -a {gtf_path} -o {base_path}/result.txt {mapping_data_path}',shell=True)
        #TPM正規化
        df = pd.read_table(f'{base_path}/result.txt',sep='\t',skiprows=1)
        df['RPK正規化'] = 1000/df.iloc[:,-2]*df.iloc[:,-1]
        df['TPM正規化'] = 1000000/df.iloc[:,-1].sum()*df.iloc[:,-1]
        df.columns = [['Gene_id','Chr','Start','End','Strand','Length',f'{count_name}_Count_data',f'{count_name}_RPK',f'{count_name}_TPM']]
        df.to_excel(f'{base_path}/{count_name}.xlsx')
        st.dataframe(df.iloc[:,[0,-3,-1]])
        st.write(f'結果ファイル名：{count_name}.xlsx）')
        #不要ファイルの除去
        os.remove(f'{base_path}/result.txt')
        os.remove(f'{base_path}/result.txt.summary')
        
def merge_data ():
    """遺伝子番号をキーとして、遺伝子発現領解析データやアノテーションデータの統合を行う
       遺伝子番号の列名は、アップロードしたxlsxデータを基にユーザーが選択
    """
    st.header('データ統合')
    count_data = st.file_uploader('カウントデータ',type='xlsx')
    merge_file = st.file_uploader('追加データ',type='xlsx',key='b')

    m_name = st.text_input('ファイル名を入力')
    
    st.write('カウントデータ')
    h1 = st.number_input('ヘッダー位置調整',0,key='a')
    df = pd.read_excel(count_data,header=h1)
    st.dataframe(df.head(5))
    
    st.write('追加データ')
    h2 = st.number_input('ヘッダー位置調整',0)
    df2 = pd.read_excel(merge_file,header=h2)
    st.dataframe(df2.head(5))
    
    st.subheader('遺伝子番号の列名の指定')
    col1,col2 = st.columns((3,10))
    with col1:
        row_num1 = st.number_input('',1,key = 'count1')
        row_name1 = df.columns.values[row_num1]
    with col2:
        name_1 = st.text_input('カウントデータ',row_name1)
        
    col2,col3 = st.columns((3,10))
    with col2:
        row_num2 = st.number_input('',1,key = 'count2')
        row_name2 = df2.columns.values[row_num2]
    with col3:
        name_2 = st.text_input('追加データ',row_name2)
        
    merge_start = st.button('データ統合')
    if merge_start:
        df3 = pd.merge(df,df2,left_on =name_1,right_on=name_2,how ='left' )
        st.write('統合データ')
        st.dataframe(df3.head(5))
        df3.to_excel(buf := BytesIO(), index=True)
        st.download_button("Download",buf.getvalue(),f"{m_name}.xlsx","application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")

def calculate_MA2 ():
    """２サンプル間の遺伝子発現量の変動解析（MAplot）を実施"""
    st.write('データアップロード')
    st.write("内部で正規化処理する為、rawデータを入力すること")
    count_data1 = st.file_uploader('実験1',type='xlsx',key='a1')
    count_data2 = st.file_uploader('実験2',type='xlsx',key='b1')
    func_data = st.file_uploader('遺伝子機能リスト',type='xlsx',key='c1')
    
    st.write('データ調整 (実験1)')
    h1 = st.number_input('ヘッダー位置調整',0,key='data1')
    ex1 = pd.read_excel(count_data1,header=h1)
    st.dataframe(ex1.head(5))    

    st.write('データ調整 (実験2)')
    h2 = st.number_input('ヘッダー位置調整',0,key='data2')
    ex2 = pd.read_excel(count_data2,header=h2)
    st.dataframe(ex2.head(5))
    
    st.write('データ調整 (遺伝子機能)')
    h3 = st.number_input('ヘッダー位置調整',0,key='z')
    ex3 = pd.read_excel(func_data,header=h3)
    st.dataframe(ex3.head(5))
    
    st.write('基本情報の入力')
    
    st.write('＜実験1＞')
    col1,col2,col3,col4,col5 = st.columns((5,5,5,5,5))
    with col1:
        ex_name1 = st.text_input('実験名(任意)','実験1',key = 'a')
    with col2:
        row_num1 = st.number_input('遺伝子番号',0,key = 'b')
        row_name1 = ex1.columns.values[row_num1]
    with col3:
        key1 = st.text_input('',row_name1,key = 'gene_id')
    with col4:
        row_num2 = st.number_input('カウント',2,key = 'c')
        row_name2 = ex1.columns.values[row_num2]
    with col5:
        count_name_1 = st.text_input('',row_name2,key = 'count1')
    
    st.write('＜実験2＞')
    col6,col7,col8,col9,col10 = st.columns((5,5,5,5,5))
    with col6:
        ex_name2 = st.text_input('実験名(任意)','実験2',key = 'd')
    with col7:
        ex2_row_num1 = st.number_input('遺伝子番号',0,key = 'gene_num')
        ex2_row_name = ex2.columns.values[ex2_row_num1]
    with col8:
        key2 = st.text_input('',ex2_row_name,key = 'gene_id2')
    with col9:
        ex2_row_num2 = st.number_input('カウント',2,key = 'count_na')
        ex2_row_name2 = ex2.columns.values[ex2_row_num2]
    with col10:
        count_name_2 = st.text_input('',ex2_row_name2,key = 'count2')
    
    st.write('＜遺伝子機能＞')
    col1,col2,col3,col4 = st.columns((5,5,5,5))
    with col1:
        row_num1 = st.number_input('',0,key = 'func1')
        row_name1 = ex3.columns.values[row_num1]
    with col2:
        key3 = st.text_input('遺伝子番号',row_name1,key='func2')
    with col3:
        row_num2 = st.number_input('',1,key = 'func3')
        row_name2 = ex3.columns.values[row_num2]
    with col4:
        func_name = st.text_input('機能情報',row_name2,key = 'func4')
    
    merge_start = st.button('実行')
    if merge_start:
        # 遺伝子番号とカウントデータ抽出、カウントデータ名をtmp1,tmp2に変更
        ex1_data = ex1[[key1,count_name_1]]
        ex2_data = ex2[[key2,count_name_2]]
        ex1_data2 = ex1_data.rename(columns={count_name_1:'tmp1'})
        ex2_data2 = ex2_data.rename(columns={count_name_2:'tmp2'})
        # 実験1、実験２の統合
        merge = pd.merge(ex1_data2,ex2_data2,left_on = key1, right_on =key2 ,how ='outer')
        # 正規化処理（実験1、実験2それぞれで、カウントデータの総数が100万リードとなるように調整）
        merge[ex_name1] = merge['tmp1']*1000000 / merge['tmp1'].sum()
        merge[ex_name2] = merge['tmp2']*1000000 / merge['tmp2'].sum()
        # カウント数が０のものを0.01に置き換え（M、A計算時のエラー回避のため)
        merge[ex_name1] = merge[ex_name1].replace(0,0.01)
        merge[ex_name2] = merge[ex_name2].replace(0,0.01)
        # M、A計算
        merge['M'] = np.log2(merge[ex_name2]) - np.log2(merge[ex_name1])
        merge['A'] = (np.log2(merge[ex_name2])+np.log2(merge[ex_name1]))/2
        # 不要列の削除
        merge2 = merge.drop(['tmp1','tmp2'],axis=1).dropna()
        # 「変化パターン」列を生成し、閾値プラス１をup、マイナス1をdownと表示
        threshold1 = 1
        threshold2 = -1
        merge2['変化パターン']='mentain'
        for i in range(0,len(merge2)):
            if threshold1 < merge2.iloc[i,-3]:
                merge2.iloc[i,-1] = 'up'
            elif merge2.iloc[i,-3] < threshold2:
                merge2.iloc[i,-1] = 'down'
                
        # 遺伝子機能データ(func_data)の調整
        func_data2 = ex3[[key3,func_name]]
        func_data3 = func_data2.rename(columns={func_name:'function'})
        # 遺伝子機能データの追加
        merge3 = pd.merge(merge2,func_data3,left_on = key1,right_on = key3,how = 'left')
        
    return merge3

        
def calculate_MA ():
    st.write('データアップロード')
    count_data1 = st.file_uploader('実験1',type='xlsx',key='a1')
    count_data2 = st.file_uploader('実験2',type='xlsx',key='b1')
    
    st.write('データ調整 (実験1)')
    h1 = st.number_input('ヘッダー位置調整',0,key='data1')
    ex1 = pd.read_excel(count_data1,header=h1)
    st.dataframe(ex1.head(5))    

    st.write('データ調整 (実験2)')
    h2 = st.number_input('ヘッダー位置調整',0,key='data2')
    ex2 = pd.read_excel(count_data2,header=h2)
    st.dataframe(ex2.head(5))

    st.write('基本情報の入力')
    st.write('＜実験1＞')
    col1,col2,col3,col4,col5 = st.columns((5,5,5,5,5))
    with col1:
        ex_name1 = st.text_input('実験名(任意)','実験1',key = 'a')
    with col2:
        row_num1 = st.number_input('遺伝子番号',0,key = 'b')
        row_name1 = ex1.columns.values[row_num1]
    with col3:
        key1 = st.text_input('',row_name1,key = 'gene_id')
    with col4:
        row_num2 = st.number_input('カウント',2,key = 'c')
        row_name2 = ex1.columns.values[row_num2]
    with col5:
        count_name_1 = st.text_input('',row_name2,key = 'count1')
    
    st.write('＜実験2＞')
    col6,col7,col8,col9,col10 = st.columns((5,5,5,5,5))
    with col6:
        ex_name2 = st.text_input('実験名(任意)','実験2',key = 'd') 
    with col7:
        ex2_row_num1 = st.number_input('遺伝子番号',0,key = 'gene_num')
        ex2_row_name = ex2.columns.values[ex2_row_num1]
    with col8:
        key2 = st.text_input('',ex2_row_name,key = 'gene_id2')
    with col9:
        ex2_row_num2 = st.number_input('カウント',2,key = 'count_na')
        ex2_row_name2 = ex2.columns.values[ex2_row_num2]
    with col10:
        count_name_2 = st.text_input('',ex2_row_name2,key = 'count2')
    
    merge_start = st.button('実行')
    if merge_start:
        # 遺伝子番号とカウントデータ抽出、カウントデータ名をtmp1,tmp2に変更
        ex1_data = ex1[[key1,count_name_1]]
        ex2_data = ex2[[key2,count_name_2]]
        ex1_data2 = ex1_data.rename(columns={count_name_1:'tmp1'})
        ex2_data2 = ex2_data.rename(columns={count_name_2:'tmp2'})
        # 実験1、実験２の統合
        merge = pd.merge(ex1_data2,ex2_data2,left_on = key1, right_on =key2 ,how ='outer')
        # 正規化処理（実験1、実験2それぞれで、カウントデータの総数が100万リードとなるように調整）
        merge[ex_name1] = merge['tmp1']*1000000 / merge['tmp1'].sum()
        merge[ex_name2] = merge['tmp2']*1000000 / merge['tmp2'].sum()
        # カウント数が０のものを0.01に置き換え（M、A計算時のエラー回避のため)
        merge[ex_name1] = merge[ex_name1].replace(0,0.01)
        merge[ex_name2] = merge[ex_name2].replace(0,0.01)
        # M、A計算
        merge['M'] = np.log2(merge[ex_name2]) - np.log2(merge[ex_name1])
        merge['A'] = (np.log2(merge[ex_name2])+np.log2(merge[ex_name1]))/2
        # 不要列の削除
        merge2 = merge.drop(['tmp1','tmp2'],axis=1).dropna()
        # 「変化パターン」列を生成し、閾値プラス１をup、マイナス1をdownと表示
        threshold1 = 1
        threshold2 = -1
        merge2['変化パターン']='mentain'
        for i in range(0,len(merge2)):
            if threshold1 < merge2.iloc[i,-3]:
                merge2.iloc[i,-1] = 'up'
            elif merge2.iloc[i,-3] < threshold2:
                merge2.iloc[i,-1] = 'down'
         # 機能情報、遺伝子タイプ列を生成し、Noneを入力
        merge2['function'] = 'None'

    return merge2
        
    
def ready_to_MAdata():
    st.markdown("""- MAデータの保存フォルダを新規作成""")
    fgroup = st.text_input('','フォルダ名',key = 'folder')
    set_folder = st.button('新規フォルダ作成')

    if set_folder:
        os.makedirs(f'{base_path}/MAplot/{fgroup}',exist_ok=True) 

    st.markdown("""- MAデータ保存フォルダの選択およびデータ保存名の入力""")

    col1, col2,col3 = st.columns((3,5,5))
    with col1:
        list_num = st.number_input('',0,key = 'list')
        lists = os.listdir(f'{base_path}/MAplot/')
    with col2:
        fgroup = st.text_input('テーマ名',lists[list_num],key = 'theme')
    with col3:
        fname = st.text_input('保存データ名','data',key = 'i')
        
    select_mode = st.radio('データ作成モードの選択',['MA計算','MA計算（＋遺伝子機能情報）'])
    if select_mode == 'MA計算':
        df = calculate_MA ()
        # 結果の表示
        st.dataframe(df.head(5))
        # ファイルの保存
        os.makedirs(f'{base_path}/MAplot/{fgroup}',exist_ok = True)
        df.to_excel(f'{base_path}/MAplot/{fgroup}/{fname}.xlsx') 
        # ダウンロード
        df.to_excel(buf := BytesIO(), index=True)
        st.download_button("Download",buf.getvalue(),f"{fname}_MAplotデータ.xlsx","application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
        
    if select_mode == 'MA計算（＋遺伝子機能情報）':
        df2 = calculate_MA2 ()
        
        # 結果の表示
        st.dataframe(df2.head(5))
        # ファイルの保存
        os.makedirs(f'{base_path}/MAplot/{fgroup}',exist_ok = True)
        df2.to_excel(f'{base_path}/MAplot/{fgroup}/{fname}.xlsx') 
        # ダウンロード
        df2.to_excel(buf := BytesIO(), index=True)
        st.download_button("Download",buf.getvalue(),f"{fname}_MAplotデータ.xlsx","application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
        
def dgeAnalysis():
    st.header('遺伝子変動解析（MAプロット）')
    ready_data = st.radio('ステップを選択',['step1_解析データ準備','step2_MAprot'])
    
    if ready_data == 'step1_解析データ準備':
        ready_to_MAdata ()
        
    if ready_data == 'step2_MAprot':
        st.subheader('MAプロット作成')
        #「MAplot」からMA計算データを読込み
        st.write('解析データの読込み')
        add = ['-選択-']
        group_list_p = os.listdir(f'{base_path}/MAplot/')
        group_list = add+group_list_p
        select_group = st.selectbox('', group_list,key='sg')
        MAdata_list_p = os.listdir(f'{base_path}/MAplot/{select_group}/')
        MAdata_list = add + MAdata_list_p
        select_MAdata = st.selectbox('', MAdata_list)
        MAdata = pd.read_excel(f'{base_path}/MAplot/{select_group}/{select_MAdata}')
        st.dataframe(MAdata)
        
        # 表示データの指定
        # 遺伝子番号
        id_name = MAdata.columns.values[1]
        value1_name = MAdata.columns.values[2]
        value2_name = MAdata.columns.values[3]
        func = 'function'
        type = 'type'
        
        # title
        st.subheader('MAplot')
        st.write(f'{value2_name}　→　{value1_name}')
        
        # MAプロット作成1
        selector = alt.selection_multi(fields=['変化パターン'],bind='legend')
        color = alt.condition(selector, alt.Color('変化パターン:N',title='変化パターン'),alt.value('lightgray'))
        base = alt.Chart(MAdata).properties(width=300, height=300)
        # プロットデータ
        prot = base.mark_point(filled=True).encode(x = alt.X('A',title='average'),y=alt.Y('M',title='log2 fold change'),color=color,tooltip=[alt.Tooltip(id_name,title='Gene_id'),alt.Tooltip(value1_name,title = f'Count:{value1_name}'),alt.Tooltip(value2_name,title = f'Count:{value2_name}'),alt.Tooltip(func,title = 'function')]).interactive().add_selection(selector)
        
        st.altair_chart(prot,use_container_width =True)
        
        # MAプロット作成2
        brush = alt.selection_interval()
        color = alt.condition(brush, alt.Color('変化パターン:N',title='変化パターン'),alt.value('lightgray'))
        base = alt.Chart(MAdata).properties(width=300, height=300)
        # プロットデータ
        prot = alt.Chart(MAdata).mark_point(filled=True).encode(x = alt.X('A',title='average'),y=alt.Y('M',title='log2 fold change'),color=color,tooltip=[alt.Tooltip(id_name,title='Gene_id'),alt.Tooltip(value1_name,title = f'Count:{value1_name}'),alt.Tooltip(value2_name,title = f'Count:{value2_name}'),alt.Tooltip(func,title = 'function')]).add_selection(brush)
        # 棒グラフ
        bar = base.mark_bar().encode(x = alt.X(id_name,title=f'{value1_name}'),y = alt.Y(value1_name),color=color,tooltip=[alt.Tooltip(id_name,title='Gene_id'),alt.Tooltip(func,title = 'function')])
        bar2 = base.mark_bar().encode(x = alt.X(id_name,title=f'{value2_name}'),y = alt.Y(value2_name),color=color,tooltip=[alt.Tooltip(id_name,title='Gene_id'),alt.Tooltip(func,title = 'function')])
        bars = alt.hconcat(bar,bar2)
        # グラフ統合
        merge_graph = alt.vconcat(prot,bars.transform_filter(brush))
        # 表示
        st.altair_chart(merge_graph,use_container_width =True)
        
        st.write('option_特定の遺伝子セットをハイライト')
        #templateのダウンローダー表示
        with st.expander("遺伝子リスト入力用シートのダウンロード", expanded=False):
           template = pd.read_excel(f'{base_path}/file/template.xlsx')
           template.to_excel(buf := BytesIO(), index=True)
           st.download_button("Download",buf.getvalue(),"template.xlsx","application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
        # gene_listのアップロード
        data = st.file_uploader('遺伝子リスト',type='xlsx',key='gene_list')
        gene_list = pd.read_excel(data)
        # keyの設定
        st.write('遺伝子番号列の設定')
        col1,col2,col3,col4 = st.columns((4,5,4,5))
        with col1:
           row_num1 = st.number_input('MAdata',0,key='a')
           row_name1 = MAdata.columns.values[row_num1]
        with col2:
           key1 = st.text_input('',row_name1,key = 'gene1')
        with col3:
           row_num2 = st.number_input('Gene_list',0,key='b')
           row_name2 = gene_list.columns.values[row_num2]
        with col4:
           key2 = st.text_input('',row_name2,key='gene2')
        
        MAdata2 = pd.merge(MAdata,gene_list,left_on =key1,right_on=key2,how = 'left')
        
        # MAプロット作成3
        selector2 = alt.selection_multi(fields=['category'],bind='legend')
        color2 = alt.condition(selector2, alt.Color('category:N',title='遺伝子タイプ'),alt.value('lightgray'))
        base2 = alt.Chart(MAdata2).properties(width=300, height=300)
        # プロットデータ
        prot2 = base2.mark_point(filled=True).encode(x = alt.X('A',title='average'),y=alt.Y('M',title='log2 fold change'),color=color2,tooltip=[alt.Tooltip(id_name,title='Gene_id'),alt.Tooltip(value1_name,title = f'Count:{value1_name}'),alt.Tooltip(value2_name,title = f'Count:{value2_name}'),alt.Tooltip(func,title = 'function')]).interactive().add_selection(selector2)

        st.altair_chart(prot2,use_container_width =True)

def change_gtf_file():
    """gffreadを用いて、GFFファイルからGTFファイルを生成"""
    genome_data_lists_ = glob.glob(f'{base_path}/*.fasta')
    genome_data_lists = [os.path.basename(i) for i in genome_data_lists_]
    genome_data = st.selectbox("ゲノムデータを選択",genome_data_lists,placeholder="＊.fastaを1種類選択")

    annotation_data_lists_ = glob.glob(f'{base_path}/*.gff')
    annotation_data_lists = [os.path.splitext(os.path.basename(ano))[0] for ano in annotation_data_lists_]
    annotation_data = st.selectbox("アノテーションデータを選択",annotation_data_lists,placeholder="＊.gffを1種類選択")

    if st.button('GTFの変換開始'):
        genome_path = f'{base_path}/{genome_data}'
        gff3_path = f'{base_path}/{annotation_data}.gff'
        subprocess.run(f'gffread {gff3_path} -g {genome_path} -E -T -o {annotation_data}.gtf' ,shell=True)

        check_data = os.path.isfile(f'{base_path}/{annotation_data}.gtf')
        if check_data:
            st.write(f'gtf生成完了！_{annotation_data}.gtf')

st.sidebar.subheader('STEP1')
select1 = st.sidebar.selectbox('readデータのQC',['-選択-','トリミング'])
st.sidebar.subheader('STEP2')
select2 = st.sidebar.selectbox('ゲノムDNAへのマッピング',['-選択-','single-end','pair-end','multi_pair_end'],key='1')
st.sidebar.subheader('STEP3')
select3 = st.sidebar.selectbox('遺伝子発現量のカウント',['-選択-','遺伝子発現量カウント'],key='2')
st.sidebar.subheader('Option1')
select4 = st.sidebar.selectbox('データ統合',['-選択-','データ統合'])
st.sidebar.subheader('Option2')
select5 = st.sidebar.selectbox('遺伝子変動解析',['-選択-','遺伝子変動解析(MAplot)'])
st.sidebar.subheader('Option3')
select6 = st.sidebar.selectbox('GFF->GTF変換',['-選択-','GTF変換ツール'])

if select1 == 'トリミング':
    read_trime()

if select2 == 'pair-end':
    p_mapping()

if select2 == 'single-end':
    s_mapping()

if select2 == 'multi_pair_end':
    multi_mapping()

if select3 == '遺伝子発現量カウント':
    p_count()

if select4 == 'データ統合':
    merge_data()

if select5 == '遺伝子変動解析(MAplot)':
    dgeAnalysis()

if select6 == 'GTF変換ツール':
    change_gtf_file()