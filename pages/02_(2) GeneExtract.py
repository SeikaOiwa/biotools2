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
base_path2 = f'{base_path_}/GeneExtract'

#表示画像の読込み
image = Image.open(f'{base_path2}/File/front2.png')
image2 = Image.open(f'{base_path2}/File/end2.png')

def main(): 
    st.image(image,use_column_width=None)   
    st.sidebar.subheader('ゲノムデータ登録・確認')
    select = st.sidebar.selectbox('選択',['-未選択-','登録情報確認','ゲノム情報登録'])
    
    st.sidebar.subheader('GeneExtract')
    select2 = st.sidebar.selectbox('抽出方法：', ['-未選択-','gene_idから抽出','位置情報から抽出','プロモーター・末端領域抽出'])         
    if select == '登録情報確認':
        genome_check()
    if select == 'ゲノム情報登録':
        genome_regist()
        
    if select2 == '位置情報から抽出':
        extract_locas()
        
    if select2 == 'gene_idから抽出' :
        extract_id ()
    
    if select2 == 'プロモーター・末端領域抽出':
        promorter_extract()     
        
# ゲノム情報の登録        
def genome_regist():
    st.header('ゲノム情報の登録')
    st.subheader('データ保管 ※データ登録済みの場合は不要') 
    st.subheader('ゲノム名を入力してください')
    strain_name = st.text_input(label='-入力-')
    make_folder = st.button(label = 'フォルダ作成')
    if make_folder:
       os.makedirs(f'{base_path2}/{strain_name}/list')
       st.write(f'「GeneExtract」/{strain_name}にゲノム情報(.fasta)、アノテーションデータ(.gtf)を保管')    
    set_genome_data = st.button('データセット')
    st.write('ゲノム情報のみでも可')
    create_g_list = st.button('遺伝子リスト構築')
    st.write('ゲノム情報＋アノテーション情報必要')
    #ゲノム基本情報を出力（ゲノム情報のみで実施）
    if set_genome_data:
        #統計データ１の作成
      subprocess.run(f'seqkit stats {base_path2}/{strain_name}/*.fasta -T > {base_path2}/{strain_name}/stat.csv',shell=True) 
        #統計データ2の作成
      subprocess.run(f'seqkit watch --fields ReadLen {base_path2}/{strain_name}/*.fasta -O {base_path2}/{strain_name}/len.png',shell=True)
      # GC含量
      subprocess.run(f'seqkit fx2tab {base_path2}/{strain_name}/*.fasta | head -n 1 | seqkit tab2fx | seqkit sliding -s 1 -W 100 | seqkit fx2tab -n -g > {base_path2}/{strain_name}/GC.csv',shell=True) 
      # コンティグリスト作成
      subprocess.run(f'seqkit seq -ni {base_path2}/{strain_name}/*.fasta > {base_path2}/{strain_name}/id_list.csv',shell=True)
      
    #遺伝子、アミノ酸リストの出力（ゲノム情報+アノテーションデータで実施）
    if create_g_list:
       #遺伝子リスト作成
       subprocess.run(f'gffread --gtf {base_path2}/{strain_name}/*.gtf -w {base_path2}/{strain_name}/list/gene_list.fasta -g {base_path2}/{strain_name}/*.fasta',shell=True)
       #gene_idリスト作成
       subprocess.run(f'seqkit seq -ni {base_path2}/{strain_name}/list/gene_list.fasta > {base_path2}/{strain_name}/list/gene_id_list.csv',shell=True)
       #アミノ酸リスト作成
       subprocess.run(f'gffread --gtf {base_path2}/{strain_name}/*.gtf -y {base_path2}/{strain_name}/list/protein_list.fasta -g {base_path2}/{strain_name}/*.fasta',shell=True)
         
# ゲノム情報の可視化 
def genome_check():
    st.header('ゲノム情報の可視化')
    st.subheader('ゲノム情報の選択')
    genome_name__ = os.listdir(f'{base_path2}/')
    add_list = ['-選択-']
    genome_name_ = add_list + genome_name__
    genome_name = [i for i in genome_name_ if not i == 'File']
    select_genome_name = st.selectbox('', genome_name)
    
    st.subheader('データ可視化')
    visual_genome = st.button('可視化')
    if visual_genome:
        #データ取込み
        df = pd.read_csv(f'{base_path2}/{select_genome_name}/stat.csv',sep = '\t')
        df2 = pd.read_csv(f'{base_path2}/{select_genome_name}/GC.csv',sep = '\t',header=None)
        df2.columns = [['id','GC含量']]
        gc_ave = df2['GC含量'].mean()
        for i in glob.glob(f'{base_path2}/{select_genome_name}/*.gtf'):       
          df3 = pd.read_table(i,header=None)
          df4 = df3.head(100)
 
        #表示
        st.subheader('統計量')
        stat_image = Image.open(f'{base_path2}/{select_genome_name}/len.png')
        st.write(df)
        st.image(stat_image,use_column_width=True) 
        st.subheader('GC含量')
        st.write(gc_ave)
        st.subheader('アノテーションデータ')
        st.write(df4)
               
#位置情報から遺伝子抽出
def extract_locas():
    st.header('ゲノム上の位置情報から遺伝子抽出')  
    st.subheader('(1)ゲノムの選択')
    genome_list_ = os.listdir(f'{base_path2}/') 
    genome_list = [i for i in genome_list_ if not i == 'README.md' and not i == 'File' and not i == '.DS_Store']  
    select_genome = st.radio('_',genome_list)  
    st.subheader('(2)位置情報の入力')
    chrom_list = pd.read_csv(f'{base_path2}/{select_genome}/id_list.csv',sep = '\t',header=None)
    chrom_list.columns = ['id'] 
    chrom = st.selectbox('コンティグ名(染色体)',chrom_list)
    start = st.text_input(label = '開始位置')
    end = st.text_input(label = '終了位置')
    
    st.subheader('(3)遺伝子領域の抽出')
    go = st.button(label='開始')
    if go:
        subprocess.run(f'seqkit subseq --chr {chrom} -r {start}:{end} {base_path2}/{select_genome}/*.fasta >{base_path2}/extract.txt', shell=True)
      
        #抽出データファイルの読み込み
        with open(f'{base_path2}/extract.txt','r') as f:
            seq = f.read()
        os.remove(f'{base_path2}/extract.txt')
        #抽出データファイルの表示            
        st.write(seq)
        st.download_button('ダウンロード',seq, file_name=f'{chrom}_start_{start}_end_{end}.txt')
        st.write('InterProScan：https://www.ebi.ac.uk/interpro/search/sequence/')
        st.image(image2,use_column_width=True)


#gene_idから抽出
def extract_id():
    st.header('Gene_idから遺伝子抽出')
    st.subheader('(1)ゲノムの選択')
    genome_list_ = os.listdir(f'{base_path2}/') 
    genome_list = [i for i in genome_list_ if not i == 'README.md' and not i == 'File' and not i == '.DS_Store']   
    select_genome = st.radio('_',genome_list)    
    st.subheader('(2)抽出領域の入力')      
    gene_list_p = pd.read_csv(f'{base_path2}/{select_genome}/list/gene_id_list.csv',sep = '\t',header=None)
    gene_list_p.columns = ['id'] 
    gene_list = gene_list_p.sort_values('id')

    gene_id = st.text_input('gene_id')
    with st.expander("gene_idリスト", expanded=False):
       gene_list2 = st.selectbox('',gene_list)
    upstream = st.number_input('上流',0)
    downstream = st.number_input('下流',0)    
    
    st.subheader('(3)遺伝子領域の抽出')
    go = st.button(label='塩基配列を抽出（イントロン含む）')
    go1 = st.button(label='塩基配列を抽出（イントロン含まない）')
    go2 = st.button(label='アミノ酸配列を抽出')
    
    if go:
     for gtf_file in glob.glob(f'{base_path2}/{select_genome}/*.gtf'):
        df = pd.read_table(gtf_file,header=None)
        df2 = df[df.iloc[:,-1].str.contains(f'"{gene_id}"')]
     chrom = df2.iloc[0,0]
     start = df2.iloc[0,3]-upstream
     end = df2.iloc[0,4]+downstream
     subprocess.run(f'seqkit subseq --chr {chrom} -r {start}:{end} {base_path2}/{select_genome}/*.fasta > {base_path2}/extract.fasta', shell=True)
     
     #相補鎖の取得
     subprocess.run(f'seqkit seq -pr {base_path2}/extract.fasta > {base_path2}/extract2.fasta', shell=True)  
     
     #配列１
     with open(f'{base_path2}/extract.fasta','r') as f:
        seq1 = f.read()
        os.remove(f'{base_path2}/extract.fasta')
        
     #配列2
     with open(f'{base_path2}/extract2.fasta','r') as f2:
        seq2 = f2.read()
        os.remove(f'{base_path2}/extract2.fasta')
    
    #抽出データファイルの表示
     st.write('配列(5→3)')
     st.write(seq1)
     st.download_button('ダウンロード',seq1, file_name=f'{gene_id}_up_{start}_down_{end}_with_intron.txt')   
     st.write('InterProScan：https://www.ebi.ac.uk/interpro/search/sequence/')     
     st.write('配列(3→5)')
     st.write(seq2)
     st.download_button('ダウンロード',seq2, file_name=f'{gene_id}_up_{start}_down_{end}_with_intron.txt') 
     st.write('InterProScan：https://www.ebi.ac.uk/interpro/search/sequence/')
     st.image(image2,use_column_width=True) 
     
    if go1:
     subprocess.run(f'seqkit grep -nrp {gene_id} {base_path2}/{select_genome}/list/gene_list.fasta > {base_path2}/extract.txt',shell=True)
     with open(f'{base_path2}/extract.txt','r') as f:
        seq = f.read()
        os.remove(f'{base_path2}/extract.txt')
                 
     #抽出データファイルの表示            
     st.write(seq)
     st.download_button('ダウンロード',seq, file_name=f'{gene_id}_without_intron.txt') 
     st.write('InterProScan：https://www.ebi.ac.uk/interpro/search/sequence/')
     st.image(image2,use_column_width=True)     
    
    if go2:
     subprocess.run(f'seqkit grep -nrp {gene_id} {base_path2}/{select_genome}/list/protein_list.fasta > {base_path2}/extract.txt',shell=True)
     with open(f'{base_path2}/extract.txt','r') as f:
        seq = f.read()
        os.remove(f'{base_path2}/extract.txt')
                 
     #抽出データファイルの表示            
     st.write(seq)
     st.download_button('ダウンロード',seq, file_name=f'{gene_id}_without_intron.txt')  
     st.write('InterProScan：https://www.ebi.ac.uk/interpro/search/sequence/')
     st.image(image2,use_column_width=True)  
     
def promorter_extract():
    st.header('プロモーター・末端領域の抽出')
    st.subheader('(1)ゲノムの選択')
    genome_list_ = os.listdir(f'{base_path2}/') 
    genome_list = [i for i in genome_list_ if not i == 'README.md' and not i == 'File' and not i == '.DS_Store']    
    select_genome = st.radio('_',genome_list)
    add_gene_num = st.number_input('抽出範囲',5)    
    
    st.subheader('(2)-1_特定遺伝子領域を対象')
    gene_list_p = pd.read_csv(f'{base_path2}/{select_genome}/list/gene_id_list.csv',sep = '\t',header=None)
    gene_list_p.columns = ['id'] 
    gene_list = gene_list_p.sort_values('id')
    gene_id = st.text_input('gene_id')
    with st.expander("gene_idリスト", expanded=False):
       gene_list2 = st.selectbox('',gene_list)

    go = st.button(label='開始')
    if go:
       #tmp.fastaに抽出範囲(add_gene_num)を含む全遺伝子を格納
       subprocess.run(f'gffread --gtf {base_path2}/{select_genome}/*.gtf -w {base_path2}/tmp.fasta -g {base_path2}/{select_genome}/*.fasta --w-add {add_gene_num}',shell=True)
       #tmp.fastaから特定遺伝子(gene_id)を抽出し、gene_id_add.fastaに格納
       subprocess.run(f'seqkit grep -nrp {gene_id} {base_path2}/tmp.fasta >  {base_path2}/gene_id_add.fasta',shell=True)
       #gene_id_add.fastaから5'領域のみを抽出し、5end.txtに格納
       subprocess.run(f'seqkit subseq -r 1:{add_gene_num} {base_path2}/gene_id_add.fasta > {base_path2}/5end.txt',shell=True)
       #gene_id_add.fastaから3'領域のみを抽出し、3end.txtに格納
       subprocess.run(f'seqkit subseq -r -{add_gene_num}:-1 {base_path2}/gene_id_add.fasta > {base_path2}/3end.txt',shell=True)
       #抽出結果の読込み
       with open(f'{base_path2}/5end.txt','r') as f:
          promoter_seq = f.read()
       with open(f'{base_path2}/3end.txt','r') as f2:
          terminater_seq = f2.read()
          os.remove(f'{base_path2}/5end.txt')       
          os.remove(f'{base_path2}/3end.txt')     
          os.remove(f'{base_path2}/tmp.fasta') 
          os.remove(f'{base_path2}/gene_id_add.fasta') 
          os.remove(f'{base_path2}/gene_id_add.fasta.seqkit.fai')
     
     #抽出結果の表示
       st.subheader('プロモーター領域')
       st.write(promoter_seq)
       st.download_button('ダウンロード',promoter_seq, file_name=f'{gene_id}_5end.txt') 
       st.subheader('ターミネーター領域')
       st.write(terminater_seq)
       st.download_button('ダウンロード',terminater_seq, file_name=f'{gene_id}_3end.txt') 
       st.image(image2,use_column_width=True)
       
    st.subheader('(2)-2_全遺伝子領域を対象')
    go2 = st.button(label='開始_')
    if go2:    
       #tmp.fastaに抽出範囲(add_gene_num)を含む全遺伝子を格納
       subprocess.run(f'gffread --gtf {base_path2}/{select_genome}/*.gtf -w {base_path2}/tmp.fasta -g {base_path2}/{select_genome}/*.fasta --w-add {add_gene_num}',shell=True)    
       #tmp.fastaから5'領域のみを抽出し、5end.txtに格納
       subprocess.run(f'seqkit subseq -r 1:{add_gene_num} {base_path2}/tmp.fasta > {base_path2}/5end.fasta',shell=True)
       #tmp.fastaから5'領域のみを抽出し、3end.txtに格納
       subprocess.run(f'seqkit subseq -r -{add_gene_num}:-1 {base_path2}/tmp.fasta > {base_path2}/3end.fasta',shell=True)
             
       #抽出結果の読込み
       upper =pd.read_table(f'{base_path2}/5end.fasta',header=None)
       end = pd.read_table(f'{base_path2}/3end.fasta',header=None)
       upper_gene_id = []
       upper_seq = []
       end_gene_id = []
       end_seq = []      
      
       for i in range(0,len(upper),2):
         id = upper.iloc[i,:]
         id2 = id.str.replace('>','')
         upper_gene_id.append(id2)
    
       for i2 in range(1,len(upper),2):
         upper_seq.append(upper.iloc[i2,:])
           
       for i3 in range(0,len(end),2):
         id3 = end.iloc[i3,:]
         id4 = id3.str.replace('>','')
         end_gene_id.append(id4)
    
       for i4 in range(1,len(end),2):
         end_seq.append(end.iloc[i4,:])
             
       df = pd.DataFrame(upper_gene_id)
       df2 = pd.DataFrame(upper_seq) 
       df3 = pd.DataFrame(end_gene_id)
       df4 = pd.DataFrame(end_seq)      
      
       dfr = df.reset_index(drop=True)
       df2r = df2.reset_index(drop=True)
       df3r = df3.reset_index(drop=True)
       df4r = df4.reset_index(drop=True)

       ex_upper = pd.concat([dfr,df2r],axis=1,ignore_index=True)
       ex_upper.columns=[['gene_id','sequence']]
       ex_end = pd.concat([df3r,df4r],axis=1,ignore_index=True)
       ex_end.columns=[['gene_id','sequence']]      
      
       os.remove(f'{base_path2}/tmp.fasta')
       os.remove(f'{base_path2}/tmp.fasta.seqkit.fai')
       os.remove(f'{base_path2}/5end.fasta')
       os.remove(f'{base_path2}/3end.fasta') 
           
    #抽出結果の表示
       st.subheader('プロモーター領域')
       st.write('top20件を表示')
       st.dataframe(ex_upper.head(20))
       ex_upper.to_excel(buf := BytesIO(), index=True)
       st.download_button("Download",buf.getvalue(),"all_genes_5end.xlsx","application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
       
       st.subheader('ターミネーター領域')
       st.write('top20件を表示')       
       st.dataframe(ex_end.head(20))
       ex_end.to_excel(buf := BytesIO(), index=True)
       st.download_button("Download",buf.getvalue(),"all_genes_3end.xlsx","application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")       
       st.image(image2,use_column_width=True)
       
if __name__ == "__main__":
    main()