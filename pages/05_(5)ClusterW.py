from Bio.Seq import Seq
import os
import shutil
import subprocess
import time
import streamlit as st
import glob

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.abspath(os.path.join(current_dir,os.pardir))
fig_path = f"{parent_dir}/File/clustero.png"

def Comparing_SeqData(data_path):
    """シーケンスデータの向きを判定
    Parameter
    ---------
    data_path: str
        clustaloの多重整列結果データ（.txt)のパス
    
    Returns
    -------
    correct_seq: str
        正しい配列（seq1もしくはseq2）
    """

    with open(data_path,"r") as f:
        data = f.readlines()

    ex_data = data[1].split("\n")[0]
    ex_data2 = ex_data.split(" ")
    ex_data3 = [i for i in ex_data2 if i != ""]

    if float(ex_data3[2])>float(ex_data3[3]):
        return "seq1"
    else:
        return "seq2"

def make_query(gene_name,gene_seq,tmp_seq_data_F_path,data_name):
    """sequenceデータ(.txt)から相補配列データを生成
       multi fasta形式で記載したtmp.fastaをseq_data_F_path下に生成
    Paramerter:
    ----------
    gene_name: str
        元配列名
    gene_seq: str
        元配列データ
    tmp_seq_data_F_path: str
        sequenceデータが保管されたフォルダパス
    data_name: str
        シーケンスファイル名（拡張子は含まない！）
    """
    
    gene_seq_data = f">{gene_name}\n{gene_seq}"
    
    data_path = f"{tmp_seq_data_F_path}/{data_name}.txt"
    
    with open(data_path,'r') as f:

        original_header_n = ">" + data_name + "\n"
        original_seq_data = f.read()

        reverse_header_n =  ">" + data_name + "_" + "\n"
        reverse_seq_data = str(Seq(original_seq_data).reverse_complement())

        body = gene_seq_data + "\n" + original_header_n + original_seq_data + "\n" + reverse_header_n + reverse_seq_data
        
        with open(f"{tmp_seq_data_F_path}/tmp.txt","w") as f2:
            f2.write(body)
    
        os.rename(f"{tmp_seq_data_F_path}/tmp.txt",f"{tmp_seq_data_F_path}/tmp.fasta")
                  
    return f"{tmp_seq_data_F_path}/tmp.fasta",original_header_n+original_seq_data,original_header_n+reverse_seq_data

def make_query2(tmp_path,seq_data):
    """マルチファスタ形式のファイル生成
    　 final_seq.txtに追記
    Paramerter:
    ----------
    clustalO_tmp: str
        final_seq.txtデータ保管先のフォルダパス
    seq_data: str
        fasta形式ファイル(>{read_name}\n{atgcc・・・})
    """
    
    with open(f'{tmp_path}/final_seq.txt',"a") as f:
        body = seq_data + "\n"
        f.write(body)
    
# path情報
seq_data_F_path = '/home/seika_oiwa_3590/notebooks/box_path/bio_team/ClusterW'

st.title("シーケンスデータのアラインメント")
st.image(fig_path)

# 解析情報の設定
folder_lists_ = os.listdir(seq_data_F_path)
folder_lists = [i for i in folder_lists_ if i != '.ipynb_checkpoints']

select_folder = st.selectbox('フォルダを選択',folder_lists,key="f1")
tmp_seq_data_F_path = f"{seq_data_F_path}/{select_folder}"

# フォルダ内データセットの表示
data_sets_ = glob.glob(f'{tmp_seq_data_F_path}/*.txt')
data_set = [os.path.basename(i) for i in data_sets_]
data_set2 = " / ".join(data_set)

st.write("フォルダ内のデータ")
st.write(data_set2)

# ベース遺伝子の情報
gene_name = st.text_input("遺伝子名を入力",key="n1")
gene_seq = st.text_area("ベース遺伝子配列を入力",key="n2")

# 解析開始
if st.button("アラインメント開始"):
    
    # 作業用フォルダ作成
    result_path = f'{tmp_seq_data_F_path}/Result'
    os.makedirs(result_path,exist_ok=True)

    # クラスター解析用のベースファイル作成(final_seq.txt),ベース配列のみ
    with open(f"{result_path}/final_seq.txt","a") as f:
        body = f">{gene_name}\n{gene_seq}\n"
        f.write(body)

    for read_path in glob.glob(f'{tmp_seq_data_F_path}/*.txt'):
        data_name = os.path.splitext(os.path.basename(read_path))[0]

        # queryデータの作成 (インプットしたseqデータの相補配列（２種類）と、ベースとなる遺伝子配列を含むmulti fastaファイルを生成)
        fasta_path,seq1,seq2 = make_query(gene_name,gene_seq,tmp_seq_data_F_path,data_name)

        # cluster解析
        time.sleep(2)
        subprocess.run(f"clustalo --full --percent-id --distmat-out=output.txt -i {fasta_path}",shell=True)

        # 相補配列のどちらが正しい向きか選択
        result_seq = Comparing_SeqData("output.txt")

        # 正しい向きの相補配列データをfinal_seq.txtに追加
        if result_seq == 'seq1':
            make_query2(f"{result_path}",seq1)
        if result_seq == 'seq2':
            make_query2(f"{result_path}",seq2)

        # 不要ファイルの削除
        os.remove("output.txt")
        os.remove(fasta_path)

    # multi-clustal
    os.rename(f'{result_path}/final_seq.txt',f'{result_path}/final_seq.fasta')
    subprocess.run(f"clustalo -t DNA -i {result_path}/final_seq.fasta -o {result_path}/multi_alingment_result_{gene_name}.fasta",shell=True)
    
    with open(f"{result_path}/multi_alingment_result_{gene_name}.fasta","r") as f:
        st.download_button("結果のダウンロード",f,f"multi_alingment_result_{gene_name}.fasta")
    
        
    
