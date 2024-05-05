# Brakerによる遺伝子領域の予測

# 参考：

- [Git_Gaius-Augustus/BRAKER](https://github.com/Gaius-Augustus/BRAKER)
- [Macでバイオインフォマティックス](https://kazumaxneo.hatenablog.com/entry/2020/08/14/133846)
- [UTM](https://envader.plus/article/66)
- [UTMdownload](https://github.com/utmapp/UTM/releases)

# BRAKERインストール

- インストール
  
`git clone https://github.com/Gaius-Augustus/BRAKER.git`

- 実行権限の付与
  
`cd [path_to_BRAKER]`

`cd scripts`　ー＞　`chmod a+x *.pl *.py`

# モジュールインストール

`conda install -c anaconda perl`

`conda install -c anaconda biopython`

`conda install -c bioconda perl-app-cpanminus`

`conda install -c bioconda perl-file-spec`

`conda install -c bioconda perl-hash-merge`

`conda install -c bioconda perl-list-util`

`conda install -c bioconda perl-module-load-conditional`

`conda install -c bioconda perl-posix`

`conda install -c bioconda perl-file-homedir`

`conda install -c bioconda perl-parallel-forkmanager`

`conda install -c bioconda perl-scalar-util-numeric`

`conda install -c bioconda perl-yaml`

`conda install -c bioconda perl-class-data-inheritable`

`conda install -c bioconda perl-exception-class`

`conda install -c bioconda perl-test-pod`

`conda install -c bioconda perl-mce`

`conda install -c bioconda perl-threaded`

`conda install -c bioconda perl-list-util`

`conda install -c bioconda perl-math-utils`

`conda install -c bioconda cdbtools`

`conda install -c bioconda perl-data-dumper`

`conda install -c eumetsat perl-yaml-xs`

# GENEMark-ETPインストール

`git clone https://github.com/gatech-genemark/GeneMark-ETP.git`

`perl change_path_in_perl_scripts.pl "/usr/bin/env perl"`

# AUGUSTUSインストール

`git clone https://github.com/Gaius-Augustus/Augustus.git`

# Bamtoolsのインストール

`git clone https://github.com/pezmaster31/bamtools.git`

`cd [path_to_bamtools]`

`mkdir build` -> `cd build` -> `cmake` -> `make` -> `make DESTDIR=/home/**`

# TSEBRAのインストール
`git clone https://github.com/Gaius-Augustus/TSEBRA`

# ProtHintのインストール
`wget https://github.com/gatech-genemark/ProtHint/releases/download/v2.6.0/ProtHint-2.6.0.tar.gz`

`tar xf ProtHint-2.6.0.tar.gz`

# UTM上の共有Fの設定コマンド

`sudo mount -t 9p -o trans=virtio share {/home/****/path_to_folder @UTM側 } -oversion=9p2000.L`

# 関連モジュールのpath設定

- BRAKER、GeneMark、AUGUSTUS、BAMTOOLS、TSEBRAのパス設定
- スクリプトに追加する

`export PATH=/path/to/BRAKER/scripts/:$PATH`

`export PATH=/path/to/GeneMark-ETP/tools/:$PATH`

`export GENEMARK_PATH=/path/to/GeneMark-ETP/bin`

`export AUGUSTUS_CONFIG_PATH=/path/to/Augustus/config`

`export AUGUSTUS_BIN_PATH=/path/to/Augustus/bin`

`export AUGUSTUS_SCRIPTS_PATH=/path/to/Augustus/scripts`

`export BAMTOOLS_PATH=/path/to/usr/bamtools/local/bin`

`export TSEBRA_PATH=/path/to/TSEBRA/bin`

`export PROTHINT_PATH=/path/to/ProtHint-2.6.0/bin`




