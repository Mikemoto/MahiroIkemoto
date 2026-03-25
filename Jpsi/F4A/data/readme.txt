

実際のデータdstファイルから解析に使うrootファイルを作るF4Aマクロとそのためのcondorの設定とかがあります。

データdstは自力で探し出す必要があります。
ここにあるのは2026年1月時点で私が見る蹴ることのできた全てのデータになります。

Fun4All_DataDST_SiliconSeedAna.C → メインとなるF4Aマクロです。
これを回すのにSIliconSeedsAnaモジュールを使用しています。
また、扱うデータはとても大きく多いので基本的にはcondorを使用するので、このマクロ単体で回すことはほぼないです。

Calo_Calib.C → EMCalのクラスタリングをする&EMCalのキャリブレーションをもとの設＝silicon seeds ana モジュールを使用する時に矛盾が発生しない設定 にするためにF4Aマクロ内で使っているものです。
詳しくは私の20251111と20251209重イオンMTスライドを見てください。

Charge_Recalc.C → 粒子の電荷を正しく計算するマクロです。F4Aマクロ内で使っています。
これを初めて使う時、F4Aマクロの中のCharge_Recalc.Cのパスを自身のものに変えてください。

localAlignmentParamsFile1.txt → INTTのジオメトリがおかしいのを治すために使っていましたが、これを入れるとF4Aマクロを回す際に衝突点データを取って来れなくなるというバグが生じるので使っていません。
使うならこのファイルの名前をlocalAlignmentParamsFile.txtに変えてください。
INTTジオメトリ問題については私の20250623重イオンスライドMTを参照してください。

run_condor.job、run_job.sh → condor

majitest.list → マジのテストをする時用のlistです。マジのテストをする時はこのlistを読み込ませてcondor投げてください。マジのテストでは50000イベントのデータを作るjobを一つ投げます。

list/ → condorでdstファイルから解析用データを作るのに使ったlistの一覧です。修論に使用した全ラン分、3種類のdstファイル全てのものがありますが、重すぎたのでpathだけ書いています。
このlistは、dstファイルのパスだけ書いてあります。
Caloのdstファイルは1ファイル50000イベント分、silicon trackerのdstファイル(seedとclus)は1ファイルに10000イベント入っているので、1つのジョブでcalo dst は1ファイル、seed dstとclus dstは5ファイルずつ読めるようにしています。このリストたちは直接的には使用しません。condor投げる時に使うlist内で使用しています。

alllist/ → condor投げる時に使うlistです。↑のlistとか出力ファイル名とかを書いています。使用したrun全てあります。

makelist.sh → 上二つのlist全てを一気に作ってくれるシェルスクリプトです。大変便利です！。run番号と、そのrunのdstファイルのpath、およびcalo dstファイルの数は自分で探して調べて書く必要があります。


F4Aマクロを回す時の環境設定
とりあえずこれらを順に打てばなんとかなる。

source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.527
cd SiliconSeedAna/build
../autogen.sh --prefix=${PWD}/../install
make
make install
source /opt/sphenix/core/bin/setup_local.sh installディレクトリのフルpath

2回目以降
source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.527
source /opt/sphenix/core/bin/setup_local.sh installディレクトリのフルpath

毎回するのがめんどくさい時
sphrnixサーバーにssh接続した時に入るメインホームディレクトリ内の
.bash_profileの中に上の二行を書く
.bash_profileはサーバーにsshでログインする時に自動的に読み込まれるものです。書き換えた後は手動で読み込み直す必要があります

$ source ./bash_profile
