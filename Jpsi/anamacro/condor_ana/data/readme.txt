解析マクロをcondorで回す時に使うものたち
データを作る時に50,000イベントずつのファイルをcondorで作っているので、そこで作ったファイルをそのまま使って解析を行う

makelist_ana.sh → 5万イベントが入ったそのrunのデータファイルすべてを一行ずつ羅列した list file を作るシェルスクリプトです。
$ ./makelist_ana.sh ラン番号

/list → ↑で作ったlistが全てぶちこんであります。今sphenixのサーバー上で見つけられる限りのデータから作ったファイルたちのpathがここにまとまっているはずです。データファイルは容量が大きすぎるのでここには置いていません、listに書いてあるpathを辿っていただければと思います

mass_ana_job.sh → ↑で作ったデータファイルのpathが書かれたlistを用いて解析マクロを回すシェルスクリプトです。condor job投げる時に使います
mass_ana_condor.job → ↑のshを使ってjobを投げます。5万イベントの解析なので1jobは基本一瞬(数十秒もかからない)で終わるはずです
$ condor_submit mass_ana_condor.job -append 'RUN=ラン番号'

merge.sh → 解析して出力したroot fileを一つにくっつけるシェルスクリプトです。愚かにもhaddコマンドで一つ一つ直列にファイルをくっつける作業をしているので時間がかかります。
$ ./merge.sh ラン番号

../DrawHist.C → mergeしたrootファイルをもとに必要なhistgramをPDFにDrawするマクロです
$ root -q -b 'DrawHist.C(ラン番号)'

Mass_reco.sh → makelist ~ 解析job投げる ~ merge ~ hist Draw までを一気にやってくれるシェルスクリプトです。
$ ./Mass_reco.sh ラン番号

ここまでは1runごとに行います

allmerge.sh → 全てのrunで同様の↑作業をした後、それぞれの root file をくっつけてくれるシェルスクリプトです。
allmerge.shで出力したroot fileを使ってDrawHist.Cを回すことで、全てのrunデータで出した結果を確認できます。修論に載せたデータの結果はこれでDrawしたものを使っています。