pythia の解析をcondor使って行うものたち
pythia のデータを作る時のcondorでは1000イベント × 100,000 jobで1億イベント作っていたので、それを一部まとめて 50000 イベント × 2000 ファイルにして、そのファイルを使ってこの解析を行っている

makelist_ana_pythia.sh → 2000ファイルのpathを一行ずつ羅列するlistを作るシェルスクリプト
input_pythia.list → ↑で作ったlist。データはこのlist内のpathを参照ください

pythia_job.sh → ↑のlistを使って解析マクロを回すシェルスクリプト。job投げる時に使う
pythia_condor.job → ↑のshを使ってjobを投げる

pythia_merge.sh → condorで回した解析結果の入ったroot fileを一つにくっつける

./DrawHist.C → 一つにまとめたroot fileの中の必要なhistをDrawする
$ root -q -b 'DrawHist.C(1)'