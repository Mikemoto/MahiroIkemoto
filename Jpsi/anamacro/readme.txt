解析マクロたち。

DrawMassDis_data.C → メイン解析マクロ。dataと名前に入っているけどsimulationでも回ります。


dphidz_fitting.C → dphiとdz分布にガウスfittingして幅を取ってくるマクロです。dphiは補正前と補正後どちらもにfittingを行います。ここで得た値をメイン解析マクロ内で使っています

dphi_alignment.C →　dphiのずれをなおすための値を計算するマクロです。具体的には、phi0 VS dphi(補正前)の2dhistをもとに、phi0を0.2 rad ずつで区切ったdphi分布にガウスfitをし、それぞれのdphiのピーク位置が中心0radにくるようにずらす関数を計算しています
。ここで得た関数をメイン解析マクロ内で使っています。


dphi_phi0_alighnment.root → ↑で求めたdphiアラインメントをする時に使う関数が入ったrootfileです。これを解析マクロに読み込ませることで使用しています。

DrawBG.C → シミュレーションデータ解析に用いたマクロです。single simulationとPythia p+p simulationの結果をもとに、規格化したsignal とpythiaの質量分布を重ね合わせるマクロです。

recalc_charge.C → まるで電荷を再計算するような名前のマクロですが、これはEvent displayを描くマクロです。
使用したデータとかはこのマクロないに書いてあるpathを見てください(ごめんなさい)

result/ → 結果ファイルをおいているpathを書いています

condor_ana/ → 解析マクロをcondor使って回す

data/ → 使用したデータのpathを書いています
