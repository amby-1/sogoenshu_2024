# このページについて
本ページは，大阪大学基礎工学研究科，総合工学演習2024 テーマC3 「数値シミュレーションによる身体性の理解」の授業用サイトです．
適宜参考にしながら課題に取り組んでもらいます．
なお，本サイトは授業終了後に予告なく閉鎖することがあります．

# 青井研究室の紹介
ヒトや動物は、冗長で複雑な筋骨格系を巧みに動かし、多様な環境で適応的に歩行する．
歩行とは、脳・身体・環境の相互作用を介して初めて形成される力学現象であり，このメカニズムを理解することができれば，生物の理解やロボットへの応用につながる．
[青井研究室](https://mechbiosys.me.es.osaka-u.ac.jp/)では，生物の運動計測から背後に潜むメカニズムについて仮説をたて，シミュレーションにより生物の動きを再現して仮説を検証，さらにはロボットの制御への応用を目指しています．

# 演習目的
人や動物は，身体の持つ力学特性をうまく使って，多様な環境で優れた運動能力を示す（身体性）．本テーマでは，特に人の走行運動を模擬できる簡単なモデルを対象に物理シミュレーションを行い，身体性についての理解を深めるとともに，数値シミュレーション・数値解析の基礎を学習する．

# 演習内容
人の走行運動の重要な特徴は，バネ（脚を近似）と質点（胴体を近似）によって構成される，アクチュエータが無い受動モデル（Spring-loaded Inverted Pendulum（SLIPモデル））で再現できる．このモデルでは，質点に適切な初期速度を与えるだけで，制御なしに脚が接地と離地を繰り返し，人の走行運動に似た運動を続けられる．これは，身体の持つ力学特性によって基本的運動を実現することができる，身体性の良い事例である．数値計算手法を復習しながら，このモデルの物理シミュレータを自作，その運動を数値解析手法で分析する．
- 数値計算法の復習：バネとダンパで釣られた質点（バネダンパ系）を対象とし，運動方程式を導出，ルンゲクッタ法等の数値計算手法を復習しながら，数値シミュレータを構築する（サンプルはC/C++, gnuplotで用意）．
- SLIPモデルの物理シミュレータ構築：1回目の復習をもとに，SLIPモデルの数値シミュレーションを構築する（サンプルはC/C++, gnuplotで用意）．
- SLIPモデルの物理シミュレータ構築・運動分析：SLIPモデルの数値シミュレーションを完成させ，その運動を分析する．適切な初期値を選べば，本当に走行運動を続けることができるのだろうか？ また，初期値を少し変えたとしても，同じような走行運動が出てくるのか？　などを基本に，各自の興味に応じて分析を行う．
以下のような走行シミュレーションを各自で作ってもらいます．
<div style="display: flex; align-items: flex-start;">
  <img src="Figs/SLIP_overview.png" alt="SLIPモデル" width="300" style="vertical-align: middle; margin-right: 20px;">
  <img src="Figs/animation.gif" alt="Animation" width="300" style="vertical-align: middle;">
</div>

# 注意点等
- Ｃ/Ｃ++言語で講義資料やサンプルプログラムは準備します．
  - 課題がこなせるのであれば，プログラミング言語は問いません．他の言語のサンプル等は用意できませんので，自己判断ください．
  - 授業当日もPCはネット環境に繋ぐ可能性があるので，odins等に接続できるように設定しておいてください．
- データはこまめにバックアップを取ってください
  - PCの不調等によりデータが失われた場合については，自己責任となります．特別対応はしない予定ですのでご注意ください．

# 授業計画
以下の計画で授業を進める予定です．

## [1日目 バネダンパでつながれた質点の運動をシミュレーションしよう](https://github.com/amby-1/sogoenshu_2024/blob/main/Day_1.md)
バネダンパ系を用いて数値計算法を復習します

## [2日目 SLIPモデルを用いたシミュレーションを作ろう１](https://github.com/amby-1/sogoenshu_2024/blob/main/Day_2-3.md)
SLIPモデルのシミュレータを作ります．

## [3日目 SLIPモデルを用いたシミュレーションを作ろう２](https://github.com/amby-1/sogoenshu_2024/blob/main/Day_2-3.md)
SLIPモデルのシミュレータを完成させ，周期的な走行運動が実現できるかを確認します．

# レポート
- 各日の授業資料には，複数の課題や応用課題（任意）が含まれ，これらの課題についての取り組みをレポートとして提出してもらう予定です
  - 応用課題については，課題がすべて終わって時間を持て余した場合に解いてください．基本的には課題を解いてもらえれば合格です．
- レポートは日にちごとに1つのPDFにまとめていただき，表紙にタイトル（何日目のレポートかを明記），日付，氏名，学籍番号を忘れずに記載してください．
- 各日程のレポートは，まず演習内容について簡潔にまとめた後，各課題について，番号順に解答を載せてください．
  - 課題の問題文もレポートには記載してください（要約可）．
  - なお，応用課題をやった人は，通常課題の後にそれも載せてください．分量に応じて加点します．
- 原則として次の授業日の前日を締め切りとします．
   - なお，3日目分のレポートについては，3日目の授業時に締め切り日を案内します（1週間は取り組み時間として確保します）．
- 提出方法はメール提出を予定しています．
  - メール内容
    - 件名：総合演習レポート?日目（授業日10/??）(名前)
    - 宛先：ambe.yuichi.es　あっと　osaka-u.ac.jp
  - 受領連絡
    - 1-2日目のレポートの場合，受領連絡は次回の授業時におこないます(メールの返信は行いません)
    - 3日目のレポートについては，受領連絡は1日以内にメール返信にておこないます．（受領連絡がない場合は再送して下さい．）
