NGT
===

Neighborhood Graph and Tree for Indexing High-dimensional Data

[トップ](/README-jp.md) / [インストール](/README-jp.md#インストール) / [コマンド](/bin/ngt/README-jp.md) / [ライセンス](/README-jp.md#ライセンス) / [関連文献](/README-jp.md#関連文献) / [About Us](http://research-lab.yahoo.co.jp/)

Command
=======



**ngt** - 高次元ベクトルデータ近傍検索



      $ ngt command [option] index [data]
        
**注：**
CygWin といった POSIXLY_CORERECT が設定されている環境では、コマンドの前にオプションを指定しなければなりません。

      $ ngt [option] command index [data]



大量（数百万から数千万データ）の高次元ベクトルデータ（数十～数千次元）に対して高速な近傍検索を提供します。

**command**:

-   *[create](#create)*
-   *[append](#append)*
-   *[search](#search)*
-   *[remove](#remove)*
-   *[prune](#prune)*
-   *[reconstruct graph](#reconstruct-graph)*

### CREATE

指定されたインデックスを生成した上で指定されたデータをインデックスに登録します。

      $ ngt create -d no_of_dimensions [-p no_of_threads] [-b no_of_batch_processes] 
          [-i index_type] [-g graph_type] [-t edge_reduction_threshold] 
          [-e search_range_coefficient] [-E no_of_edges] [-S no_of_edges_at_search_time] 
          [-o object_type] [-D distance_function] [-n no_of_registration_data] 
          index registration_data
        

*index*  
生成するインデックス名を指定します。データを登録後、本インデックス名のディレクトリが生成されてその中に複数のファイルからなるインデックスが生成されます。

*registration\_data*  
登録するベクトルデータを指定します。１行が１オブジェクト（データ）で構成され、各次元要素のデータはスペースまたはタブで区切られていなければなりません。

**-d** *no\_of\_dimensions*  
登録データの次元数を指定します。登録データファイルの各行がすべて次元要素のみで構成されている場合には指定不要ですが、次元要素に続き属性情報などが存在している場合にはこの次元数に基づき後続データを無視します。
    
**-p** *no\_of\_threads* （デフォルト=24、推奨値=コア数）   
生成時の並列処理時に利用するスレッド数を指定します。

**-b** *no\_of\_batch\_processes* (デフォルト=推奨=200）  
通常定数オブジェクトごとに一括登録処理を行います。その一括処理時のオブジェクト数を指定します。通常指定する必要はありません。

**-i** *index\_type* (__g__|__t__)  
グラフのみか、グラフに加えツリーを作成するかを選択します。
- __g__: グラフインデックスのみ生成します。検索時にグラフの探索起点ノードをランダムに決定します。
- __t__: グラフインデックスに加えてツリーインデックスを生成します。検索時にツリーを用いてグラフの探索起点ノードを決定します。（デフォルト、推奨）

**-g** *graph\_type*  
グラフインデックスのタイプを指定します。
- __a__: ANNG を生成します。 高速な登録が可能です。（デフォルト、推奨）
- __k__: KNNG を生成します。 極めて登録が低速です。（実験用なので非推奨）
- __b__: BKNNG を生成します。極めて登録が低速ですが、検索は高速です。（実験用なので非推奨）

**-t** *edge\_reduction\_threshold*（デフォルト=推奨値=0）  
過剰エッジ削減処理を実行する判断基準となる増加エッジ量を指定します。ANNGを選択した時にのみ過剰エッジ削減処理が利用可能となりますが、過剰エッジ削減処理はかなり重い処理なので、メモリ消費量を可能な限り削減したいといった場合を除き利用する必要性はありません。過剰エッジ削減処理自体を実行しない場合には0を指定します。

**-e** *search\_range\_coefficient* （デフォルト=推奨値=0.1） 
ANNGやBKNNGを指定した場合には登録データ（ノード）からエッジで接続される近傍ノードは検索によって獲得され、エッジで結合されます。その検索時の探索範囲の拡大係数です。

**-E** *no\_of\_edges* (default = recommend value = 10)  
グラフ生成時の各ノードの初期エッジ数を指定します。インデックス生成終了時にANNGやBKNNGでは指定されたエッジ数以上のエッジが付与されますが、KNNGでは指定されたエッジ数となります。

**-S** *no\_of\_edges\_at\_search\_time* （デフォルト=推奨値=40）
インデックス生成に伴う検索時及び生成後の検索時に利用するエッジ数を指定します。seachコマンドによる検索時においてエッジ数を指定しない場合にこの値が利用されます。グラフ上の各ノードの実エッジ数よりも少ないエッジ数で検索する場合に指定します。ANNGやBKNNGでは大量のエッジが生成される場合があり、エッジ数を制限することで検索性能が向上する傾向があります。エッジ数を制限しない（実エッジをすべて利用する）場合には0を指定します。0を指定した場合にはインデックスの生成が比較的遅くなりますが、検索時には最も高い性能を得られます。
    
**-o** *object\_type*  
データオブジェクトの型を指定する。
- __c__: 1バイト整数
- __f__: 4バイト浮動小数点（デフォルト）

**-D** *distance\_function*  
距離関数を指定します。
- __1__: L1距離
- __2__: L2距離（デフォルト）
- __a__: 角度
- __A__: 正規化角度。指定されたデータを正規化した上で保存します。
- __c__: コサイン類似度
- __C__: 正規化コサイン類似度。指定されたデータを正規化した上で保存します。
- __h__: ハミング距離。データオブジェクトの型は1バイト整数を指定してください。

**-n** *no\_of\_registration\_data*  
登録するデータ数を指定します。指定しない場合には指定されたファイル中のすべてのデータを登録します。

### APPEND

指定された登録データを指定されたインデックスに追加登録します。

      $ ngt append [-p no_of_threads] [-d no_of_dimensions] [-n no_of_registration_data]
          index registration_data
        

*index*  
既存のインデックスを指定します。

*registration\_data*  
登録するベクトルデータを指定します。１行が１オブジェクト（データ）で、各次元のデータはスペース又はタブで区切られていなければなりません。

**-p** *no\_of\_threads*（デフォルト=推奨値=24）  
生成時の並列処理時に利用するスレッド数を指定します。

**-d** *no\_of\_dimensions*  
登録データの次元数を指定します。登録データファイルの各行がすべて次元要素でのみ構成されている場合には指定不要ですが、次元要素に続き属性情報などが存在している場合にはこの次元数に基づき後続データを無視します。

**-n** *no\_of\_registration\_data*  
登録するデータ数を指定します。指定しない場合には指定されたファイル中のすべてのデータを登録します。

### SEARCH

指定されたクエリデータを用いてインデックスを検索します。

      $ ngt search [-i index_type] [-e search_range_coefficient] [-n no_of_searches] 
          [-E max_no_of_edges] [-r search_radius] index query_data
        

*index*  
既存のインデックス名を指定します。

*query\_data*  
クエリデータのファイル名を指定します。１行が１クエリデータであり、登録データと同様に各次元のデータはスペース又はタブで区切られていなければなりません。複数クエリを与えた場合には順次検索します。
    
**-i** *index\_type* (__g__|__t__|__s__)  
グラフのみか、グラフに加えツリーを利用するかを選択します。
- __g__: グラフインデックスのみ利用します。ツリーが存在する場合でも指定が可能です。グラフの探索起点をランダムに決定します。
- __t__: グラフインデックスに加えてツリーインデックスを利用します。ツリーが存在しない場合にはエラーとなります。ツリーを用いてグラフの探索起点を決定します。グラフのみの場合よりも高速な検索が可能です。（デフォルト、推奨）
- __s__: インデックスを利用せずに線形探索を行います。

**-e** *search\_range\_coefficient* （デフォルト=推奨値=0.1）
探索範囲の拡大係数です。大きければ精度が高くなりますが遅くなり、小さければ精度は下がりますが速くなります。0～0.3の範囲内で調整することが望ましいですが、負の値も指定可能です。

**-n** *no\_of\_searches*（デフォルト：20） 
検索結果数を指定します。

**-E** *max\_no\_of\_edges*（デフォルト=createで指定した値または40、推奨値=40） 
検索時に利用するエッジ数を指定します。グラフ上の各ノードのエッジ数よりも小さいエッジ数で検索する場合に指定します。ANNGやBKNNGでは大量のエッジが生成される場合があり、エッジ数制限を指定することで検索性能が向上する傾向があります。エッジ数制限しない（実エッジをすべて利用する）場合には0を指定します。

**-r** *search\_radius* （デフォルト=無限円）  
検索範囲を円の半径で指定する。

### REMOVE

指定されたオブジェクトをインデックスから削除します。

      $ ngt remove [-d object_id_specification_method] index object_id        

*index*  
既存のインデックス名を指定します。

*object\_id*  
削除対象のオブジェクトIDを指定します。

**-d** *object\_id\_specification\_method* (__f__|__d__) （デフォルト=f） 
削除するオブジェクトIDの指定方法を指定します。fを指定した場合には後述のオブジェクトIDの指定をファイルだとみなします。指定されたファイルには１行ごとに削除するIDが１エントリずつ指定されていなければなりません。dを指定した場合には後述のオブジェクトIDの指定はそのままオブジェクトIDの値だとみなし、削除します。

### PRUNE

指定されたインデックスのグラフ中の長いエッジを削減します。このコマンドにより検索時間が短縮されますが、性能向上には以下の reconstruct graph のパス最適化の利用をお勧めします。

      $ ngt prune -e no_of_forcedly_pruned_edges -s no_of_selectively_pruned_edge index

*index*  
既存のインデックス名を指定します。

**-e** *no\_of\_forcedly\_pruned\_edges*   
各ノードに付与できる最大エッジ数を指定します。各ノードについて指定された数以上のエッジは長い順に強制的に削除します。

**-s** *no\_of\_selectively\_pruned\_edge*   
削除対象としないエッジ数（ランク）を指定します。各ノードについて短い順にエッジを並べた時のランクが指定数を越えて、かつ、そのエッジの代替パスが存在する場合に、そのエッジを削除します。
no\_of\_forcedly\_pruned\_edgesはno\_of\_selectively\_pruned\_edgesより大きくなければなりません。

### RECONSTRUCT GRAPH

指定されたインデックスからグラフを再構成したインデックスを生成します。

      $ ngt reconstruct-graph [-m mode] -o no_of_original_edges -i no_of_reverse_edge input_index reconstructed_index

*input_index*  
既存のインデックス名を指定します。

*reconstructed_index*   
再構成されるインデックス名を指定します。

**-o** *no_of_original_edges*   
再構成されるグラフに付与する元グラフの各ノードの出エッジ数を指定します。この値は再構成されるグラフの出次数の下限値となります。

**-i** *no_of_reverse_edges*   
再構成されるグラフに付与する元グラフの各ノードの出エッジ数を指定します。ただし、出エッジの方向を反転した上で再構成されるグラフに付与されます。この値は再構成されるグラフの入次数の下限値となります。

**-m** *mode*   
グラフのパス最適化のモードを指定します。
- __s__: パス最適化を行いません。
- __S__: パス最適化を行います。




### 生成・登録

128次元１バイト整数型データのインデックスの生成

      $ cd (NGT_TOP_DIR)
      $ ngt create -d 128 -o c index ./data/sift-dataset-5k.tsv
      Data loading time=0.160748 (sec) 160.748 (msec)
      # of objects=5000
      Index creation time=0.379659 (sec) 379.659 (msec)

### 検索

ファイルで指定された3クエリによる近傍検索

      $ cd (NGT_TOP_DIR)
      $ ngt search -n 20 index ./data/sift-query-3.tsv
      Query No.1
      Rank	ID	Distance
      1	3031	239.332
      2	4079	240.002
      3	3164	244.504
      4	3718	246.763
      5	157	251.094
      6	2422	251.185
      7	1313	251.34
      8	379	252.446
      9	3521	260.158
      10	2594	261.132
      11	4627	262.381
      12	2159	263.471
      13	3519	264.909
      14	1764	265.136
      15	4400	266.156
      16	2717	266.914
      17	3168	269.637
      18	4236	270.673
      19	4700	272.725
      20	679	272.973
      Query Time= 0.000472 (sec), 0.472 (msec)
      Query No.2
      Rank	ID	Distance
      1	2726	291.983
      2	924	296.987
      3	3638	298.585
      4	858	300.376
      5	1453	306.805
      6	174	307.789
      7	2992	308.485
      8	2980	308.93
      9	1525	309.49
      10	244	309.816
      11	910	310.446
      12	3310	310.585
      13	2433	311.482
      14	1633	311.735
      15	3761	312.44
      16	407	313.252
      17	4546	313.876
      18	697	315.108
      19	34	315.563
      20	2189	316.193
      Query Time= 0.000478 (sec), 0.478 (msec)
      Query No.3
      Rank	ID	Distance
      1	762	194.286
      2	1046	212.695
      3	4906	215.244
      4	2905	216.539
      5	4142	219.479
      6	1879	219.616
      7	4398	223.352
      8	3842	223.468
      9	233	224.127
      10	2794	224.366
      11	2476	224.804
      12	1848	225.803
      13	3364	226.561
      14	4098	226.74
      15	3023	228.884
      16	4113	229.325
      17	1036	232.852
      18	1740	233.144
      19	2302	233.818
      20	2440	233.91
      Query Time= 0.00018 (sec), 0.18 (msec)
      Average Query Time= 0.000376667 (sec), 0.376667 (msec), (0.00113/3)



#### [ONNG](/README.md#onng)
```
$ ngt create -i t -g a -S 0 -e 0.0 -E no_of_edges -d dimensionality_of_data -o data_type -D distatnce_type anng-index vector-data.dat
$ ngt reconstruct-graph -m S -o outdegree -i indegree anng-index onng-index
```
e.g.  
```
$ ngt create -i t -g a -S 0 -e 0.0 -E 100 -d 128 -o c -D 2 anng-index vector-data.dat
$ ngt reconstruct-graph -m S -o 10 -i 120 anng-index onng-index
```
#### [PANNG](/README.md#panng)
```
$ ngt create -i g|t -g a -S 0 -e 0.1 -E no_of_edges -d dimensionality_of_data -o data_type -D distatnce_type panng-index vector-data.dat
$ ngt prune -e no_of_forcedly_pruned_edges -s no_of_selectively_pruned_edges panng-index
```
e.g.  
```
$ ngt create -i t -g a -S 0 -e 0.1 -E 10 -d 128 -o c -D 2 panng-index vector-data.dat
$ ngt prune -e 60 -s 30 panng-index
```
#### [ANNGT](/README.md#anngt)
```
$ ngt create -i t -g a -S 0 -e 0.1 -E no_of_edges(k) -d dimensionality_of_data -o data_type -D distance_type anngt-index vector-data.dat
```
e.g.
```  
$ ngt create -i t -g a -S 0 -e 0.1 -E 16 -d 128 -o c -D 2 anngt-index vector-data.dat
```
#### [ANNG](/README.md#anng)
```
$ ngt create -i g -g a -S 0 -e 0.1 -E no_of_edges(k) -d dimensionality_of_data -o data_type -D distance_type anng-index vector-data.dat
```
e.g.
```
$ ngt create -i g -g a -S 0 -e 0.1 -E 16 -d dimensionality_of_data -o data_type -D distance_type anng-index vector-data.dat
```
#### KNNG  
```
$ ngt create -i g -g k -S 0 -E no_of_edges(k) -d dimensionality_of_data -o data_type -D distance_type knng-index vector-data.dat
```
e.g.
```  
$ ngt create -i g -g k -S 0 -E 20 -d 128 -o c -D 2 knng-index vector-data.dat
```
