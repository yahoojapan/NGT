NGTQ
===

Neighborhood Graph and Tree for Indexing High-dimensional Data with Quantization

Command
=======



**ngtq** - 大規模高次元ベクトルデータの近傍検索



      $ ngtq command [option] index [data]
        
**注：**
CygWin といった POSIXLY_CORERECT が設定されている環境では、コマンドの前にオプションを指定しなければなりません。

      $ ngtq [option] command index [data]



十億件以上もの高次元ベクトルデータ（数十～数千次元）に対して高速な近傍検索を提供します。

**command** is one of:

-   *[create](#create)*
-   *[append](#append)*
-   *[search](#search)*

### CREATE

指定されたインデックスを生成した上で指定されたデータをインデックスに登録します。

      $ ngtq create -d no_of_dimensions [-p no_of_threads] [-o object_type] [-n no_of_registration_data] 
          [-C global_codebook_size] [-c local_codebook_size] [-N no_of_divisions] 
          [-L local_centroid_creation_mode] 
          index registration_data

*index*  
生成するインデックス名を指定します。データを登録後、本インデックス名のディレクトリが生成されてその中に複数のファイルからなるインデックスが生成されます。

*registration\_data*  
登録するベクトルデータを指定します。１行が１オブジェクト（データ）で構成され、各次元要素のデータはスペースまたはタブで区切られていなければなりません。距離関数はL2のみが利用可能です。

**-d** *no\_of\_dimensions*  
登録データの次元数を指定します。

**-p** *no\_of\_threads* (default = 24, recomended value = number of cores)   
生成時の並列処理時に利用するスレッド数を指定します。

**-o** *object\_type*  
データオブジェクトの型を指定する。
- __c__: 1バイト整数 (一部未実装)
- __f__: 4バイト浮動小数点（デフォルト、推奨）

**-n** *no\_of\_registration\_data*  
Specifies the number of data items to be registered. If not specified, all data in the specified file will be registered.

**-C** *global\_codebook\_size*  
グローバルコード（セントロイド）の数を指定します。

**-c** *local\_codebook\_size*  
ローカルコード（セントロイド）の数を指定します。

**-N** *no\_of\_divisions*  
ローカルベクトルデータ（残差データ）を生成するためのベクトルデータの分割数を指定します。

**-L** *local\_centroid\_creation\_mode*  
ローカルセントロイドを生成するモードを指定します。
- __d__: 指定された登録データの先頭をローカルセントロイドとして使用します。
- __k__: kmeans を使用してローカルセントロイドを生成します。


### APPEND

指定された登録データを指定されたインデックスに追加登録します。

      $ ngtq append [-n no_of_registration_data] index registration_data
        

*index*  
既存のインデックスを指定します。

*registration\_data*  
登録するベクトルデータを指定します。１行が１オブジェクト（データ）で、各次元のデータはスペース又はタブで区切られていなければなりません。

**-n** *no\_of\_registration\_data*  
登録するデータ数を指定します。指定しない場合には指定されたファイル中のすべてのデータを登録します。

### SEARCH

指定されたクエリデータを用いてインデックスを検索します。

      $ ngtq search [-n no_of_search_results] [-e search_range_coefficient] [-m mode] 
          [-r search_radius] [-E approximate-expansion] index query_data
        

*index*  
既存のインデックス名を指定します。

*query\_data*  
クエリデータのファイル名を指定します。１行が１クエリデータであり、登録データと同様に各次元のデータはスペース又はタブで区切られていなければなりません。複数クエリを与えた場合には順次検索します。

**-n** *no\_of\_search\_results* (default: 20)  
検索結果数を指定します。

**-e** *search\_range\_coefficient* (default = recomended value = 0.1)  
グローバルコードブックを検索する時の探索範囲の拡大係数です。大きければ精度が高くなりますが遅くなり、小さければ精度は下がりますが速くなります。0～0.3の範囲内で調整することが望ましいですが、負の値も指定可能です。

**-m** *search\_mode* (__r__|__e__|__l__|__c__|__a__)  
検索モードを
- __a__: 近似距離を用いて検索します。
- __c__: 近似距離を用いて検索します。計算済みローカル距離がキャッシュされることで検索時間が削減されます。（推奨）
- __l__: ローカル距離のルックアップテーブルによる近似距離を用いて検索します。
- __e__: 正確な距離を用いて検索します。ローカルコードブックを利用しません。
- __r__: 近似距離を用いて絞り込んだ後に正確な距離を用いて検索します。（正確な距離が必要な場合には推奨）

**-E** *approximate\_expansion*  
検索結果に対する近似検索結果の割合を指定します。例えば、割合が10で検索結果数が20の場合近似検索結果数は200となります。

