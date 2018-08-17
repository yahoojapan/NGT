ngtpy
=====

NGT python

Class Index
===========

## Member Functions

### \_\_init\_\_
指定されたインデックスをオープンし、そのインデックスのオブジェクトを生成します。

      __init__(self: ngtpy.Index, path: str, read_only: bool=False, zero_based_numbering: bool=True)

**Returns**  
なし

**path**   
オープンするインデックスのパスを指定します。

**read\_only**   
書込み不可でインデックスをオープンします。Falseは書込み可を意味します。

**zero\_based\_numbering**   
オブジェクトIDが０から始まることを指定します。Falseは１から始まることを意味します。

### close
インデックスをクローズします。

      close(self: ngtpy.Index)

**Returns**   
なし

### insert
指定されたオブジェクトをインデックスに登録します。そのオブジェクトのインデックスは生成されません。この関数でオブジェクトを登録後に以下の関数build_indexを呼び出してインデックスを生成しなければなりません。

      int insert(self: ngtpy.Index, object: numpy.ndarray[float64])

**Returns**   
登録したオブジェクトのID

**object**   
登録するオブジェクトを指定します。

### build_index
関数insertで登録されたオブジェクトのインデックスを生成します。

      build_index(self: ngtpy.Index, num_threads: int=8)

**Returns**   
なし

**num_thread**  
インデックスを生成する時に利用するスレッド数を指定します。

### batch_insert
指定した複数のオブジェクトを登録し、そのオブジェクトのインデックスを生成します。この関数はngtコマンドの"ngt append"を実行するのとほぼ同じです。この関数を呼び出す代わりにngtコマンドを使っても構いません。

      batch_insert(self: ngtpy.Index, objects: numpy.ndarray[float64], num_threads: int=8)

**Returns**  
なし

**objects**   
登録する複数のオブジェクトを指定します。

**num_thread**   
インデックスを生成する時に利用するスレッド数を指定します。


### remove
指定されたオブジェクトを削除します。

      remove(self: ngtpy.Index, object_id: int)

**Returns**  
なし

**object\_id**   
削除するオブジェクトのIDを指定します。

### save
インデックスを保存します。

      save(self: ngtpy.Index)

**Returns**  
なし

### get_object
指定されたオブジェクトを取得します。

      List[float] get_object(self: ngtpy.Index, object_id: int)

**Returns**   
指定されたオブジェクト

### search
指定されたクエリオブジェクトに対する近傍のオブジェクトを検索する。

      object search(self: ngtpy.Index, query: object, size: int=10, epsilon: float=0.1, edge_size: int=-1, with_distance: bool=True)

**Returns**   
検索結果としてタプル（ID、距離）のリスト

**query**   
クエリオブジェクトを指定します。

**size**   
検索結果として返るオブジェクトの数を指定します。

**epsilon**   
グラフの探索範囲を決定する変数イプシロンを指定します。

**edge\_size**   
グラフを探索するのに利用する各ノードのエッジ数を指定します。

**with\_distance**   
距離付きオブジェクトIDの検索結果を返すことを指定します。Falseは結果がオブジェクトIDののみのリストとなることを意味します。

FUNCTIONS
=========

### create
空のインデックスを生成します。この関数はngtコマンドの"ngt create"を実行するのとほぼ同じです。この関数を呼び出す代わりにngtコマンドを使っても構いません。

      create(path: str, dimension: int, edge_size_for_creation: int=10, edge_size_for_search: int=40, distance_type: str='L2', object_type: str='Float')

**Returns**   
なし

**path**   
インデックスのパスを指定します。

**dimension**  
登録するオブジェクトの次元数を指定します。

**edge\_size\_for\_creation**   
各ノードの初期エッジ数を指定します。

**edge\_size\_for\_search**   
検索時にグラフを探索するためのノードのエッジ数を指定します。

**distance\_type**   
オブジェクトの距離関数を指定します。
- __L1__: L1 距離
- __L2__: L2 距離（デフォルト）
- __Angle__: 角度距離
- __Normalized Angle__: 正規化角度距離。指定されたデータは自動的に正規化された上でインデックスに登録されます。
- __Cosine__: コサイン類似度
- __Normalized Cosine__: 正規化コサイン類似度。指定されたデータは自動的に正規化された上でインデックスに登録されます。
- __Hamming__: ハミング距離

**object\_type**  
オブジェクトのデータタイプを指定します。
- __Float__: 4 バイト浮動小数点
- __Byte__: 1 バイト符号なし整数

