
# python NGT


## インストール
次のコマンドでpybind11 (ngtpy) の python binding がインストールできます。
```
pip3 install ngt
```
ctypes (ngt) の python binding を利用になりたい場合には、さらに、[README-jp](../README-jp.md#build)にしたがってNGTライブラリをインストールしてください。

また、ソースコードからも pybind11 と ctypes のpython binding をインストールできます。ただし、python NGTをインストールする前にNGTライブラリを[README-jp](../README-jp.md#build)にしたがって**必ず**インストールしてください。
```
pip3 install pybind11
pip3 install numpy
cd NGT_ROOT/python
python3 setup.py sdist
pip3 install dist/ngt-x.x.x.tar.gz
```

## ドキュメント

[ngtpy (pybind11) リファレンス](README-ngtpy-jp.md)

## サンプルコード

### ngtpy (pybind11)

ngtpy(pybind11)はngt(ctypes)より処理時間を削減できます。特に短時間の検索において効果があります。

```python
  import ngtpy
  import random

  dim = 10
  objects = []
  for i in range(0, 100) :
      vector = random.sample(range(100), dim)
      objects.append(vector)

  query = objects[0]

  ngtpy.create(b"tmp", dim)
  index = ngtpy.Index(b"tmp")
  index.batch_insert(objects)
  index.save()

  result = index.search(query, 3)

  for i, o in enumerate(result) :
      print(str(i) + ": " + str(o[0]) + ", " + str(o[1]))
      object = index.get_object(o[0])
      print(object)
```

ご参考： [sample.py](sample/sample.py).

### ngt (ctypes)

```python
  from ngt import base as ngt
  import random

  dim = 10
  objects = []
  for i in range(0, 100) :
      vector = random.sample(range(100), dim)
      objects.append(vector)

  query = objects[0]
  index = ngt.Index.create(b"tmp", dim)
  index.insert(objects)
  # You can also insert objects from a file like this.
  # index.insert_from_tsv('list.tsv') 

  index.save()
  # You can load saved the index like this.
  # index = ngt.Index(b"tmp")

  result = index.search(query, 3)

  for i, o in enumerate(result) :
      print(str(i) + ": " + str(o.id) + ", " + str(o.distance))
      object = index.get_object(o.id)
      print(object)
```
