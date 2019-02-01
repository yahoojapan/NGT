
# python NGT

[日本語](README-jp.md)

## Install
Python binding with pybind11 (ngtpy) is installed as follows.
```
pip3 install ngt
```
If you would like to use python binding with ctypes (ngt), additionally you have to install the NGT library according to the [README](../README.md#build).

You can install the python bindings with pybind11 and ctypes from source code. You **MUST** install the NGT library according to the [README](../README.md#build) before installing the python bindings as follows.
```
pip3 install pybind11
pip3 install numpy
cd NGT_ROOT/python
python3 setup.py sdist
pip3 install dist/ngt-x.x.x.tar.gz
```

## Documents

[ngtpy (pybind11) reference](README-ngtpy.md)

## Simple samples

### ngtpy (pybind11)

ngtpy(pybind11) can reduce the processing times than ngt(ctypes). It is more effective especially for the short search time. 

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

See also [sample.py](sample/sample.py).

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
