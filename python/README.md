
# python NGT

## Install

You **MUST** install the NGT library according to the [README](../README.md#build) before installing the following python NGT.

```
cd NGT_ROOT/python
python setup.py sdist
pip install dist/ngt-1.0.0.tar.gz
```

## Usage

```python
  from ngt import base as ngt
  import random

  dim = 10
  objects = []
  for i in range(0, 100) :
      vector = random.sample(range(100), dim)
      objects.append(vector)

  query = objects[0]
  index = ngt.Index.create("tmp", dim)
  index.insert(objects)
  # You can also insert objects from a file like this.
  # index.insert_from_tsv('list.tsv') 

  index.save()

  result = index.search(query, 3)

  for i, o in enumerate(result) :
      print(str(i) + ": " + str(o.id) + ", " + str(o.distance))
      object = index.get_object(o.id)
      print(object)
```

