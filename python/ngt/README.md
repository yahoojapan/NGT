<h1 id="ngt">ngt</h1>


<h1 id="ngt.base.Index">Index</h1>

```python
Index(self, path)
```

NGT: Neighborhood Graph and Tree for Indexing High-dimensional Data  
  NGT provides the functionality of searching for approximate nearest neighbors in high-dimensional data.

Example:
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
  # index.insert_from_tsv('list.dat')

  index.save()

  result = index.search(query, 3)

  for i, o in enumerate(result) :
      print(str(i) + ": " + str(o.id) + ", " + str(o.distance))
      object = index.get_object(o.id)
      print(object)
```

<h2 id="ngt.base.Index.insert">insert</h2>

```python
Index.insert(self, objects, num_threads=8)
```

insert the specified objects into the index and build the index.

    objects : Inserted objects.
    return  : List of the IDs of the inserted objects.

<h2 id="ngt.base.Index.get_object">get_object</h2>

```python
Index.get_object(self, id)
```

get the specfied object by id.

    id : Object id.

<h2 id="ngt.base.Index.insert_object">insert_object</h2>

```python
Index.insert_object(self, object)
```

insert the specified object into the index.  
must call build_index after call this method.

    object : Inserted object.
    return : The ID of the inserted object.

<h2 id="ngt.base.Index.insert_from_tsv">insert_from_tsv</h2>

```python
Index.insert_from_tsv(self, path, num_threads=8, dlmt='\t')
```

insert objects in the specified file and build the index.

    path : Path of the object file.
    num_threads : Number of threads in building index.
    dlmt : Delimiter to sepalate each element in the object file.

<h2 id="ngt.base.Index.remove">remove</h2>

```python
Index.remove(self, id)
```

remove the specified object by id

    id : Object id.

<h2 id="ngt.base.Index.build_index">build_index</h2>

```python
Index.build_index(self, num_threads=8)
```

build the inserted objects into the index.

    num_threads : Number of threads in building index.

<h2 id="ngt.base.Index.search">search</h2>

```python
Index.search(self, query, k=20, epsilon=0.1)
```

search for the k nearest neighbors of the specifiecd query object.

    query   : Query object.
    k       : Number of searched objects.
    epsilon : Epsilon defines a search range.

<h2 id="ngt.base.Index.create">create</h2>

```python
Index.create(path, dimension, edge_size_for_creation=10, edge_size_for_search=40, object_type='Float', distance_type='L2')
```

create an empty index with the specified parameters.
  edge_size_for_creation : Number of edges for each node in the graph.
  edge_size_for_search   : Number of edges to search.
  object_type            : Type of the data object. (Float, Integer [Integer is 1 byte])
  distance_type          : Type of the distance function. (L1,L2,Angle,Hamming)

<h2 id="ngt.base.Index.save">save</h2>

```python
Index.save(self, path=None)
```

save the index.

    path : Path to save the index. default overwrite the files.

<h2 id="ngt.base.Index.insert_blob">insert_blob</h2>

```python
Index.insert_blob(self, objects, num_threads=8)
```

insert the specified objects into the index and build the index.
Although this is the same as the fucntion, both implementations are different. 

    objects : Inserted objects.
    num_threads : Number of threads in building index.

