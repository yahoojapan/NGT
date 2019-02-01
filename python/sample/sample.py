import csv
import ngtpy

# create an index framwork in filesystem.
ngtpy.create(path=b'index', dimension=128, distance_type="L2")

# load objects.
objects = []
with open(b'../../data/sift-dataset-5k.tsv', 'r') as fp:
    for object in csv.reader(fp, delimiter = '\t'):
        objects.append(object[0:128])

# open index.
index = ngtpy.Index(b'index')

# insert the objects.
index.batch_insert(objects)

# save the index.
index.save()

# close the index.
index.close()

# open the index.
index = ngtpy.Index(b'index')

# load query data.
with open(b'../../data/sift-query-3.tsv', 'r') as fp:
    query = list(csv.reader(fp, delimiter = '\t'))

# search for the index with the first query.
results = index.search(query[0], size=5)

print('ID\tDistance');
for result in results:
    print('{}\t{}'.format(*result))
print('# of distance computations=' + str(index.get_num_of_distance_computations()))

# get an object in the index.
object = index.get_object(4078)

# search with the object in the index.
results = index.search(object, size=5)

print('\nID\tDistance');
for result in results:
    print('{}\t{}'.format(*result))
print('# of distance computations=' + str(index.get_num_of_distance_computations()))

# insert the same objects individually. not build the index for them.
with open(b'../../data/sift-dataset-5k.tsv', 'r') as fp:
    for object in csv.reader(fp, delimiter = '\t'):
        objectID = index.insert(object[0:128])
        if objectID % 1000 == 0:
            print('Processed {} objects.'.format(objectID))

# build the index for the inserted objects to search.
index.build_index()

# search with the first query to confirm the insertion.
results = index.search(query[0], size=6)

# get the search results.
print('\nID\tDistance');
for result in results:
    print('{}\t{}'.format(*result))
print('# of distance computations=' + str(index.get_num_of_distance_computations()))

# remove an object.
index.remove(3030)

# search with the first query to confirm the removal.
results = index.search(query[0], size=6)
print('\nID\tDistance');
for result in results:
    print('{}\t{}'.format(*result))
print('# of distance computations=' + str(index.get_num_of_distance_computations()))

index.save()
index.close()
