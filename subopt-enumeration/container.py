'''
A quick and dirty container class to work as either a heap or a list
Use according to whether you care about storing structures sorted by energy,
or don't care about the order.

Note, there is no random access!  (because it's not need right now, and heaps are not RA) 
For random access, use the getList() function to get access to the store
'''
import heapq

class Container(object):
    def __init__(self, t):
        self.store = []
        self.index = 0
        
        if t == 'list':
            self.add = self.addToList
            self._next = self.nextFromList
        elif t == 'heap':
            self.add = self.addToHeap
            self._next = self.nextFromHeap
        else:
            raise Error('Unsupported type!')

    def __iter__(self):
        return self

    def addToList(self, elem):
        self.store.append(elem)

    def addToHeap(self, elem):
        heapq.heappush(self.store, elem)

    def nextFromList(self):
        if self.index < len(self.store):
            val = self.store[self.index]
            self.index += 1
            return val
        else:
            self.index = 0
            raise StopIteration

    def nextFromHeap(self):
        if self.index < len(self.store):
            return heapq.heappop(self.store)
        else:
            raise StopIteration

    def next(self):
        return self._next()

    def getList(self):
        return self.store


'''
# Some test code
container = Container('heap')
container.add(5)
container.add(101)
container.add(15)
for e in container:
    print e
'''
