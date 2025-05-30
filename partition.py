import math
from random import random, randint, choice
from functools import reduce

class Partition:

    def __init__(self, B=[[]], name=None):
        self.blocks = set([])
        self.elements = set([])
        for block in B:
            self.elements |= set(block)
            self.blocks |= set([frozenset(block)])
            if name: self.name = name

    def __str__(self):
        s = ' '.join([str([str(e) for e in block])+',' for block in self.blocks])
        s = s[0:len(s)-1]
        s = s.replace('\'','')
        s = s.replace('[', '{')
        return '{ ' + s.replace(']', '}') + ' }'

    def __len__(self):
        return len(self.blocks)



    def __add__(self,other):
        '''
            Partition addition algorithm as in [H&S]
        '''
        def compute_block_of(s):
            b0 = set([])
            b1 = self.block_containing(s) | other.block_containing(s)
            while b0 != b1:
                b0 = b1.copy()
                for block in self.blocks:
                    if len(block & b1) > 0: b1 |= block
                for block in other.blocks:
                    if len(block & b1) > 0: b1 |= block
            return b1

        new_blocks = []
        elements_done = set([])
        for e in self.elements:
            if not e in elements_done:
                b = compute_block_of(e)
                elements_done |= b
                new_blocks += [b]
        return Partition(new_blocks)

        



    def __mul__(self, other):
        '''
            Partition multiplication
        '''
        new_blocks = []
        elements_done = set([])
        for e in self.elements:
            if e not in elements_done:
                b1 = self.block_containing(e)
                b2 = other.block_containing(e)
                new_block = b1 & b2
                elements_done |= new_block
            new_blocks += [new_block]
        return Partition(new_blocks)

    def __truediv__(self, other):
        # quotient partition goes here
        return False

    def __eq__(self, other):
        if len(self.blocks) != len(other.blocks): return False
        if any([e not in other.elements for e in self.elements]): return False
        for block in self.blocks:
            corresponding = False
            for element in block:
                if corresponding is False:
                    corresponding = other.block_containing(element)
                if element not in corresponding: return False
        return True

    def __ne__(self, other):
        return not self == other

    def __ge__(self, other):
        # may be faster ways
        return self * other == other

    def __le__(self, other):
        # may be faster ways
        return self * other == self

    def __lt__(self, other):
        return self <= other and self != other

    def __gt__(self, other):
        return self >= other and self != other

    def __invert__(self):
        # identity if zero
        # else zero
        if len(self.blocks) == len(self.elements):
            return Partition([self.elements])
        return Partition([[e] for e in self.elements])

    def block_containing(self, element):
        for block in self.blocks:
            if element in block:
                return set(block)

    def order_blocks(self, order):
        def minIndex(B):
            return min([order.index(e) for e in B])
        O = list(self.blocks)
        O.sort(key=minIndex)
        return O


    def to_RGS(self, order):
        '''
        encode as a string mapping states, as ordered in the given list, to blocks, represented as space-separated numbers
        '''
        B = self.order_blocks(order)
        RGS = ""
        for e in order:
            if e not in self.elements: return False
            RGS += str(B.index(self.block_containing(e))) + " "
        return RGS[:-1]


    @staticmethod
    # return minimal partition containing the block
    def min_join(block, elements):
        E = set(elements) - set(block)
        E = [[e] for e in E]
        return Partition(E + [block])

    @staticmethod
    # return maximal partition containing the block
    def max_split(block, elements):
        return Partition([block, set(elements) - set(block)])

    @staticmethod
    # sum over list of partitions
    def sum(P):
        return reduce(lambda a,b: a + b, P)

    @staticmethod
    # product over list of partitions
    def prod(P):
        return reduce(lambda a,b: a * b, P)

    @staticmethod
    # lattice identity element, ie, set containing all elements in one set
    def identity(elements):
        return Partition([elements])

    @staticmethod
    # lattice zero element, ie, set of singletons containing each element
    def zero(elements):
        return Partition([[e] for e in elements])

    @staticmethod
    def from_RGS(rgs, order):
        '''
        convert back from RGS to Partition
        TODO move helpers out of function body if more efficient
        '''
        blocknumbers = [int(n) for n in rgs.split(" ")]
        blocks = [[] for i in range(max(blocknumbers)+1)]
        for i in range(len(order)): blocks[blocknumbers[i]] += [order[i]]

        return Partition(blocks)



    @staticmethod
    def random(elements):
        # https://pythonhosted.org/combalg-py/_modules/combalg/combalg.html

        elements = list(elements)
        
        def binomial_coefficient(n,k):
            # take care of nonsense
            if k > n or k < 0:
                return 0
            nk = 1    # (n)_k : falling factorial
            kf = 1    # k!    : k-factorial
            for i in range(0,k):
                nk *= n-i
                kf *= k-i
            return nk//kf

        g_bell = {}
        def bell_n(n):
            if n == 0:
                return 1
            if n not in g_bell:
                g_bell[n] = sum([binomial_coefficient(n-1,k) * bell_n(k) for k in range(n)])
            return g_bell[n]

        def prob(n,k):
                return float(binomial_coefficient(n-1,k-1))*bell_n(n-k)/bell_n(n)
        def label_partition(a, q):
            nc = 1 + max(q)
            eq_classes = []
            for t in range(nc):
                eq_classes.append([])
            for t in range(len(q)):
                eq_classes[q[t]].append(a[t])
            return eq_classes

        def random_permutation(elements):
            a = list(elements)
            n = len(a)
            for m in range(n-1):
                l = m + int(random() * (n - m))
                tmp = a[l]
                a[l] = a[m]
                a[m] = tmp
            return a

        n = len(elements)
        q = [0]*n
        m = n
        l = 0
        while m > 0:
            total = 0
            rho = random()
            k = 1
            while k < m:
                total += prob(m,k)
                if total >= rho:
                    break
                k += 1
            for i in range(m-k,m):
                q[i] = l
            l += 1
            m -= k
        return Partition(label_partition(elements, random_permutation(q)))

    
    @staticmethod
    def list(elements, n = False):
        elements = list(elements)
        if(len(elements) > 10): print("I would reconsider that if I were you")

        # https://codereview.stackexchange.com/a/240277
        # https://stackoverflow.com/a/30134039

        def partitions_n(n): yield from partition_n(elements,n,n)

        def partition_n(elements, m, n):
            if len(elements) == 1:
              yield [elements]
              return
        
            first = elements[0]
            for smaller in partition_n(elements[1:], m - 1, n):
              if len(smaller) > n: continue
              if len(smaller) >= m:
                for i, subset in enumerate(smaller):
                  yield smaller[:i] + [[first] + subset]  + smaller[i+1:]
              if len(smaller) < n: yield [[first]] + smaller


        if not n:
            partitions = []
            for i in range(1, len(elements)+1):
                partitions += partitions_n(i)
            return [Partition(L) for L in partitions]

        return [Partition(L) for L in partitions_n(n)]
