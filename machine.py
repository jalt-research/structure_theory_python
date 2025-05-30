from partition import *
from itertools import product
import time
from matplotlib.pyplot import plot, show, legend

class Machine:

    def __init__(self, delta, omega=None):
        self.delta = delta
        self.omega = omega
        self.S = set(delta.keys())
        self.I = set([])
        self.O = set([])

        for row in delta.values():
            self.I |= set(row.keys())

        if omega:
            for row in omega.values():
                self.O |= set(row.values())

        self.transition_inputs = set(product(self.S, self.I))

        # canonical orderings for conversion to RGS
        self.S_order = list(self.S)
        self.I_order = list(self.I)


    def __str__(self):
        out = ""
        for s in self.S:
            out += "\n" + str(s) + ":"
            for i in self.I:
                out += "\n\t" + str(i) + " ---> " + str(self.delta[s][i])
        if(self.omega):
            out += "\n"
            for s in self.S:
                out += "\n" + str(s) + ":"
                for i in self.I:
                    out += "\n\t" + str(i) + " ---> " + str(self.omega[s][i])
        return out


    @staticmethod
    def random(states, inputs, outputs=[], kind='Mealy'):
        delta = {}
        omega = {}
        for s in states:
            delta[s] = {}
            omega[s] = {}
            if kind is 'Moore': out = choice(outputs)
            for x in inputs:
                delta[s][x] = choice(states)
                if outputs:
                    if kind is 'Moore': omega[s][x] = out
                    else: omega[s][x] = choice(outputs)

        if outputs: return Machine(delta, omega)
        return Machine(delta)


    @staticmethod
    def homomorphic(M, pi_S, pi_I=None, pi_O=None):
        S = set(pi_S.blocks)
        if pi_I: I = set(pi_I.blocks)
        if pi_O: O = set(pi_O.blocks)
        # ....

    def equivalent_states(self):
        # computes output equivalence of states
        if not self.omega: return self.pI()
        def k_refine(pi): 
            update = self.p0()
            for block in pi.blocks:
                B = list(block)
                for i in range(len(B)):
                    e1 = B[i]
                    for j in range(i): 
                        e2 = B[j]
                        equiv = True
                        for x in self.I:
                            if pi.block_containing(self.delta[e1][x]) != pi.block_containing(self.delta[e2][x]):
                                equiv = False
                        if equiv: update += Partition.min_join(set([e1,e2]), self.S)
            return pi * update
        pi = self.p0()

        # states are initially equivalent if they have the same output under all inputs
        for e1 in self.S:
            for e2 in self.S:
                equiv = True
                for x in self.I:
                    if self.omega[e1][x] != self.omega[e2][x]:
                        equiv = False
                if equiv: pi += Partition.min_join(set([e1,e2]), self.S)
        k = k_refine(pi)
        # then they are refined based on their state transitions
        while k != pi:
            pi = k
            k = k_refine(pi)
        return k

    def reduced(self):
        # construct reduced machine from equivalent states

        pi = self.equivalent_states()
        delta = {}
        omega = {}
        mapping = {}
        labels = []

        for block in pi.blocks:
            label = '{' + ','.join(block) + '}'
            labels += [label]
            omega[label] = {}
            delta[label] = {}
            for element in block:
                mapping[element] = label

        for B in pi.blocks:
            e = list(B)[0]
            label = mapping[e]
            for x in self.I:
                if self.omega: omega[label][x] = self.omega[e][x]
                delta[label][x] = mapping[self.delta[e][x]]
        
        if self.omega: return Machine(delta, omega)
        return Machine(delta)



    def p0(self, on='s'):
        if(on == 's' or on == 'S'):
            return Partition([[s] for s in self.S])
        if(on == 'i' or on == 'I'):
            return Partition([[i] for i in self.I])
        # if(on == 'o' or on == 'O'):
            # return Partition([[o] for o in self.O])


    def pI(self, on='s'):
        if(on == 's' or on == 'S'):
            return Partition([self.S])
        if(on == 'i' or on == 'I'):
            return Partition([self.I])
        # if(on == 'o' or on == 'O'):
            # return Partition([self.O])


    def delta_block(self, b, i):
        return frozenset([self.delta[s][i] for s in b])


    def delta_blocks(self, p):
        B = set([])
        for b in p.blocks:
            for i in self.I:
                B |= set([self.delta_block(b,i)])
        return B


    def m(self, p, pair='ss'):
        # the m function of Hartmanis' pair algebra
        B = p.blocks
        if(pair == 'ss'): # state-state partition pair
            delta_blocks = self.delta_blocks(p)
            merging_blocks = delta_blocks | set([])
            for s in self.S:
                merge_s = set([b for b in merging_blocks if s in b])
                merged_s = frozenset([s])
                for b in merge_s:
                    merged_s |= b
                merging_blocks = (merging_blocks - merge_s) | set([merged_s])

            return Partition(merging_blocks)


    def m_trajectory(self, p, pair='ss'):
        P = [p, self.m(p)]
        while P[-1] != P[-2] and P[-1] not in P[:-1]:
            P += [self.m(P[-1])]
        return P


    def A(self, p):
        return p * self.m(p)


    def A_trajectory(self, p, n):
        # change this to be like the m_trajectory function
        for i in range(n):
            p = m(p) * p
        return p


    def is_SP(self, p):
        # for all pi >= m(p), (pi, p) is a partition pair 
        # if m(p) <= p then (p, p) is such a pair and p has s.p.
        return self.m(p) <= p


    def min_containing(self, block):
            return Partition(self.p0().blocks - set([frozenset(s) for s in block]) | set([frozenset(block)]))


    def semigroup(self):
        table = self.delta
        columns = {}
        for length in range(1,len(self.S)):
            for string in (''.join(x) for x in product(self.I_order, repeat=length)):
                print(string, len(columns.keys()))
                col = list(self.S_order)
                for symbol in string:
                    for i in range(len(col)):
                        col[i] = table[col[i]][symbol]
                try:
                    columns[str(col)] |= set([string])
                except KeyError:
                    columns[str(col)] = set([string])
        return columns


    def enumerate_SP(self, flag=False):
        '''
            Hartmanis + Stearns' algorithm for enumerating the SP lattice.
            While in the worst case the SP lattice is the same size as the full lattice of partitions (ie, grows as the Bell numbers), the worst cases are rare
        '''
        # 1: find the finest SP partitions lumping each pair of states
        order = list(self.S)
        known = {}

        def sort_by(L1, O): return [b for a,b in sorted(zip([O.index(e) for e in L1], L1))]

        def iterate(p, block):
            if p == self.min_containing(block):
                for i in self.I:
                    block_i = list(self.delta_block(block,i))
                    block_i = sort_by(block_i, order)
                    if str(block_i) in known: return known[str(block_i)]
                    min_p = self.min_containing(block_i)
                    p += min_p
            else: 
                for b in p.blocks:
                    for i in self.I:
                        block_i = list(self.delta_block(b,i))
                        block_i = sort_by(block_i, order)
                        if str(block_i) in known: return known[str(block_i)]
                        min_p = self.min_containing(block_i)
                        p += min_p
            return p


        for s,t in product(self.S, self.S):
            if order.index(s) <= order.index(t): continue
            merged = self.min_containing([s,t])
            while not self.is_SP(merged):
                merged = iterate(merged,[s,t])
            known[str([s,t])] = merged

        minimals = set([v.to_RGS(self.S_order) for v in known.values()])


        # 2: iteratively find  sums of the minimal partitions

        old = set([])
        new = minimals.copy()
        novel = new - old
        while len(novel):
            last = new
            nxt = set([])
            # only add the new combinations each step
            for p1, p2 in set(product(last,novel)) | set(product(novel, novel)):
                (p1, p2) = (Partition.from_RGS(p,self.S_order) for p in (p1,p2))
                nxt |= set([(p1 + p2).to_RGS(self.S_order)])
            old |= new
            new = nxt
            novel = new - old

        return [Partition.from_RGS(p, self.S_order) for p in old | new] + [self.p0()]

