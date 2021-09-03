# -*- coding: utf-8 -*-
"""
Simulations.
"""

import warnings
from itertools import groupby

import numpy as np

# Try to import the original SamplableSet,
# but if that doesn't work, use the mock version
# that's implemented in this package
try:
    from SamplableSet import SamplableSet
except ModuleNotFoundError as e: # pragma: no cover
    warnings.warn("Couldn't find the efficient implementation of `SamplableSet` (see github.com/gstonge/SamplableSet). Proceeding with less efficient implementation.")
    from vaccontrib.mock_samplable_set import MockSamplableSet as SamplableSet

from vaccontrib.mock_samplable_set import choice as _choice
from collections import Counter

from vaccontrib.main import get_2d_contribution_matrix


class LinearSystem():

    def __init__(self, reproduction_rate_matrix, decay_rate_vector, initial_conditions=20):

        self.set_rate_matrices(reproduction_rate_matrix, decay_rate_vector)
        self.set_initial_conditions(initial_conditions)

    def set_rate_matrices(self, reproduction_rate_matrix, decay_rate_vector):

        shape = reproduction_rate_matrix.shape
        assert(shape[0] == shape[1])
        assert(shape[0] == decay_rate_vector.shape[0])
        assert(not np.any(decay_rate_vector<=0))
        assert(np.all(reproduction_rate_matrix>=0))
        assert(np.any(reproduction_rate_matrix>0))

        nnz = reproduction_rate_matrix.nonzero()
        min_weight = np.min([np.amin(reproduction_rate_matrix[nnz]), np.amin(decay_rate_vector)])
        max_weight = np.max([np.amax(reproduction_rate_matrix), np.amax(decay_rate_vector)])

        self.A = np.array(reproduction_rate_matrix).astype(np.float64)
        self.B = np.array(decay_rate_vector).astype(np.float64)
        self.N = len(decay_rate_vector)

        self.K = A / B[None,:]
        self.C, self.y = get_2d_contribution_matrix(self.K,return_eigenvector_too=True)

        self.total_rates = reproduction_rate_matrix.sum(axis=0).flatten() + decay_rate_vector
        self.state_rates = np.vstack((self.A,self.B)).T

        print(self.state_rates)
        assert(np.all(self.state_rates.sum(axis=1)==self.total_rates))

        self.state_events = []
        for i in range(self.N):
            min_weight = self.state_rates[i][self.state_rates[i]>0].min()
            max_weight = self.state_rates[i].max()
            S = SamplableSet(min_weight, max_weight,cpp_type='int')
            for j in self.state_rates[i].nonzero()[0]:
                S[j] = self.state_rates[i][j]
            self.state_events.append(S)

        for ievent, events in enumerate(self.state_events):
            print(ievent)
            print([item for item in events])


    def set_initial_conditions(self, initial_conditions=20):

        min_weight = self.total_rates.min()
        max_weight = self.total_rates.max()

        if not hasattr(initial_conditions,'__len__'):
            self.y0 = (self.y * 20).astype(np.int32)
        else:
            self.y0 = np.array(initial_conditions).astype(np.int32)

        self.S = SamplableSet(min_weight=min_weight,max_weight=max_weight)
        i = 0

        self.states = {}

        for state, count in enumerate(self.y0):
            for j in range(int(count)):
                self.S[i] = self.total_rates[state]
                self.states[i] = state
                i += 1

        self.leftover_indices = []
        self.max_index = i - 1


    def simulate(self,t_start_measuring,t_stop_measuring,verbose=False):
        t = 0.

        active_nodes = {}


        counters = []
        ts = [t]
        total = np.sum(self.y0)
        ys = [total]

        new_tmax = 1.

        simulation_ended = False
        end_initialized = False

        while not simulation_ended:
            Lambda = self.S.total_weight()
            tau = np.random.exponential(scale=1/Lambda)
            t += tau
            node, _ = self.S.sample()
            state = self.states[node]

            event, _ = self.state_events[state].sample()

            # last event means decay
            if event == self.N:
                try:
                    counter = active_nodes.pop(node)
                    counters.append(( state, counter ))
                except KeyError as e:
                    pass
                del self.S[node]
                self.leftover_indices.append(node)
                self.states.pop(node)
                total -= 1
                #if node == self.max_index:
                #    self.max_index -= 1
            #birth
            else:
                total += 1
                birth_state = event
                if len(self.leftover_indices) > 0:
                    new_index = self.leftover_indices.pop()
                else:
                    #if len(self.S) == 0:
                    #print(node, state, event, self.S)
                    new_index = self.max_index + 1
                    self.max_index = new_index
                    #try:
                    #    new_index = max(self.S)[0] + 1
                    #except ValueError as e:
                    #    print(self.S)
                    #    print(len(self.S))
                    #    raise e

                if new_index in self.S:
                    raise ValueError("new_index is in S already")
                if not end_initialized:
                    self.S[new_index] = self.total_rates[birth_state]
                    self.states[new_index] = birth_state

                if t > t_start_measuring:
                    if t < t_stop_measuring:
                        active_nodes[new_index] = Counter()

                    try:
                        active_nodes[node][birth_state] += 1
                    except KeyError as e:
                        pass

            if t > t_stop_measuring and not end_initialized:
                end_initialized = True

            simulation_ended = t > t_stop_measuring and len(active_nodes) == 0
            simulation_ended = simulation_ended or len(self.S) == 0

            ts.append(t)
            ys.append(total)

            if verbose and t > new_tmax:
                #print(t)
                print(t, node, state, total, "len(active_nodes) =", len(active_nodes))
                new_tmax = t + 1.


        return ts, ys, counters


def get_mean_contribution_matrix_from_simulation(N,counters):
    C = np.zeros((N,N)).astype(np.float64)

    total_offspring = 0
    for state, counter in counters:
        for reproduced_state, count in counter.items():
            C[reproduced_state, state] += count
            total_offspring += count
    R = total_offspring / len(counters)
    C /= total_offspring
    C *= R

    return C

def get_mean_next_generation_matrix_from_simulation(N,counters):
    K = [ [ [] for i in range(N) ] for j in range(N) ]

    total_offspring = 0
    for state, counter in counters:
        #for reproduced_state, count in counter.items():
        for reproduced_state in range(N):
            K[reproduced_state][state].append(counter[reproduced_state])

    _K = [ [ np.mean(K[j][i]) for i in range(N) ] for j in range(N) ]
    _K = np.array(_K)

    return _K

def get_mean_eigenstate_from_simulation(N, counters):
    y = np.zeros((N,)).astype(np.float64)

    total_offspring = 0
    for state, counter in counters:
        y[state] += 1

    y /= sum(y)

    return y



if __name__=="__main__":


    A = np.arange(1,5).reshape(2,2)
    B = np.arange(1,3)
    L = LinearSystem(A, B)
    print(L.A, L.K, L.B)

    ts, ys, counters = L.simulate(1.,2.,verbose=True)


    K = A / B[None,:]
    C, y = get_2d_contribution_matrix(K,return_eigenvector_too=True)
    print("K_theory =")
    print(K)
    print("C_theory =")
    print(C)
    print("y_theory =")
    print(y)
    print("C_measured =")
    print(get_mean_contribution_matrix_from_simulation(K.shape[0],counters))
    print("y_measured =")
    print(get_mean_eigenstate_from_simulation(K.shape[0],counters))
    print("K_measured =")
    K_measured = get_mean_next_generation_matrix_from_simulation(K.shape[0],counters)
    print(K_measured)

    print()
    K = K_measured
    C, y = get_2d_contribution_matrix(K,return_eigenvector_too=True)
    print("K_false =")
    print(K)
    print("C_false =")
    print(C)
    print("y_false =")
    print(y)

    import matplotlib.pyplot as pl

    pl.plot(ts, ys)
    pl.show()
