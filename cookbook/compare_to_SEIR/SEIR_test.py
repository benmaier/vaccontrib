import numpy as np

pops = 2

I = np.eye(pops)
O = np.zeros_like(I)
latency = np.diag(np.random.rand(pops))
rec = np.diag(np.random.rand(pops))
#latency = 1

inf_rates = np.random.rand(pops,pops)


print("----SIR------")
vals, vecs = np.linalg.eig(-inf_rates.dot(np.linalg.inv(-rec)))
ndx = np.argmax(np.real(vals))
R = np.real(vals[ndx])
pop = np.real(vecs[:,ndx])

print(R, pop/pop.sum())


print("------SEIR-----")

SEIR_rec = np.block([[-latency, O],
                     [+latency,-rec],
                     ])
SEIR_inf = np.block([[O,inf_rates],
                     [O,O]
                     ])

print(-SEIR_inf.dot(np.linalg.inv(SEIR_rec)))

vals, vecs = np.linalg.eig(-SEIR_inf.dot(np.linalg.inv(SEIR_rec)))
ndx = np.argmax(np.real(vals))
R = np.real(vals[ndx])
pop = np.real(vecs[:,ndx])

print(R, pop/pop.sum())
print(R, (pop[:pops]+pop[pops:])/(pop[:pops]+pop[pops:]).sum())
print(R, (pop[:pops])/(pop[:pops]).sum())
print(R, (pop[pops:])/(pop[pops:]).sum())


print(vals)
print(vecs)
