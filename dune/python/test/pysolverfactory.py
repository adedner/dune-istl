import dune.istl
import numpy as np

def buildMat(N):
    A = dune.istl.BCRSMatrix11(N,N,3,1.5, dune.istl.BuildMode.implicit)
    #A.setImplicitBuildModeParameters(3,1)
    for i in range(N):
        if i>0:
            A[i,i-1] = -1
        A[i,i] = 2
        if i<N-1:
            A[i,i+1] = -1
    A.compress()
    return A

N = 1000
A = buildMat(N)
b = dune.istl.BlockVector1(np.ones([N,1]))

op = dune.istl.matrixAdapter(A, b)

prec_config = {"type": "ssor"}
defaults = {"reduction": 1e-8,
          "verbose": 1,
          "maxit": 1000,
          "preconditioner": prec_config}
types = ["cgsolver", "minressolver", "cholmod"]

for solver_type in types:
    print(solver_type)
    solver = dune.istl.getSolverFromFactory(op, {"type":solver_type, **defaults},
                                            {"dune/istl/solvers.hh", "dune/istl/cholmod.hh"})

    x = b.copy()
    res = solver(x,b.copy())
    print(res)
    r = b.copy()
    op.applyscaleadd(-1., x, r)
    print("recomputed residual:", r.two_norm)
