import dune.common

from ._istl import BCRSMatrix as BCRSMatrix11
from ._istl import BlockVector as BlockVector1
from ._istl import *

from dune.generator.generator import SimpleGenerator
from dune.common.hashit import hashIt
from dune.generator import Constructor
generator = SimpleGenerator("BCRSMatrix","Dune::Python")
generatorvec = SimpleGenerator("BlockVector","Dune::Python")
generatormatrixindexset = SimpleGenerator("MatrixIndexSet","Dune::Python")

def load(includes ,typeName ,constructors=None, methods=None):
    includes = includes + ["dune/python/istl/bcrsmatrix.hh"]
    typeHash = "istlbcrsmatrix_" + hashIt(typeName)
    return generator.load(includes ,typeName ,typeHash ,constructors ,methods)

def loadvec(includes ,typeName ,constructors=None, methods=None):
    includes = includes + ["dune/python/common/fvector.hh",
                           "dune/python/istl/bvector.hh"]
    typeHash = "istlbvector_" + hashIt(typeName)
    return generatorvec.load(includes ,typeName ,typeHash ,constructors ,methods)

def loadmatrixindexset(includes ,typeName ,constructors=None, methods=None):
    includes = includes + ["dune/python/istl/matrixindexset.hh"]
    typeHash = "matrixindexset_" + hashIt(typeName)
    return generatormatrixindexset.load(includes ,typeName ,typeHash ,constructors ,methods)

def MatrixIndexSet(rows, cols):
    includes = ["dune/istl/matrixindexset.hh"]
    typeName = "Dune::MatrixIndexSet"
    return loadmatrixindexset(includes, typeName).MatrixIndexSet(rows,cols)

def BCRSMatrix(blockType):
    includes = ["dune/istl/bcrsmatrix.hh"]
    try:
        typeName = "Dune::BCRSMatrix<" + blockType.cppTypeName + " >"
        return load(includes,typeName).BCRSMatrix
    except AttributeError:
        #check whether blocksize is 1,1
        if blockType[0] == blockType[1] == 1:
            return BCRSMatrix11

        typeName = "Dune::BCRSMatrix< Dune::FieldMatrix< double,"\
                + str(blockType[0]) + "," + str(blockType[1]) \
                + " > >"
        return load(includes, typeName).BCRSMatrix

def bcrsMatrix(size=0, *args, **kwargs):
    blockType = kwargs.get("blockType",[1,1])
    if size != 0 :
        return BCRSMatrix(blockType)(size,*args)
    return BCRSMatrix(blockType)()


def BlockVector(blockSize):
    if blockSize == 1:
        return BlockVector1
    typeName = "Dune::BlockVector< Dune::FieldVector< double ," + str(blockSize) + " > >"
    includes = ["dune/python/istl/bvector.hh"]
    # todo: provide other constructors
    return loadvec(includes, typeName).BlockVector

def blockVector(size, blockSize=1):
    if blockSize == 1:
        return BlockVector1(size)
    typeName = "Dune::BlockVector< Dune::FieldVector< double ," + str(blockSize) + " > >"
    includes = ["dune/istl/bvector.hh"]
    # todo: provide other constructors
    return loadvec(includes, typeName).BlockVector(size)

def matrixAdapter(mat, u, v = None):
    matType = mat.cppTypeName
    domType = u.cppTypeName
    if v is None:
        ranType = domType
    else:
        ranType = v.cppTypeName
    typeName = "Dune::MatrixAdapter<"+matType+","+domType+","+ranType+">"
    includes = {"dune/istl/operators.hh", *mat.cppIncludes, *u.cppIncludes}
    if v is not None:
        includes += v.cppIncludes
    ctor = Constructor(["const "+matType+"& mat"], ["return new DuneType(mat);"], ['"matrix"_a'])
    gen = SimpleGenerator("LinearOperator", "Dune::Python")
    typeHash = "matrixadapter_"+hashIt(typeName)
    return gen.load(list(includes), typeName, typeHash, ctor, holder="std::shared_ptr").LinearOperator(mat)

def getSolverFromFactory(op, config, includes = {}):
    operatorType = op.cppTypeName
    typeName = "Dune::SolverFactory< " + operatorType + " >"
    includes = ["dune/python/istl/solverfactory.hh", *includes]
    typeHash = "solverfactory_"+hashIt([typeName, *includes]) # add includes to the hash to enforce rebuild if custom includes are added
    gen = SimpleGenerator("SolverFactory", "Dune::Python")
    return gen.load(includes, typeName, typeHash).SolverFactory.get(op, config)
