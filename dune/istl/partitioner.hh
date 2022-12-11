// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ISTL_PARTITIONER_HH
#define DUNE_ISTL_PARTITIONER_HH

#include <vector>

namespace Dune {

template<class IndexType, class WeightType>
class Graph
{
  using index_type = IndexType;
  using weight_type = WeightType;

public:
  /**
   * \param n  number of vertices.
   * \param m  number of edges.
   **/
  Graph (index_type n, index_type m, int nvwxt = 1, int nadjwgt = 0, int nvsize = 0)
    : xadj_(n+1)
    , adjncy_(2*m)
    , vwxt_(n*nvwxt, 1.0)
    , adjwgt_(2*m*nadjwgt, 1.0)
    , vsize_(n*nvsize, 1.0)
  {}

  index_type numVertices () const
  {
    return xadj_.size() - 1;
  }

  index_type numEdges () const
  {
    return  adjncy_.size() / 2;
  }

  int numVertexConstraints () const
  {
    return vwxt_.size() / numVertices();
  }

  int numEdgeConstraints () const
  {
    return adjwgt_.size() / (2*numEdges());
  }

private:
  std::vector<index_type> xadj_;
  std::vector<index_type> adjncy_;
  std::vector<weight_type> vwxt_;
  std::vector<weight_type> adjwgt_;
  std::vector<weight_type> vsize_;
};


template<class ImplType>
class GraphPartitioner
{
  using Impl = ImplType;

public:
  /**
   * \param nparts  The number of parts to partition the graph.
   **/
  std::vector<int> part (int nparts) const
  {
    return impl().part(nparts);
  }

private:
  Impl const& impl () const
  {
    return static_cast<Impl const&>(*this);
  }
};


template<class GraphType>
class SimpleGraphPartitioner
    : public GraphPartitioner<SimpleGraphPartitioner<GraphType>>
{
  using graph_type = GraphType;

public:
  /**
   * \param graph   The graph data structure.
   **/
  SimpleGraphPartitioner (graph_type const& graph)
    : graph_(graph)
  {}

public:
  /**
   * \param nparts  The number of parts to partition the graph.
   **/
  std::vector<int> part (int nparts) const
  {
    std::vector<int> p(graph_.numVertices(), 0);
    std::size_t blocksize = p.size() / nparts;
    int b = p.size() - blocksize * nparts;
    int a = nparts - b;

    for (int i = 0; i < a; ++i)
      for (std::size_t j = i*blocksize; j < (i+1)*blocksize; ++j)
        p[j] = i;

    std::size_t shift = a*blocksize;
    for (int i = 0; i < b; ++i)
      for (std::size_t j = i*(blocksize+1); j < (i+1)*(blocksize+1); ++j)
        p[shift+j] = a+i;

    return p;
  }

private:
  const graph_type& graph_;
};

#if HAVE_METIS

template<class GraphType>
class MetisGraphPartitioner
    : public GraphPartitioner<MetisGraphPartitioner<GraphType>>
{
  using graph_type = GraphType;

#if HAVE_PARMETIS && defined(IDXTYPEWIDTH)
  using idx_t = ::idx_t;
#elif HAVE_PARMETIS && defined(HAVE_SCOTCH_NUM_TYPE)
  using idx_t = SCOTCH_Num;
#elif HAVE_PARMETIS
  using idx_t = int;
#else
  using idx_t = std::size_t;
#endif

public:
  enum PartType
  {
    RECURSIVE = 1,
    KWAY = 2,
    UNKNOWN = 3,
  };

  struct MetisError
      : public Dune::Exception
  {
  public:
    MetisError (int errcode) noexcept
      : errcode_(errcode)
    {}

    const char* what() const noexcept override
    {
      switch (errcode_) {
        case METIS_ERROR_INPUT:
          return "Input error.";
        case METIS_ERROR_MEMORY:
          return "Could not allocate the required memory.";
        case METIS_ERROR:
          return "Some other type of error.";
        default:
          return "Unknown METIS error."
      }
    }

  public:
    int errcode_;
  };

  /**
   * \param graph   The graph data structure.
   **/
  MetisGraphPartitioner (graph_type const& graph, PartType partType)
    : graph_(graph)
    , partType_(partType)
  {}

public:
  /**
   * \param np  The number of parts to partition the graph.
   **/
  std::vector<int> part (int np) const
  {
    std::vector<int> p(graph_.numVertices(), 0);
    idx_t nvtxs = graph_.numVertices(), ncon = graph_.numVertexConstraints(), nparts = np;

    int ret = 0;
    switch (partType_) {
    case RECURSIVE:
      ret = METIS_PartGraphRecursive(&nvtxs, &ncon,
        graph_.xadj().data(), graph_.adjncy().data(), graph_.vwgt().data(), graph_.vsize().data(),
        graph_.adjwgt().data(), nparts, tpwgts, options_, &objval_, part.data());
      break;
    case KWAY:
      ret = METIS_PartGraphKway(&nvtxs, &ncon,
        graph_.xadj().data(), graph_.adjncy().data(), graph_.vwgt().data(), graph_.vsize().data(),
        graph_.adjwgt().data(), nparts, tpwgts, options_, &objval_, part.data());
      break;
    }

    if (ret != METIS_OK)
      throw new MetisError(ret);
    return p;
  }

  idx_t objval () const
  {
    return objval_;
  }

private:
  const graph_type& graph_;
  PartType partType_;

  mutable idx_t objval_;
};

#endif // HAVE_METIS

} // end namespace Dune

#endif // DUNE_ISTL_PARTITIONER_HH