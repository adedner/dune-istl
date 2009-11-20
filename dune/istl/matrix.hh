// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MATRIX_HH
#define DUNE_MATRIX_HH

/** \file
    \brief A dynamic dense block matrix class
 */

#include <vector>
#include <dune/istl/bvector.hh>

namespace Dune {

  /** \addtogroup ISTL_SPMV
     \{
   */
  /** \brief A generic dynamic matrix

      This matrix is currently implemented as a BlockVector of BlockVectors.
      That makes the code fairly simple, as we get all iterators for free.
      However, it is not the best way as far as efficiency is concerned.
   */
  template<class T, class A=ISTLAllocator>
  class Matrix
  {
  public:

    /** \brief Export the type representing the underlying field */
    typedef typename T::field_type field_type;

    /** \brief Export the type representing the components */
    typedef T block_type;

    /** \brief Export the allocator */
    typedef A allocator_type;

    /** \brief The type implementing a matrix row */
    typedef BlockVector<T> row_type;

    /** \brief Type for indices and sizes */
    typedef typename A::size_type size_type;

    /** \brief Iterator over the matrix rows */
    typedef typename BlockVector<row_type>::iterator RowIterator;

    /** \brief Iterator for the entries of each row */
    typedef typename row_type::iterator ColIterator;

    /** \brief Const iterator over the matrix rows */
    typedef typename BlockVector<row_type>::const_iterator ConstRowIterator;

    /** \brief Const iterator for the entries of each row */
    typedef typename row_type::const_iterator ConstColIterator;

    enum {
      //! The number of nesting levels the matrix contains.
      blocklevel = T::blocklevel+1
    };

    /** \brief Create empty matrix */
    Matrix() : data_(0), cols_(0)
    {}

    /** \brief Create uninitialized matrix of size rows x cols
     */
    Matrix(size_type rows, size_type cols) : data_(rows), cols_(cols)
    {
      for (size_type i=0; i<rows; i++)
        data_[i].resize(cols);
    }

    /** \brief Change the matrix size
     *
     * The way the data is handled is unpredictable.
     */
    void setSize(size_type rows, size_type cols) {
      data_.resize(rows);
      for (size_type i=0; i<rows; i++)
        data_[i].resize(cols);
      cols_ = cols;
    }

    /** \brief Get iterator to first row */
    RowIterator begin()
    {
      return data_.begin();
    }

    /** \brief Get iterator to one beyond last row */
    RowIterator end()
    {
      return data_.end();
    }

    /** \brief Get iterator to last row */
    RowIterator rbegin()
    {
      return data_.rbegin();
    }

    /** \brief Get iterator to one before first row */
    RowIterator rend()
    {
      return data_.rend();
    }

    /** \brief Get const iterator to first row */
    ConstRowIterator begin() const
    {
      return data_.begin();
    }

    /** \brief Get const iterator to one beyond last row */
    ConstRowIterator end() const
    {
      return data_.end();
    }

    /** \brief Get const iterator to last row */
    ConstRowIterator rbegin() const
    {
      return data_.rbegin();
    }

    /** \brief Get const iterator to one before first row */
    ConstRowIterator rend() const
    {
      return data_.rend();
    }

    /** \brief Assignment from scalar */
    Matrix& operator= (const field_type& t)
    {
      typedef typename BlockVector<row_type,allocator_type>::size_type
      size_type;

      for (size_type i=0; i<data_.size(); i++)
        data_[i] = t;
      return *this;
    }

    /** \brief The index operator */
    row_type& operator[](size_type row) {
#ifdef DUNE_ISTL_WITH_CHECKING
      if (row<0)
        DUNE_THROW(ISTLError, "Can't access negative rows!");
      if (row>=N())
        DUNE_THROW(ISTLError, "Row index out of range!");
#endif
      return data_[row];
    }

    /** \brief The const index operator */
    const row_type& operator[](size_type row) const {
#ifdef DUNE_ISTL_WITH_CHECKING
      if (row<0)
        DUNE_THROW(ISTLError, "Can't access negative rows!");
      if (row>=N())
        DUNE_THROW(ISTLError, "Row index out of range!");
#endif
      return data_[row];
    }

    /** \brief Return the number of rows */
    size_type N() const {
      return data_.size();
    }

    /** \brief Return the number of columns */
    size_type M() const {
      return cols_;
    }

    /** \brief The number of scalar rows */
    size_type rowdim() const {
#ifdef DUNE_ISTL_WITH_CHECKING
      if (M()==0)
        DUNE_THROW(ISTLError, "Can't compute rowdim() when there are no columns!");
#endif
      size_type dim = 0;
      for (size_type i=0; i<data_.size(); i++)
        dim += data_[i][0].rowdim();

      return dim;
    }

    /** \brief The number of scalar columns */
    size_type coldim() const {
#ifdef DUNE_ISTL_WITH_CHECKING
      if (N()==0)
        DUNE_THROW(ISTLError, "Can't compute coldim() when there are no rows!");
#endif
      size_type dim = 0;
      for (size_type i=0; i<data_[0].size(); i++)
        dim += data_[0][i].coldim();

      return dim;
    }

    /** \brief The number of scalar rows */
    size_type rowdim(size_type r) const {
#ifdef DUNE_ISTL_WITH_CHECKING
      if (r<0 || r>=N())
        DUNE_THROW(ISTLError, "Rowdim for nonexisting row " << r << " requested!");
      if (M()==0)
        DUNE_THROW(ISTLError, "Can't compute rowdim() when there are no columns!");
#endif
      return data_[r][0].rowdim();
    }

    /** \brief The number of scalar columns */
    size_type coldim(size_type c) const {
#ifdef DUNE_ISTL_WITH_CHECKING
      if (c<0 || c>=M())
        DUNE_THROW(ISTLError, "Coldim for nonexisting column " << c << " requested!");
      if (N()==0)
        DUNE_THROW(ISTLError, "Can't compute coldim() when there are no rows!");
#endif
      return data_[0][c].coldim();
    }

    /** \brief Multiplication with a scalar */
    Matrix<T>& operator*=(const field_type& scalar) {
      for (size_type row=0; row<data_.size(); row++)
        for (size_type col=0; col<cols_; col++)
          (*this)[row][col] *= scalar;

      return (*this);
    }

    /** \brief Multiplication with a scalar */
    Matrix<T>& operator/=(const field_type& scalar) {
      for (size_type row=0; row<data_.size(); row++)
        for (size_type col=0; col<cols_; col++)
          (*this)[row][col] /= scalar;

      return (*this);
    }

    /*! \brief Add the entries of another matrix to this one.
     *
     * \param b The matrix to add to this one. Its size has to
     * be the same as the size of this matrix.
     */
    Matrix& operator+= (const Matrix& b) {
#ifdef DUNE_ISTL_WITH_CHECKING
      if(N()!=b.N() || M() != b.M())
        DUNE_THROW(RangeError, "Matrix sizes do not match!");
#endif
      for (size_type row=0; row<data_.size(); row++)
        (*this)[row] += b[row];

      return (*this);
    }

    /*! \brief Subtract the entries of another matrix from this one.
     *
     * \param b The matrix to add to this one. Its size has to
     * be the same as the size of this matrix.
     */
    Matrix& operator-= (const Matrix& b) {
#ifdef DUNE_ISTL_WITH_CHECKING
      if(N()!=b.N() || M() != b.M())
        DUNE_THROW(RangeError, "Matrix sizes do not match!");
#endif
      for (size_type row=0; row<data_.size(); row++)
        (*this)[row] -= b[row];

      return (*this);
    }

    /** \brief Return the transpose of the matrix */
    Matrix transpose() const {
      Matrix out(N(), M());
      for (size_type i=0; i<M(); i++)
        for (size_type j=0; j<N(); j++)
          out[j][i] = (*this)[i][j];

      return out;
    }

    //! Multiplication of the transposed matrix times a vector
    template <class X, class Y>
    Y transposedMult(const X& vec) {
#ifdef DUNE_ISTL_WITH_CHECKING
      if (N()!=vec.size())
        DUNE_THROW(ISTLError, "Vector size doesn't match the number of matrix rows!");
#endif
      Y out(M());
      out = 0;

      for (size_type i=0; i<out.size(); i++ ) {
        for ( size_type j=0; j<vec.size(); j++ )
          out[i] += (*this)[j][i]*vec[j];
      }

      return out;
    }

    /// Generic matrix multiplication.
    friend Matrix<T> operator*(const Matrix<T>& m1, const Matrix<T>& m2) {
      Matrix<T> out(m1.N(), m2.M());
      out.clear();

      for (size_type i=0; i<out.N(); i++ ) {
        for ( size_type j=0; j<out.M(); j++ )
          for (size_type k=0; k<m1.M(); k++)
            out[i][j] += m1[i][k]*m2[k][j];
      }

      return out;
    }

    /// Generic matrix-vector multiplication.
    template <class X, class Y>
    friend Y operator*(const Matrix<T>& m, const X& vec) {
#ifdef DUNE_ISTL_WITH_CHECKING
      if (M()!=vec.size())
        DUNE_THROW(ISTLError, "Vector size doesn't match the number of matrix columns!");
#endif
      Y out(m.N());
      out = 0;

      for (size_type i=0; i<out.size(); i++ ) {
        for ( size_type j=0; j<vec.size(); j++ )
          out[i] += m[i][j]*vec[j];
      }

      return out;
    }

    //! y = A x
    template <class X, class Y>
    void mv(const X& x, Y& y) const
    {
#ifdef DUNE_ISTL_WITH_CHECKING
      if (x.N()!=M()) DUNE_THROW(ISTLError,"vector/matrix size mismatch!");
      if (y.N()!=N()) DUNE_THROW(ISTLError,"vector/matrix size mismatch!");
#endif

      for (size_type i=0; i<data_.size(); i++) {
        y[i]=0;
        for (size_type j=0; j<cols_; j++)
          (*this)[i][j].umv(x[j], y[i]);

      }

    }

    //! y = A^T x
    template<class X, class Y>
    void mtv (const X& x, Y& y) const
    {
#ifdef DUNE_ISTL_WITH_CHECKING
      if (x.N()!=N()) DUNE_THROW(ISTLError,"index out of range");
      if (y.N()!=M()) DUNE_THROW(ISTLError,"index out of range");
#endif
      for(int i=0; i<y.N(); ++i)
        y[i]=0;
      umtv(x,y);
    }

    //! y += A x
    template <class X, class Y>
    void umv(const X& x, Y& y) const
    {
#ifdef DUNE_ISTL_WITH_CHECKING
      if (x.N()!=M()) DUNE_THROW(ISTLError,"vector/matrix size mismatch!");
      if (y.N()!=N()) DUNE_THROW(ISTLError,"vector/matrix size mismatch!");
#endif

      for (size_type i=0; i<data_.size(); i++) {

        for (size_type j=0; j<cols_; j++)
          (*this)[i][j].umv(x[j], y[i]);

      }

    }

    //! y -= A x
    template<class X, class Y>
    void mmv (const X& x, Y& y) const
    {
#ifdef DUNE_ISTL_WITH_CHECKING
      if (x.N()!=M()) DUNE_THROW(ISTLError,"vector/matrix size mismatch!");
      if (y.N()!=N()) DUNE_THROW(ISTLError,"vector/matrix size mismatch!");
#endif
      ConstRowIterator endi=end();
      for (ConstRowIterator i=begin(); i!=endi; ++i)
      {
        ConstColIterator endj = (*i).end();
        for (ConstColIterator j=(*i).begin(); j!=endj; ++j)
          (*j).mmv(x[j.index()],y[i.index()]);
      }
    }

    /** \brief \f$ y += \alpha A x \f$ */
    template <class X, class Y>
    void usmv(const field_type& alpha, const X& x, Y& y) const
    {
#ifdef DUNE_ISTL_WITH_CHECKING
      if (x.N()!=M()) DUNE_THROW(ISTLError,"vector/matrix size mismatch!");
      if (y.N()!=N()) DUNE_THROW(ISTLError,"vector/matrix size mismatch!");
#endif

      for (size_type i=0; i<data_.size(); i++) {

        for (size_type j=0; j<cols_; j++)
          (*this)[i][j].usmv(alpha, x[j], y[i]);

      }

    }

    //! y += A^T x
    template<class X, class Y>
    void umtv (const X& x, Y& y) const
    {
#ifdef DUNE_ISTL_WITH_CHECKING
      if (x.N()!=N()) DUNE_THROW(ISTLError,"vector/matrix size mismatch!");
      if (y.N()!=M()) DUNE_THROW(ISTLError,"vector/matrix size mismatch!");
#endif
      ConstRowIterator endi=end();
      for (ConstRowIterator i=begin(); i!=endi; ++i)
      {
        ConstColIterator endj = (*i).end();
        for (ConstColIterator j=(*i).begin(); j!=endj; ++j)
          (*j).umtv(x[i.index()],y[j.index()]);
      }
    }

    //! y -= A^T x
    template<class X, class Y>
    void mmtv (const X& x, Y& y) const
    {
#ifdef DUNE_ISTL_WITH_CHECKING
      if (x.N()!=N()) DUNE_THROW(ISTLError,"vector/matrix size mismatch!");
      if (y.N()!=M()) DUNE_THROW(ISTLError,"vector/matrix size mismatch!");
#endif
      ConstRowIterator endi=end();
      for (ConstRowIterator i=begin(); i!=endi; ++i)
      {
        ConstColIterator endj = (*i).end();
        for (ConstColIterator j=(*i).begin(); j!=endj; ++j)
          (*j).mmtv(x[i.index()],y[j.index()]);
      }
    }

    //! y += alpha A^T x
    template<class X, class Y>
    void usmtv (const field_type& alpha, const X& x, Y& y) const
    {
#ifdef DUNE_ISTL_WITH_CHECKING
      if (x.N()!=N()) DUNE_THROW(ISTLError,"vector/matrix size mismatch!");
      if (y.N()!=M()) DUNE_THROW(ISTLError,"vector/matrix size mismatch!");
#endif
      ConstRowIterator endi=end();
      for (ConstRowIterator i=begin(); i!=endi; ++i)
      {
        ConstColIterator endj = (*i).end();
        for (ConstColIterator j=(*i).begin(); j!=endj; ++j)
          (*j).usmtv(alpha,x[i.index()],y[j.index()]);
      }
    }

    //! y += A^H x
    template<class X, class Y>
    void umhv (const X& x, Y& y) const
    {
#ifdef DUNE_ISTL_WITH_CHECKING
      if (x.N()!=N()) DUNE_THROW(ISTLError,"vector/matrix size mismatch!");
      if (y.N()!=M()) DUNE_THROW(ISTLError,"vector/matrix size mismatch!");
#endif
      ConstRowIterator endi=end();
      for (ConstRowIterator i=begin(); i!=endi; ++i)
      {
        ConstColIterator endj = (*i).end();
        for (ConstColIterator j=(*i).begin(); j!=endj; ++j)
          (*j).umhv(x[i.index()],y[j.index()]);
      }
    }

    //! y -= A^H x
    template<class X, class Y>
    void mmhv (const X& x, Y& y) const
    {
#ifdef DUNE_ISTL_WITH_CHECKING
      if (x.N()!=N()) DUNE_THROW(ISTLError,"vector/matrix size mismatch!");
      if (y.N()!=M()) DUNE_THROW(ISTLError,"vector/matrix size mismatch!");
#endif
      ConstRowIterator endi=end();
      for (ConstRowIterator i=begin(); i!=endi; ++i)
      {
        ConstColIterator endj = (*i).end();
        for (ConstColIterator j=(*i).begin(); j!=endj; ++j)
          (*j).mmhv(x[i.index()],y[j.index()]);
      }
    }

    //! y += alpha A^H x
    template<class X, class Y>
    void usmhv (const field_type& alpha, const X& x, Y& y) const
    {
#ifdef DUNE_ISTL_WITH_CHECKING
      if (x.N()!=N()) DUNE_THROW(ISTLError,"vector/matrix size mismatch!");
      if (y.N()!=M()) DUNE_THROW(ISTLError,"vector/matrix size mismatch!");
#endif
      ConstRowIterator endi=end();
      for (ConstRowIterator i=begin(); i!=endi; ++i)
      {
        ConstColIterator endj = (*i).end();
        for (ConstColIterator j=(*i).begin(); j!=endj; ++j)
          (*j).usmhv(alpha,x[i.index()],y[j.index()]);
      }
    }

    //===== norms

    //! frobenius norm: sqrt(sum over squared values of entries)
    double frobenius_norm () const
    {
      return std::sqrt(frobenius_norm2());
    }

    //! square of frobenius norm, need for block recursion
    double frobenius_norm2 () const
    {
      double sum=0;
      for (size_type i=0; i<N(); ++i)
        for (size_type j=0; j<M(); ++j)
          sum += data_[i][j].frobenius_norm2();
      return sum;
    }

    //! infinity norm (row sum norm, how to generalize for blocks?)
    double infinity_norm () const
    {
      double max=0;
      for (size_type i=0; i<N(); ++i) {
        double sum=0;
        for (size_type j=0; j<M(); j++)
          sum += data_[i][j].infinity_norm();
        max = std::max(max,sum);
      }
      return max;
    }

    //! simplified infinity norm (uses Manhattan norm for complex values)
    double infinity_norm_real () const
    {
      double max=0;
      for (size_type i=0; i<N(); ++i) {
        double sum=0;
        for (size_type j=0; j<M(); j++)
          sum += data_[i][j].infinity_norm_real();
        max = std::max(max,sum);
      }
      return max;
    }

    //===== query

    //! return true if (i,j) is in pattern
    bool exists (size_type i, size_type j) const
    {
#ifdef DUNE_ISTL_WITH_CHECKING
      if (i<0 || i>=N()) DUNE_THROW(ISTLError,"row index out of range");
      if (j<0 || i>=M()) DUNE_THROW(ISTLError,"column index out of range");
#endif
      return true;
    }

  protected:

    BlockVector<row_type, allocator_type> data_;

    size_type cols_;
  };
  /** \} */
} // end namespace Dune

#endif