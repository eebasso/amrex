#ifndef AMREX_REALVECT_H_
#define AMREX_REALVECT_H_
#include <AMReX_Config.H>

#include <AMReX_Box.H>
#include <AMReX_REAL.H>
#include <AMReX_SPACE.H>
#include <AMReX_IntVect.H>
#include <AMReX_Utility.H>
#include <AMReX_Math.H>

#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <iosfwd>
#include <vector>

namespace amrex
{

/**
 * \brief A Real vector in SpaceDim-dimensional space
 *
 * The class RealVect is an implementation of a Real vector in a
 * SpaceDim-dimensional space.
 *
 * RealVect values are accessed using the operator[] function, as for a normal
 * C++ array.  In addition, the basic arithmetic operators have been overloaded
 * to implement scaling and translation operations.
 */
class RealVect
{
public:

  /**
     \name Constructors and Accessors
  */
  /*@{*/

  ///
  /**
     Construct a RealVect whose components are zero.
  */
  constexpr RealVect () noexcept {} // NOLINT

  explicit RealVect (const std::vector<Real>& vr) noexcept : vect{AMREX_D_DECL(vr[0],vr[1],vr[2])} {
      BL_ASSERT(vr.size() == AMREX_SPACEDIM);
  }

  ///
  /**
     Construct a RealVect given the specific values for its
     coordinates.  AMREX_D_DECL is a macro that sets the constructor to
     take AMREX_SPACEDIM arguments.
  */
#if (AMREX_SPACEDIM > 1)
  AMREX_GPU_HOST_DEVICE
  constexpr RealVect (AMREX_D_DECL(Real i, Real j, Real k)) noexcept : vect{AMREX_D_DECL(i,j,k)} {}
#endif

  AMREX_GPU_HOST_DEVICE
  explicit constexpr RealVect (Real i) noexcept : vect{AMREX_D_DECL(i,i,i)} {}

  ///
  /**
     Construct a RealVect setting the coordinates to the corresponding values in
     the Real array a.
  */
  AMREX_GPU_HOST_DEVICE
  explicit RealVect (const Real* a) noexcept : vect{AMREX_D_DECL(a[0],a[1],a[2])} {}

  ///
  /**
     Construct a RealVect from an IntVect by coercing each component
     from <tt>int</tt> to Real.
  */
  AMREX_GPU_HOST_DEVICE
  RealVect (const IntVect & iv) noexcept : vect{AMREX_D_DECL((Real)iv[0],(Real)iv[1],(Real)iv[2])} {}

  ///
  /**
     Prevent rvalue assignment to the <i>i</i>'th coordinate of the
     RealVect.
  */
  Real& operator[] (int i) && = delete;

  ///
  /**
     Returns a modifiable lvalue reference to the <i>i</i>'th coordinate of the
     RealVect.
  */
  AMREX_GPU_HOST_DEVICE
  inline
  Real& operator[] (int i) & noexcept
  { AMREX_ASSERT(i>=0 && i < SpaceDim); return vect[i]; }

  ///
  /**
     Returns the <i>i</i>'th coordinate of the RealVect.
  */
  AMREX_GPU_HOST_DEVICE
  inline
  const Real& operator[] (int i) const& noexcept
  { AMREX_ASSERT(i>=0 && i < SpaceDim); return vect[i]; }

  /*@}*/

  /**
     \name Iterators
  */
  /*@{*/

  ///
  /**
     Return pointer to the first component.
  */
  [[nodiscard]] AMREX_GPU_HOST_DEVICE
  inline
  Real* begin () noexcept { return vect; }

  ///
  /**
     Return pointer to the first component.
  */
  [[nodiscard]] AMREX_GPU_HOST_DEVICE
  inline
  const Real* begin () const noexcept { return vect; }

  ///
  /**
     Return pointer to the one past last component
  */
  [[nodiscard]] AMREX_GPU_HOST_DEVICE
  inline
  Real* end () noexcept { return vect + AMREX_SPACEDIM; }

  ///
  /**
     Return pointer to the one past last component
  */
  [[nodiscard]] AMREX_GPU_HOST_DEVICE
  inline
  const Real* end () const noexcept { return vect + AMREX_SPACEDIM; }

  /*@}*/

  /**
     \name Comparison Operators
  */
  /*@{*/

  ///
  /**
     Returns true if this RealVect is equivalent to argument RealVect.  All
     comparisons between analogous components must be satisfied.
  */
  AMREX_GPU_HOST_DEVICE
  bool operator== (const RealVect& p) const noexcept {
      return AMREX_D_TERM(vect[0] == p[0], && vect[1] == p[1], && vect[2] == p[2]);
  }

  ///
  /**
     Returns true if this RealVect is different from argument RealVect.
     All comparisons between analogous components must be satisfied.
  */
  AMREX_GPU_HOST_DEVICE
  bool operator!= (const RealVect& p) const noexcept {
      return AMREX_D_TERM(vect[0] != p[0], || vect[1] != p[1], || vect[2] != p[2]);
  }

  ///
  /**
     Returns true if this RealVect is less than argument RealVect.  All
     comparisons between analogous components must be satisfied.  Note
     that, since the comparison is component-wise, it is possible for
     an RealVect to be neither greater than, less than, nor equal to
     another.
  */
  AMREX_GPU_HOST_DEVICE inline
  bool operator< (const RealVect& p) const noexcept;

  ///
  /**
     Returns true if this RealVect is less than or equal to argument
     RealVect.  All comparisons between analogous components must be
     satisfied.  Note that, since the comparison is component-wise, it
     is possible for an RealVect to be neither greater than or equal
     to, less than or equal to, nor equal to another.
  */
  AMREX_GPU_HOST_DEVICE inline
  bool operator<= (const RealVect& p) const noexcept;

  ///
  /**
     Returns true if this RealVect is greater than argument RealVect.
     All comparisons between analogous components must be satisfied.
     Note that, since the comparison is component-wise, it is possible
     for an RealVect to be neither greater than, less than, nor equal
     to another.
  */
  AMREX_GPU_HOST_DEVICE inline
  bool operator> (const RealVect& p) const noexcept;

  ///
  /**
     Returns true if this RealVect is greater than or equal to argument
     RealVect.  All comparisons between analogous components must be
     satisfied.  Note that, since the comparison is component-wise, it
     is possible for an RealVect to be neither greater than or equal
     to, less than or equal to, nor equal to another.
  */
  AMREX_GPU_HOST_DEVICE inline
  bool operator>= (const RealVect& p) const noexcept;

  /*@}*/

  /**
     \name Arithmetic Operators
  */
  /*@{*/

  ///
  /**
     Modifies this RealVect by addition of a scalar to each component.
  */
  AMREX_GPU_HOST_DEVICE inline
  RealVect& operator+= (Real s) noexcept;

  ///
  /**
     Returns a RealVect that is this RealVect with a scalar s added
     to each component.
  */
  AMREX_GPU_HOST_DEVICE inline
  RealVect operator+ (Real s) const noexcept;

  ///
  /**
     Modifies this RealVect by component-wise addition by argument.
  */
  AMREX_GPU_HOST_DEVICE inline
  RealVect& operator+= (const RealVect& p) noexcept;

  ///
  /**
     Modifies this RealVect by subtraction of a scalar from each
     component.
  */
  AMREX_GPU_HOST_DEVICE inline
  RealVect& operator-= (Real s) noexcept;

  ///
  /**
     Modifies this RealVect by component-wise subtraction by argument.
  */
  AMREX_GPU_HOST_DEVICE inline
  RealVect& operator-= (const RealVect& p) noexcept;

  ///
  /**
     Returns a RealVect that is this RealVect with a scalar s subtracted
     from each component.
  */
  AMREX_GPU_HOST_DEVICE inline
  RealVect operator- (Real s) const noexcept;

  ///
  /**
     Modifies this RealVect by multiplying each component by a scalar.
  */
  AMREX_GPU_HOST_DEVICE inline
  RealVect& operator*= (Real s) noexcept;

  ///
  /**
     Returns the result of the scalar product with another RealVect
  */
  [[nodiscard]] AMREX_GPU_HOST_DEVICE inline
  Real dotProduct(const RealVect& a_rhs) const noexcept;

#if (AMREX_SPACEDIM == 3)
  ///
  /**
     Returns the result of the cross product with another RealVect
  */
  [[nodiscard]] AMREX_GPU_HOST_DEVICE inline
  RealVect crossProduct(const RealVect& a_rhs) const noexcept;
#endif

  ///
  /**
     Modifies this RealVect by component-wise multiplication by argument.
  */
  AMREX_GPU_HOST_DEVICE inline
  RealVect& operator*= (const RealVect& p) noexcept;

//XXX  ///
//XXX  /**
//XXX     Returns component-wise product of this RealVect and argument.
//XXX  */
//XXX  RealVect operator* (const RealVect& p) const;

  ///
  /**
     Returns a RealVect that is this RealVect with each component
     multiplied by a scalar.
  */
  AMREX_GPU_HOST_DEVICE inline
  RealVect operator* (Real s) const noexcept;

  ///
  /**
     Modifies this RealVect by dividing each component by a scalar.
  */
  AMREX_GPU_HOST_DEVICE inline
  RealVect& operator/= (Real s) noexcept;

  ///
  /**
     Modifies this RealVect by component-wise division by argument.
  */
  AMREX_GPU_HOST_DEVICE inline
  RealVect& operator/= (const RealVect& p) noexcept;

//XXX  ///
//XXX  /**
//XXX     Returns component-wise quotient of this RealVect by argument.
//XXX  */
//XXX  RealVect operator/ (const RealVect& p) const;

  ///
  /**
     Returns a RealVect that is this RealVect with each component
     divided by a scalar.
  */
  AMREX_GPU_HOST_DEVICE inline
  RealVect operator/ (Real s) const noexcept;

  ///
  /**
     Modifies this RealVect by multiplying each component by a scalar.
  */
  AMREX_GPU_HOST_DEVICE inline
  RealVect& scale (Real s) noexcept;

  ///
  /**
     Returns an IntVect whose components are the std::floor of the RealVect
     components.
  */
  [[nodiscard]] AMREX_GPU_HOST_DEVICE inline
  IntVect floor () const noexcept;

  ///
  /**
     Returns an IntVect whose components are the std::floor of the RealVect
     components.
  */
  [[nodiscard]] AMREX_GPU_HOST_DEVICE inline
  IntVect ceil () const noexcept;

  ///
  /**
     Returns an IntVect whose components are the std::floor of the RealVect
     components.
  */
  [[nodiscard]] AMREX_GPU_HOST_DEVICE inline
  IntVect round () const noexcept;

  /*@}*/

  /**
     \name Other arithmetic operators
  */
  /*@{*/

  ///
  /**
     Modifies this RealVect by taking component-wise min with RealVect
     argument.
  */
  AMREX_GPU_HOST_DEVICE inline
  RealVect& min (const RealVect& p) noexcept;

  ///
  /**
     Returns the RealVect that is the component-wise minimum of two
     argument RealVects.
  */
  friend AMREX_GPU_HOST_DEVICE inline
  RealVect min (const RealVect& p1, const RealVect& p2) noexcept;

  ///
  /**
     Modifies this RealVect by taking component-wise max with RealVect
     argument.
  */
  AMREX_GPU_HOST_DEVICE inline
  RealVect& max (const RealVect& p) noexcept;

  ///
  /**
     Returns the RealVect that is the component-wise maximum of two
     argument RealVects.
  */
  friend AMREX_GPU_HOST_DEVICE inline
  RealVect max (const RealVect& p1, const RealVect& p2) noexcept;

  /*@}*/

  /**
     \name Unary operators
  */
  /*@{*/

  ///
  /**
     Unary plus -- for completeness.
  */
  AMREX_GPU_HOST_DEVICE inline
  RealVect operator+ () const noexcept;

  ///
  /**
     Unary minus -- negates all components of this RealVect.
  */
  AMREX_GPU_HOST_DEVICE inline
  RealVect operator- () const noexcept;

  ///
  /**
     Sum of all components of this RealVect.
  */
  [[nodiscard]] AMREX_GPU_HOST_DEVICE inline
  Real sum () const noexcept;

  ///
  /**
     sqrt(sum squares)
  */
  [[nodiscard]] AMREX_GPU_HOST_DEVICE inline
  Real vectorLength() const noexcept;

  ///
  /**
     sum squares--no square root
  */
  [[nodiscard]] AMREX_GPU_HOST_DEVICE inline
  Real radSquared() const noexcept;

  ///
  /**
     Product of all components of this RealVect.
  */
  [[nodiscard]] AMREX_GPU_HOST_DEVICE inline
  Real product () const noexcept;

  ///
  /**
     Component with the minimum value of this RealVect (returns 0 if they are all the same).
     a_doAbs : if true then take the absolute value before comparing
  */
  [[nodiscard]] AMREX_GPU_HOST_DEVICE inline
  int minDir(const bool& a_doAbs) const noexcept;

  ///
  /**
     Component with the maximum value of this RealVect (returns 0 if they are all the same).
     a_doAbs : if true then take the absolute value before comparing
  */
  [[nodiscard]] AMREX_GPU_HOST_DEVICE inline
  int maxDir(const bool& a_doAbs) const noexcept;

  /*@}*/

  /**
     \name Data pointer functions
  */
  /*@{*/

  ///
  /**
     Only for sending stuff to Fortran
   */
  [[nodiscard]] AMREX_GPU_HOST_DEVICE
  const Real* dataPtr() const noexcept { return vect; }

  ///
  /**
     Only for sending stuff to Fortran
   */
  AMREX_GPU_HOST_DEVICE
  Real* dataPtr() noexcept { return vect; }

  /*@}*/

  /**
     \name Constants
  */
  /*@{*/

  ///
  /**
     Returns a basis vector in the given coordinate direction.<br>
     In 2-D:<br>
     BASISREALV(0) == (1.,0.);
     BASISREALV(1) == (0.,1.).<br>
     In 3-D:<br>
     BASISREALV(0) == (1.,0.,0.);
     BASISREALV(1) == (0.,1.,0.);
     BASISREALV(2) == (0.,0.,1.).<br>
     Note that the coordinate directions are based at zero.
  */
  friend AMREX_GPU_HOST_DEVICE
  RealVect BASISREALV(int dir) noexcept;

  ///
  /**
     This is a RealVect all of whose components are equal to zero.
  */
  static AMREX_EXPORT const RealVect Zero;

  ///
  /**
     This is a RealVect all of whose components are equal to one.
  */
  static AMREX_EXPORT const RealVect Unit;

  /*@}*/

  /**
     \name Static public members functions
  */
  /*@{*/

  ///
  /**
     The zero vector in AMREX_SPACEDIM-dimensional space.
  */
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  static constexpr RealVect TheZeroVector () noexcept {
      return RealVect(0.);
  }

  /**
     The unit vector in AMREX_SPACEDIM-dimensional space.
  */
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  static constexpr RealVect TheUnitVector () noexcept {
      return RealVect(1.);
  }

  /*@}*/

  /**
     \name Arithmetic friend functions
  */
  /*@{*/

  ///
  /**
     Returns a RealVect that is a RealVect <i>p</i> with
     a scalar <i>s</i> added to each component.
  */
  friend AMREX_GPU_HOST_DEVICE
  RealVect operator+ (Real s, const RealVect& p) noexcept;

  ///
  /**
     Returns <i>s - p</i>.
  */
  friend AMREX_GPU_HOST_DEVICE
  RealVect operator- (Real s, const RealVect& p) noexcept;

  ///
  /**
     Returns a RealVect that is a RealVect <i>p</i> with each component
     multiplied by a scalar <i>s</i>.
  */
  friend AMREX_GPU_HOST_DEVICE
  RealVect operator* (Real s, const RealVect& p) noexcept;
  ///
  /**
     Returns a RealVect that is a RealVect <i>p</i> with each component
     divided by a scalar <i>s</i>.
  */
  friend AMREX_GPU_HOST_DEVICE
  RealVect operator/ (Real s, const RealVect& p) noexcept;

  ///
  /**
     Returns component-wise sum of RealVects <i>s</i> and <i>p</i>.
  */
  friend AMREX_GPU_HOST_DEVICE
  RealVect operator+ (const RealVect& s, const RealVect& p) noexcept;

  ///
  /**
     Returns <i>s - p</i>.
  */
  friend AMREX_GPU_HOST_DEVICE
  RealVect operator- (const RealVect& s, const RealVect& p) noexcept;

  ///
  /**
     Returns component-wise product of <i>s</i> and <i>p</i>.
  */
  friend AMREX_GPU_HOST_DEVICE
  RealVect operator* (const RealVect& s, const RealVect& p) noexcept;
  ///
  /**
     Returns component-wise quotient <i>p / s</i>.
  */
  friend AMREX_GPU_HOST_DEVICE
  RealVect operator/ (const RealVect& s, const RealVect& p) noexcept;

  ///
  /**
     Returns a RealVect obtained by multiplying each of the components
     of the given RealVect by a scalar.
  */
  friend inline AMREX_GPU_HOST_DEVICE
  RealVect scale (const RealVect& p, Real s) noexcept;

  /*@}*/

  ///
  /**
     Print to the given output stream in ASCII.
  */
  friend std::ostream& operator<< (std::ostream& ostr, const RealVect& p);

  friend std::istream& operator>> (std::istream& is, RealVect& iv);

protected:
    /**
        The individual components of this RealVectND.
    */
    Real vect[dim] = {};
};

// Template deduction guide for RealVectND
template <class...Args,
    std::enable_if_t<
        IsConvertible_v<Real, Args...>,
    int> = 0>
AMREX_GPU_HOST_DEVICE // __device__ for HIP
RealVectND(Real, Real, Args...) -> RealVectND<sizeof...(Args)+2>;

// Template deduction guide for RealVectND
template<int dim>
AMREX_GPU_HOST_DEVICE // __device__ for HIP
RealVectND(const IntVectND<dim>&) -> RealVectND<dim>;

// Template deduction guide for RealVectND
template<int dim>
AMREX_GPU_HOST_DEVICE // __device__ for HIP
RealVectND(const GpuArray<Real, dim>&) -> RealVectND<dim>;

template <int dim>
inline constexpr const RealVectND<dim> RealVectND<dim>::Zero{Real(0)};

template <int dim>
inline constexpr const RealVectND<dim> RealVectND<dim>::Unit{Real(1)};

using RealVect = RealVectND<AMREX_SPACEDIM>;

template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim>&
RealVectND<dim>::operator-= (Real s) noexcept {
    for (int i=0; i<dim; ++i) {
        vect[i] -= s;
    }
    return *this;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim>&
RealVectND<dim>::operator*= (Real s) noexcept {
    for (int i=0; i<dim; ++i) {
        vect[i] *= s;
    }
    return *this;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim>&
RealVectND<dim>::operator-= (const RealVectND<dim>& p) noexcept {
    for (int i=0; i<dim; ++i) {
        vect[i] -= p[i];
    }
    return *this;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim>
RealVectND<dim>::operator+ () const noexcept {
    return RealVectND<dim>(*this);
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim>
RealVectND<dim>::operator- () const noexcept {
    RealVectND<dim> retval;
    for (int i=0; i<dim; ++i) {
        retval[i] = -vect[i];
    }
    return retval;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim>&
RealVectND<dim>::scale (Real s) noexcept {
    for (int i=0; i<dim; ++i) {
        vect[i] *= s;
    }
    return *this;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
IntVectND<dim>
RealVectND<dim>::floor () const noexcept {
    IntVectND<dim> retval;
    for (int i=0; i<dim; ++i) {
        retval[i] = static_cast<int>(std::floor(vect[i]));
    }
    return retval;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
IntVectND<dim>
RealVectND<dim>::ceil () const noexcept {
    IntVectND<dim> retval;
    for (int i=0; i<dim; ++i) {
        retval[i] = static_cast<int>(std::ceil(vect[i]));
    }
    return retval;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
IntVectND<dim>
RealVectND<dim>::round () const noexcept {
    IntVectND<dim> retval;
    for (int i=0; i<dim; ++i) {
        retval[i] = static_cast<int>(std::round(vect[i]));
    }
    return retval;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
Real
RealVectND<dim>::sum () const noexcept {
    Real retval = vect[0];
    for (int i=1; i<dim; ++i) {
        retval += vect[i];
    }
    return retval;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
Real
RealVectND<dim>::vectorLength () const noexcept {
    Real len = this->radSquared();
    len = std::sqrt(len);
    return len;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
Real
RealVectND<dim>::radSquared () const noexcept {
    Real retval = vect[0] * vect[0];
    for (int i=1; i<dim; ++i) {
        retval += vect[i] * vect[i];
    }
    return retval;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
Real
RealVectND<dim>::product () const noexcept {
    Real retval = vect[0];
    for (int i=1; i<dim; ++i) {
        retval *= vect[i];
    }
    return retval;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
bool
RealVectND<dim>::operator< (const RealVectND<dim>& p) const noexcept {
    bool retval = vect[0] < p[0];
    for (int i=1; i<dim; ++i) {
        retval = retval && vect[i] < p[i];
    }
    return retval;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
bool
RealVectND<dim>::operator<= (const RealVectND<dim>& p) const noexcept {
    bool retval = vect[0] <= p[0];
    for (int i=1; i<dim; ++i) {
        retval = retval && vect[i] <= p[i];
    }
    return retval;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
bool
RealVectND<dim>::operator> (const RealVectND<dim>& p) const noexcept {
    bool retval = vect[0] > p[0];
    for (int i=1; i<dim; ++i) {
        retval = retval && vect[i] > p[i];
    }
    return retval;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
bool
RealVectND<dim>::operator>= (const RealVectND<dim>& p) const noexcept {
    bool retval = vect[0] >= p[0];
    for (int i=1; i<dim; ++i) {
        retval = retval && vect[i] >= p[i];
    }
    return retval;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim>&
RealVectND<dim>::min (const RealVectND<dim>& p) noexcept {
    for (int i=0; i<dim; ++i) {
        vect[i] = amrex::min(vect[i], p.vect[i]);
    }
    return *this;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim>&
RealVectND<dim>::max (const RealVectND<dim>& p) noexcept {
    for (int i=0; i<dim; ++i) {
        vect[i] = amrex::max(vect[i], p.vect[i]);
    }
    return *this;
}

/**
    Returns the RealVectND that is the component-wise minimum of two
    argument RealVectNDs.
*/
template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim> min (const RealVectND<dim>& p1, const RealVectND<dim>& p2) noexcept {
    RealVectND<dim> p(p1);
    return p.min(p2);
}

/**
    Returns the RealVectND that is the component-wise maximum of two
    argument RealVectNDs.
*/
template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim> max (const RealVectND<dim>& p1, const RealVectND<dim>& p2) noexcept {
    RealVectND<dim> p(p1);
    return p.max(p2);
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
Real
RealVectND<dim>::dotProduct (const RealVectND<dim>& a_rhs) const noexcept {
    Real retval  = vect[0] * a_rhs.vect[0];
    for (int i=1; i<dim; ++i) {
        retval += vect[i] * a_rhs.vect[i];
    }
    return retval;
}

template <int dim>
template <int N, std::enable_if_t<( N==3 ), int>>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim>
RealVectND<dim>::crossProduct (const RealVectND<dim>& a_rhs) const noexcept {
    RealVectND<dim> tmp(vect[1] * a_rhs[2] - vect[2] * a_rhs[1],
                        vect[2] * a_rhs[0] - vect[0] * a_rhs[2],
                        vect[0] * a_rhs[1] - vect[1] * a_rhs[0]);
    return tmp;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim>&
RealVectND<dim>::operator+= (Real s) noexcept {
    for (int i=0; i<dim; ++i) {
        vect[i] += s;
    }
    return *this;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim>&
RealVectND<dim>::operator+= (const RealVectND<dim>& p) noexcept {
    for (int i=0; i<dim; ++i) {
        vect[i] += p[i];
    }
    return *this;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim>&
RealVectND<dim>::operator*= (const RealVectND<dim>& p) noexcept {
    for (int i=0; i<dim; ++i) {
        vect[i] *= p[i];
    }
    return *this;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim>
RealVectND<dim>::operator* (Real s) const noexcept {
    RealVectND<dim> retval;
    for (int i=0; i<dim; ++i) {
        retval[i] = vect[i] * s;
    }
    return retval;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim>
RealVectND<dim>::operator- (Real s) const noexcept {
    RealVectND<dim> retval;
    for (int i=0; i<dim; ++i) {
        retval[i] = vect[i] - s;
    }
    return retval;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim>
RealVectND<dim>::operator+ (Real s) const noexcept {
    RealVectND<dim> retval;
    for (int i=0; i<dim; ++i) {
        retval[i] = vect[i] + s;
    }
    return retval;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim>&
RealVectND<dim>::operator/= (Real s) noexcept {
    for (int i=0; i<dim; ++i) {
        vect[i] /= s;
    }
    return *this;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim>&
RealVectND<dim>::operator/= (const RealVectND<dim>& p) noexcept {
    for (int i=0; i<dim; ++i) {
        vect[i] /= p[i];
    }
    return *this;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim>
RealVectND<dim>::operator/ (Real s) const noexcept {
    RealVectND<dim> retval;
    for (int i=0; i<dim; ++i) {
        retval[i] = vect[i] / s;
    }
    return retval;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
int
RealVectND<dim>::minDir (const bool& a_doAbs) const noexcept {
    int mDir = 0;
    for (int idir = 0; idir < dim; idir++) {
        if (a_doAbs) {
            if (std::abs(vect[idir]) < std::abs(vect[mDir])) {
                mDir = idir;
            }
        } else {
            if (vect[idir] < vect[mDir]) {
                mDir = idir;
            }
        }
    }
    return mDir;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
int
RealVectND<dim>::maxDir (const bool& a_doAbs) const noexcept {
    int mDir = 0;
    for (int idir = 0; idir < dim; idir++) {
        if (a_doAbs) {
            if (std::abs(vect[idir]) > std::abs(vect[mDir])) {
                mDir = idir;
            }
        } else {
            if (vect[idir] > vect[mDir]) {
                mDir = idir;
            }
        }
    }
    return mDir;
}

/**
    Returns a basis vector in the given coordinate direction.<br>
    In 2-D:<br>
    BASISREALV(0) == (1.,0.);
    BASISREALV(1) == (0.,1.).<br>
    In 3-D:<br>
    BASISREALV(0) == (1.,0.,0.);
    BASISREALV(1) == (0.,1.,0.);
    BASISREALV(2) == (0.,0.,1.).<br>
    Note that the coordinate directions are based at zero.
*/
template <int dim=AMREX_SPACEDIM>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim> BASISREALV (int dir) noexcept {
    AMREX_ASSERT(dir >= 0 && dir < dim);
    RealVectND<dim> tmp{0.};
    tmp[dir] = 1;
    return tmp;
}

/**
    \name Arithmetic functions
*/
/*@{*/

/**
     Returns a RealVectND that is a RealVectND <i>p</i> with each component
    divided by a scalar <i>s</i>.
*/
template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim> operator/ (Real s, const RealVectND<dim>& p) noexcept {
    RealVectND<dim> retval;
    for (int i=0; i<dim; ++i) {
        retval[i] = s / p[i];
    }
    return retval;
}

/**
    Returns a RealVectND that is a RealVectND <i>p</i> with
    a scalar <i>s</i> added to each component.
*/
template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim> operator+ (Real s, const RealVectND<dim>& p) noexcept {
    RealVectND<dim> retval;
    for (int i=0; i<dim; ++i) {
        retval[i] = p[i] + s;
    }
    return retval;
}

/**
     Returns <i>s - p</i>.
*/
template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim> operator- (Real s, const RealVectND<dim>& p) noexcept {
    RealVectND<dim> retval;
    for (int i=0; i<dim; ++i) {
        retval[i] = s - p[i];
    }
    return retval;
}

/**
     Returns a RealVectND that is a RealVectND <i>p</i> with each component
    multiplied by a scalar <i>s</i>.
*/
template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim> operator* (Real s, const RealVectND<dim>& p) noexcept {
    RealVectND<dim> retval;
    for (int i=0; i<dim; ++i) {
        retval[i] = s * p[i];
    }
    return retval;
}

/**
     Returns component-wise quotient <i>p / s</i>.
*/
template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim> operator/ (const RealVectND<dim>& s, const RealVectND<dim>& p) noexcept {
    RealVectND<dim> retval;
    for (int i=0; i<dim; ++i) {
        retval[i] = s[i] / p[i];
    }
    return retval;
}

/**
     Returns component-wise sum of RealVectNDs <i>s</i> and <i>p</i>.
*/
template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim> operator+ (const RealVectND<dim>& s, const RealVectND<dim>& p) noexcept {
    RealVectND<dim> retval;
    for (int i=0; i<dim; ++i) {
        retval[i] = p[i] + s[i];
    }
    return retval;
}

/**
     Returns <i>s - p</i>.
*/
template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim> operator- (const RealVectND<dim>& s, const RealVectND<dim>& p) noexcept {
    RealVectND<dim> retval;
    for (int i=0; i<dim; ++i) {
        retval[i] = s[i] - p[i];
    }
    return retval;
}

/**
     Returns component-wise product of <i>s</i> and <i>p</i>.
*/
template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim> operator* (const RealVectND<dim>& s, const RealVectND<dim>& p) noexcept {
    RealVectND<dim> retval;
    for (int i=0; i<dim; ++i) {
        retval[i] = p[i] * s[i];
    }
    return retval;
}

/**
    Returns a RealVectND obtained by multiplying each of the components
    of the given RealVectND by a scalar.
*/
template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim> scale (const RealVectND<dim>& p, Real s) noexcept {
    RealVectND<dim> retval;
    for (int i=0; i<dim; ++i) {
        retval[i] = s * p[i];
    }
    return retval;
}

/// \cond DOXYGEN_IGNORE
namespace detail {
    std::ostream& real_vector_write (std::ostream& os, const Real* p, int dim);
    std::istream& real_vector_read (std::istream& is, Real* p, int dim);

    template<int dim>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr
    void RealVectCat_imp (Real*& dst, const RealVectND<dim>& src) noexcept {
        for (int i=0; i<dim; ++i) {
            dst[i] = src[i];
        }
        dst += dim;
    }

    template<int dim>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr
    void RealVectSplit_imp2 (RealVectND<dim>& dst, const Real*& src) noexcept {
        for (int i=0; i<dim; ++i) {
            dst[i] = src[i];
        }
        src += dim;
    }

    template<class T, std::size_t...Ns>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr
    T RealVectSplit_imp (T& retval, std::index_sequence<Ns...>, const Real * src) noexcept {
        (RealVectSplit_imp2(amrex::get<Ns>(retval), src), ...);
        return retval;
    }
}
/// \endcond

template<int dim>
std::ostream&
operator<< (std::ostream& os, const RealVectND<dim>& p)
{
    return detail::real_vector_write(os, p.begin(), dim);
}

template<int dim>
std::istream&
operator>> (std::istream& is, RealVectND<dim>& p)
{
    return detail::real_vector_read(is, p.begin(), dim);
}

/**
* \brief Returns a RealVectND obtained by concatenating the input RealVectNDs.
* The dimension of the return value equals the sum of the dimensions of the inputted RealVectNDs.
*/
template<int d, int...dims>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
constexpr RealVectND<detail::get_sum<d, dims...>()>
RealVectCat (const RealVectND<d>& v, const RealVectND<dims>&...vects) noexcept {
    RealVectND<detail::get_sum<d, dims...>()> retval {0.};
    Real* dst = retval.begin();
    detail::RealVectCat_imp(dst, v);
    (detail::RealVectCat_imp(dst, vects), ...);
    return retval;
}

/**
* \brief Returns a tuple of RealVectND obtained by splitting the input RealVectND
* according to the dimensions specified by the template arguments.
*/
template<int d, int...dims>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
constexpr GpuTuple<RealVectND<d>, RealVectND<dims>...>
RealVectSplit (const RealVectND<detail::get_sum<d, dims...>()>& v) noexcept {
    GpuTuple<RealVectND<d>, RealVectND<dims>...> retval;
    return detail::RealVectSplit_imp(retval,
                                    std::make_index_sequence<1 + sizeof...(dims)>(),
                                    v.begin());
}

/**
* \brief Returns a new RealVectND of size new_dim and
* assigns the first new_dim values of iv to it.
*/
template<int new_dim, int old_dim>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
constexpr RealVectND<new_dim>
RealVectShrink (const RealVectND<old_dim>& iv) noexcept {
    return iv.template shrink<new_dim>();
}

/**
* \brief Returns a new RealVectND of size new_dim and
* assigns all values of iv to it and fill_extra to the remaining elements.
*/
template<int new_dim, int old_dim>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
constexpr RealVectND<new_dim>
RealVectExpand (const RealVectND<old_dim>& iv, Real fill_extra=0) noexcept {
    return iv.template expand<new_dim>(fill_extra);
}

/**
* \brief Returns a new RealVectND of size new_dim
* by either shrinking or expanding iv.
*/
template<int new_dim, int old_dim>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
constexpr RealVectND<new_dim>
RealVectResize (const RealVectND<old_dim>& iv, Real fill_extra=0) noexcept {
    return iv.template resize<new_dim>(fill_extra);
}

/*@}*/

}  // namespace amrex

// Spcialize std::tuple_size for RealVectND. Used by structured bindings.
template<int dim>
struct std::tuple_size<amrex::RealVectND<dim>> {
    static constexpr std::size_t value = dim;
};

// Spcialize std::tuple_element for RealVectND. Used by structured bindings.
template<std::size_t s, int dim>
struct std::tuple_element<s, amrex::RealVectND<dim>> {
    using type = amrex::Real;
};

#endif
