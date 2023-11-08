#ifndef AMREX_DIM3_H_
#define AMREX_DIM3_H_
#include <AMReX_Config.H>

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IntVect.H>

#include <iostream>

namespace amrex {

struct Dim3 { int x; int y; int z; };
struct XDim3 { Real x; Real y; Real z; };

template <typename T,
          typename std::enable_if <std::is_same<T,  Dim3>::value ||
                                   std::is_same<T, XDim3>::value>::type* = nullptr>
std::ostream& operator<< (std::ostream& os, const T& d)
{
    os << '(' << d.x << ',' << d.y << ',' << d.z << ')';
    return os;
}

// template <class T>
// struct Dim3T
// {
//     AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
//     Dim3T () noexcept = default;

//     AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
//     Dim3T (T a_x, T a_y, T a_z) noexcept { define(a_x,a_y,a_z); }

//     AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
//     void define (T a_x, T a_y, T a_z) noexcept { x=a_x; y=a_y; z=a_z; }

//     // AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
//     // T& operator[] (int d) noexcept { return d==0 ? x : d==1 ? y : z; }

//     // AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
//     // const T& operator[] (int d) const noexcept { return d==0 ? x : d==1 ? y : z; }

//     AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
//     bool operator== (T val) const noexcept { return x == val && y == val && z == val; }

//     AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
//     bool operator!= (T val) const noexcept { return x != val || y != val || z != val; }

//     AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
//     bool operator== (const Dim3T<T>& rhs) const noexcept { return x == rhs.x && y == rhs.y && z == rhs.z; }

//     AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
//     bool operator!= (const Dim3T<T>& rhs) const noexcept { return x != rhs.x || y != rhs.y || z != rhs.z; }

//     //! Unary plus -- for completeness.
//     AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
//     Dim3T<T> operator+ () const noexcept { return *this; }
//     //! Unary Minus -- negates all components.
//     AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
//     Dim3T<T> operator- () const noexcept { return {-x, -y, -z}; }
//     //! Modifies Dim3 by addition of a scalar to each component.
//     AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
//     Dim3T<T>& operator+= (T s) noexcept
// {
//         x += s; y += s, z += s; return *this;
//     }
//     //! Modifies Dim3 by component-wise addition with argument.
//     AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
//     Dim3T<T>& operator+= (const Dim3T<T>& rhs) noexcept
//     {
//         x += rhs.x; y += rhs.y; z += rhs.z; return *this;
//     }
//     //! Modifies IntVect by multiplication of a scalar to each component.
//     AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
//     Dim3T<T>& operator*= (T s) noexcept
//     {
//         x *= s; y *= s; z *= s; return *this;
//     }
//     //! Modifies IntVect by component-wise multiplication with argument.
//     AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
//     Dim3T<T>& operator*= (const Dim3T<T>& rhs) noexcept
//     {
//         x *= rhs.x, y *= rhs.y, z *= rhs.z; return *this;
//     }
//     //! Modifies IntVect by division by a scalar to each component.
//     AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
//     Dim3T<T>& operator/= (T s) noexcept
//     {
//         x /= s, y /= s, z /= s; return *this;
//     }
//     //! Modifies IntVect by component-wise division with argument.
//     AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
//     Dim3T<T>& operator/= (const Dim3T<T>& rhs) noexcept
//     {
//         x /= rhs.x, y /= rhs.y, z /= rhs.z; return *this;
//     }
//     //! Modifies IntVect by subtraction of a scalar to each component.
//     AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
//     Dim3T<T>& operator-= (T s) noexcept
//     {
//         x -= s, y -= s, z -= s; return *this;
//     }
//     //! Modifies IntVect by component-wise subtraction with argument.
//     AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
//     Dim3T<T>& operator-= (const Dim3T<T>& rhs) noexcept
//     {
//         x -= rhs.x, y -= rhs.y, z -= rhs.z; return *this;
//     }
//     //! Returns component-wise sum of IntVect and argument.
//     AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
//     Dim3T<T> operator+ (const Dim3T<T>& rhs) const noexcept
//     {
//         return Dim3T<T>(x + rhs.x, y + rhs.y, z + rhs.z);
//     }
//     //! Return an IntVect that is this IntVect + s.
//     AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
//     Dim3T<T> operator+ (T s) const noexcept
//     {
//         return Dim3T<T>(x + s, y + s, z + s);
//     }
//     //! Returns component-wise difference of IntVect and argument.
//     AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
//     Dim3T<T> operator- (const Dim3T<T>& rhs) const noexcept
//     {
//         return Dim3T<T>(x - rhs.x, y - rhs.y, z - rhs.z);
//     }
//     //! Return an IntVect that is this IntVect - s.
//     AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
//     Dim3T<T> operator- (T s) const noexcept
//     {
//         return Dim3T<T>(x - s, y - s, z - s);
//     }
//     //! Returns component-wise product of IntVect and argument.
//     AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
//     Dim3T<T> operator* (const Dim3T<T>& rhs) const noexcept
//     {
//         return Dim3T<T>(x * rhs.x, y * rhs.y, z * rhs.z);
//     }
//     //! Returns component-wise product of IntVect and s.
//     AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
//     Dim3T<T> operator* (T s) const noexcept
//     {
//         return Dim3T<T>(x * s, y * s, z * s);
//     }
//     //! Returns component-wise division of IntVect by argument.
//     AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
//     Dim3T<T> operator/ (const Dim3T<T>& rhs) const noexcept
//     {
//         return Dim3T<T>(x / rhs.x, y / rhs.y, z / rhs.z);
//     }
//     //! Returns component-wise division of IntVect by s.
//     AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
//     Dim3T<T> operator/ (T s) const noexcept
//     {
//         return Dim3T<T>(x / s, y / s, z / s);
//     }

//     AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
//     T dot(const Dim3T<T>& rhs) const noexcept
//     {
//         return x*rhs.x + y*rhs.y + z*rhs.z;
//     }

//     AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
//     Dim3T<T> cross(const Dim3T<T>& rhs) const noexcept
//     {
//         return Dim3T<T>(y*rhs.z - z*rhs.y, z*rhs.x - x*rhs.z, x*rhs.y - y*rhs.x);
//     }

//     T x; T y; T z;

// };

// struct Dim3
//     : public Dim3T<int>
// {
//     // AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
//     // Dim3 () { define(0,0,0); }

//     AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
//     Dim3 (int a_x, int a_y, int a_z) { define(a_x,a_y,a_z); };

//     AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
//     Dim3 (const IntVect& iv)
//     {
//         AMREX_D_PICK(define(iv[0],0    ,0    );,
//                      define(iv[0],iv[1],0    );,
//                      define(iv[0],iv[1],iv[2]););
//     }

//     IntVect toIntVect () { return IntVect(*this); }

//     bool operator== (const IntVect& iv) { return iv == *this; }

//     bool operator!= (const IntVect& iv) { return iv != *this; }
// };

// struct XDim3
//     : public Dim3T<Real>
// {

// };
// template <class T>
// std::ostream& operator<< (std::ostream& os, const Dim3T<T>& d)
// {
//     os << '(' << d.x << ',' << d.y << ',' << d.z << ')';
//     return os;
// }

}

#endif
