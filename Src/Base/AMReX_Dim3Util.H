#ifndef AMREX_DIM3_UTIL_H_
#define AMREX_DIM3_UTIL_H_
#include <AMReX_Config.H>

#include <AMReX_Dim3.H>
#include <AMReX_IndexType.H>
#include <AMReX_IntVect.H>

namespace amrex {

template<class D>
D Dim3_plus (const D& a, const D& b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }

template<class D>
D Dim3_minus (const D& a, const D& b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }


// template<class D>
// D Dim3_times (const XDim3& a, const XDim3& b) { return {a.x * b.x, a.y * b.y, a.z * b.z}; }
// template<class D>
// D Dim3_divide (const XDim3& num, const XDim3& den) { return {num.x / den.x, num.y / den.y, num.z / den.z}; }


XDim3 Dim3_times (Real s, const XDim3& v) { return {s * v.x, s * v.y, s * v.z}; }

template<class D>
D Dim3_negate (const D& v) { return { -v.x, -v.y, -v.z}; }

int Dim3_dot (const Dim3& a, const Dim3& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
Real Dim3_dot (const XDim3& a, const XDim3& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }

Dim3 Dim3_cross (const Dim3& a, const Dim3& b) { return {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x}; }
XDim3 Dim3_cross (const XDim3& a, const XDim3& b) { return {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x}; }

Real Dim3_norm (const XDim3& a) { return Real(std::sqrt(Real(Dim3_dot(a,a)))); }

template<class D>
void Dim3_plus_equals (D& lhs, const D& rhs) { lhs.x += rhs.x; lhs.y += rhs.y; lhs.z += rhs.z; }

}

#endif /* AMREX_DIM3_UTIL_H_ */