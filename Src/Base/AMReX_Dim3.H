#ifndef AMREX_DIM3_H_
#define AMREX_DIM3_H_
#include <AMReX_Config.H>

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

#include <iostream>

namespace amrex {

template <class VT>
struct Dim3Template
{
    VT x; VT y; VT z;
    VT& operator[] (int d) { return d==0 ? x : d==1 ? y : z; }
    const VT& operator[] (int d) const { return d==0 ? x : d==1 ? y : z; }
};

struct Dim3 : public Dim3Template<int> {};

struct XDim3 : public Dim3Template<Real> {};

template <class VT>
std::ostream& operator<< (std::ostream& os, const Dim3Template<VT>& d)
{
    os << '(' << d.x << ',' << d.y << ',' << d.z << ')';
    return os;
}

}

#endif
