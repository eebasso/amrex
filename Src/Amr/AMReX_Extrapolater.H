#ifndef AMREX_EXTRAPOLATER_H_
#define AMREX_EXTRAPOLATER_H_
#include <AMReX_Config.H>

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>

namespace amrex::Extrapolater
{
    // finebnd: boundary cells covered by fine cells (including periodically shifted fine cells)
    const int finebnd = 1;
    // crsebnd: boundary cells not covered by fine cells
    const int crsebnd = 0;
    // physbnd: boundary cells outside the domain (excluding periodic boundaries)
    const int physbnd = 0;
    // interior: interior cells
    const int interior = 1;

    //! It is expected that FillBoundary (w/ periodicity) has been called on mf.
    void FirstOrderExtrap (MultiFab& mf, const Geometry& geom, int scomp, int ncomp);

    //! It is expected that FillBoundary (w/ periodicity) has been called on mf.
    void FirstOrderExtrap (MultiFab& mf, const Geometry& geom, int scomp, int ncomp, int nGrow);

}

#endif
