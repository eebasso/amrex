#ifndef MY_TEST_H_
#define MY_TEST_H_

#include <AMReX_EBFabFactory.H>
#include <AMReX_MLMG.H>
#include <AMReX_Array.H>

using namespace amrex;

class MyTest
{
public:

    MyTest ();

    void write ();

private:

    void init_arrays ();

    void test_set_eb_data ();

    amrex::Geometry geom;
    amrex::BoxArray box_arr;
    amrex::DistributionMapping dmap;
    std::unique_ptr<amrex::EBFArrayBoxFactory> factory;

    GpuArray<Real,AMREX_SPACEDIM> dx_iso;
    GpuArray<Real,AMREX_SPACEDIM> dx_aniso;

    MultiFab m_levelset;
    FabArray<EBCellFlagFab> m_cellflag;

    MultiFab m_volfrac;
    MultiFab m_volcent;

    MultiFab m_bndryarea;
    MultiFab m_bndrycent;
    MultiFab m_bndrynorm;

    Array<MultiFab,AMREX_SPACEDIM> m_areafrac;
    Array<MultiFab,AMREX_SPACEDIM> m_facecent;
    Array<MultiFab,AMREX_SPACEDIM> m_edgecent;

    iMultiFab m_cutcellmask;

    // amrex::MultiFab m_apx;
    // amrex::MultiFab m_apy;
    // amrex::MultiFab m_vfrac;
    // amrex::MultiFab m_barea;



};

#endif
