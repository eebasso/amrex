#ifndef MY_TEST_H_
#define MY_TEST_H_

#include <AMReX_Array.H>
#include <AMReX_EB2_MultiGFab.H>
#include <AMReX_EBFabFactory.H>

class MyTest
{
public:

    MyTest ();

    void write ();

private:

    void define_multifabs ();

    void test_eb2 (const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> &dx);

    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> m_dx_iso{AMREX_D_DECL(0.2, 0.2, 0.2)};
    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> m_dx_aniso{AMREX_D_DECL(7.0, 11.0, 13.0)};

    amrex::BoxArray m_grid;
    amrex::DistributionMapping m_dmap;

    // amrex::EB2::MultiGFab m_mgf;

    amrex::MultiFab m_levelset;
    amrex::FabArray<amrex::EBCellFlagFab> m_cellflag;

    amrex::MultiFab m_volfrac;
    amrex::MultiFab m_volcent;

    amrex::MultiFab m_bndryarea;
    amrex::MultiFab m_bndrycent;
    amrex::MultiFab m_bndrynorm;

    amrex::Array<amrex::MultiFab,AMREX_SPACEDIM> m_areafrac;
    amrex::Array<amrex::MultiFab,AMREX_SPACEDIM> m_facecent;
    amrex::Array<amrex::MultiFab,AMREX_SPACEDIM> m_edgecent;

    // iMultiFab m_cutcellmask;
};

#endif
