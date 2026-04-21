#ifndef MY_TEST_H_
#define MY_TEST_H_

#include <AMReX_EBFabFactory.H>
#include <AMReX_MLMG.H>
#include <AMReX_Array.H>
#include <AMReX_EB2_MultiGFab.H>

using namespace amrex;

class MyTest
{
public:

    MyTest ();

    void write ();

private:

    void define_multifabs ();

    // void init_fabs ();

    void test_eb2 (const GpuArray<Real,AMREX_SPACEDIM> &dx);

    GpuArray<Real,AMREX_SPACEDIM> m_dx_iso{AMREX_D_DECL(0.2, 0.2, 0.2)};
    GpuArray<Real,AMREX_SPACEDIM> m_dx_aniso{AMREX_D_DECL(7.0, 11.0, 13.0)};

    amrex::Geometry m_geom;
    amrex::BoxArray m_box_arr;
    amrex::DistributionMapping m_dmap;

    amrex::EB2::MultiGFab m_mgf;

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

    // Box domain{IntVect(AMREX_D_DECL(0,0,0)), IntVect(AMREX_D_DECL(0,0,0))};

    // Array<BaseFab<EB2::Type_t>, AMREX_SPACEDIM> m_facetype;

    // iMultiFab m_cutcellmask;

    // FArrayBox m_levelset;
    // EBCellFlagFab m_cellflag;

    // FArrayBox m_volfrac;
    // FArrayBox m_volcent;

    // FArrayBox m_bndryarea;
    // FArrayBox m_bndrycent;
    // FArrayBox m_bndrynorm;

    // Array<FArrayBox,AMREX_SPACEDIM> m_areafrac;
    // Array<FArrayBox,AMREX_SPACEDIM> m_facecent;
    // Array<FArrayBox,AMREX_SPACEDIM> m_edgecent;

    // iMultiFab m_cutcellmask;
};

#endif
