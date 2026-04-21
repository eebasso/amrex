#include "MyTest.H"

#include <AMReX_MLEBABecLap.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_EB2.H>

using namespace amrex;

MyTest::MyTest ()
{
    init_arrays();
    test_eb2();
}

void
MyTest::init_arrays ()
{
    Box domain(IntVect(AMREX_D_DECL(0,0,0)), IntVect(AMREX_D_DECL(0,0,0)));

    m_geom.define(domain);
    dx_iso = GpuArray<Real,AMREX_SPACEDIM>{AMREX_D_DECL(0.2, 0.2, 0.2)};
    dx_aniso = GpuArray<Real,AMREX_SPACEDIM>{AMREX_D_DECL(0.3, 0.5, 0.7)};

    m_box_arr.define(domain);
    m_dmap.define(m_box_arr);

    MFInfo mfinfo{};
    IntVect ng{0};
    // mfinfo.SetTag("Tests::EB2::MyTest");

    m_cellflag.define(m_box_arr, m_dmap, 1, ng, mfinfo);
    m_volfrac.define(m_box_arr, m_dmap, 1, ng, mfinfo);
    m_volcent.define(m_box_arr, m_dmap, AMREX_SPACEDIM, ng, mfinfo);

    m_bndryarea.define(m_box_arr, m_dmap, 1, ng, mfinfo);
    m_bndrycent.define(m_box_arr, m_dmap, AMREX_SPACEDIM, ng, mfinfo);
    m_bndrynorm.define(m_box_arr, m_dmap, AMREX_SPACEDIM, ng, mfinfo);

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {

        BoxArray box_arr_fcent = amrex::convert(m_box_arr, IntVect::TheDimensionVector(idim));

        m_areafrac[idim].define(box_arr_fcent, m_dmap, 1, ng, mfinfo);
        m_facecent[idim].define(box_arr_fcent, m_dmap, AMREX_SPACEDIM-1, ng, mfinfo);

        IntVect edge_type{1};
        edge_type[idim] = 0;
        BoxArray box_arr_ecent = amrex::convert(m_box_arr, edge_type);

        m_edgecent[idim].define(box_arr_ecent, m_dmap, 1, ng, mfinfo);
    }
}

void
MyTest::test_eb2 ()
{
    Real apxm = Real(0.3);
    Real apxp = Real(1.0);
    Real apym = Real(0.4);
    Real apyp = Real(1.0);

    for (MFIter mfi(m_box_arr, m_dmap); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        const Box& nbx = amrex::surroundingNodes(bx);
        Array4<Real> const& apx_arr = m_apx.array(mfi);
        Array4<Real> const& apy_arr = m_apy.array(mfi);

        auto& gfab = m_mgf[mfi];
        const Box& vbx = gfab.validbox();

        auto& levelset = gfab.getLevelSet();
        if (iter == 0) {
            gshop.fillFab(levelset, m_geom, gshop_run_on, bounding_box);
#ifdef AMREX_USE_GPU
            if (hybrid) {
                levelset.prefetchToDevice();
            }
#endif
        }

        auto& cellflag = m_cellflag[mfi];

        gfab.buildTypes(cellflag);

        const Array4<const Real> & clst = levelset.const_array();
        const Array4<Real> &  lst = levelset.array();
        const Array4<EBCellFlag> &cfg = m_cellflag.array(mfi);
        const Array4<Real> &vfr = m_volfrac.array(mfi);
        const Array4<Real> &ctr = m_volcent.array(mfi);
        const Array4<Real> &bar = m_bndryarea.array(mfi);
        const Array4<Real> &bct = m_bndrycent.array(mfi);
        const Array4<Real> &bnm = m_bndrynorm.array(mfi);

        auto& facetype = gfab.getFaceType();
        AMREX_D_TERM(
            Array4<Real> const& apx = m_areafrac[0].array(mfi);
            Array4<Real> const& fcx = m_facecent[0].array(mfi);
            Array4<EB2::Type_t> const& ftx = facetype[0].array();
            ,
            Array4<Real> const& apy = m_areafrac[1].array(mfi);
            Array4<Real> const& fcy = m_facecent[1].array(mfi);

            ,
            Array4<Real> const& apz = m_areafrac[2].array(mfi);
            Array4<EB2::Type_t> const& ftz = facetype[2].array();
        );
        AMREX_D_TERM(
            ,
            Array4<Real> const& fcz = m_facecent[2].array(mfi);

        );

        AMREX_D_TERM(
            ,
            Array4<EB2::Type_t> const& fty = facetype[1].array();,

        );

#if AMREX_SPACEDIM == 2
        EB2::build_faces


        // EB2::set_eb_data(
        //     0, 0, apx_arr, apy_arr, dx_aniso,
        //     vcent, bcent, levset
        // )
#else
        EB2::set_eb_data(
            0, 0, apx_arr, apy_arr, dx_aniso,
            vcent, bcent, levset
        )
#endif

    }


}

void
MyTest::write ()
{

}