#include "MyTest.H"

#include <AMReX_MLEBABecLap.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_EB2.H>

MyTest::MyTest ()
{
    init_arrays();
    test_set_eb_data();
}

void
MyTest::init_arrays ()
{
    Box domain(IntVect(AMREX_D_DECL(0,0,0)), IntVect(AMREX_D_DECL(0,0,0)));

    geom.define(domain);
    dx_iso = GpuArray<Real,AMREX_SPACEDIM>{AMREX_D_DECL(0.2, 0.2, 0.2)};
    dx_aniso = GpuArray<Real,AMREX_SPACEDIM>{AMREX_D_DECL(0.3, 0.5, 0.7)};

    box_arr.define(domain);
    dmap.define(box_arr);

    MFInfo mfinfo{};
    IntVect ng{0};
    // mfinfo.SetTag("Tests::EB2::MyTest");

    m_cellflag.define(box_arr, dmap, 1, ng, mfinfo);
    m_volfrac.define(box_arr, dmap, 1, ng, mfinfo);
    m_volcent.define(box_arr, dmap, AMREX_SPACEDIM, ng, mfinfo);

    m_bndryarea.define(box_arr, dmap, 1, ng, mfinfo);
    m_bndrycent.define(box_arr, dmap, AMREX_SPACEDIM, ng, mfinfo);
    m_bndrynorm.define(box_arr, dmap, AMREX_SPACEDIM, ng, mfinfo);

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {

        BoxArray box_arr_fcent = amrex::convert(box_arr, IntVect::TheDimensionVector(idim));

        m_areafrac[idim].define(box_arr_fcent, dmap, 1, ng, mfinfo);
        m_facecent[idim].define(box_arr_fcent, dmap, AMREX_SPACEDIM-1, ng, mfinfo);

        IntVect edge_type{1};
        edge_type[idim] = 0;
        BoxArray box_arr_ecent = amrex::convert(box_arr, edge_type);

        m_edgecent[idim].define(box_arr_ecent, dmap, 1, ng, mfinfo);
    }

    // const EB2::IndexSpace& eb_is = EB2::IndexSpace::top();
    // const EB2::Level& eb_level = eb_is.getLevel(geom);
    // factory = std::make_unique<EBFArrayBoxFactory>
    //     (eb_level, geom, box_arr, dmap, Vector<int>{2,2,2}, EBSupport::full);
    // factory = std::make_unique<EBFArrayBoxFactory>();

    // m_levelset.define(box_arr, dmap, 1, 0, MFInfo());
    // m_cellflag.define(box_arr, dmap, 1, 0, MFInfo());

    // m_volfrac.define(box_arr, dmap, 1, 0, MFInfo(), *factory);
    // m_volcent.define(box_arr, dmap, 1, 0, MFInfo(), *factory);

    // m_bndryarea.define(box_arr, dmap, 1, 0, MFInfo(), *factory);
    // m_bndrycent.define(box_arr, dmap, 1, 0, MFInfo(), *factory);
    // m_bndrynorm.define(box_arr, dmap, 1, 0, MFInfo(), *factory);

    // for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    //     m_areafrac[d].define(box_arr, dmap, 1, 0, MFInfo(), *factory);
    //     m_facecent[d].define(box_arr, dmap, 1, 0, MFInfo(), *factory);
    //     m_edgecent[d].define(box_arr, dmap, 1, 0, MFInfo(), *factory);
    // }
}

void
MyTest::test_set_eb_data ()
{
    Real apxm = Real(0.3);
    Real apxp = Real(1.0);
    Real apym = Real(0.4);
    Real apyp = Real(1.0);

    for (MFIter mfi(box_arr, dmap); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        const Box& nbx = amrex::surroundingNodes(bx);
        Array4<Real> const& apx_arr = m_apx.array(mfi);
        Array4<Real> const& apy_arr = m_apy.array(mfi);


        amrex::EB2::set_eb_data(
            0, 0, apx_arr, apy_arr, dx_aniso,
            vfrac,
        )
    }
}

void
MyTest::write ()
{

}