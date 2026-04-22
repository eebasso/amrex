#include "MyTest.H"

#include <AMReX_EB2.H>
#include <AMReX_EB2_C.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

MyTest::MyTest ()
{
    define_multifabs();
    amrex::Print() << "Tests::EB2::MyTest: test isotropic cell sizes (dx = dy)";
    test_eb2(m_dx_iso);
    amrex::Print() << "Tests::EB2::MyTest: test anisotropic cell sizes (dx =/= dy)";
    test_eb2(m_dx_aniso);
}

void
MyTest::define_multifabs ()
{
    Box domain(IntVect(AMREX_D_DECL(0,0,0)), IntVect(AMREX_D_DECL(0,0,0)));

    m_geom.define(domain);

    m_grid.define(domain);
    m_dmap.define(m_grid);

    // FArrayBox m_levelset = FArrayBox{domain, 1};

    // MFInfo mfinfo{};
    // m_mgf.define(m_grid, m_dmap);
    // const int ng = amrex::EB2::GFab::ng;
    const int ng = 0;
    // IntVect ng{0};
    MFInfo mfinfo;
    // mfinfo.SetTag("Tests::EB2::MyTest");

    IntVect node_type = IntVect::TheNodeVector();
    m_levelset.define(amrex::convert(m_grid, node_type), m_dmap, 1, ng, mfinfo);

    m_cellflag.define(m_grid, m_dmap, 1, ng, mfinfo);
    m_volfrac.define(m_grid, m_dmap, 1, ng, mfinfo);
    m_volcent.define(m_grid, m_dmap, AMREX_SPACEDIM, ng, mfinfo);

    m_bndryarea.define(m_grid, m_dmap, 1, ng, mfinfo);
    m_bndrycent.define(m_grid, m_dmap, AMREX_SPACEDIM, ng, mfinfo);
    m_bndrynorm.define(m_grid, m_dmap, AMREX_SPACEDIM, ng, mfinfo);

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        IntVect face_type = IntVect::TheDimensionVector(idim);
        IntVect edge_type{1}; edge_type[idim] = 0;
        m_areafrac[idim].define(amrex::convert(m_grid, face_type), m_dmap, 1, ng, mfinfo);
        m_facecent[idim].define(amrex::convert(m_grid, face_type), m_dmap, AMREX_SPACEDIM-1, ng, mfinfo);
        m_edgecent[idim].define(amrex::convert(m_grid, edge_type), m_dmap, 1, ng, mfinfo);
    }
}

void
MyTest::test_eb2 (
    const GpuArray<Real,AMREX_SPACEDIM> &dx
)
{
    for (MFIter mfi(m_grid, m_dmap); mfi.isValid(); ++mfi)
    {
        // const Box& bx = mfi.validbox();
        // const Box& nbx = amrex::surroundingNodes(bx);
        // Array4<Real> const& apx_arr = m_apx.array(mfi);
        // Array4<Real> const& apy_arr = m_apy.array(mfi);
        // auto& gfab = m_mgf[mfi];
        // const Box& vbx = gfab.validbox();
        // auto& levelset = gfab.getLevelSet();

        // auto& cellflag = m_cellflag[mfi];
        // gfab.buildTypes(cellflag);
        // const Array4<EBCellFlag> &cflag_arr = m_cellflag.array(mfi);

        const Array4<Real> &vfrac_arr = m_volfrac.array(mfi);
        const Array4<Real> &vcent_arr = m_volcent.array(mfi);
        const Array4<Real> &barea_arr = m_bndryarea.array(mfi);
        const Array4<Real> &bcent_arr = m_bndrycent.array(mfi);
        const Array4<Real> &bnorm_arr = m_bndrynorm.array(mfi);
        const Array4<Real> &levset_arr = m_levelset.array(mfi);

        // auto& facetype = gfab.getFaceType();
        AMREX_D_TERM(
            const Array4<Real> &apx_arr = m_areafrac[0].array(mfi);
            const Array4<Real> &fcentx_arr = m_facecent[0].array(mfi);
            // const Array4<EB2::Type_t> &ftypex_arr = facetype[0].array();
            ,
            const Array4<Real> &apy_arr = m_areafrac[1].array(mfi);
            const Array4<Real> &fcenty_arr = m_facecent[1].array(mfi);
            // const Array4<EB2::Type_t> &ftypey_arr = facetype[0].array();

            ,
            const Array4<Real> &apz_arr = m_areafrac[2].array(mfi);
            const Array4<Real> &fcentz_arr = m_facecent[1].array(mfi);
            // const Array4<EB2::Type_t> &ftypez_arr = facetype[2].array();
        );

#if AMREX_SPACEDIM == 2

        Real apXm = Real(0.3);
        // Real apxp = Real(0.0);
        Real apYm = Real(0.4);
        // Real apyp = Real(0.0);

        apx_arr(0,0,0) = apXm;
        apx_arr(1,0,0) = Real(0.0);
        apy_arr(0,0,0) = apYm;
        apy_arr(0,1,0) = Real(0.0);

        levset_arr(0,0,0) = Real(1.0);
        levset_arr(1,0,0) = Real(-1.0);
        levset_arr(0,1,0) = Real(-1.0);
        levset_arr(1,1,0) = Real(-1.0);


        Real apnorm_exact =  hypot(apXm*dx[1], apYm*dx[0]);

        Real vfrac_exact = 0.5*apXm*apYm;
        Real vcentx_exact = (1./ vfrac_exact)*(1./12.)*apXm*apYm*(-3.0 + 2.0*apYm);
        Real vcenty_exact = (1./ vfrac_exact)*(1./12.)*apYm*apXm*(-3.0 + 2.0*apXm);
        Real barea_exact = (apnorm_exact * apnorm_exact) / hypot(apXm*dx[1]*dx[1], apYm*dx[0]*dx[0]);
        Real bcentx_exact = 0.5 - 0.5*apYm;
        Real bcenty_exact = 0.5 - 0.5*apXm;
        Real bnormx_exact = apXm*dx[1] / apnorm_exact;
        Real bnormy_exact = apYm*dx[0] / apnorm_exact;

        EB2::Test::set_eb_data_wrapper(
            0, 0,
            apx_arr, apy_arr, dx,
            vfrac_arr, vcent_arr,
            barea_arr, bcent_arr,
            bnorm_arr, levset_arr
        );

        Real vfrac = vfrac_arr(0, 0, 0);
        Real vcentx = vcent_arr(0, 0, 0, 0);
        Real vcenty = vcent_arr(0, 0, 0, 1);
        Real barea = barea_arr(0, 0, 0);
        Real bcentx = bcent_arr(0, 0, 0, 0);
        Real bcenty = bcent_arr(0, 0, 0, 1);
        Real bnormx = bnorm_arr(0, 0, 0, 0);
        Real bnormy = bnorm_arr(0, 0, 0, 1);

        Real tiny = 1.e-9;

        AMREX_ASSERT(amrex::abs(vfrac - vfrac_exact) < tiny);
        AMREX_ASSERT(amrex::abs(vcentx - vcentx_exact) < tiny);
        AMREX_ASSERT(amrex::abs(vcenty - vcenty_exact) < tiny);
        AMREX_ASSERT(amrex::abs(barea - barea_exact) < tiny);
        AMREX_ASSERT(amrex::abs(bcentx - bcentx_exact) < tiny);
        AMREX_ASSERT(amrex::abs(bcenty - bcenty_exact) < tiny);
        AMREX_ASSERT(amrex::abs(bnormx - bnormx_exact) < tiny);
        AMREX_ASSERT(amrex::abs(bnormy - bnormy_exact) < tiny);

#else
        // amrex::EB2::Testing::set_eb_data_wrapper(
        //     0, 0, apx_arr, apy_arr, m_dx_aniso, vcent_arr, bcent_arr, levset
        // )
#endif

    }


}

void
MyTest::write ()
{

}