#ifndef AMREX_ML_EDGECURLCURL_H_
#define AMREX_ML_EDGECURLCURL_H_

#include <AMReX_Config.H>
#include <AMReX_MLLinOp.H>

namespace amrex {

class MultiFabSet;

/**
 * \brief Edge-centered embedded-boundary curl curl operator
 *
 * Solves `omega A + curl (nu curl A) = rhs`
 *
 * where `A`, `rhs`, and `omega` are edge-centered fields,
 * and `nu` is a face-centered field. Works both with and
 * without embedded-boundaries enabled.
 */
class MLEdgeCurlCurl
    : public MLLinOpT<MultiFabSet>

{
public:
    void Fapply (int amrlev, int mglev, MultiFabSet& out, const MultiFabSet& in) const;
protected:
    MultiFabSet m_omega;
    MultiFabSet m_nu;
    MultiFabSet m_flux;
#ifdef AMREX_USE_EB
    MultiFabSet m_A_eb_line_integral;
#endif

};

void
MLEdgeCurlCurl::Fapply (int amrlev, int mglev, MultiFabSet& out, const MultiFabSet& in) const
{
    BL_PROFILE("MLEdgeCurlCurl::Fapply()");

    const auto dxinv = m_geom[amrlev][mglev].InvCellSizeArray();

#ifdef AMREX_USE_EB
    const auto *factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][mglev].get());
    const auto &edgecent = factory->getEdgeCent();
    const auto &levset_mf = factory->getLevelSet();
#endif

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(out[0],TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box &box = mfi.tilebox();
        const Array4<const Real> &inX = in[0].const_array(mfi);
        const Array4<const Real> &inY = in[1].const_array(mfi);
        const Array4<const Real> &inZ = in[2].const_array(mfi);
        const Array4<Real> &outX = out[0].array(mfi);
        const Array4<Real> &outY = out[1].array(mfi);
        const Array4<Real> &outZ = out[2].array(mfi);
        // Array4<int const> const& dmarr = dmask.const_array(mfi);
#ifdef AMREX_USE_EB
        bool cutfab = edgecent[0]->ok(mfi);
        if (cutfab) {
            AMREX_D_TERM(const Array4<const Real> &ecx = edgecent[0]->const_array(mfi);,
                         const Array4<const Real> &ecy = edgecent[1]->const_array(mfi);,
                         const Array4<const Real> &ecz = edgecent[2]->const_array(mfi));
            const auto &levset = levset_mf.const_array(mfi);
            AMREX_D_TERM(const Array4<const Real> &AebLineInt_XY = m_A_eb_line_integral[2].const_array(mfi);,
                         const Array4<const Real> &AebLineInt_XZ = m_A_eb_line_integral[1].const_array(mfi),
                         const Array4<const Real> &AebLineInt_XY = m_A_eb_line_integral[0].const_array(mfi)
                        );
            if (phieb == std::numeric_limits<Real>::lowest()) {
                auto const& phiebarr = m_phi_eb[amrlev].const_array(mfi);
                    AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
                    {
                        mlebndfdlap_adotx_eb(i,j,k,yarr,xarr,levset,dmarr,AMREX_D_DECL(ecx,ecy,ecz),
                                             phiebarr, AMREX_D_DECL(bx,by,bz));
                    });
            } else {
#if (AMREX_SPACEDIM == 2)
                if (m_rz) {
                    AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
                    {
                        mlebndfdlap_adotx_rz_eb(i,j,k,yarr,xarr,levset,dmarr,ecx,ecy,
                                                phieb, sig0, dx0, dx1, xlo, alpha);
                    });
                } else
#endif
                {
                    AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
                    {
                        mlebndfdlap_adotx_eb(i,j,k,yarr,xarr,levset,dmarr,AMREX_D_DECL(ecx,ecy,ecz),
                                             phieb, AMREX_D_DECL(bx,by,bz));
                    });
                }
            }
        } else
#endif // AMREX_USE_EB
        {
#if (AMREX_SPACEDIM == 2)
            if (m_rz) {
                AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
                {
                    mlebndfdlap_adotx_rz(i,j,k,yarr,xarr,dmarr,sig0,dx0,dx1,xlo,alpha);
                });
            } else
#endif
            {
                AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
                {
                    mlebndfdlap_adotx(i,j,k,yarr,xarr,dmarr,AMREX_D_DECL(bx,by,bz));
                });
            }
        }
    }
}

}

#endif // AMREX_ML_EDGECURLCURL_H_