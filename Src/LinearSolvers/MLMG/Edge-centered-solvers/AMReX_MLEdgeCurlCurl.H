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
    MultiFabSet m_Aeb_dr_line_integral;
#endif

};

void
MLEdgeCurlCurl::Fapply (int amrlev, int mglev, MultiFabSet& out, const MultiFabSet& in) const
{
    BL_PROFILE("MLEdgeCurlCurl::Fapply()");

    AMREX_D_TERM(const MultiFab& in_mf_x = in[0];,
                 const MultiFab& in_mf_y = in[1];,
                 const MultiFab& in_mf_z = in[2];);
    AMREX_D_TERM(MultiFab& out_mf_x = out[0];,
                 MultiFab& out_mf_y = out[1];,
                 MultiFab& out_mf_z = out[2];);
    const auto dr = m_geom[amrlev][mglev].CellSizeArray();
    const auto drinv = m_geom[amrlev][mglev].InvCellSizeArray();

#ifdef AMREX_USE_EB
    const auto *factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][mglev].get());
    const auto &edgecent = factory->getEdgeCent();
    const auto &areafrac = factory->getAreaFrac();
    const auto &levset_mf = factory->getLevelSet();
    const MultiFab& Aeb_dr_mf_xy = m_Aeb_dr_line_integral[0];
    AMREX_3D_ONLY(const MultiFab& Aeb_dr_mf_xz = m_Aeb_dr_line_integral[1];
                  const MultiFab& Aeb_dr_mf_yz = m_Aeb_dr_line_integral[2];);
#endif

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(out_mf_x,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        AMREX_D_TERM(const Box box_x = mfi.tilebox(out_mf_x.ixType().toIntVect());,
                     const Box box_y = mfi.tilebox(out_mf_y.ixType().toIntVect());,
                     const Box box_z = mfi.tilebox(out_mf_z.ixType().toIntVect()););
        AMREX_D_TERM(const Array4<const Real> &inX = in_mf_x.const_array(mfi);,
                     const Array4<const Real> &inY = in_mf_y.const_array(mfi);,
                     const Array4<const Real> &inZ = in_mf_z.const_array(mfi););
        AMREX_D_TERM(const Array4<      Real> &outX = out_mf_x.array(mfi);,
                     const Array4<      Real> &outY = out_mf_x.array(mfi);,
                     const Array4<      Real> &outZ = out_mf_x.array(mfi););
        // Array4<int const> const& dmarr = dmask.const_array(mfi);
#ifdef AMREX_USE_EB
        bool cutfab = edgecent[0]->ok(mfi);
        if (cutfab) {
            AMREX_D_TERM(const Array4<const Real> &ecx = edgecent[0]->const_array(mfi);,
                         const Array4<const Real> &ecy = edgecent[1]->const_array(mfi);,
                         const Array4<const Real> &ecz = edgecent[2]->const_array(mfi););
            AMREX_2D_ONLY(


            );
            AMREX_3D_ONLY(const Array4<const Real> &aFrac_yz = areafrac[0]->const_array(mfi);
                          const Array4<const Real> &aFrac_zx = areafrac[1]->const_array(mfi);
                          const Array4<const Real> &aFrac_xy = areafrac[2]->const_array(mfi);
            );
            const auto &levset = levset_mf.const_array(mfi);
#endif

            if (Aeb == std::numeric_limits<Real>::lowest()) {
                const auto &Aeb_dr_xy = Aeb_dr_mf_xy.const_array(mfi);
                const auto &Aeb_dr_xz = Aeb_dr_mf_xz.const_array(mfi);
                const auto &Aeb_dr_yz = Aeb_dr_mf_yz.const_array(mfi);
                AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM (
                    box_x, tbox_x,
                    {
                        mledgecurlcurl_apply_eb_x(tbox_x,outX,
                            inX,inY,inZ,
                            Aeb_dr_xy,Aeb_dr_xz,Aeb_dr_yz,
                            ecx,ecy,ecz,
                            apxy,apxz,apyz,
                            dr,drinv);
                    },
                    box_y, tbox_y,
                    {
                        mledgecurlcurl_apply_eb_y();
                    },
                    box_z, tbox_z,
                    {
                        mledgecurlcurl_apply_eb_z();
                    }
                );
            } else {
                {
                    AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
                    {
                        mledgecurlcurl_apply_eb(i,j,k,yarr,xarr,levset,dmarr,AMREX_D_DECL(ecx,ecy,ecz),
                                             phieb, AMREX_D_DECL(bx,by,bz));
                    });
                }
            }
        } else
#endif // AMREX_USE_EB
        {
            mledgecurlcurl_apply(
                box_x,box_y,box_z,
                outX,outY,outZ,
                inX,inY,inZ,
                nu12,nu13,nu23,
                omega,dr,drinv);
        }
    }
}

AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM (
                xbx, txbx,
                {
                  mlebtensor_flux_x(txbx, Ax, fx, apx, fcx, bscalar, msk, face_only, xbx);
                }
                , ybx, tybx,
                {
                  mlebtensor_flux_y(tybx, Ay, fy, apy, fcy, bscalar, msk, face_only, ybx);
                }
                , zbx, tzbx,
                {
                  mlebtensor_flux_z(tzbx, Az, fz, apz, fcz, bscalar, msk, face_only, zbx);
                }
            );

}

#endif // AMREX_ML_EDGECURLCURL_H_