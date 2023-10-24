#ifndef AMREX_MLABECLAP_3D_K_H_
#define AMREX_MLABECLAP_3D_K_H_
#include <AMReX_Config.H>
#include <AMReX_Array4.H>
#include <AMReX_MultiFab.H>

namespace amrex {

// AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE

void
mledgecurlcurl_apply_array4 (
    Array4<Real> Ax, Array4<Real> Ay, Array4<Real> Az,
    Array4<Real> Jx, Array4<Real> Jy, Array4<Real> Jz,
    const Array4<Real> alpha, const Array4<Real> nu,
    const GpuArray<Real,AMREX_SPACEDIM>& dxinv)
{
    
}

void
edge_loop_integral (
    MultiFab Ax, MultiFab Ay, MultiFab Az,
    MultiFab fluxBxy, MultiFab fluxBzx, MultiFab fluxByz
#ifdef AMREX_USE_EB
    , MultiFab Aeblineintegral
#endif
)
{
    Ax.const_arrays();
}


void
edge_loop_integral_on_xy_face (
    IntVect i_j_k_, IntVect iij_k_, IntVect i_jjk_, IntVect iijjk_,
    Array4<Real const> const& Ax, Array4<Real const> const& Ay,
    Array4<Real> const& fluxBxy,
#ifdef AMREX_USE_EB
    Array4<Real const> const& Aeb_lineint_xy,
    Array4<Real const> const& ecx, Array4<Real const> const& ecy,
#endif
    Array4<Real const> const& omega, Array4<Real const> const& nu

)
{
    
}

Real
EdgeLoopIntegralOnFace (
    Real A1m, Real A1p, Real A2m, Real A2p,
#ifdef AMREX_USE_EB
    Real Aeb_lineint_xy, Real ec1m, Real ec1p, Real ec2m, Real ec2p,
#endif
    Real dx1m, Real dx1p, Real dx2m, Real dx2p, bool face_is_cut)
{
#ifdef AMREX_USE_EB
    if (face_is_cut) {
        Real loop_integral = Aeb_lineint_xy;
        loop_integral += A1m * dx1m * GetLengthFromEdgeCentroid(ec1m);
        loop_integral += A2p * dx2p * GetLengthFromEdgeCentroid(ec2p);
        loop_integral -= A1p * dx1p * GetLengthFromEdgeCentroid(ec1p);
        loop_integral -= A2m * dx2m * GetLengthFromEdgeCentroid(ec2m);
        return loop_integral;
    }
#endif
    return A1m * dx1m + A2p * dx2p - A1p * dx1p - A2m * dx2m;
}



Real
GetLengthFromEdgeCentroid (Real edgeCentroid) {
    if (edgeCentroid == Real(1.0)) {
        return Real(1.0);
    } else if (edgeCentroid == Real(-1.0)) {
        return Real(0.0);
    } else {
        return Real(1.0) - Real(2.0) * std::abs(edgeCentroid);
    }
}

}
#endif