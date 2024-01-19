#ifndef AMREX_MLEDGECURLCURL_3D_K_H_
#define AMREX_MLEDGECURLCURL_3D_K_H_
#include <AMReX_Config.H>

#include <AMReX_Array4.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Array.H>
#include <AMReX_Box.H>
#include <AMReX_REAL.H>

namespace amrex {

// AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE

void
mledgecurlcurl_apply (
    const Box& box_x, const Box& box_y, const Box& box_z,
    const Array4<      Real>& out_x, const Array4<      Real>& out_y, const Array4<      Real>& out_z,
    const Array4<const Real>& in_x , const Array4<const Real>& in_y , const Array4<const Real>& in_z ,
    const Array4<const Real>& nu_xy, const Array4<const Real>& nu_zx, const Array4<const Real>& nu_yz,
    const Array4<const Real>& alpha,
    const GpuArray<Real,AMREX_SPACEDIM>& drinv)
{
    
}

Real
mledgecurlcurl_apply_out_x (
    int i, int j, int k, int n,
    const Array4<const Real>& in_x , const Array4<const Real>& in_y , const Array4<const Real>& in_z ,
    const Array4<const Real>& nu_xy, const Array4<const Real>& nu_zx,
    const Array4<const Real>& alpha, const GpuArray<Real,AMREX_SPACEDIM>& drinv)
{
    // Real Hxy_py = mledgecurlcurl_get_Hxy(i,j  ,k,n,in_x,in_y,nu_xy,drinv);
    // Real Hxy_my = mledgecurlcurl_get_Hxy(i,j-1,k,n,in_x,in_y,nu_xy,drinv);

    // Real Hzx_pz = mledgecurlcurl_get_Hzx();
    // Real Hzx_mz = mledgecurlcurl_get_Hzx();

    Real dxinv = drinv[0];
    Real dyinv = drinv[1];
    Real dzinv = drinv[2];

    Real Hxy_py = mledgecurlcurl_beta_times_curl_on_face(i,j  ,k,n,in_x,in_y,nu_xy,dxinv,dyinv);
    Real Hxy_my = mledgecurlcurl_beta_times_curl_on_face(i,j-1,k,n,in_x,in_y,nu_xy,dxinv,dyinv);

    Real Hzx_pz = mledgecurlcurl_beta_times_curl_on_face(i,j,k  ,n,in_z,in_x,nu_zx,dz,dx);
    Real Hzx_mz = mledgecurlcurl_beta_times_curl_on_face();


    Real nuCurlxy_py = mledgecurlcurl_curl_on_face(i,j  ,k,n,in_x,in_y,nu_xy,drinv);

    return alpha(i,j,k)*in_x(i,j,k,n) + drinv[1]*(Hxy_py - Hxy_my) - drinv[2]*(Hzx_pz - Hzx_mz);
}

// Real
// mledgecurlcurl_beta_times_curl_on_face (
//     int i, int j, int k, int n,
//     const Array4<const Real>& in_1,  const Array4<const Real>& in_2,
//     const Array4<const Real>& beta_12, Real drinv1, Real drinv2)
// {
//     return beta_12(i,j,k)*(drinv1*(in_2(i+1,j,k,n) - in_2(i,j,k,n)) - drinv2*(in_1(i,j+1,k,n) - in_1(i,j,k,n)));
// }

Real
mledgecurlcurl_get_Hxy (
    int i, int j, int k, int n,
    const Array4<const Real>& in_x,  const Array4<const Real>& in_y,
    const Array4<const Real>& nu_xy, const GpuArray<Real,AMREX_SPACEDIM>& drinv)
{
    return nu_xy(i,j,k)*(drinv[0]*(in_y(i+1,j,k,n) - in_y(i,j,k,n)) - drinv[1]*(in_x(i,j+1,k,n) - in_x(i,j,k,n)));
}



Real
mledgecurlcurl_curl_on_face (
    int i, int j, int k, int n,
    const Array4<const Real>& in_1,  const Array4<const Real>& in_2,
    Real drinv1, Real drinv2)
{
    return drinv1*(in_2(i+1,j,k,n) - in_2(i,j,k,n)) - drinv2*(in_1(i,j+1,k,n) - in_1(i,j,k,n));
}

void
mledgecurlcurl_smooth (
    const Box& box_x, const Box& box_y, const Box& box_z,
    const Array4<      Real>& sol_x, const Array4<      Real>& sol_y, const Array4<      Real>& sol_z,
    const Array4<const Real>& rhs_x, const Array4<const Real>& rhs_y, const Array4<const Real>& rhs_z ,
    const Array4<const Real>& nu_xy, const Array4<const Real>& nu_zx, const Array4<const Real>& nu_yz,
    const Array4<const Real>& alpha,
    const GpuArray<Real,AMREX_SPACEDIM>& drinv)
{

}

void
mledgecurlcurl_smooth_x (
    int i, int j, int k, int n,
    const Array4<      Real>& sol_x, const Array4<      Real>& sol_y, const Array4<      Real>& sol_z,
    const Array4<const Real>& rhs_x, const Array4<const Real>& rhs_y, const Array4<const Real>& rhs_z ,
    const Array4<const Real>& nu_xy, const Array4<const Real>& nu_zx, const Array4<const Real>& nu_yz,
    const Array4<const Real>& alpha,
    const GpuArray<Real,AMREX_SPACEDIM>& drinv)
{
    const Real Lsol_x = mledgecurlcurl_apply_out_x(i,j,k,n,sol_x,sol_y,sol_z,nu_xy,nu_zx,alpha,drinv);
    const Real den = alpha(i,j,k,n)
                   + (nu_xy(i,j,k) + nu_xy(i,j-1,k))*drinv[1]*drinv[1]
                   + (nu_zx(i,j,k) + nu_zx(i,j,k-1))*drinv[2]*drinv[2];
    sol_x(i,j,k,n) += (rhs_x(i,j,k,n) - Lsol_x) / den;
}

Real
mledgecurlcurl_loop_on_face (
    Real A1, Real A1_p2, Real A2, Real A2_p1,
    Real dr1, Real dr1_p2, Real dr2, Real dr2_p1)
{
    return A1*dr1 + A2_p1*dr2_p1 - A1_p2*dr1_p2 - A2*dr2;
}

}

// void
// edge_loop_integral_on_xy_face (
//     Array4<Real>& fluxBxy,
//     IntVect i_j_k_, IntVect iij_k_, IntVect i_jjk_,
//     const Array4<const Real>& in_x, const Array4<const Real>& in_y,
// #ifdef AMREX_USE_EB
//     Array4<Real const> const& Aeb_lineint_xy,
//     Array4<Real const> const& ecx, Array4<Real const> const& ecy,
//     const Array4<const Real>& areaFrac,
// #endif
//     const GpuArray<Real,AMREX_SPACEDIM>& dx
// )
// {
//     fluxBxy(i_j_k_) = EdgeLoopIntegralOnFace(
//         in_x(i_j_k_), in_x(iij_k_), in_y(i_j_k_), in_y(i_jjk_),
// #ifdef AMREX_USE_EB
//         Aeb_lineint_xy(i_j_k_), ecx(i_j_k_), ecx(i_jjk_), ecy(i_j_k_), ecy(iij_k_), areaFrac(i_j_k_) != Real(1.0),
// #endif
//         dx[0], dx[0], dx[1], dx[1]
//     );
// }

#endif /* AMREX_MLEDGECURLCURL_3D_K_H_ */
