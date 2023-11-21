#ifndef AMREX_MLABECLAP_3D_K_H_
#define AMREX_MLABECLAP_3D_K_H_
#include <AMReX_Config.H>

#include <AMReX_Array4.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Array.H>



namespace amrex {

// AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE

void
mledgecurlcurl_apply (
    const Box& box_x, const Box& box_y, const Box& box_z,
    const Array4<      Real>& out_x, const Array4<      Real>& out_y, const Array4<      Real>& out_z,
    const Array4<const Real>& in_x , const Array4<const Real>& in_y , const Array4<const Real>& in_z ,
    const Array4<const Real>& nu_xy, const Array4<const Real>& nu_zx, const Array4<const Real>& nu_yz,
    const Array4<const Real>& omega,
    const GpuArray<Real,AMREX_SPACEDIM>& drinv)
{
    
}

Real
mledgecurlcurl_apply_out_x (
    int i, int j, int k, int n,
    const Array4<const Real>& in_x , const Array4<const Real>& in_y , const Array4<const Real>& in_z ,
    const Array4<const Real>& nu_xy, const Array4<const Real>& nu_zx,
    const Array4<const Real>& omega, const GpuArray<Real,AMREX_SPACEDIM>& drinv)
{
    Real Hxy_py, Hxy_my, Hzx_pz, Hzx_mz;

    Hxy_py = mledgecurlcurl_get_Hxy(i,j  ,k,n,in_x,in_y,nu_xy,drinv);
    Hxy_my = mledgecurlcurl_get_Hxy(i,j-1,k,n,in_x,in_y,nu_xy,drinv);

    Hzx_pz = mledgecurlcurl_get_Hzx();
    Hzx_mz = mledgecurlcurl_get_Hzx();

    return omega(i,j,k)*in_x(i,j,k,n) + drinv[1]*(Hxy_py - Hxy_my) - drinv[2]*(Hzx_pz - Hzx_mz);
}

Real
mledgecurlcurl_get_Hxy (
    int i, int j, int k, int n,
    const Array4<const Real>& Ax,  const Array4<const Real>& Ay,
    const Array4<const Real>& nu_xy, const GpuArray<Real,AMREX_SPACEDIM>& drinv)
{
    return nu_xy(i,j,k)*((Ay(i+1,j,k,n) - Ay(i,j,k,n))*drinv[0] - (Ax(i,j+1,k,n) - Ax(i,j,k,n))*drinv[1]);
}

void
mledgecurlcurl_smooth (
    const Box& box_x, const Box& box_y, const Box& box_z,
    const Array4<      Real>& sol_x, const Array4<      Real>& sol_y, const Array4<      Real>& sol_z,
    const Array4<const Real>& rhs_x, const Array4<const Real>& rhs_y, const Array4<const Real>& rhs_z ,
    const Array4<const Real>& nu_xy, const Array4<const Real>& nu_zx, const Array4<const Real>& nu_yz,
    const Array4<const Real>& omega,
    const GpuArray<Real,AMREX_SPACEDIM>& drinv)
{

}

void
mledgecurlcurl_smooth_x (
    int i, int j, int k, int n,
    const Array4<      Real>& sol_x, const Array4<      Real>& sol_y, const Array4<      Real>& sol_z,
    const Array4<const Real>& rhs_x, const Array4<const Real>& rhs_y, const Array4<const Real>& rhs_z ,
    const Array4<const Real>& nu_xy, const Array4<const Real>& nu_zx, const Array4<const Real>& nu_yz,
    const Array4<const Real>& omega,
    const GpuArray<Real,AMREX_SPACEDIM>& drinv)
{
    const Real Lsol_x = mledgecurlcurl_apply_out_x(i,j,k,n,sol_x,sol_y,sol_z,nu_xy,nu_zx,omega,drinv);
    const Real den = omega(i,j,k,n)
                   + (nu_xy(i,j,k) + nu_xy(i,j-1,k))*drinv[1]*drinv[1]
                   + (nu_zx(i,j,k) + nu_zx(i,j,k-1))*drinv[2]*drinv[2];
    sol_x(i,j,k,n) += (rhs_x(i,j,k,n) - Lsol_x) / den;
}


















void
mledgecurlcurl_eb_apply (
    const Box& box_x, const Box& box_y, const Box& box_z,
    const Array4<      Real>& out_x, const Array4<      Real>& out_y, const Array4<      Real>& out_z,
    const Array4<const Real>& in_x , const Array4<const Real>& in_y , const Array4<const Real>& in_z ,
    const Array4<const Real>& nu_xy , const Array4<const Real>& nu_zx , const Array4<const Real>& nu_yz ,
    const Array4<const Real>& omega,
    const GpuArray<Real,AMREX_SPACEDIM>& dr, const GpuArray<Real,AMREX_SPACEDIM>& drinv
)
{
    
}

void
mledgecurlcurl_eb_apply_out_x_ijkn (
    int i, int j, int k, int n, const Array4<Real>& out_x,
    const Array4<const Real>& Ax,  const Array4<const Real>& Ay,  const Array4<const Real>& Az,
    const Array4<const Real>& ecx, const Array4<const Real>& ecy, const Array4<const Real>& ecz,
    const Array4<const Real>& Aeb_dr_xy, const Array4<const Real>& Aeb_dr_zx, const Array4<const Real>& Aeb_dr_yz,
    const Array4<const Real>& aFrac_xy,  const Array4<const Real>& aFrac_zx,  const Array4<const Real>& aFrac_yz,
    const Array4<const Real>& omega,
    const Array4<const Real>& nu_xy, const Array4<const Real>& nu_zx, const Array4<const Real>& nu_yz,
    const GpuArray<Real,AMREX_SPACEDIM>& dr, const GpuArray<Real,AMREX_SPACEDIM>& drinv)
{
    Real Hxy_py, Hxy_my, Hzx_pz, Hzx_mz;

    Hxy_py = mledgecurlcurl_eb_get_Hxy(i,j  ,k,n,Ax,Ay,Aeb_dr_xy,ecx,ecy,aFrac_xy,nu_xy,dr,drinv);
    Hxy_my = mledgecurlcurl_eb_get_Hxy(i,j-1,k,n,Ax,Ay,Aeb_dr_xy,ecx,ecy,aFrac_xy,nu_xy,dr,drinv);

    Hzx_pz = mledgecurlcurl_get_Hzx(i,j,k  ,n,Az,Ax,AMREX_EB_ONLY_ARGS(Aeb_dr_zx,ecz,ecx,aFrac_zx,)nu_zx,dr,drinv);
    Hzx_mz = mledgecurlcurl_get_Hzx(i,j,k-1,n,Az,Ax,AMREX_EB_ONLY_ARGS(Aeb_dr_zx,ecz,ecx,aFrac_zx,)nu_zx,dr,drinv);

    out_x(i,j,k,n) = omega(i,j,k)*Ax(i,j,k,n) + drinv[1]*(Hxy_py - Hxy_my) + drinv[2]*(Hzx_mz - Hzx_pz);
}

Real
mledgecurlcurl_eb_get_Hxy (
    int i, int j, int k, int n,
    const Array4<const Real>& Ax,  const Array4<const Real>& Ay,  const Array4<const Real>& Aeb_dr_xy,
    const Array4<const Real>& ecx, const Array4<const Real>& ecy, const Array4<const Real>& aFrac_xy,
    const Array4<const Real>& nu_xy,
    const GpuArray<Real,AMREX_SPACEDIM>& dr, const GpuArray<Real,AMREX_SPACEDIM>& drinv)
{
    if (aFrac_xy(i,j,k) == Real(1.0)) {
        // return nu_xy(i,j,k)*drinv[0]*drinv[1]*mledgecurlcurl_loop_on_face(
        //             Ax(i,j,k,n), Ax(i,j+1,k,n), Ay(i,j,k,n), Ay(i+1,j,k,n),
        //             dr[0], dr[0], dr[1], dr[1]);
        return mledgecurlcurl_get_Hxy(i,j,k,n,Ax,Ay,nu_xy,drinv);
    } else {
        return nu_xy(i,j,k)*drinv[0]*drinv[1]*mledgecurlcurl_eb_loop_on_face(
            Ax(i,j,k,n), Ax(i,j+1,k,n), Ay(i,j,k,n), Ay(i+1,j,k,n), Aeb_dr_xy(i,j,k,n),
            ecx(i,j,k),  ecx(i,j+1,k),  ecy(i,j,k),  ecy(i+1,j,k),
            dr[0], dr[0], dr[1], dr[1]);
    }
}

Real
mledgecurlcurl_loop_on_face (
    Real A1, Real A1_p2, Real A2, Real A2_p1,
    Real dr1, Real dr1_p2, Real dr2, Real dr2_p1)
{
    return A1*dr1 + A2_p1*dr2_p1 - A1_p2*dr1_p2 - A2*dr2;
}

/**
 * \brief Loop integral of edge centered field on face
 */
Real
mledgecurlcurl_eb_loop_on_face (
    Real A1, Real A1_p2, Real A2, Real A2_p1, Real Aeb_dr_on_face,
    Real ec1, Real ec1_p2, Real ec2, Real ec2_p1,
    Real dr1, Real dr1_p2, Real dr2, Real dr2_p1)
{
    Real loop_integral = Aeb_dr_on_face;
    loop_integral += A1    * dr1    * LengthFractionFromEdgeCentroid(ec1);
    loop_integral += A2_p1 * dr2_p1 * LengthFractionFromEdgeCentroid(ec2_p1);
    loop_integral -= A1_p2 * dr1_p2 * LengthFractionFromEdgeCentroid(ec1_p2);
    loop_integral -= A2    * dr2    * LengthFractionFromEdgeCentroid(ec2);
    return loop_integral;
}



Real
LengthFractionFromEdgeCentroid (Real edgeCentroid) {
    if (edgeCentroid == Real(1.0)) {
        return Real(1.0);
    } else if (edgeCentroid == Real(-1.0)) {
        return Real(0.0);
    } else {
        return Real(1.0) - Real(2.0) * Math::abs(edgeCentroid);
    }
}

}
#endif

// void
// edge_loop_integral_on_xy_face (
//     Array4<Real>& fluxBxy,
//     IntVect i_j_k_, IntVect iij_k_, IntVect i_jjk_,
//     const Array4<const Real>& Ax, const Array4<const Real>& Ay,
// #ifdef AMREX_USE_EB
//     Array4<Real const> const& Aeb_lineint_xy,
//     Array4<Real const> const& ecx, Array4<Real const> const& ecy,
//     const Array4<const Real>& areaFrac,
// #endif
//     const GpuArray<Real,AMREX_SPACEDIM>& dx
// )
// {
//     fluxBxy(i_j_k_) = EdgeLoopIntegralOnFace(
//         Ax(i_j_k_), Ax(iij_k_), Ay(i_j_k_), Ay(i_jjk_),
// #ifdef AMREX_USE_EB
//         Aeb_lineint_xy(i_j_k_), ecx(i_j_k_), ecx(i_jjk_), ecy(i_j_k_), ecy(iij_k_), areaFrac(i_j_k_) != Real(1.0),
// #endif
//         dx[0], dx[0], dx[1], dx[1]
//     );
// }