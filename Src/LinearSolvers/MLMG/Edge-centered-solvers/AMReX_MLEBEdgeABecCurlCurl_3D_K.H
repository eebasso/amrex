#ifndef AMREX_MLEBEDGECURLCURL_3D_K_H_
#define AMREX_MLEBEDGECURLCURL_3D_K_H_
#include <AMReX_Config.H>

#include <AMReX_Array4.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Array.H>
#include <AMReX_Box.H>
#include <AMReX_REAL.H>


namespace amrex {

void
mlebabedgecurlcurl_apply (
    const Box& box_x, const Box& box_y, const Box& box_z,
    const Array4<      Real>& out_x, const Array4<      Real>& out_y, const Array4<      Real>& out_z,
    const Array4<const Real>& in_x , const Array4<const Real>& in_y , const Array4<const Real>& in_z ,
    const Array4<const Real>& beta_xy, const Array4<const Real>& beta_zx, const Array4<const Real>& beta_yz,
    const Array4<const Real>& alpha,
    const GpuArray<Real,AMREX_SPACEDIM>& dr, const GpuArray<Real,AMREX_SPACEDIM>& drinv)
{
    
}

void
mleb_AB_curlcurl_apply_out_x (
    int i, int j, int k, int n, const Array4<Real>& out_x,
    const Array4<const Real>& in_x,  const Array4<const Real>& in_y,  const Array4<const Real>& in_z,
    const Array4<const Real>& ecx, const Array4<const Real>& ecy, const Array4<const Real>& ecz,
    const Array4<const Real>& in_eb_dr_xy, const Array4<const Real>& in_eb_dr_zx, const Array4<const Real>& in_eb_dr_yz,
    const Array4<const Real>& aFrac_xy,  const Array4<const Real>& aFrac_zx,  const Array4<const Real>& aFrac_yz,
    const Array4<const Real>& alpha,
    const Array4<const Real>& beta_xy, const Array4<const Real>& beta_zx, const Array4<const Real>& beta_yz,
    const GpuArray<Real,AMREX_SPACEDIM>& dr, const GpuArray<Real,AMREX_SPACEDIM>& drinv)
{
    Real Hxy_py, Hxy_my, Hzx_pz, Hzx_mz;

    Real Hxy_py = mlebedgecurlcurl_get_Hxy(i,j  ,k,n,in_x,in_y,in_eb_dr_xy,ecx,ecy,aFrac_xy,beta_xy,dr,drinv);
    Real Hxy_my = mlebedgecurlcurl_get_Hxy(i,j-1,k,n,in_x,in_y,in_eb_dr_xy,ecx,ecy,aFrac_xy,beta_xy,dr,drinv);

    Real Hzx_pz = mledgecurlcurl_get_Hzx(i,j,k  ,n,in_z,in_x,AMREX_EB_ONLY_ARGS(in_eb_dr_zx,ecz,ecx,aFrac_zx,)beta_zx,dr,drinv);
    Real Hzx_mz = mledgecurlcurl_get_Hzx(i,j,k-1,n,in_z,in_x,AMREX_EB_ONLY_ARGS(in_eb_dr_zx,ecz,ecx,aFrac_zx,)beta_zx,dr,drinv);

    out_x(i,j,k,n) = alpha(i,j,k)*in_x(i,j,k,n) + drinv[1]*(Hxy_py - Hxy_my) + drinv[2]*(Hzx_mz - Hzx_pz);
}

Real
mledgecurlcurl_apply_out_x (
    int i, int j, int k, int n,
    const Array4<const Real>& in_x,  const Array4<const Real>& in_y,  const Array4<const Real>& in_z,
    const Array4<const Real>& ecx, const Array4<const Real>& ecy, const Array4<const Real>& ecz,
    const Array4<const Real>& in_eb_dr_xy, const Array4<const Real>& in_eb_dr_zx, const Array4<const Real>& in_eb_dr_yz,
    const Array4<const Real>& aFrac_xy,  const Array4<const Real>& aFrac_zx,  const Array4<const Real>& aFrac_yz,
    const Array4<const Real>& alpha,
    const Array4<const Real>& beta_xy, const Array4<const Real>& beta_zx, const Array4<const Real>& beta_yz,
    const GpuArray<Real,AMREX_SPACEDIM>& dr, const GpuArray<Real,AMREX_SPACEDIM>& drinv)
{
    Real Hxy_py = mlebedgecurlcurl_get_Hxy(i,j  ,k,n,in_x,in_y,in_eb_dr_xy,ecx,ecy,aFrac_xy,beta_xy,dr,drinv);
    Real Hxy_my = mlebedgecurlcurl_get_Hxy(i,j-1,k,n,in_x,in_y,in_eb_dr_xy,ecx,ecy,aFrac_xy,beta_xy,dr,drinv);

    Real Hzx_pz = mledgecurlcurl_get_Hzx(i,j,k  ,n,in_z,in_x,in_eb_dr_zx,ecz,ecx,aFrac_zx,beta_zx,dr,drinv);
    Real Hzx_mz = mledgecurlcurl_get_Hzx(i,j,k-1,n,in_z,in_x,in_eb_dr_zx,ecz,ecx,aFrac_zx,beta_zx,dr,drinv);

    return alpha(i,j,k)*in_x(i,j,k,n) + drinv[1]*(Hxy_py - Hxy_my) + drinv[2]*(Hzx_mz - Hzx_pz);
}

Real
mlebedgecurlcurl_get_Hxy (
    int i, int j, int k, int n,
    const Array4<const Real>& in_x,  const Array4<const Real>& in_y,  const Array4<const Real>& in_eb_dr_xy,
    const Array4<const Real>& ecx, const Array4<const Real>& ecy, const Array4<const Real>& aFrac_xy,
    const Array4<const Real>& beta_xy, Real dxinv, Real dyinv)
{
    Real result;
    if (aFrac_xy(i,j,k) == Real(1.0)) {
        result = dxinv*(in_y(i+1,j,k,n) - in_y(i,j,k,n)) - dyinv*(in_x(i,j+1,k,n) - in_x(i,j,k,n));
    } else {
        result = in_eb_dr_xy
        return beta_xy(i,j,k)*drinv[0]*drinv[1]*
            mledgecurlcurl_eb_loop_on_face(
                in_x(i,j,k,n), in_x(i,j+1,k,n), in_y(i,j,k,n), in_y(i+1,j,k,n), in_eb_dr_xy(i,j,k,n),
                ecx(i,j,k),  ecx(i,j+1,k),  ecy(i,j,k),  ecy(i+1,j,k),
                dr[0], dr[0], dr[1], dr[1]) 
    }
    return beta_xy(i,j,k)*result;
}


// Real
// mlebabedgecurlcurl_get_Hxy (
//     int i, int j, int k, int n,
//     const Array4<const Real>& in_x,  const Array4<const Real>& in_y,  const Array4<const Real>& in_eb_dr_xy,
//     const Array4<const Real>& ecx, const Array4<const Real>& ecy, const Array4<const Real>& aFrac_xy,
//     const Array4<const Real>& beta_xy,
//     const GpuArray<Real,AMREX_SPACEDIM>& dr, const GpuArray<Real,AMREX_SPACEDIM>& drinv)
// {
//     if (aFrac_xy(i,j,k) == Real(1.0)) {
//         // return beta_xy(i,j,k)*drinv[0]*drinv[1]*mledgecurlcurl_loop_on_face(
//         //             in_x(i,j,k,n), in_x(i,j+1,k,n), in_y(i,j,k,n), in_y(i+1,j,k,n),
//         //             dr[0], dr[0], dr[1], dr[1]);
//         return mledgecurlcurl_get_Hxy(i,j,k,n,in_x,in_y,beta_xy,drinv);
//     } else {
//         return beta_xy(i,j,k)*drinv[0]*drinv[1]*
//             mledgecurlcurl_eb_loop_on_face(
//                 in_x(i,j,k,n), in_x(i,j+1,k,n), in_y(i,j,k,n), in_y(i+1,j,k,n), in_eb_dr_xy(i,j,k,n),
//                 ecx(i,j,k),  ecx(i,j+1,k),  ecy(i,j,k),  ecy(i+1,j,k),
//                 dr[0], dr[0], dr[1], dr[1]) 
//     }
// }


Real
mlebedgecurlcurl_curl_on_face (
    Real A1, Real A1_p2, Real A2, Real A2_p1, Real Aeb_dr_on_face,
    Real eFrac1, Real eFrac1_p2, Real eFrac2, Real eFrac2_p1,
    Real drinv1, Real drinv2)
{
    Real curl_on_face = drinv1*drinv2*Aeb_dr_on_face; // 2 ops
    curl_on_face += drinv1 * (A2_p1 * eFrac2_p1 - A2 * eFrac2); // 5 ops
    curl_on_face -+ drinv2 * (A1_p2 * eFrac1_p2 - A1 * eFrac1); // 5 ops
    return curl_on_face;
}

Real
mlebedgecurlcurl_loop_on_face (
    Real A1, Real A1_p2, Real A2, Real A2_p1, Real Aeb_dr_on_face,
    Real eFrac1, Real eFrac1_p2, Real eFrac2, Real eFrac2_p1,
    Real dr1, Real dr2)
{
    Real loop_on_face = Aeb_dr_on_face; // 0 ops
    loop_on_face += dr2 * (A2_p1 * eFrac2_p1 - A2 * eFrac2); // 5 ops
    loop_on_face -+ dr1 * (A1_p2 * eFrac1_p2 - A1 * eFrac1); // 5 ops
    return loop_on_face;
}


/**
 * \brief Loop integral of edge centered field on face
 */
Real
mlebabedgecurlcurl_loop_on_face (
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


#endif /* AMREX_MLEBEDGECURLCURL_3D_K_H_ */