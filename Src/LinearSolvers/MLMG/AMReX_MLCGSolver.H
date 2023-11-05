
#ifndef AMREX_MLCGSOLVER_H_
#define AMREX_MLCGSOLVER_H_
#include <AMReX_Config.H>

#include <AMReX_MLLinOp.H>

namespace amrex {

template <typename MF>
class MLCGSolverT
{
public:

    using FAB = typename MF::fab_type;
    using RT  = typename MF::value_type;

    enum struct Type { BiCGStab, CG };

    MLCGSolverT (MLLinOpT<MF>& _lp, Type _typ = Type::BiCGStab);
    ~MLCGSolverT ();

    MLCGSolverT (const MLCGSolverT<MF>& rhs) = delete;
    MLCGSolverT (MLCGSolverT<MF>&& rhs) = delete;
    MLCGSolverT<MF>& operator= (const MLCGSolverT<MF>& rhs) = delete;
    MLCGSolverT<MF>& operator= (MLCGSolverT<MF>&& rhs) = delete;

    void setSolver (Type _typ) noexcept { solver_type = _typ; }

    /**
    * solve the system, Lp(solnL)=rhsL to relative err, tolerance
    * Returns an int indicating success or failure.
    * 0 means success
    * 1 means failed for loss of precision
    * 2 means iterations exceeded
    */
    int solve (MF& solnL, const MF& rhsL, RT eps_rel, RT eps_abs);

    void setVerbose (int _verbose) { verbose = _verbose; }
    [[nodiscard]] int getVerbose () const { return verbose; }

    void setMaxIter (int _maxiter) { maxiter = _maxiter; }
    [[nodiscard]] int getMaxIter () const { return maxiter; }

    void setNGhost(int _nghost) {nghost = IntVect(_nghost);}
    [[nodiscard]] int getNGhost() {return nghost[0];}

    [[nodiscard]] RT dotxy (const MF& r, const MF& z, bool local = false);
    [[nodiscard]] RT norm_inf (const MF& res, bool local = false);
    int solve_bicgstab (MF& solnL, const MF& rhsL, RT eps_rel, RT eps_abs);
    int solve_cg (MF& solnL, const MF& rhsL, RT eps_rel, RT eps_abs);

    [[nodiscard]] int getNumIters () const noexcept { return iter; }

    void PrintDots (const std::string& title, const char* names, const MF& x, const MF& y);

private:

    MLLinOpT<MF>& Lp;
    Type solver_type;
    const int amrlev = 0;
    const int mglev;
    int verbose   = 0;
    int maxiter   = 100;
    IntVect nghost = IntVect(0);
    IntVect ngrow = IntVect(0);
    int iter = -1;
};

static bool printedHeader = false;

template <typename MF>
MLCGSolverT<MF>::MLCGSolverT (MLLinOpT<MF>& _lp, Type _typ)
    : Lp(_lp), solver_type(_typ), mglev(_lp.NMGLevels(0)-1)
{}

template <typename MF> MLCGSolverT<MF>::~MLCGSolverT () = default;

template <typename MF>
int
MLCGSolverT<MF>::solve (MF& sol, const MF& rhs, RT eps_rel, RT eps_abs)
{
    if (solver_type == Type::BiCGStab) {
        return solve_bicgstab(sol,rhs,eps_rel,eps_abs);
    } else {
        return solve_cg(sol,rhs,eps_rel,eps_abs);
    }
}

template <typename MF>
int
MLCGSolverT<MF>::solve_bicgstab (MF& sol, const MF& rhs, RT eps_rel, RT eps_abs)
{
    BL_PROFILE("MLCGSolver::bicgstab");
    if (!printedHeader) {
        amrex::Print() << "MLCGSolver_BiCGStab: print-test branch" << "\n";
        printedHeader = true;
    }

    std::string MLCGSolverTitle = "MLCGSolver_BiCGStab: ";
    std::string TitleHalf_Iter = "MLCGSolver_BiCGStab: Half Iter           ";
    std::string TitleIteration = "MLCGSolver_BiCGStab: Iteration           ";
    std::string TitlePost_Calc = "MLCGSolver_BiCGStab: Post_Calc           ";

    std::string Half_Iter = TitleHalf_Iter + std::to_string(0);
    std::string Iteration = TitleIteration + std::to_string(0);
    std::string Post_Calc = TitlePost_Calc + std::to_string(0);

    amrex::Print() << "MLCGSolver_BiCGStab: nghost = " << nghost << "\n";
    amrex::Print() << "MLCGSolver_BiCGStab: ngrow  = " << sol.nGrowVect() << "\n";

    ngrow = sol.nGrowVect();
    nghost = sol.nGrowVect();
    verbose = 3;

    const int ncomp = sol.nComp();

    const BoxArray& ba = sol.boxArray();
    const DistributionMapping& dm = sol.DistributionMap();
    const auto& factory = sol.Factory();

    // MF ph(ba, dm, ncomp, sol.nGrowVect(), MFInfo(), factory);
    // MF sh(ba, dm, ncomp, sol.nGrowVect(), MFInfo(), factory);
    // ph.setVal(RT(0.0));
    // sh.setVal(RT(0.0));
    // PrintDots(Iteration,"Dots(ph,ph) initialized",ph,ph);
    // PrintDots(Iteration,"Dots(sh,sh) initialized",sh,sh);

    MF sorig(ba, dm, ncomp, nghost, MFInfo(), factory);
    MF p    (ba, dm, ncomp, nghost, MFInfo(), factory);
    MF r    (ba, dm, ncomp, nghost, MFInfo(), factory);
    // MF s    (ba, dm, ncomp, nghost, MFInfo(), factory);
    MF rh   (ba, dm, ncomp, nghost, MFInfo(), factory);
    MF v    (ba, dm, ncomp, nghost, MFInfo(), factory);
    MF t    (ba, dm, ncomp, nghost, MFInfo(), factory);
    sorig.setVal(RT(0.0));
    p.setVal(RT(0.0));
    r.setVal(RT(0.0));
    // s.setVal(RT(0.0));
    rh.setVal(RT(0.0));
    v.setVal(RT(0.0));
    t.setVal(RT(0.0));


    PrintDots(Iteration,"Dots(sol,sol) before Lp.correctionResidual",sol,sol);
    Lp.correctionResidual(amrlev, mglev, r, sol, rhs, MLLinOpT<MF>::BCMode::Homogeneous);
    PrintDots(Iteration,"Dots(sol,sol) after  Lp.correctionResidual",sol,sol);

    // Then normalize
    PrintDots(Iteration,"Dots(r,r) before Lp.normalize",r,r);
    Lp.normalize(amrlev, mglev, r);
    PrintDots(Iteration,"Dots(r,r) after  Lp.normalize",r,r);

    sorig.LocalCopy(sol,0,0,ncomp,nghost);
    PrintDots(Iteration,"Dots(sorig,sorig)",sorig,sorig);

    rh.LocalCopy   (r  ,0,0,ncomp,nghost);
    PrintDots(Iteration,"Dots(rh,rh)",rh,rh);

    sol.setVal(RT(0.0));
    PrintDots(Iteration,"Dots(sol,sol) after sol.setVal(0)",sol,sol);

    RT rnorm = norm_inf(r);
    const RT rnorm0 = rnorm;

    if ( verbose > 0 )
    {
        amrex::Print() << "MLCGSolver_BiCGStab: Initial error (error0) =        " << rnorm0 << '\n';
    }
    int ret = 0;
    iter = 1;
    RT alpha, beta, omega, rho, rhTv;

    if ( rnorm0 == 0 || rnorm0 < eps_abs )
    {
        if ( verbose > 0 )
        {
            amrex::Print() << "MLCGSolver_BiCGStab: niter = 0,"
                           << ", rnorm = " << rnorm
                           << ", eps_abs = " << eps_abs << std::endl;
        }
        return ret;
    }

    rho = dotxy(rh,r);
    PrintDots(Post_Calc,"Dots(rh,r)",rh,r);
    p.LocalCopy(r,0,0,ncomp,nghost);
    PrintDots(Post_Calc, "Dots(p,p)",p,p);

    for (; iter <= maxiter; ++iter)
    {
        Half_Iter = TitleHalf_Iter + std::to_string(iter);
        Iteration = TitleIteration + std::to_string(iter);
        Post_Calc = TitlePost_Calc + std::to_string(iter);
        if ( rho == 0 )
        {
            ret = 1; break;
        }
        // ph.LocalCopy(p,0,0,ncomp,nghost);
        // ph.LocalCopy(p,0,0,ncomp,sol.nGrowVect());
        PrintDots(Half_Iter,"Dots(p,p) before Lp.apply",p,p);
        Lp.apply(amrlev, mglev, v, p, MLLinOpT<MF>::BCMode::Homogeneous, MLLinOpT<MF>::StateMode::Correction);
        PrintDots(Half_Iter,"Dots(p,p) after  Lp.apply",p,p);

        PrintDots(Half_Iter,"Dots(v,v) before Lp.normalize",v,v);
        Lp.normalize(amrlev, mglev, v);
        PrintDots(Half_Iter,"Dots(v,v) after  Lp.normalize",v,v);

        rhTv = dotxy(rh,v);
        PrintDots(Half_Iter,"Dots(rh,v)",rh,v);

        if ( rhTv == RT(0.0) )
        {
            ret = 2; break;
        }
        alpha = rho/rhTv;
        // amrex::Print() << Half_Iter << " alpha: " << alpha << "\n";
        MF::Saxpy(sol, alpha, p, 0, 0, ncomp, nghost); // sol += alpha * p
        // MF::LinComb(s, RT(1.0), r, 0, -alpha, v, 0, 0, ncomp, nghost); // s = r - alpha * v
        MF::Saxpy(r,  -alpha, v, 0, 0, ncomp, nghost); // r -= alpha * v
        PrintDots(Half_Iter,"Dots(sol,sol)",sol,sol);
        PrintDots(Half_Iter,"Dots(r  ,r  )",r  ,r  );

        rnorm = norm_inf(r);

        if ( verbose > 2 && ParallelDescriptor::IOProcessor() )
        {
            amrex::Print() << "MLCGSolver_BiCGStab: Half Iter "
                           << std::setw(11) << iter
                           << " rel. err. "
                           << rnorm/(rnorm0) << '\n';
        }

        if ( rnorm < eps_rel*rnorm0 || rnorm < eps_abs ) { break; }

        // sh.LocalCopy(s,0,0,ncomp,nghost);
        PrintDots(Iteration,"Dots(r,r) before Lp.apply",r,r);
        Lp.apply(amrlev, mglev, t, r, MLLinOpT<MF>::BCMode::Homogeneous, MLLinOpT<MF>::StateMode::Correction);
        PrintDots(Iteration,"Dots(r,r) after  Lp.apply",r,r);

        PrintDots(Iteration,"Dots(t,t) before Lp.normalize",t,t);
        Lp.normalize(amrlev, mglev, t);
        PrintDots(Iteration,"Dots(t,t) after  Lp.normalize",t,t);
        //
        // This is a little funky.  I want to elide one of the reductions
        // in the following two dotxy()s.  We do that by calculating the "local"
        // values and then reducing the two local values at the same time.
        //
        RT tvals[2] = { dotxy(t,t,true), dotxy(t,r,true) };

        BL_PROFILE_VAR("MLCGSolver::ParallelAllReduce", blp_par);
        ParallelAllReduce::Sum(tvals,2,Lp.BottomCommunicator());
        BL_PROFILE_VAR_STOP(blp_par);

        if ( tvals[0] == RT(0.0) )
        {
            ret = 3; break;
        }
        omega = tvals[1]/tvals[0];
        PrintDots(Iteration,"Dots(r,t)",t,r);

        MF::Saxpy(sol, omega, r, 0, 0, ncomp, nghost); // sol += omega * r
        // MF::LinComb(r, RT(1.0), s, 0, -omega, t, 0, 0, ncomp, nghost); // r = s - omega * t
        MF::Saxpy(r,  -omega, t,  0, 0, ncomp, nghost); // r -= omega * t
        PrintDots(Iteration,"Dots(sol,sol)",sol,sol);
        PrintDots(Iteration,"Dots(r  ,r  )",r,r);

        rnorm = norm_inf(r);

        if ( verbose > 2 )
        {
            amrex::Print() << "MLCGSolver_BiCGStab: Iteration "
                           << std::setw(11) << iter
                           << " rel. err. "
                           << rnorm/(rnorm0) << '\n';
        }

        if ( rnorm < eps_rel*rnorm0 || rnorm < eps_abs ) { break; }

        if ( omega == RT(0.0) )
        {
            ret = 4; break;
        }
        RT rho_1 = rho;
        rho = dotxy(rh,r);
        PrintDots(Post_Calc,"Dots(rh,r)",rh,r);

        beta = (rho/rho_1)*(alpha/omega);
        // beta = rho/(rhTv*omega);
        // amrex::Print() << Iteration << " beta: " << beta << '\n';

        PrintDots(Post_Calc, "Dots(p,p) before     ",p,p);
        MF::Saxpy(p, -omega, v, 0, 0, ncomp, nghost); // p += -omega*v
        PrintDots(Post_Calc, "Dots(p,p) after Saxpy",p,p);
        MF::Xpay(p, beta, r, 0, 0, ncomp, nghost); // p = r + beta*p
        PrintDots(Post_Calc, "Dots(p,p) after  Xpay",p,p);
    }

    if ( verbose > 0 )
    {
        amrex::Print() << "MLCGSolver_BiCGStab: Final: Iteration "
                       << std::setw(4) << iter
                       << " rel. err. "
                       << rnorm/(rnorm0) << '\n';
    }

    if ( ret == 0 && rnorm > eps_rel*rnorm0 && rnorm > eps_abs)
    {
        if ( verbose > 0 && ParallelDescriptor::IOProcessor() ) {
            amrex::Warning("MLCGSolver_BiCGStab:: failed to converge!");
        }
        ret = 8;
    }

    std::string IterationFinal = "MLCGSolver_BiCGStab: Final: Iteration    " + std::to_string(iter);

    if ( ( ret == 0 || ret == 8 ) && (rnorm < rnorm0) )
    {
        // keep updated solution
        amrex::Print() << "MLCGSolver_BiCGstab: Keep new solution.\n";
        sol.LocalAdd(sorig, 0, 0, ncomp, nghost);
        PrintDots(IterationFinal,"Dots(sol,sol)",sol,sol);
    }
    else
    {
        // reset solution to original copy
        amrex::Print() << "MLCGSolver_BiCGstab: Reset solution.\n";
        sol.setVal(RT(0.0));
        sol.LocalAdd(sorig, 0, 0, ncomp, nghost);
        PrintDots(IterationFinal,"Dots(sol,sol)",sol,sol);
    }

    return ret;
}

template <typename MF>
int
MLCGSolverT<MF>::solve_cg (MF& sol, const MF& rhs, RT eps_rel, RT eps_abs)
{
    BL_PROFILE("MLCGSolver::cg");

    const int ncomp = sol.nComp();

    const BoxArray& ba = sol.boxArray();
    const DistributionMapping& dm = sol.DistributionMap();
    const auto& factory = sol.Factory();

    MF p(ba, dm, ncomp, sol.nGrowVect(), MFInfo(), factory);
    p.setVal(RT(0.0));

    MF sorig(ba, dm, ncomp, nghost, MFInfo(), factory);
    MF r    (ba, dm, ncomp, nghost, MFInfo(), factory);
    MF z    (ba, dm, ncomp, nghost, MFInfo(), factory);
    MF q    (ba, dm, ncomp, nghost, MFInfo(), factory);

    sorig.LocalCopy(sol,0,0,ncomp,nghost);

    Lp.correctionResidual(amrlev, mglev, r, sol, rhs, MLLinOpT<MF>::BCMode::Homogeneous);

    sol.setVal(RT(0.0));

    RT       rnorm    = norm_inf(r);
    const RT rnorm0   = rnorm;

    if ( verbose > 0 )
    {
        amrex::Print() << "MLCGSolver_CG: Initial error (error0) :        " << rnorm0 << '\n';
    }

    RT rho_1 = 0;
    int  ret = 0;
    iter = 1;

    if ( rnorm0 == 0 || rnorm0 < eps_abs )
    {
        if ( verbose > 0 ) {
            amrex::Print() << "MLCGSolver_CG: niter = 0,"
                           << ", rnorm = " << rnorm
                           << ", eps_abs = " << eps_abs << std::endl;
        }
        return ret;
    }

    for (; iter <= maxiter; ++iter)
    {
        z.LocalCopy(r,0,0,ncomp,nghost);

        RT rho = dotxy(z,r);

        if ( rho == 0 )
        {
            ret = 1; break;
        }
        if (iter == 1)
        {
            p.LocalCopy(z,0,0,ncomp,nghost);
        }
        else
        {
            RT beta = rho/rho_1;
            MF::Xpay(p, beta, z, 0, 0, ncomp, nghost); // p = z + beta * p
        }
        Lp.apply(amrlev, mglev, q, p, MLLinOpT<MF>::BCMode::Homogeneous, MLLinOpT<MF>::StateMode::Correction);

        RT alpha;
        RT pw = dotxy(p,q);
        if ( pw != RT(0.0))
        {
            alpha = rho/pw;
        }
        else
        {
            ret = 1; break;
        }

        if ( verbose > 2 )
        {
            amrex::Print() << "MLCGSolver_cg:"
                           << " iter " << iter
                           << " rho " << rho
                           << " alpha " << alpha << '\n';
        }
        MF::Saxpy(sol, alpha, p, 0, 0, ncomp, nghost); // sol += alpha * p
        MF::Saxpy(r, -alpha, q, 0, 0, ncomp, nghost); // r += -alpha * q
        rnorm = norm_inf(r);

        if ( verbose > 2 )
        {
            amrex::Print() << "MLCGSolver_cg:       Iteration"
                           << std::setw(4) << iter
                           << " rel. err. "
                           << rnorm/(rnorm0) << '\n';
        }

        if ( rnorm < eps_rel*rnorm0 || rnorm < eps_abs ) { break; }

        rho_1 = rho;
    }

    if ( verbose > 0 )
    {
        amrex::Print() << "MLCGSolver_cg: Final Iteration"
                       << std::setw(4) << iter
                       << " rel. err. "
                       << rnorm/(rnorm0) << '\n';
    }

    if ( ret == 0 &&  rnorm > eps_rel*rnorm0 && rnorm > eps_abs )
    {
        if ( verbose > 0 && ParallelDescriptor::IOProcessor() ) {
            amrex::Warning("MLCGSolver_cg: failed to converge!");
        }
        ret = 8;
    }

    if ( ( ret == 0 || ret == 8 ) && (rnorm < rnorm0) )
    {
        sol.LocalAdd(sorig, 0, 0, ncomp, nghost);
    }
    else
    {
        sol.setVal(RT(0.0));
        sol.LocalAdd(sorig, 0, 0, ncomp, nghost);
    }

    return ret;
}

template <typename MF>
auto
MLCGSolverT<MF>::dotxy (const MF& r, const MF& z, bool local) -> RT
{
    BL_PROFILE_VAR_NS("MLCGSolver::ParallelAllReduce", blp_par);
    if (!local) { BL_PROFILE_VAR_START(blp_par); }
    RT result = Lp.xdoty(amrlev, mglev, r, z, local);
    if (!local) { BL_PROFILE_VAR_STOP(blp_par); }
    return result;
}

template <typename MF>
auto
MLCGSolverT<MF>::norm_inf (const MF& res, bool local) -> RT
{
    int ncomp = res.nComp();
    RT result = res.norminf(0,ncomp,IntVect(0),true);
    if (!local) {
        BL_PROFILE("MLCGSolver::ParallelAllReduce");
        ParallelAllReduce::Max(result, Lp.BottomCommunicator());
    }
    return result;
}

// template <typename MF>
// auto
// MLCGSolverT<MF>::Dots (const MF& x, const MF& y) -> Array<RT, 3>
// {
//     RT dot1 = dotxy(x,y);
//     RT dot2 = Dot(x,0,y,0,x.nComp(),nghost);
//     RT dot3 = Dot(x,0,y,0,x.nComp(),ngrow);
//     return {dot1, dot2, dot3};
// }

template <typename MF>
void
MLCGSolverT<MF>::PrintDots (const std::string& title, const char* name, const MF& x, const MF& y)
{
    amrex::Print() << title << " " << name << " -> ";
    amrex::Print() << "dotxy: " ;
    amrex::Print() << dotxy(x,y);
    amrex::Print() << ", Dot(...,nghost): ";
    amrex::Print() << Dot(x,0,y,0,x.nComp(),nghost);
    amrex::Print() << ", Dot(...,ngrow): ";
    if (x.nGrowVect().allGE(ngrow) && y.nGrowVect().allGE(ngrow)) {
        amrex::Print() << Dot(x,0,y,0,x.nComp(),ngrow) << "\n";
    }
    else {
        amrex::Print() << "N/A" << "\n";
    }
}


using MLCGSolver = MLCGSolverT<MultiFab>;

}

#endif /*_CGSOLVER_H_*/
