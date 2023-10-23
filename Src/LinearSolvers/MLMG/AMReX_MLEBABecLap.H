#ifndef AMREX_MLEBABECLAP_H_
#define AMREX_MLEBABECLAP_H_
#include <AMReX_Config.H>

#include <AMReX_EBFabFactory.H>
#include <AMReX_MLCellABecLap.H>
#include <AMReX_Array.H>
#include <limits>

namespace amrex {

/**
 * Multi-Level cell-centered embedded boundary AB coefficients Laplacian solver. 
 *
 * (A * alpha - B * (del dot beta grad)) phi = rhs
 *
 * where `A` and `B` are scalar constants, and `alpha` and `beta` are scalar fields.
 *
 * Inherits from MLCellABecLapT<MultiFab>
 */
class MLEBABecLap
    : public MLCellABecLap
{
public:

    MLEBABecLap () = default;
    MLEBABecLap (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info,
                 const Vector<EBFArrayBoxFactory const*>& a_factory,
                 int a_ncomp = 1);

    ~MLEBABecLap () override;

    MLEBABecLap (const MLEBABecLap&) = delete;
    MLEBABecLap (MLEBABecLap&&) = delete;
    MLEBABecLap& operator= (const MLEBABecLap&) = delete;
    MLEBABecLap& operator= (MLEBABecLap&&) = delete;

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info,
                 const Vector<EBFArrayBoxFactory const*>& a_factory,
                 int a_ncomp = 1);

    void setPhiOnCentroid ();

    /**
     * Set scalar constants A and B in the equation:
     * (A \alpha - B \nabla \cdot \beta \nabla ) \phi = f
     * for the Multi-Level Embedded Boundary AB Laplacian Solver.
     */
    void setScalars (Real a, Real b);

    /**
     * \brief Sets alpha as a scalar field to values from a single component
     * mutlifab.
     *
     * \param[in] amrlev The level of the multifab for the solver, with
     *                    \p amrlev = 0 always being the lowest level in the
     *                    AMR hierarchy represented in the solve.
     * \param[in] alpha  Multifab of alpha values.
     */
    void setAlphaCoeffs (int amrlev, const MultiFab& alpha);

    /**
     * \brief Sets alpha as a single scalar constant value across
     * the mutlifab.
     *
     * \param[in] amrlev The level of the multifab for the solver, with
     *                    \p amrlev = 0 always being the lowest level in the
     *                    AMR hierarchy represented in the solve.
     * \param[in] alpha  Single scalar value to populate across mutlifab.
     */
    void setAlphaCoeffs (int amrlev, Real alpha);

    /**
     * \brief Sets beta as a field with location specified by \p a_beta_loc to be
     * the values defined in the supplied multifabs (one for each space dimension).
     *
     * \param amrlev     The level of the multifab for the solver, with
     *                    \p amrlev = 0 always being the lowest level in the
     *                    AMR hierarchy represented in the solve.
     * \param beta       Array of Multifabs of beta values.
     * \param a_beta_loc The grid location of beta
     */
    void setBetaCoeffs (int amrlev, const Array<MultiFab const*,AMREX_SPACEDIM>& beta,
                     Location a_beta_loc);

    /**
     * \brief Sets beta as a face-centered field to be the values defined
     * in the supplied multifabs (one for each space dimension).
     *
     * \param[in] amrlev The level of the multifab for the solver, with
     *                    \p amrlev = 0 always being the lowest level in the
     *                    AMR hierarchy represented in the solve.
     * \param[in] beta   Array of Multifabs of beta values.
     */
    void setBetaCoeffs (int amrlev, const Array<MultiFab const*,AMREX_SPACEDIM>& beta)
        {setBetaCoeffs (amrlev, beta, Location::FaceCenter);}

    /**
     * \brief Sets beta as a single scalar constant value across
     * the mutlifabs (one for each dimension).
     *
     * \param[in] amrlev The level of the multifab for the solver, with
     *                    \p amrlev = 0 always being the lowest level in the
     *                    AMR hierarchy represented in the solve.
     * \param[in] beta   Single scalar value to populate across mutlifabs.
     */
    void setBetaCoeffs (int amrlev, Real beta);

    /**
     * \brief Set each beta component to a constant value corresponding to the
     * respective component of the supplied vector.
     *
     * \param[in] amrlev The level of the multifab for the solver, with
     *                    \p amrlev = 0 always being the lowest level in the
     *                    AMR hierarchy represented in the solve.
     * \param[in] beta   Vector of constant values.
     */
    void setBetaCoeffs (int amrlev, Vector<Real> const& beta);

    // Equivalent to setAlphaCoeffs. Use of this method is discouraged to avoid confusion between A and alpha
    [[deprecated("Use setAlphaCoeffs instead.")]]
    void setACoeffs (int amrlev, const MultiFab& alpha) { setAlphaCoeffs(amrlev,alpha); }
    // Equivalent to setAlphaCoeffs. Use of this method is discouraged to avoid confusion between A and alpha
    [[deprecated("Use setAlphaCoeffs instead.")]]
    void setACoeffs (int amrlev, Real alpha) { setAlphaCoeffs(amrlev,alpha); }
    // Equivalent to setBetaCoeffs. Use of this method is discouraged to avoid confusion between B and beta
    [[deprecated("Use setBetaCoeffs instead.")]]
    void setBCoeffs (int amrlev, const Array<MultiFab const*,AMREX_SPACEDIM>& beta,
                     Location a_beta_loc) { setBetaCoeffs(amrlev,beta,a_beta_loc); }
    // Equivalent to setBetaCoeffs. Use of this method is discouraged to avoid confusion between B and beta
    [[deprecated("Use setBetaCoeffs instead.")]]
    void setBCoeffs (int amrlev, const Array<MultiFab const*,AMREX_SPACEDIM>& beta)
        {setBetaCoeffs (amrlev, beta);}
    // Equivalent to setBetaCoeffs. Use of this method is discouraged to avoid confusion between B and beta
    [[deprecated("Use setBetaCoeffs instead.")]]
    void setBCoeffs (int amrlev, Real beta) { setBetaCoeffs(amrlev,beta); }
    // Equivalent to setBetaCoeffs. Use of this method is discouraged to avoid confusion between B and beta
    [[deprecated("Use setBetaCoeffs instead.")]]
    void setBCoeffs (int amrlev, Vector<Real> const& beta) { setBetaCoeffs(amrlev,beta); }

    // Tells the solver that EB boundaries have Dirichlet bc's specified by "phi"
    void setEBDirichlet (int amrlev, const MultiFab& phi, const MultiFab& beta);
    void setEBDirichlet (int amrlev, const MultiFab& phi, Real beta);
    void setEBDirichlet (int amrlev, const MultiFab& phi, Vector<Real> const& beta);

    // Tells the solver that EB boundaries have homogeneous Dirichlet bc's
    void setEBHomogDirichlet (int amrlev,                      const MultiFab& beta);
    void setEBHomogDirichlet (int amrlev,                      Real beta);
    void setEBHomogDirichlet (int amrlev,                      Vector<Real> const& beta);

    int getNComp () const override { return m_ncomp; }

    bool needsUpdate () const override {
        return (m_needs_update || MLCellABecLap::needsUpdate());
    }
    void update () override;

    std::unique_ptr<FabFactory<FArrayBox> > makeFactory (int amrlev, int mglev) const final;

    bool isCrossStencil () const override { return false; }

    void applyBC (int amrlev, int mglev, MultiFab& in, BCMode bc_mode, StateMode s_mode,
                          const MLMGBndry* bndry=nullptr, bool skip_fillboundary=false) const final;
    void apply (int amrlev, int mglev, MultiFab& out, MultiFab& in, BCMode bc_mode,
                        StateMode s_mode, const MLMGBndry* bndry=nullptr) const override;
    void compGrad (int amrlev, const Array<MultiFab*,AMREX_SPACEDIM>& grad,
                           MultiFab& sol, Location loc) const final;

    void prepareForSolve () override;
    bool isSingular (int amrlev) const override { return m_is_singular[amrlev]; }
    bool isBottomSingular () const override { return m_is_singular[0]; }
    void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const final;
    void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs, int redblack) const final;
    void FFlux (int amrlev, const MFIter& mfi,
                        const Array<FArrayBox*,AMREX_SPACEDIM>& flux,
                        const FArrayBox& sol, Location loc,
                        int face_only=0) const final;

    void normalize (int amrlev, int mglev, MultiFab& mf) const final;

    Real getAScalar () const final { return m_a_scalar; }
    Real getBScalar () const final { return m_b_scalar; }
    MultiFab const* getACoeffs (int amrlev, int mglev) const final
        { return &(m_alpha_coeffs[amrlev][mglev]); }
    Array<MultiFab const*,AMREX_SPACEDIM> getBCoeffs (int amrlev, int mglev) const final
        { return amrex::GetArrOfConstPtrs(m_beta_coeffs[amrlev][mglev]); }

    std::unique_ptr<MLLinOp> makeNLinOp (int /*grid_size*/) const final {
        amrex::Abort("MLABecLaplacian::makeNLinOp: Not implemented");
        return std::unique_ptr<MLLinOp>{};
    }

    void restriction (int, int, MultiFab& crse, MultiFab& fine) const final;

    void interpolation (int amrlev, int fmglev, MultiFab& fine, const MultiFab& crse) const final;

    void averageDownSolutionRHS (int camrlev, MultiFab& crse_sol, MultiFab& crse_rhs,
                                         const MultiFab& fine_sol, const MultiFab& fine_rhs) final;

    void getEBFluxes (const Vector<MultiFab*>& a_flux,
                              const Vector<MultiFab*>& a_sol) const override;

#if defined(AMREX_USE_HYPRE) && (AMREX_SPACEDIM > 1)
    [[nodiscard]] std::unique_ptr<Hypre> makeHypre (Hypre::Interface hypre_interface) const override;
#endif

#ifdef AMREX_USE_PETSC
    [[nodiscard]] std::unique_ptr<PETScABecLap> makePETSc () const override;
#endif

protected:

    int m_ncomp = 1;

    bool m_needs_update = true;

    Location m_beta_loc; // Location of coefficients: face centers or face centroids
    Location m_phi_loc;  // Location of solution variable: cell centers or cell centroids

    Real m_a_scalar = std::numeric_limits<Real>::quiet_NaN();
    Real m_b_scalar = std::numeric_limits<Real>::quiet_NaN();
    Vector<Vector<MultiFab> > m_alpha_coeffs;
    Vector<Vector<Array<MultiFab,AMREX_SPACEDIM> > > m_beta_coeffs;
    Vector<Vector<iMultiFab> > m_cc_mask;

    Vector<std::unique_ptr<MultiFab> > m_eb_phi;
    Vector<Vector<std::unique_ptr<MultiFab> > > m_eb_beta_coeffs;

    Vector<int> m_is_singular;

    mutable int m_is_eb_inhomog;

    //
    // functions
    //
    bool isEBDirichlet   () const noexcept { return m_eb_phi[0] != nullptr; }

    void averageDownCoeffsSameAmrLevel (int amrlev, Vector<MultiFab>& a,
                                        Vector<Array<MultiFab,AMREX_SPACEDIM> >& b,
                                        const Vector<MultiFab*>& b_eb);
    void averageDownCoeffs ();
    void averageDownCoeffsToCoarseAmrLevel (int flev);
};

}

#endif
