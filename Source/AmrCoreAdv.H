#ifndef AmrCoreAdv_H_
#define AmrCoreAdv_H_

#include <AMReX_Config.H>

#include <Kernels.H>
#include <AmrCoreAdvUtil.H>

#ifdef AMREX_PARTICLES
#   include <AMReX_AmrParticles.H>
#   include <AMReX_Particles.H>
#endif

#include <AMReX_AmrCore.H>
#include <AMReX_BCRec.H>
#include <AMReX_ErrorList.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>

#ifdef AMREX_USE_OMP
#   include <omp.h>
#endif

#include <limits>
#include <memory>
#include <string>
#include <variant>

class AmrCoreAdv
    : public amrex::AmrCore
{
public:

    ////////////////
    // public member functions

    // constructor - reads in parameters from inputs file
    //             - sizes multilevel arrays and data structures
    AmrCoreAdv ();
    ~AmrCoreAdv () override;

    AmrCoreAdv (AmrCoreAdv&& rhs) noexcept = default;
    AmrCoreAdv& operator= (AmrCoreAdv&& rhs) noexcept = default;

    AmrCoreAdv (const AmrCoreAdv& rhs) = delete;
    AmrCoreAdv& operator= (const AmrCoreAdv& rhs) = delete;

    // advance solution to final time
    void Evolve ();

    // initializes multilevel data
    void InitData ();

    // Make a new level using provided BoxArray and DistributionMapping and
    // fill with interpolated coarse level data.
    // overrides the pure virtual function in AmrCore
    void MakeNewLevelFromCoarse (int lev, amrex::Real time, const amrex::BoxArray& ba,
                                 const amrex::DistributionMapping& dm) override;

    // Remake an existing level using provided BoxArray and DistributionMapping and
    // fill with existing fine and coarse data.
    // overrides the pure virtual function in AmrCore
    void RemakeLevel (int lev, amrex::Real time, const amrex::BoxArray& ba,
                      const amrex::DistributionMapping& dm) override;

    // Delete level data
    // overrides the pure virtual function in AmrCore
    void ClearLevel (int lev) override;

    // Make a new level from scratch using provided BoxArray and DistributionMapping.
    // Only used during initialization.
    // overrides the pure virtual function in AmrCore
    void MakeNewLevelFromScratch (int lev, amrex::Real time, const amrex::BoxArray& ba,
                                  const amrex::DistributionMapping& dm) override;

    // tag all cells for refinement
    // overrides the pure virtual function in AmrCore
    void ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;

    // Advance phi at a single level for a single time step, update flux registers
    void AdvancePhiAtLevel (int lev, amrex::Real time, amrex::Real dt_lev, int iteration, int ncycle);

    // Define the advection velocity as the curl of a scalar field
    void DefineVelocityAtLevel (int lev, amrex::Real time);

    // compute dt from CFL considerations
    amrex::Real EstTimeStep (int lev, amrex::Real time);

#ifdef AMREX_PARTICLES
    amrex::AmrTracerParticleContainer* theTracerPC () { return TracerPC.get(); }
#endif

    template <typename MF>
    void MakeNewLevelFromCoarse (int lev, amrex::Real time,
				 const amrex::BoxArray& ba,
				 const amrex::DistributionMapping& dm);

    template <typename MF>
    void RemakeLevel (int lev, amrex::Real time,
		      const amrex::BoxArray& ba,
		      const amrex::DistributionMapping& dm);

    template <typename MF>
    void MakeNewLevelFromScratch (int lev, amrex::Real time,
				  const amrex::BoxArray& ba,
				  const amrex::DistributionMapping& dm);

    template <typename MF>
    void ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real time, int /*ngrow*/);

    template <typename MF>
    void AdvancePhiAtLevel (int lev, amrex::Real time, amrex::Real dt_lev, int iteration, int ncycle);

    template <typename MF>
    void DefineVelocityAtLevel (int lev, amrex::Real time);

private:

    ////////////////
    // private member functions

    // read in some parameters from inputs file
    void ReadParameters();

    // set covered coarse cells to be the average of overlying fine cells
    void AverageDown ();

    // more flexible version of AverageDown() that lets you average down across multiple levels
    void AverageDownTo (int crse_lev);

    // compute a new multifab by coping in phi from valid region and filling ghost cells
    // works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
    template <typename MF>
    void FillPatch (int lev, amrex::Real time, MF& mf, int icomp, int ncomp);

    // fill an entire multifab by interpolating from the coarser level
    // this comes into play when a new level of refinement appears
    template <typename MF>
    void FillCoarsePatch (int lev, amrex::Real time, MF& mf, int icomp, int ncomp);

    template <typename MF>
    // Pack pointers to phi_old and/or phi_new and associated times.
    void GetData (int lev, amrex::Real time, amrex::Vector<MF>& data,
                  amrex::Vector<amrex::Real>& datatime);

    // Advance a level by dt - includes a recursive call for finer levels
    void timeStepWithSubcycling (int lev, amrex::Real time, int iteration);

    // a wrapper for EstTimeStep
    void ComputeDt ();

    // get plotfile name
    [[nodiscard]] std::string PlotFileName (int lev) const;

    // put together an array of multifabs for writing
    [[nodiscard]] amrex::Vector<amrex::MultiFab> PlotFileMF () const;

    // set plotfile variables names
    [[nodiscard]] static amrex::Vector<std::string> PlotFileVarNames ();

    // write plotfile to disk
    void WritePlotFile () const;

    // write checkpoint file to disk
    void WriteCheckpointFile () const;

    // read checkpoint file from disk
    void ReadCheckpointFile ();

    [[nodiscard]] bool DoublePrecisionOnLevel (int lev) const {
	return (lev < 0 && lev >= double_precision.size()) || double_precision[lev];
    }

    template <typename MF>
    bool SamePrecision (int lev) const {
	if (DoublePrecisionOnLevel(lev)) {
	    return std::is_same_v<MF,amrex::MultiFab>;
	} else {
	    return std::is_same_v<MF,amrex::fMultiFab>;
	}
    }

    ////////////////
    // private data members

    amrex::Vector<int> istep;      // which step?
    amrex::Vector<int> nsubsteps;  // how many substeps on each level?

    // keep track of old time, new time, and time step at each level
    amrex::Vector<amrex::Real> t_new;
    amrex::Vector<amrex::Real> t_old;
    amrex::Vector<amrex::Real> dt;

    // array of multifabs to store the solution at each level of refinement
    // after advancing a level we use "swap".
    amrex::Vector<std::variant<amrex::MultiFab,amrex::fMultiFab>> phi_new;
    amrex::Vector<std::variant<amrex::MultiFab,amrex::fMultiFab>> phi_old;

    // this is essentially a 2*DIM integer array storing the physical boundary
    // condition types at the lo/hi walls in each direction
    amrex::Vector<amrex::BCRec> bcs;  // 1-component

    // stores fluxes at coarse-fine interface for synchronization
    // this will be sized "nlevs_max+1"
    // NOTE: the flux register associated with flux_reg[lev] is associated
    // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
    // therefore flux_reg[0] and flux_reg[nlevs_max] are never actually
    // used in the reflux operation
    amrex::Vector<std::unique_ptr<amrex::FluxRegister> > flux_reg;

    // Velocity on all faces at all levels
    amrex::Vector< amrex::Array<std::variant<amrex::MultiFab,amrex::fMultiFab>,
				AMREX_SPACEDIM> > facevel;

    ////////////////
    // runtime parameters

    // maximum number of steps and stop time
    int max_step = std::numeric_limits<int>::max();
    amrex::Real stop_time = std::numeric_limits<amrex::Real>::max();

    // if >= 0 we restart from a checkpoint
    std::string restart_chkfile;

    // advective cfl number - dt = cfl*dx/umax
    amrex::Real cfl = 0.7;

    // how often each level regrids the higher levels of refinement
    // (after a level advances that many time steps)
    int regrid_int = 2;

    // hyperbolic refluxing as part of multilevel synchronization
    int do_reflux = 1;

    // do we subcycle in time?
    int do_subcycle = 1;

    // plotfile prefix and frequency
    std::string plot_file {"plt"};
    int plot_int = -1;

    // checkpoint prefix and frequency
    std::string chk_file {"chk"};
    int chk_int = -1;

    //number of ghost cells on facevel
    int nghost = 2;

    amrex::Vector<int> double_precision;

    std::string error_fn;
    amrex::Vector<amrex::AMRErrorTag> error_tag;

#ifdef AMREX_PARTICLES
    void init_particles ();
    int       do_tracers = 0;
    std::unique_ptr<amrex::AmrTracerParticleContainer> TracerPC;
#endif
};

template <typename MF>
void AmrCoreAdv::MakeNewLevelFromCoarse (int lev, amrex::Real time,
					 const amrex::BoxArray& ba,
					 const amrex::DistributionMapping& dm)
{
    using namespace amrex;

    const int ncomp = std::get<MF>(phi_new[lev-1]).nComp();
    const int ng = std::get<MF>(phi_new[lev-1]).nGrow();

    phi_new[lev] = MF(ba, dm, ncomp, ng);
    phi_old[lev] = MF(ba, dm, ncomp, ng);

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    // This clears the old MultiFab and allocates the new one
    for (int idim = 0; idim < AMREX_SPACEDIM; idim++)
    {
        facevel[lev][idim] = MF(amrex::convert(ba,IntVect::TheDimensionVector(idim)), dm, 1, nghost);
    }

    if (lev > 0 && do_reflux) {
        flux_reg[lev] = std::make_unique<FluxRegister>(ba, dm, refRatio(lev-1), lev, ncomp);
    }

    if (DoublePrecisionOnLevel(lev)) {
	FillCoarsePatch(lev, time, std::get<MultiFab>(phi_new[lev]), 0, ncomp);
    } else {
	FillCoarsePatch(lev, time, std::get<fMultiFab>(phi_new[lev]), 0, ncomp);
    }
}

template <typename MF>
void AmrCoreAdv::RemakeLevel (int lev, amrex::Real time,
			      const amrex::BoxArray& ba,
			      const amrex::DistributionMapping& dm)
{
    using namespace amrex;

    const int ncomp = std::get<MF>(phi_new[lev]).nComp();
    const int ng = std::get<MF>(phi_new[lev]).nGrow();

    MF new_state(ba, dm, ncomp, ng);
    MF old_state(ba, dm, ncomp, ng);

    FillPatch(lev, time, new_state, 0, ncomp);

    phi_new[lev] = std::move(new_state);
    phi_old[lev] = std::move(old_state);

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    // This clears the old MultiFab and allocates the new one
    for (int idim = 0; idim < AMREX_SPACEDIM; idim++)
    {
        facevel[lev][idim] = MF(amrex::convert(ba,IntVect::TheDimensionVector(idim)), dm, 1, nghost);
    }

    if (lev > 0 && do_reflux) {
        flux_reg[lev] = std::make_unique<FluxRegister>(ba, dm, refRatio(lev-1), lev, ncomp);
    }
}

template <typename MF>
void AmrCoreAdv::MakeNewLevelFromScratch (int lev, amrex::Real time,
					  const amrex::BoxArray& ba,
					  const amrex::DistributionMapping& dm)
{
    using namespace amrex;

    const int ncomp = 1;
    const int ng = 0;

    phi_new[lev] = MF(ba, dm, ncomp, ng);
    phi_old[lev] = MF(ba, dm, ncomp, ng);

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    // This clears the old MultiFab and allocates the new one
    for (int idim = 0; idim < AMREX_SPACEDIM; idim++)
    {
	facevel[lev][idim] = MF(amrex::convert(ba,IntVect::TheDimensionVector(idim)), dm, 1, nghost);
    }

    if (lev > 0 && do_reflux) {
        flux_reg[lev] = std::make_unique<FluxRegister>(ba, dm, refRatio(lev-1), lev, ncomp);
    }

    MF& state = std::get<MF>(phi_new[lev]);

    using RT = typename MF::value_type;

    GpuArray<RT,AMREX_SPACEDIM> problo{AMREX_D_DECL(RT(geom[lev].ProbLo(0)),
						    RT(geom[lev].ProbLo(1)),
						    RT(geom[lev].ProbLo(2)))};
    GpuArray<RT,AMREX_SPACEDIM> dx{AMREX_D_DECL(RT(geom[lev].CellSize(0)),
						RT(geom[lev].CellSize(1)),
					        RT(geom[lev].CellSize(2)))};

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(state,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        auto const& fab = state[mfi].array();
        const Box& box = mfi.tilebox();

        amrex::ParallelFor(box,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            initdata(i, j, k, fab, problo, dx);
        });
    }
}

template <typename MF>
void AmrCoreAdv::ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real time, int)
{
    using namespace amrex;
    using RT = typename MF::value_type;

    static bool first = true;
    static Vector<RT> phierr;

    // only do this during the first call to ErrorEst
    if (first)
    {
        first = false;
        // read in an array of "phierr", which is the tagging threshold
        // in this example, we tag values of "phi" which are greater than phierr
        // for that particular level
        // in subroutine state_error, you could use more elaborate tagging, such
        // as more advanced logical expressions, or gradients, etc.
        ParmParse pp("adv");
        int n = pp.countval("phierr");
        if (n > 0) {
            pp.getarr("phierr", phierr, 0, n);
        }
    }

    const int clearval = TagBox::CLEAR;
    const int   tagval = TagBox::SET;

    if (lev < phierr.size())
    {
        const auto& state = std::get<MF>(phi_new[lev]);

#ifdef AMREX_USE_OMP
#pragma omp parallel if(Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(state,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx  = mfi.tilebox();
            const auto statefab = state.array(mfi);
            const auto tagfab  = tags.array(mfi);
            auto phierror = phierr[lev];

            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                state_error(i, j, k, tagfab, statefab, phierror, tagval);
            });
        }
    }

    for (auto const& etag : error_tag) {
        etag(tags, nullptr, clearval, tagval, time, lev, geom[lev]);
    }
}

template <typename MF>
void
AmrCoreAdv::GetData (int lev, amrex::Real time, amrex::Vector<MF>& data,
		     amrex::Vector<amrex::Real>& datatime)
{
    using namespace amrex;

    data.clear();
    datatime.clear();

    using OMF = std::conditional_t<std::is_same_v<MF,MultiFab>, fMultiFab, MultiFab>;

    if (amrex::almostEqual(time, t_new[lev], 5))
    {
	if (SamePrecision<MF>(lev)) {
	    auto const& mf = std::get<MF>(phi_new[lev]);
	    data.push_back(MF(mf, amrex::make_alias, 0, mf.nComp()));
	} else {
	    auto const& mf = std::get<OMF>(phi_new[lev]);
	    data.push_back(amrex::cast<MF>(mf));
	}
        datatime.push_back(t_new[lev]);
    }
    else if (amrex::almostEqual(time, t_old[lev], 5))
    {
	if (SamePrecision<MF>(lev)) {
	    auto const& mf = std::get<MF>(phi_old[lev]);
	    data.push_back(MF(mf, amrex::make_alias, 0, mf.nComp()));
	} else {
	    auto const& mf = std::get<OMF>(phi_old[lev]);
	    data.push_back(amrex::cast<MF>(mf));
	}
        datatime.push_back(t_old[lev]);
    }
    else
    {
	if (SamePrecision<MF>(lev)) {
	    auto const& mf_old = std::get<MF>(phi_old[lev]);
	    auto const& mf_new = std::get<MF>(phi_new[lev]);
	    data.push_back(MF(mf_old, amrex::make_alias, 0, mf_old.nComp()));
	    data.push_back(MF(mf_new, amrex::make_alias, 0, mf_old.nComp()));
	} else {
	    auto const& mf_old = std::get<OMF>(phi_old[lev]);
	    auto const& mf_new = std::get<OMF>(phi_new[lev]);
	    data.push_back(amrex::cast<MF>(mf_old));
	    data.push_back(amrex::cast<MF>(mf_new));
	}
        datatime.push_back(t_old[lev]);
        datatime.push_back(t_new[lev]);
    }
}

// compute a new multifab by coping in phi from valid region and filling ghost cells
// works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
template <typename MF>
void AmrCoreAdv::FillPatch (int lev, amrex::Real time, MF& mf, int icomp, int ncomp)
{
    using namespace amrex;

    if (lev == 0)
    {
        Vector<MF> smf0;
        Vector<Real> stime;
        GetData(0, time, smf0, stime);
	auto smf = GetVecOfPtrs(smf0);

	Util::PhysBCFunctNoOp<MF> physbc{};
	amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, icomp, ncomp,
				    geom[lev], physbc, 0);
    }
    else
    {
        Vector<MultiFab> cmf0, fmf0;
        Vector<Real> ctime, ftime;
        GetData(lev-1, time, cmf0, ctime);
        GetData(lev  , time, fmf0, ftime);
	Vector<MultiFab*> cmf = GetVecOfPtrs(cmf0);
	Vector<MultiFab*> fmf = GetVecOfPtrs(fmf0);

        Interpolater* mapper = &cell_cons_interp;
	Util::PhysBCFunctNoOp<MultiFab> cphysbc{};
	Util::PhysBCFunctNoOp<MultiFab> fphysbc{};

	MultiFab mftmp;
	if constexpr (std::is_same_v<MultiFab,MF>) {
	    mftmp = MultiFab(mf, amrex::make_alias, 0, mf.nComp());
	} else {
	    mftmp.define(mf.boxArray(), mf.DistributionMap(), mf.nComp(), mf.nGrowVect());
	}

	amrex::FillPatchTwoLevels(mftmp, time, cmf, ctime, fmf, ftime,
				  0, icomp, ncomp, geom[lev-1], geom[lev],
				  cphysbc, 0, fphysbc, 0, refRatio(lev-1),
				  mapper, bcs, 0);

	if constexpr (!std::is_same_v<MultiFab,MF>) {
	    amrex::LocalCopy(mf, mftmp, 0, 0, mf.nComp(), mf.nGrowVect());
	}
    }
}

// fill an entire multifab by interpolating from the coarser level
// this comes into play when a new level of refinement appears
template <typename MF>
void
AmrCoreAdv::FillCoarsePatch (int lev, amrex::Real time, MF& mf, int icomp, int ncomp)
{
    using namespace amrex;

    BL_ASSERT(lev > 0);

    Vector<MultiFab> cmf;
    Vector<Real> ctime;
    GetData(lev-1, time, cmf, ctime);

    if (cmf.size() != 1) {
        amrex::Abort("FillCoarsePatch: how did this happen?");
    }

    Interpolater* mapper = &cell_cons_interp;
    Util::PhysBCFunctNoOp<MultiFab> cphysbc{};
    Util::PhysBCFunctNoOp<MultiFab> fphysbc{};

    MultiFab mftmp;
    if constexpr (std::is_same_v<MultiFab,MF>) {
	mftmp = MultiFab(mf, amrex::make_alias, 0, mf.nComp());
    } else {
	mftmp.define(mf.boxArray(), mf.DistributionMap(), mf.nComp(), mf.nGrowVect());
    }

    amrex::InterpFromCoarseLevel(mftmp, time, cmf[0], 0, icomp, ncomp, geom[lev-1], geom[lev],
				 cphysbc, 0, fphysbc, 0, refRatio(lev-1),
				 mapper, bcs, 0);

    if constexpr (!std::is_same_v<MultiFab,MF>) {
	amrex::LocalCopy(mf, mftmp, 0, 0, mf.nComp(), mf.nGrowVect());
    }
}

#endif
