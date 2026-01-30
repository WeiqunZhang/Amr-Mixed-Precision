
#include <AmrCoreAdv.H>
#include <AmrCoreAdvUtil.H>

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>

#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

using namespace amrex;

// constructor - reads in parameters from inputs file
//             - sizes multilevel arrays and data structures
//             - initializes BCRec boundary condition object
AmrCoreAdv::AmrCoreAdv ()
{
    ReadParameters();

    // Geometry on all levels has been defined already.

    // No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

    int nlevs_max = max_level + 1;

    istep.resize(nlevs_max, 0);
    nsubsteps.resize(nlevs_max, 1);
    if (do_subcycle) {
        for (int lev = 1; lev <= max_level; ++lev) {
            nsubsteps[lev] = MaxRefRatio(lev-1);
        }
    }

    t_new.resize(nlevs_max, 0.0);
    t_old.resize(nlevs_max, -1.e100);
    dt.resize(nlevs_max, 1.e100);

    phi_new.resize(nlevs_max);
    phi_old.resize(nlevs_max);

    facevel.resize(nlevs_max);

    // periodic boundaries
    int bc_lo[] = {BCType::int_dir, BCType::int_dir, BCType::int_dir};
    int bc_hi[] = {BCType::int_dir, BCType::int_dir, BCType::int_dir};

/*
    // walls (Neumann)
    int bc_lo[] = {amrex::BCType::foextrap, amrex::BCType::foextrap, amrex::BCType::foextrap};
    int bc_hi[] = {amrex::BCType::foextrap, amrex::BCType::foextrap, amrex::BCType::foextrap};
*/

    bcs.resize(1);     // Setup 1-component
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        // lo-side BCs
        if (bc_lo[idim] == BCType::int_dir  ||  // periodic uses "internal Dirichlet"
            bc_lo[idim] == BCType::foextrap ||  // first-order extrapolation
            bc_lo[idim] == BCType::ext_dir ) {  // external Dirichlet
            bcs[0].setLo(idim, bc_lo[idim]);
        }
        else {
            amrex::Abort("Invalid bc_lo");
        }

        // hi-side BCSs
        if (bc_hi[idim] == BCType::int_dir  ||  // periodic uses "internal Dirichlet"
            bc_hi[idim] == BCType::foextrap ||  // first-order extrapolation
            bc_hi[idim] == BCType::ext_dir ) {  // external Dirichlet
            bcs[0].setHi(idim, bc_hi[idim]);
        }
        else {
            amrex::Abort("Invalid bc_hi");
        }
    }

    // stores fluxes at coarse-fine interface for synchronization
    // this will be sized "nlevs_max+1"
    // NOTE: the flux register associated with flux_reg[lev] is associated
    // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
    // therefore flux_reg[0] is never actually used in the reflux operation
    flux_reg.resize(nlevs_max+1);

    // amrex::AMRErrorTag supports a number of common tagging
    // approaches. Here we set it up to use amrex::Parser.
    if (! error_fn.empty()) {
        amrex::Parser parser(error_fn);
        parser.registerVariables({"x","y","z","t"});
        error_tag.emplace_back(std::move(parser));
    }
}

AmrCoreAdv::~AmrCoreAdv () = default;

// advance solution to final time
void
AmrCoreAdv::Evolve ()
{
    Real cur_time = t_new[0];
    int last_plot_file_step = 0;

    for (int step = istep[0]; step < max_step && cur_time < stop_time; ++step)
    {
        amrex::Print() << "\nCoarse STEP " << step+1 << " starts ..." << '\n';

        ComputeDt();

        int lev = 0;
        int iteration = 1;
	timeStepWithSubcycling(lev, cur_time, iteration);

        cur_time += dt[0];

        // sum phi to check conservation
        Real sum_phi;
	if (DoublePrecisionOnLevel(0)) {
	    sum_phi = std::get<MultiFab>(phi_new[0]).sum(0,IntVect(0));
	} else {
	    sum_phi = std::get<fMultiFab>(phi_new[0]).sum(0,IntVect(0));
	}

        amrex::Print() << "Coarse STEP " << step+1 << " ends." << " TIME = " << cur_time
                       << " DT = " << dt[0] << " Sum(Phi) = " << sum_phi << '\n';

        // sync up time
        for (lev = 0; lev <= finest_level; ++lev) {
            t_new[lev] = cur_time;
        }

        if (plot_int > 0 && (step+1) % plot_int == 0) {
            last_plot_file_step = step+1;
            WritePlotFile();
        }

        if (chk_int > 0 && (step+1) % chk_int == 0) {
            WriteCheckpointFile();
        }

#ifdef AMREX_MEM_PROFILING
        {
            std::ostringstream ss;
            ss << "[STEP " << step+1 << "]";
            MemProfiler::report(ss.str());
        }
#endif

        if (cur_time >= stop_time - 1.e-6*dt[0]) { break; }
    }

    if (plot_int > 0 && istep[0] > last_plot_file_step) {
        WritePlotFile();
    }
}

// initializes multilevel data
void
AmrCoreAdv::InitData ()
{
    if (restart_chkfile.empty()) {
        // start simulation from the beginning
        const Real time = 0.0;
        InitFromScratch(time);
        AverageDown();

#ifdef AMREX_PARTICLES
        if (do_tracers) {
            init_particles();
        }
#endif

        if (chk_int > 0) {
            WriteCheckpointFile();
        }
    }
    else {
        // restart from a checkpoint
        ReadCheckpointFile();
    }
    if (plot_int > 0) {
        WritePlotFile();
    }
}

// Make a new level using provided BoxArray and DistributionMapping and
// fill with interpolated coarse level data.
// overrides the pure virtual function in AmrCore
void
AmrCoreAdv::MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba,
                                    const DistributionMapping& dm)
{
    if (DoublePrecisionOnLevel(lev)) {
	MakeNewLevelFromCoarse<MultiFab>(lev, time, ba, dm);
    } else {
	MakeNewLevelFromCoarse<fMultiFab>(lev, time, ba, dm);
    }
}

// Remake an existing level using provided BoxArray and DistributionMapping and
// fill with existing fine and coarse data.
// overrides the pure virtual function in AmrCore
void
AmrCoreAdv::RemakeLevel (int lev, Real time, const BoxArray& ba,
                         const DistributionMapping& dm)
{
    if (DoublePrecisionOnLevel(lev)) {
	RemakeLevel<MultiFab>(lev, time, ba, dm);
    } else {
	RemakeLevel<fMultiFab>(lev, time, ba, dm);
    }
}

// Delete level data
// overrides the pure virtual function in AmrCore
void
AmrCoreAdv::ClearLevel (int lev)
{
    phi_new[lev] = std::variant<MultiFab,fMultiFab>{};
    phi_old[lev] = std::variant<MultiFab,fMultiFab>{};
    flux_reg[lev].reset(nullptr);
}

// Make a new level from scratch using provided BoxArray and DistributionMapping.
// Only used during initialization.
// overrides the pure virtual function in AmrCore
void AmrCoreAdv::MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba,
                                          const DistributionMapping& dm)
{
    if (DoublePrecisionOnLevel(lev)) {
	MakeNewLevelFromScratch<MultiFab>(lev, time, ba, dm);
    } else {
	MakeNewLevelFromScratch<fMultiFab>(lev, time, ba, dm);
    }
}

// tag all cells for refinement
// overrides the pure virtual function in AmrCore
void
AmrCoreAdv::ErrorEst (int lev, TagBoxArray& tags, Real time, int ngrow)
{
    if (DoublePrecisionOnLevel(lev)) {
	ErrorEst<MultiFab>(lev, tags, time, ngrow);
    } else {
	ErrorEst<fMultiFab>(lev, tags, time, ngrow);
    }
}

// read in some parameters from inputs file
void
AmrCoreAdv::ReadParameters ()
{
    {
        ParmParse pp;  // Traditionally, max_step and stop_time do not have prefix.
        pp.query("max_step", max_step);
        pp.query("stop_time", stop_time);
    }

    {
        ParmParse pp("amr"); // Traditionally, these have prefix, amr.

        pp.query("regrid_int", regrid_int);
        pp.query("plot_file", plot_file);
        pp.query("plot_int", plot_int);
        pp.query("chk_file", chk_file);
        pp.query("chk_int", chk_int);
        pp.query("restart",restart_chkfile);
    }

    {
        ParmParse pp("adv");

        pp.query("cfl", cfl);
        pp.query("do_reflux", do_reflux);
        pp.query("do_subcycle", do_subcycle);
        pp.query("errfn", error_fn);
	pp.queryarr("double_precision", double_precision);
    }

#ifdef AMREX_PARTICLES
    {
        ParmParse pp("amr");
        pp.query("do_tracers", do_tracers);
    }
#endif
}

// set covered coarse cells to be the average of overlying fine cells
void
AmrCoreAdv::AverageDown ()
{
    for (int lev = finest_level-1; lev >= 0; --lev)
    {
	AverageDownTo(lev);
    }
}

// more flexible version of AverageDown() that lets you average down across multiple levels
void
AmrCoreAdv::AverageDownTo (int crse_lev)
{
    int fine_lev = crse_lev+1;
    if (DoublePrecisionOnLevel(fine_lev) && DoublePrecisionOnLevel(crse_lev)) {
	Util::average_down(std::get<MultiFab>(phi_new[fine_lev]),
			   std::get<MultiFab>(phi_new[crse_lev]),
			   refRatio(crse_lev));
    } else if (!DoublePrecisionOnLevel(fine_lev) && DoublePrecisionOnLevel(crse_lev)) {
	Util::average_down(std::get<fMultiFab>(phi_new[fine_lev]),
			   std::get<MultiFab>(phi_new[crse_lev]),
			   refRatio(crse_lev));
    } else if (DoublePrecisionOnLevel(fine_lev) && !DoublePrecisionOnLevel(crse_lev)) {
	Util::average_down(std::get<MultiFab>(phi_new[fine_lev]),
			   std::get<fMultiFab>(phi_new[crse_lev]),
			   refRatio(crse_lev));
    } else {
	Util::average_down(std::get<fMultiFab>(phi_new[fine_lev]),
			   std::get<fMultiFab>(phi_new[crse_lev]),
			   refRatio(crse_lev));
    }
}

// Advance a level by dt
// (includes a recursive call for finer levels)
void
AmrCoreAdv::timeStepWithSubcycling (int lev, Real time, int iteration)
{
    if (regrid_int > 0)  // We may need to regrid
    {

        // help keep track of whether a level was already regridded
        // from a coarser level call to regrid
        static Vector<int> last_regrid_step(max_level+1, 0);

        // regrid changes level "lev+1" so we don't regrid on max_level
        // also make sure we don't regrid fine levels again if
        // it was taken care of during a coarser regrid
        if (lev < max_level && istep[lev] > last_regrid_step[lev])
        {
            if (istep[lev] % regrid_int == 0)
            {
                // regrid could add newly refine levels (if finest_level < max_level)
                // so we save the previous finest level index
                int old_finest = finest_level;
                regrid(lev, time);

                // mark that we have regridded this level already
                for (int k = lev; k <= finest_level; ++k) {
                    last_regrid_step[k] = istep[k];
                }

                // if there are newly created levels, set the time step
                for (int k = old_finest+1; k <= finest_level; ++k) {
                    dt[k] = dt[k-1] / MaxRefRatio(k-1);
                }

#ifdef AMREX_PARTICLES
                if (do_tracers) {
                    TracerPC->Redistribute(lev);
                }
#endif
            }
        }
    }

    if (Verbose()) {
        amrex::Print() << "[Level " << lev << " step " << istep[lev]+1 << "] ";
        amrex::Print() << "ADVANCE with time = " << t_new[lev]
                       << " dt = " << dt[lev] << '\n';
    }

    // Advance a single level for a single time step, and update flux registers

    t_old[lev] = t_new[lev];
    t_new[lev] += dt[lev];

    Real t_nph = t_old[lev] + 0.5*dt[lev];

    DefineVelocityAtLevel(lev, t_nph);
    AdvancePhiAtLevel(lev, time, dt[lev], iteration, nsubsteps[lev]);


#ifdef AMREX_PARTICLES
    if (do_tracers) {
        TracerPC->AdvectWithUmac(facevel[lev].data(),lev,dt[lev]);
    }
#endif

    ++istep[lev];

    if (Verbose())
    {
        amrex::Print() << "[Level " << lev << " step " << istep[lev] << "] ";
        amrex::Print() << "Advanced " << CountCells(lev) << " cells" << '\n';
    }

    if (lev < finest_level)
    {
        // recursive call for next-finer level
        for (int i = 1; i <= nsubsteps[lev+1]; ++i)
        {
            timeStepWithSubcycling(lev+1, time+(i-1)*dt[lev+1], i);
        }

        if (do_reflux)
        {
            // update lev based on coarse-fine flux mismatch
            MultiFab mftmp;
            if (SamePrecision<MultiFab>(lev)) {
                auto& mf = std::get<MultiFab>(phi_new[lev]);
                mftmp = MultiFab(mf, amrex::make_alias, 0, mf.nComp());
            } else {
                auto& mf = std::get<fMultiFab>(phi_new[lev]);
                mftmp = amrex::cast<MultiFab>(mf);
            } 
            flux_reg[lev+1]->Reflux(mftmp, 1.0, 0, 0, mftmp.nComp(), geom[lev]);
            if (!SamePrecision<MultiFab>(lev)) {
                auto& mf = std::get<fMultiFab>(phi_new[lev]);
                amrex::LocalCopy(mf, mftmp, 0, 0, mf.nComp(), mf.nGrowVect());
            }
        }

        AverageDownTo(lev); // average lev+1 down to lev
    }


#ifdef AMREX_PARTICLES
    if (do_tracers) {
        int redistribute_ngrow = 0;
        if ((iteration < nsubsteps[lev]) || (lev == 0)){
            if (lev == 0){
                redistribute_ngrow = 0;
            } else {
                redistribute_ngrow = iteration;
            }
            TracerPC->Redistribute(lev, TracerPC->finestLevel(), redistribute_ngrow);
        }
    }
#endif
}

// a wrapper for EstTimeStep
void
AmrCoreAdv::ComputeDt ()
{
    Vector<Real> dt_tmp(finest_level+1);

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        dt_tmp[lev] = EstTimeStep(lev, t_new[lev]);
    }
    ParallelDescriptor::ReduceRealMin(dt_tmp.data(), int(dt_tmp.size()));

    constexpr Real change_max = 1.1;
    Real dt_0 = dt_tmp[0];
    int n_factor = 1;

    for (int lev = 0; lev <= finest_level; ++lev) {
        dt_tmp[lev] = std::min(dt_tmp[lev], change_max*dt[lev]);
        n_factor *= nsubsteps[lev];
        dt_0 = std::min(dt_0, n_factor*dt_tmp[lev]);
    }

    // Limit dt's by the value of stop_time.
    const Real eps = 1.e-3*dt_0;

    if (t_new[0] + dt_0 > stop_time - eps) {
        dt_0 = stop_time - t_new[0];
    }

    dt[0] = dt_0;

    for (int lev = 1; lev <= finest_level; ++lev) {
        dt[lev] = dt[lev-1] / nsubsteps[lev];
    }
}

// compute dt from CFL considerations
Real
AmrCoreAdv::EstTimeStep (int lev, Real time)
{
    BL_PROFILE("AmrCoreAdv::EstTimeStep()");

    Real dt_est = std::numeric_limits<Real>::max();

    const Real* dx  =  geom[lev].CellSize();

    if (time == Real(0.0)) {
       DefineVelocityAtLevel(lev,time);
    } else {
       Real t_nph_predicted = time + 0.5 * dt[lev];
       DefineVelocityAtLevel(lev,t_nph_predicted);
    }

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        Real est;
	if (DoublePrecisionOnLevel(lev)) {
	    est = std::get<MultiFab>(facevel[lev][idim]).norminf(0,1,IntVect(0),true);
	} else {
	    est = std::get<fMultiFab>(facevel[lev][idim]).norminf(0,1,IntVect(0),true);
	}
        dt_est = amrex::min(dt_est, dx[idim]/est);
    }

    dt_est *= cfl;

    return dt_est;
}

// get plotfile name
std::string
AmrCoreAdv::PlotFileName (int lev) const
{
    return amrex::Concatenate(plot_file, lev, 5);
}

// put together an array of multifabs for writing
Vector<MultiFab>
AmrCoreAdv::PlotFileMF () const
{
    Vector<MultiFab> r;
    for (int i = 0; i <= finest_level; ++i) {
	if (SamePrecision<MultiFab>(i)) {
	    auto const& mf = std::get<MultiFab>(phi_new[i]);
	    r.push_back(MultiFab(mf, amrex::make_alias, 0, mf.nComp()));
	} else {
	    auto const& mf = std::get<fMultiFab>(phi_new[i]);
	    r.push_back(amrex::cast<MultiFab>(mf));
	}
    }
    return r;
}

// set plotfile variable names
Vector<std::string>
AmrCoreAdv::PlotFileVarNames ()
{
    return {"phi"};
}

// write plotfile to disk
void
AmrCoreAdv::WritePlotFile () const
{
    const std::string& plotfilename = PlotFileName(istep[0]);
    const auto mf0 = PlotFileMF();
    auto mf = GetVecOfConstPtrs(mf0);
    const auto& varnames = PlotFileVarNames();

    amrex::Print() << "Writing plotfile " << plotfilename << "\n";

    amrex::WriteMultiLevelPlotfile(plotfilename, finest_level+1, mf, varnames,
                                   Geom(), t_new[0], istep, refRatio());

#ifdef AMREX_PARTICLES
        if (do_tracers) {
            TracerPC->WritePlotFile(plotfilename, "particles");
        }
#endif
}

void
AmrCoreAdv::WriteCheckpointFile () const
{
    AMREX_ALWAYS_ASSERT("xxxxx TODO: WriteCheckpointFile");
#if 0 // WriteCheckpointFile

    // chk00010            write a checkpoint file with this root directory
    // chk00010/Header     this contains information you need to save (e.g., finest_level, t_new, etc.) and also
    //                     the BoxArrays at each level
    // chk00010/Level_0/
    // chk00010/Level_1/
    // etc.                these subdirectories will hold the MultiFab data at each level of refinement

    // checkpoint file name, e.g., chk00010
    const std::string& checkpointname = amrex::Concatenate(chk_file,istep[0]);

    amrex::Print() << "Writing checkpoint " << checkpointname << "\n";

    const int nlevels = finest_level+1;

    // ---- prebuild a hierarchy of directories
    // ---- dirName is built first.  if dirName exists, it is renamed.  then build
    // ---- dirName/subDirPrefix_0 .. dirName/subDirPrefix_nlevels-1
    // ---- if callBarrier is true, call ParallelDescriptor::Barrier()
    // ---- after all directories are built
    // ---- ParallelDescriptor::IOProcessor() creates the directories
    amrex::PreBuildDirectorHierarchy(checkpointname, "Level_", nlevels, true);

    // write Header file
   if (ParallelDescriptor::IOProcessor()) {

       std::string HeaderFileName(checkpointname + "/Header");
       VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
       std::ofstream HeaderFile;
       HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
       HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                                               std::ofstream::trunc |
                                               std::ofstream::binary);
       if( ! HeaderFile.good()) {
           amrex::FileOpenFailed(HeaderFileName);
       }

       HeaderFile.precision(17);

       // write out title line
       HeaderFile << "Checkpoint file for AmrCoreAdv\n";

       // write out finest_level
       HeaderFile << finest_level << "\n";

       // write out array of istep
       for (auto s : istep) {
           HeaderFile << s << " ";
       }
       HeaderFile << "\n";

       // write out array of dt
       for (auto dti : dt) {
           HeaderFile << dti << " ";
       }
       HeaderFile << "\n";

       // write out array of t_new
       for (auto t_new_i : t_new) {
           HeaderFile << t_new_i << " ";
       }
       HeaderFile << "\n";

       // write the BoxArray at each level
       for (int lev = 0; lev <= finest_level; ++lev) {
           boxArray(lev).writeOn(HeaderFile);
           HeaderFile << '\n';
       }
   }

   // write the MultiFab data to, e.g., chk00010/Level_0/
   for (int lev = 0; lev <= finest_level; ++lev) {
       VisMF::Write(phi_new[lev],
                    amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "phi"));
   }

#ifdef AMREX_PARTICLES
            if (do_tracers) {
                TracerPC->Checkpoint(checkpointname, "particles", true);
            }
#endif

#endif
}

#ifdef AMREX_PARTICLES
void
AmrCoreAdv::init_particles ()
{
  if (do_tracers)
    {
      BL_ASSERT(TracerPC == nullptr);

      TracerPC = std::make_unique<AmrTracerParticleContainer>(this);

      AmrTracerParticleContainer::ParticleInitData pdata = {{AMREX_D_DECL(0.0, 0.0, 0.0)},{},{},{}};

      TracerPC->SetVerbose(0);
      TracerPC->InitOnePerCell(0.5, 0.5, 0.5, pdata);
      TracerPC->Redistribute();
    }
}
#endif


namespace {
// utility to skip to next line in Header
void GotoNextLine (std::istream& is)
{
    constexpr std::streamsize bl_ignore_max { 100000 };
    is.ignore(bl_ignore_max, '\n');
}
}

void
AmrCoreAdv::ReadCheckpointFile ()
{
    AMREX_ALWAYS_ASSERT("xxxxx TODO: ReadCheckpointFile");
#if 0 // ReadCheckpointFile

    amrex::Print() << "Restart from checkpoint " << restart_chkfile << "\n";

    // Header
    std::string File(restart_chkfile + "/Header");

    VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);

    std::string line, word;

    // read in title line
    std::getline(is, line);

    // read in finest_level
    is >> finest_level;
    GotoNextLine(is);

    // read in array of istep
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            istep[i++] = std::stoi(word);
        }
    }

    // read in array of dt
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            dt[i++] = std::stod(word);
        }
    }

    // read in array of t_new
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            t_new[i++] = std::stod(word);
        }
    }

    for (int lev = 0; lev <= finest_level; ++lev) {

        // read in level 'lev' BoxArray from Header
        BoxArray ba;
        ba.readFrom(is);
        GotoNextLine(is);

        // create a distribution mapping
        DistributionMapping dm { ba, ParallelDescriptor::NProcs() };

        // set BoxArray grids and DistributionMapping dmap in AMReX_AmrMesh.H class
        SetBoxArray(lev, ba);
        SetDistributionMap(lev, dm);

        // build MultiFab and FluxRegister data
        int ncomp = 1;
        int ng = 0;
        phi_old[lev].define(grids[lev], dmap[lev], ncomp, ng);
        phi_new[lev].define(grids[lev], dmap[lev], ncomp, ng);

        if (lev > 0 && do_reflux) {
            flux_reg[lev] = std::make_unique<FluxRegister>(grids[lev], dmap[lev], refRatio(lev-1), lev, ncomp);
        }

        // build face velocity MultiFabs
        for (int idim = 0; idim < AMREX_SPACEDIM; idim++)
        {
            facevel[lev][idim] = MultiFab(amrex::convert(ba,IntVect::TheDimensionVector(idim)), dm, 1, 1);
        }
    }

    // read in the MultiFab data
    for (int lev = 0; lev <= finest_level; ++lev) {
        VisMF::Read(phi_new[lev],
                    amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "phi"));
    }

#ifdef AMREX_PARTICLES
    if (do_tracers) {
        BL_ASSERT(TracerPC == nullptr);
        TracerPC = std::make_unique<AmrTracerParticleContainer>(this);
        TracerPC->Restart(this->restart_chkfile, "particles");
    }
#endif

#endif
}
