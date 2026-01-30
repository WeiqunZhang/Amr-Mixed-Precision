#include <AdvancePhiAtLevel.H>

// Advance a single level for a single time step, updates flux registers
void
AmrCoreAdv::AdvancePhiAtLevel (int lev, amrex::Real time, amrex::Real dt_lev,
			       int iteration, int ncycle)
{
    if (DoublePrecisionOnLevel(lev)) {
	AdvancePhiAtLevel<amrex::MultiFab>(lev, time, dt_lev, iteration, ncycle);
    } else {
	AdvancePhiAtLevel<amrex::fMultiFab>(lev, time, dt_lev, iteration, ncycle);
    }
}
