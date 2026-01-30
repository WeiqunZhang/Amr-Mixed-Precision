#include <DefineVelocity.H>

void
AmrCoreAdv::DefineVelocityAtLevel (int lev, amrex::Real time)
{
    if (DoublePrecisionOnLevel(lev)) {
	DefineVelocityAtLevel<amrex::MultiFab>(lev, time);
    } else {
	DefineVelocityAtLevel<amrex::fMultiFab>(lev, time);
    }
}
