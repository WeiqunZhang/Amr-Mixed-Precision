#ifndef AMRCORE_ADV_UTIL_H_
#define AMRCORE_ADV_UTIL_H_

#include <AMReX_FabArray.H>
#include <AMReX_MFParallelFor.H>

namespace Util {

template <typename FMF, typename CMF,
	  std::enable_if_t<amrex::IsFabArray_v<FMF> && amrex::IsFabArray_v<CMF>, int> = 0>
void average_down (FMF const& S_fine, CMF& S_crse, amrex::IntVect const& ratio)
{
    using namespace amrex;
    using CRT = typename CMF::value_type;

    AMREX_ALWAYS_ASSERT(ratio == 2 && S_fine.nComp() == 1 && AMREX_SPACEDIM == 2);

    CMF cfine(amrex::coarsen(S_fine.boxArray(),ratio), S_fine.DistributionMap(),
              S_fine.nComp(), 0);

    auto const& fma = S_fine.const_arrays();
    auto const& cma = cfine.arrays();

    ParallelFor(cfine, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
    {
	auto const& fa = fma[b];
	auto const& ca = cma[b];
	int ii = 2*i;
	int jj = 2*j;
	int kk = 2*k;
	ca(i,j,k) = CRT(fa(ii,jj  ,kk) + fa(ii+1,jj  ,kk) +
			fa(ii,jj+1,kk) + fa(ii+1,jj+1,kk)) * CRT(0.25);
    });
    Gpu::streamSynchronize();

    S_crse.ParallelCopy(cfine, 0, 0, S_fine.nComp());
}

template <typename MF>
class PhysBCFunctNoOp
{
public:
    void operator() (MF& /*mf*/, int /*dcomp*/, int /*ncomp*/, amrex::IntVect const& /*nghost*/,
                     amrex::Real /*time*/, int /*bccomp*/) {}
};

}

#endif
