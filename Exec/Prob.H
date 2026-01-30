#ifndef PROB_H_
#define PROB_H_

#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>

template <typename T>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
initdata(int i, int j, int k, amrex::Array4<T> const& phi,
         amrex::GpuArray<T,AMREX_SPACEDIM> const& prob_lo,
         amrex::GpuArray<T,AMREX_SPACEDIM> const& dx)
{
    T y = prob_lo[1] + (T(0.5)+j) * dx[1];
    T x = prob_lo[0] + (T(0.5)+i) * dx[0];
    T r2 = (std::pow(x-T(0.5), 2) + std::pow((y-T(0.75)),2)) / T(0.01);
    phi(i,j,k) = T(1.0) + std::exp(-r2);
}

#endif
