#ifndef slope_K_H_
#define slope_K_H_

#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>

template <typename T>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void slopex2(int i, int j, int k,
             amrex::Array4<T> const& dq,
             amrex::Array4<T const> const& q)
{
    T dlft = q(i,j,k) - q(i-1,j,k);
    T drgt = q(i+1,j,k) - q(i,j,k);
    T dcen = T(0.5)*(dlft+drgt);
    T dsgn = std::copysign(T(1.0), dcen);
    T dslop = T(2.0) * ((std::abs(dlft) < std::abs(drgt)) ?
                                std::abs(dlft) : std::abs(drgt));
    T dlim = (dlft*drgt >= T(0.0)) ? dslop : T(0.0);
    dq(i,j,k) = dsgn*amrex::min(dlim, std::abs(dcen));
}

template <typename T>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void slopex4(int i, int j, int k,
             amrex::Array4<T> const& dq4,
             amrex::Array4<T const> const& q,
             amrex::Array4<T const> const& dq)
{
    T dlft = q(i,j,k) - q(i-1,j,k);
    T drgt = q(i+1,j,k) - q(i,j,k);
    T dcen = T(0.5)*(dlft+drgt);
    T dsgn = std::copysign(T(1.0), dcen);
    T dslop = T(2.0) * ((std::abs(dlft) < std::abs(drgt)) ?
                                std::abs(dlft) : std::abs(drgt));
    T dlim = (dlft*drgt >= T(0.0)) ? dslop : T(0.0);
    T dq1 = T(4.0/3.0)*dcen - T(1.0/6.0)*(dq(i+1,j,k) + dq(i-1,j,k));
    dq4(i,j,k) = dsgn*amrex::min(dlim, std::abs(dq1));
}

// ***********************************************************

template <typename T>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void slopey2(int i, int j, int k,
             amrex::Array4<T> const& dq,
             amrex::Array4<T const> const& q)
{
    T dlft = q(i,j,k) - q(i,j-1,k);
    T drgt = q(i,j+1,k) - q(i,j,k);
    T dcen = T(0.5)*(dlft+drgt);
    T dsgn = std::copysign(T(1.0), dcen);
    T dslop = T(2.0) * ((std::abs(dlft) < std::abs(drgt)) ?
                                std::abs(dlft) : std::abs(drgt));
    T dlim = (dlft*drgt >= T(0.0)) ? dslop : T(0.0);
    dq(i,j,k) = dsgn*amrex::min(dlim, std::abs(dcen));
}

template <typename T>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void slopey4(int i, int j, int k,
             amrex::Array4<T> const& dq4,
             amrex::Array4<T const> const& q,
             amrex::Array4<T const> const& dq)
{
    T dlft = q(i,j,k) - q(i,j-1,k);
    T drgt = q(i,j+1,k) - q(i,j,k);
    T dcen = T(0.5)*(dlft+drgt);
    T dsgn = std::copysign(T(1.0), dcen);
    T dslop = T(2.0) * ((std::abs(dlft) < std::abs(drgt)) ?
                                std::abs(dlft) : std::abs(drgt));
    T dlim = (dlft*drgt >= T(0.0)) ? dslop : T(0.0);
    T dq1 = T(4.0/3.0)*dcen - T(1.0/6.0)*(dq(i,j+1,k) + dq(i,j-1,k));
    dq4(i,j,k) = dsgn*amrex::min(dlim, std::abs(dq1));
}

// ***********************************************************

template <typename T>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void slopez2(int i, int j, int k,
             amrex::Array4<T> const& dq,
             amrex::Array4<T const> const& q)
{
    T dlft = q(i,j,k) - q(i,j,k-1);
    T drgt = q(i,j,k+1) - q(i,j,k);
    T dcen = T(0.5)*(dlft+drgt);
    T dsgn = std::copysign(T(1.0), dcen);
    T dslop = T(2.0) * ((std::abs(dlft) < std::abs(drgt)) ?
                                std::abs(dlft) : std::abs(drgt));
    T dlim = (dlft*drgt >= T(0.0)) ? dslop : T(0.0);
    dq(i,j,k) = dsgn*amrex::min(dlim, std::abs(dcen));
}

template <typename T>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void slopez4(int i, int j, int k,
             amrex::Array4<T> const& dq4,
             amrex::Array4<T const> const& q,
             amrex::Array4<T const> const& dq)
{
    T dlft = q(i,j,k) - q(i,j,k-1);
    T drgt = q(i,j,k+1) - q(i,j,k);
    T dcen = T(0.5)*(dlft+drgt);
    T dsgn = std::copysign(T(1.0), dcen);
    T dslop = T(2.0) * ((std::abs(dlft) < std::abs(drgt)) ?
                                std::abs(dlft) : std::abs(drgt));
    T dlim = (dlft*drgt >= T(0.0)) ? dslop : T(0.0);
    T dq1 = T(4.0/3.0)*dcen - T(1.0/6.0)*(dq(i,j,k+1) + dq(i,j,k-1));
    dq4(i,j,k) = dsgn*amrex::min(dlim, std::abs(dq1));
}

#endif
