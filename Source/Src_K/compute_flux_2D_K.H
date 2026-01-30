#ifndef compute_flux_2d_H_
#define compute_flux_2d_H_

#include <AMReX_Array4.H>
#include <AMReX_REAL.H>

template <typename T>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void flux_x(int i, int j, int k,
            amrex::Array4<T> const& px,
            amrex::Array4<T const> const& phi,
            amrex::Array4<T const> const& vx,
            amrex::Array4<T const> const& slope,
            T dtdx)
{
    px(i,j,k) = ( (vx(i,j,k) < 0) ?
                phi(i  ,j,k) - slope(i  ,j,k)*(T(0.5) + T(0.5)*dtdx*vx(i,j,k)) :
                phi(i-1,j,k) + slope(i-1,j,k)*(T(0.5) - T(0.5)*dtdx*vx(i,j,k)) );
}

template <typename T>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void flux_y(int i, int j, int k,
            amrex::Array4<T> const& py,
            amrex::Array4<T const> const& phi,
            amrex::Array4<T const> const& vy,
            amrex::Array4<T const> const& slope,
            T dtdy)
{
    py(i,j,k) = ( (vy(i,j,k) < 0) ?
                phi(i,j  ,k) - slope(i,j  ,k)*(T(0.5) + T(0.5)*dtdy*vy(i,j,k)) :
                phi(i,j-1,k) + slope(i,j-1,k)*(T(0.5) - T(0.5)*dtdy*vy(i,j,k)) );
}

template <typename T>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void create_flux_x(int i, int j, int k,
                   amrex::Array4<T> const& fx,
                   amrex::Array4<T const> const& vx,
                   amrex::Array4<T const> const& vy,
                   amrex::Array4<T const> const& px,
                   amrex::Array4<T const> const& py,
                   T dtdy)
{
    fx(i,j,k) = ( (vx(i,j,k) < 0) ?
                (px(i,j,k) - T(0.5)*dtdy * ( T(0.5)*(vy(i  ,j+1,k  ) + vy(i  ,j,k)) * (py(i  ,j+1,k  )-py(i  ,j,k))))*vx(i,j,k) :
                (px(i,j,k) - T(0.5)*dtdy * ( T(0.5)*(vy(i-1,j+1,k  ) + vy(i-1,j,k)) * (py(i-1,j+1,k  )-py(i-1,j,k))))*vx(i,j,k) );
}

template <typename T>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void create_flux_y(int i, int j, int k,
                   amrex::Array4<T> const& fy,
                   amrex::Array4<T const> const& vx,
                   amrex::Array4<T const> const& vy,
                   amrex::Array4<T const> const& py,
                   amrex::Array4<T const> const& px,
                   T dtdx)
{
    fy(i,j,k) = ( (vy(i,j,k) < 0) ?
                (py(i,j,k) - T(0.5)*dtdx * ( T(0.5)*(vx(i+1,j  ,k  ) + vx(i,j  ,k)) * (px(i+1,j  ,k  )-px(i,j  ,k))))*vy(i,j,k) :
                (py(i,j,k) - T(0.5)*dtdx * ( T(0.5)*(vx(i+1,j-1,k  ) + vx(i,j-1,k)) * (px(i+1,j-1,k  )-px(i,j-1,k))))*vy(i,j,k) );
}

#endif
