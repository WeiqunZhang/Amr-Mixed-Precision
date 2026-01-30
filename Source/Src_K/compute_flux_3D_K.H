#ifndef compute_flux_3d_H_
#define compute_flux_3d_H_

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
void flux_z(int i, int j, int k,
            amrex::Array4<T> const& pz,
            amrex::Array4<T const> const& phi,
            amrex::Array4<T const> const& vz,
            amrex::Array4<T const> const& slope,
            T dtdz)
{
    pz(i,j,k) = ( (vz(i,j,k) < 0) ?
                phi(i,j,k  ) - slope(i,j,k  )*(T(0.5) + T(0.5)*dtdz*vz(i,j,k)) :
                phi(i,j,k-1) + slope(i,j,k-1)*(T(0.5) - T(0.5)*dtdz*vz(i,j,k)) );
}

template <typename T>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void flux_xy(int i, int j, int k,
             amrex::Array4<T> const& pxy,
             amrex::Array4<T const> const& vx,
             amrex::Array4<T const> const& vy,
             amrex::Array4<T const> const& px,
             amrex::Array4<T const> const& py,
             T dtdy)
{
    pxy(i,j,k) = ( (vx(i,j,k) < 0) ?
                 px(i,j,k) - dtdy/T(3.0) * ( T(0.5)*(vy(i,  j+1,k) + vy(i  ,j,k)) * (py(i  ,j+1,k) - py(i  ,j,k))) :
                 px(i,j,k) - dtdy/T(3.0) * ( T(0.5)*(vy(i-1,j+1,k) + vy(i-1,j,k)) * (py(i-1,j+1,k) - py(i-1,j,k))) );
}

template <typename T>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void flux_xz(int i, int j, int k,
             amrex::Array4<T> const& pxz,
             amrex::Array4<T const> const& vx,
             amrex::Array4<T const> const& vz,
             amrex::Array4<T const> const& px,
             amrex::Array4<T const> const& pz,
             T dtdz)
{
    pxz(i,j,k) = ( (vx(i,j,k) < 0) ?
                 px(i,j,k) - dtdz/T(3.0) * ( T(0.5)*(vz(i,  j,k+1) + vz(i  ,j,k)) * (pz(i  ,j,k+1) - pz(i  ,j,k))) :
                 px(i,j,k) - dtdz/T(3.0) * ( T(0.5)*(vz(i-1,j,k+1) + vz(i-1,j,k)) * (pz(i-1,j,k+1) - pz(i-1,j,k))) );
}

template <typename T>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void flux_yx(int i, int j, int k,
             amrex::Array4<T> const& pyx,
             amrex::Array4<T const> const& vx,
             amrex::Array4<T const> const& vy,
             amrex::Array4<T const> const& px,
             amrex::Array4<T const> const& py,
             T dtdx)
{
    pyx(i,j,k) = ( (vy(i,j,k) < 0) ?
                 py(i,j,k) - dtdx/T(3.0) * ( T(0.5)*(vx(i+1,j  ,k) + vx(i,j  ,k)) * (px(i+1,j  ,k) - px(i,j  ,k))) :
                 py(i,j,k) - dtdx/T(3.0) * ( T(0.5)*(vx(i+1,j-1,k) + vx(i,j-1,k)) * (px(i+1,j-1,k) - px(i,j-1,k))) );
}

template <typename T>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void flux_yz(int i, int j, int k,
             amrex::Array4<T> const& pyz,
             amrex::Array4<T const> const& vy,
             amrex::Array4<T const> const& vz,
             amrex::Array4<T const> const& py,
             amrex::Array4<T const> const& pz,
             T dtdz)
{
    pyz(i,j,k) = ( (vy(i,j,k) < 0) ?
                 py(i,j,k) - dtdz/T(3.0) * ( T(0.5)*(vz(i,  j,k+1) + vz(i,j  ,k)) * (pz(i,j  ,k+1) - pz(i,j  ,k))) :
                 py(i,j,k) - dtdz/T(3.0) * ( T(0.5)*(vz(i,j-1,k+1) + vz(i,j-1,k)) * (pz(i,j-1,k+1) - pz(i,j-1,k))) );
}

template <typename T>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void flux_zx(int i, int j, int k,
             amrex::Array4<T> const& pzx,
             amrex::Array4<T const> const& vx,
             amrex::Array4<T const> const& vz,
             amrex::Array4<T const> const& px,
             amrex::Array4<T const> const& pz,
             T dtdx)
{
    pzx(i,j,k) = ( (vz(i,j,k) < 0) ?
                 pz(i,j,k) - dtdx/T(3.0) * ( T(0.5)*(vx(i+1,j,k  ) + vx(i,j,k  )) * (px(i+1,j,k  ) - px(i,j,k  ))) :
                 pz(i,j,k) - dtdx/T(3.0) * ( T(0.5)*(vx(i+1,j,k-1) + vx(i,j,k-1)) * (px(i+1,j,k-1) - px(i,j,k-1))) );
}

template <typename T>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void flux_zy(int i, int j, int k,
             amrex::Array4<T> const& pzy,
             amrex::Array4<T const> const& vy,
             amrex::Array4<T const> const& vz,
             amrex::Array4<T const> const& py,
             amrex::Array4<T const> const& pz,
             T dtdy)
{
    pzy(i,j,k) = ( (vz(i,j,k) < 0) ?
                 pz(i,j,k) - dtdy/T(3.0) * ( T(0.5)*(vy(i,j+1,k  ) + vy(i,j,k  )) * (py(i,j+1,k  ) - py(i,j,k  ))) :
                 pz(i,j,k) - dtdy/T(3.0) * ( T(0.5)*(vy(i,j+1,k-1) + vy(i,j,k-1)) * (py(i,j+1,k-1) - py(i,j,k-1))) );
}

template <typename T>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void create_flux_x(int i, int j, int k,
                   amrex::Array4<T> const& fx,
                   amrex::Array4<T const> const& vx,
                   amrex::Array4<T const> const& vy,
                   amrex::Array4<T const> const& vz,
                   amrex::Array4<T const> const& px,
                   amrex::Array4<T const> const& pyz,
                   amrex::Array4<T const> const& pzy,
                   T dtdy, T dtdz)
{
    T f = ( (vx(i,j,k) < 0) ?
                px(i,j,k) - T(0.5)*dtdy * ( T(0.5)*(vy(i  ,j+1,k  ) + vy(i  ,j,k)) * (pyz(i  ,j+1,k  )-pyz(i  ,j,k)))
                          - T(0.5)*dtdz * ( T(0.5)*(vz(i  ,j  ,k+1) + vz(i  ,j,k)) * (pzy(i  ,j  ,k+1)-pzy(i  ,j,k))) :
                px(i,j,k) - T(0.5)*dtdy * ( T(0.5)*(vy(i-1,j+1,k  ) + vy(i-1,j,k)) * (pyz(i-1,j+1,k  )-pyz(i-1,j,k)))
                          - T(0.5)*dtdz * ( T(0.5)*(vz(i-1,j  ,k+1) + vz(i-1,j,k)) * (pzy(i-1,j  ,k+1)-pzy(i-1,j,k))) );

    fx(i,j,k) = vx(i,j,k)*f;
}

template <typename T>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void create_flux_y(int i, int j, int k,
                   amrex::Array4<T> const& fy,
                   amrex::Array4<T const> const& vx,
                   amrex::Array4<T const> const& vy,
                   amrex::Array4<T const> const& vz,
                   amrex::Array4<T const> const& py,
                   amrex::Array4<T const> const& pxz,
                   amrex::Array4<T const> const& pzx,
                   T dtdx, T dtdz)
{
    T f = ( (vy(i,j,k) < 0) ?
                py(i,j,k) - T(0.5)*dtdx * ( T(0.5)*(vx(i+1,j  ,k  ) + vx(i,j  ,k)) * (pxz(i+1,j  ,k  )-pxz(i,j  ,k)))
                          - T(0.5)*dtdz * ( T(0.5)*(vz(i,  j  ,k+1) + vz(i,j  ,k)) * (pzx(i,  j  ,k+1)-pzx(i,j  ,k))) :
                py(i,j,k) - T(0.5)*dtdx * ( T(0.5)*(vx(i+1,j-1,k  ) + vx(i,j-1,k)) * (pxz(i+1,j-1,k  )-pxz(i,j-1,k)))
                          - T(0.5)*dtdz * ( T(0.5)*(vz(i  ,j-1,k+1) + vz(i,j-1,k)) * (pzx(i  ,j-1,k+1)-pzx(i,j-1,k))) );

    fy(i,j,k) = vy(i,j,k)*f;
}

template <typename T>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void create_flux_z(int i, int j, int k,
                   amrex::Array4<T> const& fz,
                   amrex::Array4<T const> const& vx,
                   amrex::Array4<T const> const& vy,
                   amrex::Array4<T const> const& vz,
                   amrex::Array4<T const> const& pz,
                   amrex::Array4<T const> const& pxy,
                   amrex::Array4<T const> const& pyx,
                   T dtdx, T dtdy)
{
    T f = ( (vz(i,j,k) < 0) ?
                pz(i,j,k) - T(0.5)*dtdx * ( T(0.5)*(vx(i+1,j  ,k  ) + vx(i,j,k  )) * (pxy(i+1,j  ,k  )-pxy(i,j,k  )))
                          - T(0.5)*dtdy * ( T(0.5)*(vy(i,  j+1,k  ) + vy(i,j,k  )) * (pyx(i,  j+1,k  )-pyx(i,j,k  ))) :
                pz(i,j,k) - T(0.5)*dtdx * ( T(0.5)*(vx(i+1,j  ,k-1) + vx(i,j,k-1)) * (pxy(i+1,j  ,k-1)-pxy(i,j,k-1)))
                          - T(0.5)*dtdy * ( T(0.5)*(vy(i  ,j+1,k-1) + vy(i,j,k-1)) * (pyx(i  ,j+1,k-1)-pyx(i,j,k-1))) );

    fz(i,j,k) = vz(i,j,k)*f;
}

#endif
