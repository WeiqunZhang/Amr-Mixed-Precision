#ifndef FACE_VELOCITY_H_
#define FACE_VELOCITY_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>

template <typename T>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void get_face_velocity_psi(int i, int j,
                           const T time,
                           amrex::Array4<T> const& psi,
                           amrex::GeometryData const& geomdata)
{
    using namespace amrex;
    constexpr T PI = T(3.1415926535897932384626);

    const T prob_lo_x = T(geomdata.ProbLo(0));
    const T prob_lo_y = T(geomdata.ProbLo(1));
    const T dx        = T(geomdata.CellSize(0));
    const T dy        = T(geomdata.CellSize(1));

    T y = dx*(T(0.5)+j) + prob_lo_y;
    T x = dy*(T(0.5)+i) + prob_lo_x;
    psi(i,j,0) = std::pow(std::sin(PI*x), 2) * std::pow(std::sin(PI*y), 2)
                * std::cos(PI*time/T(2.0)) * (T(1.0)/PI);
}

template <typename T>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void get_face_velocity_x(int i, int j, int k,
                         amrex::Array4<T> const& vx,
                         amrex::Array4<T> const& psi,
                         T dx, T dy, amrex::Array4<T const> const& phi)
{
    vx(i,j,k) = -( (psi(i,j+1,0)+psi(i-1,j+1,0)) - (psi(i,j-1,0)+psi(i-1,j-1,0)) ) * (T(0.25)/dy);
    T gyp = (phi(i  ,j+1,k)-phi(i  ,j-1,k))/dy*T(0.5);
    T gym = (phi(i-1,j+1,k)-phi(i-1,j-1,k))/dy*T(0.5);
    T cga = (gyp - gym)/dx;
    vx(i,j,k) += cga;
    T gxp = (phi(i,j+1,k)-phi(i-1,j+1,k))/dx;
    T gx0 = (phi(i,j  ,k)-phi(i-1,j  ,k))/dx;
    T gxm = (phi(i,j-1,k)-phi(i-1,j-1,k))/dx;
    T cgb = ((gxp - gx0)*(T(1)+T(1.5)*std::numeric_limits<T>::epsilon()) +
             (gx0 - gxm)*(T(1)-T(1.5)*std::numeric_limits<T>::epsilon())) / (T(2)*dy);
    vx(i,j,k) -= cgb;
}

template <typename T>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void get_face_velocity_y(int i, int j, int k,
                         amrex::Array4<T> const& vy,
                         amrex::Array4<T> const& psi,
                         T dx, T dy, amrex::Array4<T const> const& phi)
{
    vy(i,j,k) =  ( (psi(i+1,j,0)+psi(i+1,j-1,0)) - (psi(i-1,j,0)+psi(i-1,j-1,0)) ) * (T(0.25)/dx);
    T cga = ((phi(i+1,j,k)-phi(i-1,j,k))/dx*T(0.5) - (phi(i+1,j-1,k)-phi(i-1,j-1,k))/dx*T(0.5))/dy;
    vy(i,j,k) += cga;
    T gyp = (phi(i+1,j,k)-phi(i+1,j-1,k))/dy;
    T gy0 = (phi(i  ,j,k)-phi(i  ,j-1,k))/dy;
    T gym = (phi(i-1,j,k)-phi(i-1,j-1,k))/dy;
    T cgb = ((gyp - gy0)*(T(1)+T(1.5)*std::numeric_limits<T>::epsilon()) +
             (gy0 - gym)*(T(1)-T(1.5)*std::numeric_limits<T>::epsilon())) / (T(2)*dx);
    vy(i,j,k) -= cgb;
}

template <typename T>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void get_face_velocity_z(int i, int j, int k,
                         amrex::Array4<T> const& vz)
{
    vz(i,j,k) =  T(0.0);
}

#endif
