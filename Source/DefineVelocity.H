#include <AmrCoreAdv.H>
#include <Kernels.H>

#include <AMReX_MultiFabUtil.H>

template <typename MF>
void AmrCoreAdv::DefineVelocityAtLevel (int lev, amrex::Real a_time)
{
    using namespace amrex;
    using FAB = typename MF::fab_type;
    using RT = typename MF::value_type;

    RT time = RT(a_time);

    GpuArray<RT,AMREX_SPACEDIM> dx{AMREX_D_DECL(RT(geom[lev].CellSize(0)),
						RT(geom[lev].CellSize(1)),
					        RT(geom[lev].CellSize(2)))};

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
	FAB psifab;
        for (MFIter mfi(std::get<MF>(phi_new[lev]),TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {

        // ======== GET FACE VELOCITY =========
            GpuArray<Box, AMREX_SPACEDIM> nbx;
            AMREX_D_TERM(nbx[0] = mfi.nodaltilebox(0);,
                         nbx[1] = mfi.nodaltilebox(1);,
                         nbx[2] = mfi.nodaltilebox(2););

            AMREX_D_TERM(const Box& ngbxx = amrex::grow(mfi.nodaltilebox(0),nghost);,
                         const Box& ngbxy = amrex::grow(mfi.nodaltilebox(1),nghost);,
                         const Box& ngbxz = amrex::grow(mfi.nodaltilebox(2),nghost););

            GpuArray<Array4<RT>, AMREX_SPACEDIM> vel{ AMREX_D_DECL( std::get<MF>(facevel[lev][0]).array(mfi),
                                                                    std::get<MF>(facevel[lev][1]).array(mfi),
                                                                    std::get<MF>(facevel[lev][2]).array(mfi)) };

            const Box& psibox = Box(IntVect(AMREX_D_DECL(std::min(ngbxx.smallEnd(0)-1, ngbxy.smallEnd(0)-1),
                                                         std::min(ngbxx.smallEnd(1)-1, ngbxy.smallEnd(1)-1),
                                                         0)),
                                    IntVect(AMREX_D_DECL(std::max(ngbxx.bigEnd(0),   ngbxy.bigEnd(0)+1),
                                                         std::max(ngbxx.bigEnd(1)+1, ngbxy.bigEnd(1)),
                                                         0)));

	    psifab.resize(psibox, 1, The_Async_Arena());
            Array4<RT> psi = psifab.array();
            GeometryData geomdata = geom[lev].data();

            amrex::ParallelFor(psibox,
            [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                get_face_velocity_psi(i, j, time, psi, geomdata);
            });

            amrex::ParallelFor
                (AMREX_D_DECL(ngbxx,ngbxy,ngbxz),
                 AMREX_D_DECL(
                     [=] AMREX_GPU_DEVICE (int i, int j, int k)
                     {
                         get_face_velocity_x(i, j, k, vel[0], psi, dx[1]);
                     },
                     [=] AMREX_GPU_DEVICE (int i, int j, int k)
                     {
                         get_face_velocity_y(i, j, k, vel[1], psi, dx[0]);
                     },
                     [=] AMREX_GPU_DEVICE (int i, int j, int k)
                     {
                         get_face_velocity_z(i, j, k, vel[2]);
                     }));
        }
    }
}
