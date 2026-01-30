#ifndef TAGGING_H
#define TAGGING_H

#include <AMReX_Array4.H>

template <typename T>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
state_error (int i, int j, int k,
             amrex::Array4<char> const& tag,
             amrex::Array4<T const> const& state,
             T phierr, char tagval)
{
    if (state(i,j,k) > phierr) {
        tag(i,j,k) = tagval;
    }
}

#endif
