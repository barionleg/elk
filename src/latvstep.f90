
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine latvstep
use modmain
use modmpi
implicit none
! local variables
integer i
real(8) t1
do i=1,nstrain
! compute product of current and previous stress tensor components
  t1=stress(i)*stressp(i)
! if component is in the same direction then increase step size parameter
  if (t1.gt.0.d0) then
    taulatv(i)=taulatv(i)+tau0latv
  else
    taulatv(i)=tau0latv
  end if
  t1=taulatv(i)*(stress(i)+stressp(i))
  avec(:,:)=avec(:,:)-t1*strain(:,:,i)
end do
! compute the new unit cell volume
call reciplat(avec,bvec,omega,omegabz)
! preserve the volume if required (added by Ronald Cohen)
if (latvopt.eq.2) then
  t1=(omega0/omega)**(1.d0/3.d0)
  avec(:,:)=t1*avec(:,:)
  omega=omega0
end if
! each MPI process should have numerically identical lattice vectors
call mpi_bcast(avec,9,mpi_double_precision,0,mpicom,ierror)
end subroutine

