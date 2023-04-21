
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine jtotk(ik,pmat,evecsv,evecsvt)
use modmain
implicit none
! arguments
integer, intent(in) :: ik
complex(8), intent(in) :: pmat(nstsv,nstsv,3)
complex(8), intent(in) :: evecsv(nstsv,nstsv),evecsvt(nstsv,nstsv)
! local variables
integer ist,l
real(8) wo
complex(8) z1
! allocatable arrays
complex(8), allocatable :: a(:,:),b(:,:)
! external functions
complex(8), external :: zdotc
allocate(a(nstsv,nstsv),b(nstsv,nstsv))
do l=1,3
! form the momentum operator matrix elements in the first-variational basis
  call zgemm('N','C',nstsv,nstsv,nstsv,zone,pmat(:,:,l),nstsv,evecsv,nstsv, &
   zzero,a,nstsv)
  call zgemm('N','N',nstsv,nstsv,nstsv,zone,evecsv,nstsv,a,nstsv,zzero,b,nstsv)
! add to the total current
  call zgemm('N','N',nstsv,nstsv,nstsv,zone,b,nstsv,evecsvt,nstsv,zzero,a,nstsv)
  do ist=1,nstsv
    wo=occsv(ist,ik)
    if (abs(wo).lt.epsocc) cycle
    wo=wo*wkpt(ik)
    z1=zdotc(nstsv,evecsvt(:,ist),1,a(:,ist),1)
!$OMP ATOMIC
    jtot(l)=jtot(l)+wo*dble(z1)
  end do
end do
deallocate(a,b)
end subroutine

