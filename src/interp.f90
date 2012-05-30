module interp
use Precision
use hm_init
implicit none


contains

real(dl) function interpf(x,y,x_0)
! returns interpolated value y_0 at point x_0 from arrays x and y
real(dl), dimension(:), intent(in) :: x,y 
real(dl), intent(in) :: x_0
real(dl), dimension(size(x)) :: y2

call spline(x, y, size(x), 1d40, 1d40, y2)
call splint(x, y, y2, size(x), x_0, interpf)
end function interpf

   SUBROUTINE spline(x, y, n, yp1, ypn, y2)
!   use nrtype
!
! Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e.
! y(i)=f(x(i)), with x(1)<x(2)<...<x(n), and given values yp1 and ypn for
! the first derivative of the interpolating function at points 1 and n,
! respectively, this routine returns an array y2(1:n) of length n which
! contains the second derivatives of the interpolating function at the
! tabulated points x(i).  If yp1 and/or ypn are equal to 1.e30 or larger,
! the routine is signaled to set the corresponding boundary condition for a
! natural spline with zero second derivative on that boundary.
! Parameter: nmax is the largest anticipiated value of n
! (adopted from Numerical Recipes in FORTRAN 77)
!
!   INTEGER, PARAMETER :: DP=KIND(1.0D0)
   INTEGER:: n
   INTEGER, PARAMETER:: nmax=500
   REAL(dl):: yp1, ypn, x(n), y(n), y2(n)
   INTEGER:: i, k
   REAL(dl):: p, qn, sig, un, u(nmax)

     if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
     else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
     endif

     do i=2, n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/&
             & (x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
     enddo

     if (ypn.gt..99e30) then
        qn=0.
        un=0.
     else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
     endif

     y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)

     do k=n-1, 1, -1
        y2(k)=y2(k)*y2(k+1)+u(k)
     enddo

     return
     END SUBROUTINE spline

   SUBROUTINE splint(xa, ya, y2a, n, x, y)
!   USE nrtype
!
! Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function
! (with the xa(i) in order), and given the array y2a(1:n), which is the output
! from the subroutine spline, and given a value of x, this routine returns a
! cubic spline interpolated value y.
! (adopted from Numerical Recipes in FORTRAN 77)
!
!   INTEGER, PARAMETER :: DP = KIND(1.0D0)
   INTEGER:: n
   REAL(dl):: x, y, xa(n), y2a(n), ya(n)
   INTEGER:: k, khi, klo
   REAL(dl):: a, b, h

     klo=1
     khi=n
1   if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if (xa(k).gt.x) then
           khi=k
        else
           klo=k
        endif
        goto 1
     endif

     h=xa(khi)-xa(klo)
     !if (h.eq.0.) write(*,*) 'bad xa input in splint'

     a=(xa(khi)-x)/h
     b=(x-xa(klo))/h
     y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.

     return
     END SUBROUTINE splint

    SUBROUTINE splie2(x1a,x2a,ya,m,n,y2a)
    INTEGER m,n,NN
    REAL(dl) x1a(m),x2a(n),y2a(m,n),ya(m,n)
    PARAMETER (NN=200)
    !USES spline
    INTEGER j,k
    REAL(dl) y2tmp(NN),ytmp(NN)
    do j=1,m
     do k=1,n
        ytmp(k)=ya(j,k)
     end do
      call spline(x2a,ytmp,n,1.d30,1.d30,y2tmp)
     do k=1,n
        y2a(j,k)=y2tmp(k)
     end do
     end do
     END SUBROUTINE splie2

     SUBROUTINE splin2(x1a,x2a,ya,y2a,m,n,x1,x2,y)
    INTEGER m,n,NN
    REAL(dl) x1,x2,y,x1a(m),x2a(n),y2a(m,n),ya(m,n)
    PARAMETER (NN=200)
    !USES spline,splint
    INTEGER j,k
    REAL(dl) y2tmp(NN),ytmp(NN),yytmp(NN)
    do j=1,m
     do k=1,n
        ytmp(k)=ya(j,k)
       y2tmp(k)=y2a(j,k)
     end do
      call splint(x2a,ytmp,y2tmp,n,x2,yytmp(j))
    end do
    call spline(x1a,yytmp,m,1.d30,1.d30,y2tmp)
    call splint(x1a,yytmp,y2tmp,m,x1,y)
    END SUBROUTINE splin2

end module interp
