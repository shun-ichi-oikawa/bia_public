subroutine advance_hmt ( Y, eps, dt, N, Nhmt )
	implicit none
	integer :: i, N, Nhmt, sw = 0
	real(8) :: eps, dt
	real(8), dimension(6*N) :: Y
	real(8), dimension(3,N) :: r0, v0
	real(8), dimension(3,N) :: r1, v1
!	if (sw == 0) then
	r0(1,1:N) = Y(1:6*N-5:6);	r0(2,1:N) = Y(2:6*N-4:6); r0(3,1:N) = Y(3:6*N-3:6)
	v0(1,1:N) = Y(4:6*N-2:6);	v0(2,1:N) = Y(5:6*N-1:6); v0(3,1:N) = Y(6:6*N-0:6)
	r1 = r0; v1 = v0;
!	else
		sw = 1
!	end if
!print '(a,es11.3)', "x   [0]=", r0(1,1)
	do i = 1, Nhmt
		call Hermite(N,eps,dt/Nhmt,r1,v1)
	end do
	Y(1:6*N-5:6) = r1(1,1:N); Y(2:6*N-4:6) = r1(2,1:N); Y(3:6*N-3:6) = r1(3,1:N)
	Y(4:6*N-2:6) = v1(1,1:N); Y(5:6*N-1:6) = v1(2,1:N); Y(6:6*N-0:6) = v1(3,1:N)
	return
	contains
	subroutine Hermite(N,eps,dt,r,v)
		implicit none
		integer N
		real(8), dimension(3,N) :: r , v
		real(8), dimension(3,N) :: r0, v0, a0, j0
		real(8), dimension(3,N) :: rp, vp, ap, jp
		real(8)                 :: dt, eps
		r0 = r; v0 = v
		call make_aj(N,eps,r0,v0,a0,j0)
		! 予測子
		rp = r0 + ( v0 + ( a0 + j0 * dt / 3 ) * dt / 2 ) * dt
		vp =        v0 + ( a0 + j0 * dt / 2 ) * dt
		call make_aj(N,eps,rp,vp,ap,jp)
		! 修正子
		v  = v0 + ( a0 + ap + ( j0 - jp ) * dt / 6 ) * dt / 2
		r  = r0 + ( v0 + v  + ( a0 - ap ) * dt / 6 ) * dt / 2
	end subroutine Hermite
	subroutine make_aj(N,eps,r,v,a,j)
		integer :: N, i, k
		real(8), dimension(3,N) :: r, v, a, j
		real(8), dimension(3  ) :: r_ik, g_ik, daik, djik
		real(8)                 :: dist2, dist3, dt, eps, g_r=0
		a = 0; j = 0
		do i = 1, N - 1
			do k = i + 1, N
				r_ik(:) = r(:,i) - r(:,k); dist2 = sum(r_ik**2); dist3 = dist2 * sqrt(dist2)
				daik(:) = eps * r_ik(:) / dist3
				a (:,i) = a(:,i) + daik
				a (:,k) = a(:,k) - daik
				g_ik(:) = v(:,i) - v (:,k)
				g_r = sum( g_ik * r_ik )
				djik(:) = eps * ( g_ik(:) - 3 * sum( g_ik * r_ik ) * r_ik(:) / dist2 ) / dist3
				j (:,i) = j(:,i) + djik(:)
				j (:,k) = j(:,k) - djik(:)
			end do
		end do
	end subroutine make_aj
end subroutine advance_hmt
