! gxtb_overlap.f90 -- U2a: the overlap integral engine, verbatim port of
! prototype/overlap.py (Obara-Saika 1-D recursion, PySCF primitive-norm convention,
! contracted self-overlap normalization). Spherical = cartesian for l <= 1, which
! covers the v0.1 elements (s/p); d/f arrive with the spherical transform in a later
! rung. Shell inputs carry RAW adapted contraction coefficients; normalization
! happens here, exactly as the python engine does it.
module gxtb_overlap
   implicit none
   private
   public :: shell_t, overlap_matrix, nao_of

   integer, parameter :: dp = selected_real_kind(15)
   real(dp), parameter :: pi = 3.14159265358979323846_dp

   type :: shell_t
      integer :: at, l, nprim
      real(dp), allocatable :: exp(:), coef(:)
   end type

contains

   integer function nao_of(shells)
      type(shell_t), intent(in) :: shells(:)
      integer :: i
      nao_of = 0
      do i = 1, size(shells)
         nao_of = nao_of + 2*shells(i)%l + 1
      end do
   end function

   pure function prim_norm(a, l) result(nrm)
      real(dp), intent(in) :: a
      integer, intent(in) :: l
      real(dp) :: nrm, dfac
      integer :: k
      dfac = 1.0_dp
      do k = 2*l - 1, 1, -2
         dfac = dfac*real(k, dp)
      end do
      nrm = (2.0_dp*a/pi)**0.75_dp*(4.0_dp*a)**(real(l, dp)/2.0_dp)/sqrt(dfac)
   end function

   subroutine normalize(sh, cn)
      ! normalized-primitive coefficients, then unit contracted self-overlap
      type(shell_t), intent(in) :: sh
      real(dp), intent(out) :: cn(:)
      real(dp) :: dfac, s_rad, norm2
      integer :: i, j, k
      do i = 1, sh%nprim
         cn(i) = sh%coef(i)*prim_norm(sh%exp(i), sh%l)
      end do
      dfac = 1.0_dp
      do k = 2*sh%l - 1, 1, -2
         dfac = dfac*real(k, dp)
      end do
      norm2 = 0.0_dp
      do i = 1, sh%nprim
         do j = 1, sh%nprim
            s_rad = pi**1.5_dp*dfac/2.0_dp**sh%l &
                    /(sh%exp(i) + sh%exp(j))**(real(sh%l, dp) + 1.5_dp)
            norm2 = norm2 + cn(i)*cn(j)*s_rad
         end do
      end do
      cn = cn/sqrt(norm2)
   end subroutine

   subroutine os_1d(la, lb, pa, pb, mu, s)
      ! Obara-Saika 1-D table s(0:la, 0:lb)
      integer, intent(in) :: la, lb
      real(dp), intent(in) :: pa, pb, mu
      real(dp), intent(out) :: s(0:, 0:)
      real(dp) :: inv2mu
      integer :: i, j
      s = 0.0_dp
      s(0, 0) = 1.0_dp
      inv2mu = 0.5_dp/mu
      do i = 1, la
         s(i, 0) = pa*s(i - 1, 0)
         if (i > 1) s(i, 0) = s(i, 0) + real(i - 1, dp)*inv2mu*s(i - 2, 0)
      end do
      do j = 1, lb
         do i = 0, la
            s(i, j) = pb*s(i, j - 1)
            if (i > 0) s(i, j) = s(i, j) + real(i, dp)*inv2mu*s(i - 1, j - 1)
            if (j > 1) s(i, j) = s(i, j) + real(j - 1, dp)*inv2mu*s(i, j - 2)
         end do
      end do
   end subroutine

   subroutine cart_prims(a, b, ra, rb, la, lb, blk)
      ! cartesian primitive-pair overlap block; component order for l=1 is (x,y,z)
      ! (PySCF's p ordering, which the oracle shares)
      real(dp), intent(in) :: a, b, ra(3), rb(3)
      integer, intent(in) :: la, lb
      real(dp), intent(out) :: blk(:, :)
      real(dp) :: mu, p(3), ab2, pref
      real(dp) :: t1(0:la, 0:lb), t2(0:la, 0:lb), t3(0:la, 0:lb)
      integer :: ia, ib, ax(3), bx(3)
      mu = a + b
      p = (a*ra + b*rb)/mu
      ab2 = sum((ra - rb)**2)
      pref = exp(-a*b/mu*ab2)*(pi/mu)**1.5_dp
      call os_1d(la, lb, p(1) - ra(1), p(1) - rb(1), mu, t1)
      call os_1d(la, lb, p(2) - ra(2), p(2) - rb(2), mu, t2)
      call os_1d(la, lb, p(3) - ra(3), p(3) - rb(3), mu, t3)
      do ia = 1, ncart(la)
         ax = comp(la, ia)
         do ib = 1, ncart(lb)
            bx = comp(lb, ib)
            blk(ia, ib) = pref*t1(ax(1), bx(1))*t2(ax(2), bx(2))*t3(ax(3), bx(3))
         end do
      end do
   end subroutine

   pure integer function ncart(l)
      integer, intent(in) :: l
      ncart = (l + 1)*(l + 2)/2
   end function

   pure function comp(l, i) result(c)
      ! l = 0: (0,0,0); l = 1: x, y, z. (d/f need the full lexicographic list +
      ! spherical transform -- out of U2a's s/p scope, guarded by the caller.)
      integer, intent(in) :: l, i
      integer :: c(3)
      c = 0
      if (l == 1) c(i) = 1
   end function

   subroutine overlap_matrix(shells, xyz, s)
      type(shell_t), intent(in) :: shells(:)
      real(dp), intent(in) :: xyz(:, :)          ! (3, natoms), Bohr
      real(dp), intent(out) :: s(:, :)
      integer :: na, nb, ia, ib, ip, jp, oa, ob
      real(dp) :: blk(3, 3), acc(3, 3)
      real(dp), allocatable :: cna(:), cnb(:)
      s = 0.0_dp
      oa = 0
      do ia = 1, size(shells)
         na = 2*shells(ia)%l + 1
         allocate (cna(shells(ia)%nprim))
         call normalize(shells(ia), cna)
         ob = 0
         do ib = 1, ia
            nb = 2*shells(ib)%l + 1
            allocate (cnb(shells(ib)%nprim))
            call normalize(shells(ib), cnb)
            acc = 0.0_dp
            do ip = 1, shells(ia)%nprim
               do jp = 1, shells(ib)%nprim
                  call cart_prims(shells(ia)%exp(ip), shells(ib)%exp(jp), &
                                  xyz(:, shells(ia)%at), xyz(:, shells(ib)%at), &
                                  shells(ia)%l, shells(ib)%l, blk(1:na, 1:nb))
                  acc(1:na, 1:nb) = acc(1:na, 1:nb) + cna(ip)*cnb(jp)*blk(1:na, 1:nb)
               end do
            end do
            s(oa + 1:oa + na, ob + 1:ob + nb) = acc(1:na, 1:nb)
            s(ob + 1:ob + nb, oa + 1:oa + na) = transpose(acc(1:na, 1:nb))
            deallocate (cnb)
            ob = ob + nb
         end do
         deallocate (cna)
         oa = oa + na
      end do
   end subroutine

end module gxtb_overlap
