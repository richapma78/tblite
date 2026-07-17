! test_elem.f90 -- the U1 gate: every element-law value must match the python
! engine's elem() (elem-reference.dat) to 1e-10. Bit-level agreement is expected
! (same inputs, same arithmetic); any miss is a decode or parse defect, never
! something to tune.
program test_elem
   use gxtb_param
   use gxtb_elem
   implicit none
   integer, parameter :: dp = selected_real_kind(15)
   type(params_t) :: p
   type(derived_t) :: d
   type(elem_t) :: e
   character(len=512) :: pfile, dfile, rfile, line
   integer :: u, ios, z, nsh, l, nfail, ntot
   real(dp) :: rcx, rkb, rsg, rpi, rk, rl2, rmu, ru, rl5

   call get_command_argument(1, pfile)
   call get_command_argument(2, dfile)
   call get_command_argument(3, rfile)
   call load_params(trim(pfile), p)
   call load_derived(trim(dfile), d)

   nfail = 0; ntot = 0
   open (newunit=u, file=trim(rfile), status='old', action='read')
   read (u, '(a)') line              ! header comment
   do
      read (u, *, iostat=ios) z, nsh, rcx, rkb, rsg, rpi
      if (ios /= 0) exit
      e = make_elem(z, nsh, p, d)
      call check(e%cx, rcx, z, 'cx')
      call check(e%kb, rkb, z, 'kb')
      call check(e%kd_sg, rsg, z, 'kd_sg')
      call check(e%kd_pi, rpi, z, 'kd_pi')
      do l = 1, nsh
         read (u, *) rk, rl2, rmu, ru, rl5
         call check(e%k(l), rk, z, 'k')
         call check(e%l2(l), rl2, z, 'L2')
         call check(e%mu(l), rmu, z, 'MU')
         call check(e%u(l), ru, z, 'U')
         call check(e%l5(l), rl5, z, 'L5')
      end do
   end do
   close (u)
   if (nfail == 0) then
      print '(a,i0,a)', 'U1 GATE PASSED: ', ntot, ' values match to 1e-10'
   else
      print '(a,i0,a,i0)', 'U1 GATE FAILED: ', nfail, ' of ', ntot
      stop 1
   end if

contains

   subroutine check(got, want, zz, tag)
      real(dp), intent(in) :: got, want
      integer, intent(in) :: zz
      character(*), intent(in) :: tag
      ntot = ntot + 1
      if (abs(got - want) > 1.0e-10_dp*max(1.0_dp, abs(want))) then
         nfail = nfail + 1
         print '(a,i0,1x,a,2es22.14)', '  MISMATCH z=', zz, tag, got, want
      end if
   end subroutine

end program test_elem
