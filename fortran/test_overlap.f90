! test_overlap.f90 -- the U2a gate: the Fortran integral engine must reproduce the
! python engine's S (adapted basis) and S-tilde (metric-scaled basis) for all six
! gate systems to 1e-10. Shell specs and reference matrices come from
! export_fortran.py; geometry is implicit in the exported shells' atom indices plus
! the coordinate table below (same constants as the exporter).
program test_overlap
   use gxtb_overlap
   implicit none
   integer, parameter :: dp = selected_real_kind(15)
   real(dp), parameter :: BOHR = 1.8897261254578281_dp
   real(dp), parameter :: pi = 3.14159265358979323846_dp
   character(len=8) :: names(6) = [character(len=8) :: &
                                   'h2', 'f2', 'hf', 'h2o', 'ch4', 'nh3']
   integer :: isys, ib, nfail
   real(dp) :: dmax

   nfail = 0
   do isys = 1, 6
      do ib = 1, 2
         call gate(trim(names(isys)), ib, dmax)
         if (dmax > 1.0e-10_dp) then
            nfail = nfail + 1
            print '(a,a,a,i0,a,es10.2,a)', '  ', trim(names(isys)), ' basis ', ib, &
               '  max|dS| ', dmax, '  MISMATCH'
         else
            print '(a,a,a,i0,a,es10.2)', '  ', trim(names(isys)), ' basis ', ib, &
               '  max|dS| ', dmax
         end if
      end do
   end do
   if (nfail == 0) then
      print '(a)', 'U2a GATE PASSED: S and S-tilde match python on all six systems'
   else
      print '(a,i0)', 'U2a GATE FAILED: ', nfail
      stop 1
   end if

contains

   subroutine geom(name, xyz, nat)
      character(*), intent(in) :: name
      real(dp), allocatable, intent(out) :: xyz(:, :)
      integer, intent(out) :: nat
      real(dp) :: ang, roh, a4, rnh, st, ct
      integer :: k
      ang = 104.5_dp*pi/180.0_dp
      roh = 0.9572_dp*BOHR
      a4 = 1.087_dp*BOHR/sqrt(3.0_dp)
      rnh = 1.012_dp*BOHR
      st = 0.9262_dp; ct = -0.3770_dp
      select case (name)
      case ('h2'); nat = 2; allocate (xyz(3, 2)); xyz = 0.0_dp
         xyz(3, 2) = 1.4_dp
      case ('f2'); nat = 2; allocate (xyz(3, 2)); xyz = 0.0_dp
         xyz(3, 2) = 2.668_dp
      case ('hf'); nat = 2; allocate (xyz(3, 2)); xyz = 0.0_dp
         xyz(3, 2) = 1.733_dp
      case ('h2o'); nat = 3; allocate (xyz(3, 3)); xyz = 0.0_dp
         xyz(1, 2) = roh*sin(ang/2); xyz(3, 2) = roh*cos(ang/2)
         xyz(1, 3) = -roh*sin(ang/2); xyz(3, 3) = roh*cos(ang/2)
      case ('ch4'); nat = 5; allocate (xyz(3, 5)); xyz = 0.0_dp
         xyz(:, 2) = [a4, a4, a4]; xyz(:, 3) = [a4, -a4, -a4]
         xyz(:, 4) = [-a4, a4, -a4]; xyz(:, 5) = [-a4, -a4, a4]
      case ('nh3'); nat = 4; allocate (xyz(3, 4)); xyz = 0.0_dp
         do k = 0, 2
            xyz(1, k + 2) = rnh*st*cos(2.0_dp*pi*k/3.0_dp)
            xyz(2, k + 2) = rnh*st*sin(2.0_dp*pi*k/3.0_dp)
            xyz(3, k + 2) = rnh*ct
         end do
      end select
   end subroutine

   subroutine gate(name, ibas, dmax)
      character(*), intent(in) :: name
      integer, intent(in) :: ibas
      real(dp), intent(out) :: dmax
      type(shell_t), allocatable :: shells(:)
      real(dp), allocatable :: xyz(:, :), s(:, :), sref(:, :)
      character(len=256) :: fs, fm
      integer :: u, nsh, i, j, nat, nao, nref
      write (fs, '(a,i0,a,a,a)') 'data/shells', ibas, '-', name, '.dat'
      write (fm, '(a,i0,a,a,a)') 'data/S', ibas, '-', name, '.dat'
      call geom(name, xyz, nat)
      open (newunit=u, file=trim(fs), status='old', action='read')
      read (u, *) nsh
      allocate (shells(nsh))
      do i = 1, nsh
         read (u, *) shells(i)%at, shells(i)%l, shells(i)%nprim
         allocate (shells(i)%exp(shells(i)%nprim), shells(i)%coef(shells(i)%nprim))
         do j = 1, shells(i)%nprim
            read (u, *) shells(i)%exp(j), shells(i)%coef(j)
         end do
      end do
      close (u)
      nao = nao_of(shells)
      allocate (s(nao, nao), sref(nao, nao))
      call overlap_matrix(shells, xyz, s)
      open (newunit=u, file=trim(fm), status='old', action='read')
      read (u, *) nref
      if (nref /= nao) then
         print *, 'nao mismatch for ', name, nref, nao
         stop 1
      end if
      do i = 1, nao
         read (u, *) (sref(i, j), j=1, nao)
      end do
      close (u)
      dmax = maxval(abs(s - sref))
      deallocate (shells, xyz, s, sref)
   end subroutine

end program test_overlap
