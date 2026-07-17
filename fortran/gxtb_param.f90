! gxtb_param.f90 -- U1: parse the g-xTB parameter file (two global rows, then per
! element: a single-integer header line followed by nine rows L1..L9 of reals) and
! the derived-constants export (key-value). No method constants live in this code:
! everything loads from the two files at start (the house rule).
module gxtb_param
   implicit none
   private
   public :: params_t, load_params, derived_t, load_derived

   integer, parameter :: dp = selected_real_kind(15)
   integer, parameter :: maxcol = 16, maxz = 103

   type :: params_t
      real(dp) :: g1(maxcol) = 0.0_dp, g2(maxcol) = 0.0_dp
      integer :: ng1 = 0, ng2 = 0
      real(dp) :: row(maxz, 9, maxcol) = 0.0_dp
      integer :: ncol(maxz, 9) = 0
      logical :: has(maxz) = .false.
   end type

   type :: derived_t
      real(dp) :: srule_p2_a, srule_p2_b, srule_h, srule_he
      real(dp) :: kb_scale, cx_denom
   end type

contains

   subroutine load_params(path, p)
      character(*), intent(in) :: path
      type(params_t), intent(out) :: p
      character(len=4096) :: line
      integer :: u, ios, z, nrow, nt
      real(dp) :: vals(maxcol)
      z = 0; nrow = 0
      open (newunit=u, file=path, status='old', action='read')
      do
         read (u, '(a)', iostat=ios) line
         if (ios /= 0) exit
         if (len_trim(line) == 0) cycle
         call tokens(line, vals, nt)
         if (nt == 1 .and. is_header(line)) then
            z = nint(vals(1)); nrow = 0
            p%has(z) = .true.
         else if (z == 0) then
            if (p%ng1 == 0) then
               p%g1(1:nt) = vals(1:nt); p%ng1 = nt
            else
               p%g2(1:nt) = vals(1:nt); p%ng2 = nt
            end if
         else
            nrow = nrow + 1
            if (nrow <= 9) then
               p%row(z, nrow, 1:nt) = vals(1:nt)
               p%ncol(z, nrow) = nt
            end if
         end if
      end do
      close (u)
   end subroutine

   logical function is_header(line)
      character(*), intent(in) :: line
      character(len=:), allocatable :: t
      integer :: i
      t = trim(adjustl(line))
      is_header = len(t) > 0
      do i = 1, len(t)
         if (index('0123456789', t(i:i)) == 0) is_header = .false.
      end do
   end function

   subroutine tokens(line, vals, nt)
      character(*), intent(in) :: line
      real(dp), intent(out) :: vals(maxcol)
      integer, intent(out) :: nt
      integer :: ios
      vals = 0.0_dp
      do nt = maxcol, 1, -1
         read (line, *, iostat=ios) vals(1:nt)
         if (ios == 0) return
      end do
      nt = 0
   end subroutine

   subroutine load_derived(path, d)
      character(*), intent(in) :: path
      type(derived_t), intent(out) :: d
      character(len=256) :: line, key
      real(dp) :: val
      integer :: u, ios
      open (newunit=u, file=path, status='old', action='read')
      do
         read (u, '(a)', iostat=ios) line
         if (ios /= 0) exit
         if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
         read (line, *, iostat=ios) key, val
         if (ios /= 0) cycle
         select case (trim(key))
         case ('srule_p2_a'); d%srule_p2_a = val
         case ('srule_p2_b'); d%srule_p2_b = val
         case ('srule_h'); d%srule_h = val
         case ('srule_he'); d%srule_he = val
         case ('kb_scale'); d%kb_scale = val
         case ('cx_denom'); d%cx_denom = val
         end select
      end do
      close (u)
   end subroutine

end module gxtb_param
