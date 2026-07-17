! gxtb_elem.f90 -- U1: the decoded element laws, verbatim from the python engine's
! elem() (gxtb_engine.py). Row map (decoded across the campaign): L2 = levels (row 2),
! k = L4 direct (row 4), L5 = MFX row (row 5), U = row 6, mu = row 7, L8 = row 8
! (kdiat sigma/pi at cols 1/2, the kb knob at col 4). s-rule: measured anchors for
! H/He, the period-2 line a*(z-2)+b for C/N/O/F. c_x = s_rule/cx_denom;
! kb = kb_scale*L8(4). nsh arrives as input until U2 ports the basis parser (labeled).
module gxtb_elem
   use gxtb_param, only: params_t, derived_t
   implicit none
   private
   public :: elem_t, make_elem

   integer, parameter :: dp = selected_real_kind(15)

   type :: elem_t
      integer :: z, nsh
      real(dp) :: cx, kb, kd_sg, kd_pi
      real(dp) :: k(3), l2(3), mu(3), u(3), l5(3)
   end type

contains

   function srule(z, d) result(s)
      integer, intent(in) :: z
      type(derived_t), intent(in) :: d
      real(dp) :: s
      select case (z)
      case (1); s = d%srule_h
      case (2); s = d%srule_he
      case default; s = d%srule_p2_a*real(z - 2, dp) + d%srule_p2_b
      end select
   end function

   function make_elem(z, nsh, p, d) result(e)
      integer, intent(in) :: z, nsh
      type(params_t), intent(in) :: p
      type(derived_t), intent(in) :: d
      type(elem_t) :: e
      integer :: l
      e%z = z; e%nsh = nsh
      e%cx = srule(z, d)/d%cx_denom
      e%kb = d%kb_scale*p%row(z, 8, 4)
      e%kd_sg = p%row(z, 8, 1)
      e%kd_pi = p%row(z, 8, 2)
      do l = 1, nsh
         e%k(l) = p%row(z, 4, l)
         e%l2(l) = p%row(z, 2, l)
         e%mu(l) = p%row(z, 7, l)
         e%u(l) = p%row(z, 6, l)
         e%l5(l) = p%row(z, 5, l)
      end do
   end function

end module gxtb_elem
