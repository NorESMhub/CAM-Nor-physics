module oslo_aero_share

   ! Module to replace the need for OSLO_AERO ifdef in NorESM physics

   implicit none
   public

   integer :: nbmodes = 0
   logical :: use_oslo_aero = .false.

end module oslo_aero_share
