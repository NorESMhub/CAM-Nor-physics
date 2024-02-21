module oslo_aero_ocean

   ! Module to replace the need for OSLO_AERO ifdef in NorESM physics

   implicit none
   public

contains

   subroutine oslo_aero_ocean_adv(state, pbuf2d)
      use ppgrid,         only : begchunk, endchunk
      use physics_types,  only : physics_state
      use physics_buffer, only : physics_buffer_desc

      type(physics_state), intent(in)    :: state(begchunk:endchunk)
      type(physics_buffer_desc), pointer :: pbuf2d(:,:)
   end subroutine oslo_aero_ocean_adv

end module oslo_aero_ocean
