module oslo_aero_microp

   ! Module to replace the need for OSLO_AERO ifdef in NorESM physics

   implicit none
   public

contains

   subroutine oslo_aero_microp_register()
   end subroutine oslo_aero_microp_register

   subroutine oslo_aero_microp_init()
   end subroutine oslo_aero_microp_init

  subroutine oslo_aero_microp_run (state, ptend_all, deltatin, pbuf)
     use shr_kind_mod,   only: r8=>shr_kind_r8
     use physics_types,  only: physics_state, physics_ptend
     use physics_buffer, only: physics_buffer_desc

     type(physics_state),         intent(in)    :: state
     type(physics_ptend),         intent(out)   :: ptend_all
     real(r8),                    intent(in)    :: deltatin
     type(physics_buffer_desc),   pointer       :: pbuf(:)
  end subroutine oslo_aero_microp_run

end module oslo_aero_microp
