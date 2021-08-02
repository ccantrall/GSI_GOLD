module fuv_setup
    use kinds, only: i_kind
    implicit none
    private
    public:: setup 
          interface setup; module procedure setupfuv; end interface

contains
subroutine setupfuv(obsLL,odiagLL,lunin,mype,bwork,awork,nele,nobs,is,fuv_diagsave,init_pass)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    setupfuv     compute rhs of oi for far ultraviolet airglow
!   prgmmr: cantrall
!   date: 2021-07-21
!
! abstract:  For assimilation of far ultraviolet airglow observations 
!            (GOLD)
!            this routine:
!              a) reads obs assigned to given mpi task (geographic region),
!              b) simulates obs from guess,
!              c) apply some quality control to obs,
!              d) load weight and innovation arrays used in minimization
!              e) collects statistics for runtime diagnostic output
!              f) writes additional diagnostic information to output file
!
! program history log:
!   2021-07-21  cantrall  -  first version of setupfuv: 
!---
!
!   input argument list:
!     lunin    - unit from which to read observations
!     mype     - mpi task id
!     nele     - number of data elements per observation
!     nobs     - number of observations
!
!   output argument list:
!     bwork    - array containing information about obs-ges statistics
!     awork    - array containing information for data counts and gross checks
!
! attributes:
!   language: Fortran 90 and/or above
!   machine:  
!
!$$$                          
    use mpeu_util, only: die,perr
    use kinds, only: r_kind,r_single,r_double,i_kind
    
    use guess_grids, only: nfldsig
    !use guess_grids, only: hrdifsig,
    !use gridmod, only: dx_gfs
    !use gridmod, only: region_dx,region_dy      ! dx, dy (:,:)
    !use gridmod, only: wrf_mass_regional
  !--
    use gridmod, only: nsig
    !use gridmod, only: lat2,lon2,get_ij,nlat_sfc,nlon_sfc
    !use gridmod, only: regional,nsig, &
    !                   eta1_ll,pt_ll,aeta1_ll
    !use gridmod, only: latlon11
  !--
    use m_obsdiagNode, only: obs_diag
    use m_obsdiagNode, only: obs_diags
    !use m_obsdiagNode, only: obsdiagLList_nextNode
    !use m_obsdiagNode, only: obsdiagNode_set
    !use m_obsdiagNode, only: obsdiagNode_get
    !use m_obsdiagNode, only: obsdiagNode_assert
  
    use obsmod, only: rmiss_single,lobsdiagsave
    !use obsmod, only: nobskeep,lobsdiag_allocated
    use obsmod, only: netcdf_diag, binary_diag, dirname, ianldate, time_offset
    
    use nc_diag_write_mod, only: nc_diag_init, nc_diag_header, nc_diag_metadata, &
         nc_diag_write, nc_diag_data2d
    *********use nc_diag_read_mod, only: nc_diag_read_init, nc_diag_read_get_dim, nc_diag_read_close
    
    use obsmod, only: luse_obsdiag
    
    use m_obsNode, only: obsNode
    
    use m_fuvNode, only: fuvNode
    use m_fuvNode, only: fuvNode_appendto
    
    use m_obsLList , only: obsLList
    
    use gsi_4dvar, only: nobs_bins,hr_obsbin
    
    use constants, only: zero,one,r1000, &
         tiny_r_kind,three,half,two,cg_term,huge_single,&
         wgtlim, qcmin, rad2deg
    use constants, only: one_tenth,qmin,ten,t0c,five,r0_05,r60
    
    use jfunc, only: jiter,last,miter
    
    use qcmod, only: dfact,dfact1,npres_print
    
    use goldinfo, only: ctype_gold,ntype_gold,itype_gold, &
                        isubtype_gold,typeuse_gold,grosserr_gold
    
    use m_dtime, only: dtime_setup, dtime_check
  !--
    use gsi_bundlemod, only: gsi_bundlegetpointer
    
    use gsi_metguess_mod, only: gsi_metguess_get,gsi_metguess_bundle
  
    use mpimod, only: ierror,mpi_comm_world,mpi_rtype,mpi_itype,mpi_sum

    use state_vectors, only: nsdim
  !--
  !--
  
    implicit none
  
  ! Declare passed variables
    type(obsLList ),target,dimension(:),intent(in):: obsLL
    type(obs_diags),target,dimension(:),intent(in):: odiagLL
  
    logical                                           ,intent(in   ) :: fuv_diagsave
    integer(i_kind)                                   ,intent(in   ) :: lunin,mype,nele,nobs
    real(r_kind),dimension(100+7*nsig)                ,intent(inout) :: awork
    real(r_kind),dimension(npres_print,ntype_gold,5,3),intent(inout) :: bwork
    integer(i_kind)                                   ,intent(in   ) :: is ! ndat index
    logical                                           ,intent(in   ) :: init_pass
  
  ! Declare local parameter
    character(len=*),parameter:: myname="setupfuv"
  
  ! Declare local variables
    real(r_kind):: fuvges0,fuvges,grsmlt,dlat,dlon,dtime,obserror, &
                   obserrlm,residual,ratio,dfuv
    real(r_kind) error,ddiff 
    real(r_kind) ressw2,ress,scale,val2,val,valqc
    real(r_kind) rat_err2,exp_arg,term,ratio_errors,rwgt
!    real(r_kind) cg_light,wgross,wnotgross,wgt,arg
    real(r_kind) errinv_input,errinv_adjst,errinv_final
    real(r_kind) err_input,err_adjst,err_final,tfact
    real(r_kind) rfuv,presw
    real(r_kind),dimension(nele,nobs):: data
    real(r_single),allocatable,dimension(:)::cdiagbuf
    real(r_single),allocatable,dimension(:,:)::rdiagbuf

    ! Local variables
!    integer(i_kind)                   :: it,k,istatus,ier,nsig_read
!    real(r_kind), pointer             :: flashrate  (:,:,:)  ! lightning flash rate
!    real(r_kind), pointer             :: flashrate_h(:,:,:)  ! lightning flash rate
!    real(r_kind), pointer             :: htot_h    (:,:,:)  ! lightning flash rate, non-h, cloud-res
!    real(r_kind), pointer             :: dx  (:,:)  ! 
!    real(r_kind), pointer             :: dy  (:,:)  ! 
!    real(r_kind),allocatable          :: sigmadot(:,:,:,:)  !! vert. vel in sigma
  !----
  
    !integer(i_kind),allocatable       :: kbot(:)    
!    real(r_kind),allocatable          :: kbot(:,:,:)    
  
!    real(r_kind),allocatable          :: kvert(:,:,:)
!    real(r_kind)                      :: sum_loc,sum_gbl
!    real(r_kind)                      :: r0,w0
!    real(r_kind)                      :: eps
!    real(r_kind)                      :: eps0
!    real(r_kind),dimension(lat2,lon2,nsig,nfldsig)      :: cwgues
  
!    integer(i_kind),dimension(12)     :: light_ij
!    integer(i_kind)                   :: ix,ixp,iy,iyp
!    integer(i_kind)                   :: jtime,jtimep
  !---
!    integer(i_kind) ikxx,nn,ibin,ioff
    integer(i_kind) i,nchar,nreal,j,jj,ii,l,mm1,im,jm,km,ioff
    integer(i_kind) ilon,ilat,ihgt,ifuvob,ioza,isza
    integer(i_kind) ier2,iuse,ilate,ilone,ikx
    integer(i_kind) satellite_id
!    integer(i_kind) nobs_loc,nobs_gbl
  
!    logical,allocatable               :: wmaxflag(:,:,:)
    logical,dimension(nobs):: luse,muse
    integer(i_kind),dimension(nobs):: ioid ! initial (pre-distribution) obs ID
    logical proceed
  
  ! Declare external calls for code analysis
    external:: mpi_barrier
    external:: mpi_allreduce
    external:: mpi_finalize 
    external:: mpi_reduce
!    external:: sumslightbias
     
      
  ! File(s) for postprocessing
!   character :: post_file*40
  
!    logical:: in_curbin, in_anybin
    integer(i_kind) :: istat
    type(fuvNode),pointer:: my_head
    type(obs_diag ),pointer:: my_diag
    type(obs_diags),pointer:: my_diagLL
  
  ! Guess fields
    real(r_kind),allocatable,dimension(:,:,:  ) :: ges_ps
    real(r_kind),allocatable,dimension(:,:,:  ) :: ges_z
    real(r_kind),allocatable,dimension(:,:,:,:) :: ges_o
    real(r_kind),allocatable,dimension(:,:,:,:) :: ges_o2
    real(r_kind),allocatable,dimension(:,:,:,:) :: ges_tv

  ! Derived fields
    real(r_kind),allocatable,dimension(:,:,:,:) :: der_n2
  
 !   type(obsLList),pointer,dimension(:):: lighthead
 !   lighthead => obsLL(:)
  !--
  
!    grsmlt=three  ! multiplier factor for gross check, an appropriate magnitude
                  ! is yet to be determined.
!    mm1=mype+1
!    scale=one
  
  ! Check to see if required guess fields are available
    call check_vars_(proceed)
    if (.not.proceed) return  ! not all vars available, simply return
  
  ! If require guess vars available, extract from bundle ...
    call init_vars_

! ***** code ******

! ***** code ******
! End of routine
  return
  contains

  subroutine check_vars_ (proceed)
  logical,intent(inout) :: proceed
  integer(i_kind) ivar, istatus
! Check to see if required guess fields are available
  call gsi_metguess_get ('var::z' , ivar, istatus )
  proceed=proceed.and.ivar>0
  call gsi_metguess_get ('var::tv', ivar, istatus )
  proceed=proceed.and.ivar>0
  call gsi_metguess_get ('var::o' , ivar, istatus )
  proceed=proceed.and.ivar>0
  call gsi_metguess_get ('var::o2' , ivar, istatus )
  proceed=proceed.and.ivar>0
  
  end subroutine check_vars_

  subroutine init_vars_

    real(r_kind),dimension(:,:  ),pointer:: rank2=>NULL()
    real(r_kind),dimension(:,:,:),pointer:: rank3=>NULL()
    character(len=5) :: varname
    integer(i_kind) ifld, istatus
  
  ! If require guess vars available, extract from bundle ...
    if(size(gsi_metguess_bundle)==nfldsig) then
  !    get z ...
       varname='z'
       call gsi_bundlegetpointer(gsi_metguess_bundle(1),trim(varname),rank2,istatus)
       if (istatus==0) then
           if(allocated(ges_z))then
              write(6,*) trim(myname), ': ', trim(varname), ' already incorrectly alloc '
              call stop2(999)
           endif
           allocate(ges_z(size(rank2,1),size(rank2,2),nfldsig))
           ges_z(:,:,1)=rank2
           do ifld=2,nfldsig
              call gsi_bundlegetpointer(gsi_metguess_bundle(ifld),trim(varname),rank2,istatus)
              ges_z(:,:,ifld)=rank2
           enddo
       else
           write(6,*) trim(myname),': ', trim(varname), ' not found in met bundle, ier= ',istatus
           call stop2(999)
       endif
  !    get tv ...
       varname='tv'
       call gsi_bundlegetpointer(gsi_metguess_bundle(1),trim(varname),rank3,istatus)
       if (istatus==0) then
           if(allocated(ges_tv))then
              write(6,*) trim(myname), ': ', trim(varname), ' already incorrectly alloc '
              call stop2(999)
           endif
           allocate(ges_tv(size(rank3,1),size(rank3,2),size(rank3,3),nfldsig))
           ges_tv(:,:,:,1)=rank3
           do ifld=2,nfldsig
              call gsi_bundlegetpointer(gsi_metguess_bundle(ifld),trim(varname),rank3,istatus)
              ges_tv(:,:,:,ifld)=rank3
           enddo
       else
           write(6,*) trim(myname),': ', trim(varname), ' not found in met bundle, ier= ',istatus
           call stop2(999)
       endif
  !   get o
    varname='o'
       call gsi_bundlegetpointer(gsi_metguess_bundle(1),trim(varname),rank3,istatus)
       if (istatus==0) then
           if(allocated(ges_o))then
              write(6,*) trim(myname), ': ', trim(varname), ' already incorrectly alloc '
              call stop2(999)
           endif
           allocate(ges_o(size(rank3,1),size(rank3,2),size(rank3,3),nfldsig))
           ges_o(:,:,:,1)=rank3
           do ifld=2,nfldsig
              call gsi_bundlegetpointer(gsi_metguess_bundle(ifld),trim(varname),rank3,istatus)
              ges_o(:,:,:,ifld)=rank3
           enddo
       else
           write(6,*) trim(myname),': ', trim(varname), ' not found in met bundle, ier= ',istatus
           call stop2(999)
       endif
  !   get o2
       varname='o2'
       call gsi_bundlegetpointer(gsi_metguess_bundle(1),trim(varname),rank3,istatus)
       if (istatus==0) then
           if(allocated(ges_o2))then
              write(6,*) trim(myname), ': ', trim(varname), ' already incorrectly alloc '
              call stop2(999)
           endif
           allocate(ges_o2(size(rank3,1),size(rank3,2),size(rank3,3),nfldsig))
           ges_o2(:,:,:,1)=rank3
           do ifld=2,nfldsig
              call gsi_bundlegetpointer(gsi_metguess_bundle(ifld),trim(varname),rank3,istatus)
              ges_o2(:,:,:,ifld)=rank3
           enddo
       else
           write(6,*) trim(myname),': ', trim(varname), ' not found in met bundle, ier= ',istatus
           call stop2(999)
       endif
    endif
!    get ps ...
    varname='ps'
    call gsi_bundlegetpointer(gsi_metguess_bundle(1),trim(varname),rank2,istatus)
    if (istatus==0) then
        if(allocated(ges_ps))then
           write(6,*) trim(myname), ': ', trim(varname), ' already incorrectly alloc '
           call stop2(999)
        endif
        allocate(ges_ps(size(rank2,1),size(rank2,2),nfldsig))
        ges_ps(:,:,1)=rank2
        do ifld=2,nfldsig
           call gsi_bundlegetpointer(gsi_metguess_bundle(ifld),trim(varname),rank2,istatus)
           ges_ps(:,:,ifld)=rank2
        enddo
    else
        write(6,*) trim(myname),': ', trim(varname), ' not found in met bundle, ier= ',istatus
        call stop2(999)
    endif
  end subroutine init_vars_
  
  subroutine init_netcdf_diag_
  character(len=80) string
  character(len=128) diag_fuv_file
  integer(i_kind) ncd_fileid,ncd_nobs
  logical append_diag
  logical,parameter::verbose=.false.
     write(string,900) jiter
900  format('gold_fuv_',i2.2,'.nc4')
     diag_fuv_file=trim(dirname) // trim(string)

     inquire(file=diag_fuv_file, exist=append_diag)

     if (append_diag) then
        call nc_diag_read_init(diag_fuv_file,ncd_fileid)
        ncd_nobs = nc_diag_read_get_dim(ncd_fileid,'nobs')
        call nc_diag_read_close(diag_fuv_file)

        if (ncd_nobs > 0) then
           if(verbose) print *,'file ' // trim(diag_fuv_file) // ' exists. Appending.  nobs,mype=',ncd_nobs,mype
        else
           if(verbose) print *,'file ' // trim(diag_fuv_file) // ' exists but contains no obs.  Not appending. nobs,mype=',ncd_nobs,mype
           append_diag = .false. ! if there are no obs in existing file, then do not try to append
        endif
     end if

     call nc_diag_init(diag_fuv_file, append=append_diag)

     if (.not. append_diag) then ! don't write headers on append - the module will break?
        call nc_diag_header("date_time",ianldate )
        call nc_diag_header("Number_of_state_vars", nsdim          )
     endif
  end subroutine init_netcdf_diag_

  subroutine contents_binary_diag_(odiag)
  type(obs_diag),pointer,intent(in):: odiag

        cdiagbuf(ii)    = satellite_id       ! satellite id

        rdiagbuf(1,ii)  = 1                  ! observation type
        rdiagbuf(2,ii)  = 1                  ! observation subtype

        rdiagbuf(3,ii)  = data(ilate,i)      ! observation latitude (degrees)
        rdiagbuf(4,ii)  = data(ilone,i)      ! observation longitude (degrees)
        rdiagbuf(5,ii)  = presw              ! observation pressure (hPa)
        rdiagbuf(6,ii)  = data(ihgt,i)       ! observation height (meters)
        rdiagbuf(7,ii)  = (dtime*r60)-time_offset  ! obs time (sec relative to analysis time)
        rdiagbuf(8,ii)  = rmiss_single       ! input prepbufr qc or event mark
        rdiagbuf(9,ii) = rmiss_single       ! setup qc or event mark
        rdiagbuf(10,ii) = data(iuse,i)       ! read_prepbufr data usage flag
        if(muse(i)) then
           rdiagbuf(11,ii) = one             ! analysis usage flag (1=use,-1=not used)
        else
           rdiagbuf(11,ii) = -one
        endif

        rdiagbuf(12,ii) = errinv_input         ! prepbufr inverse obs error (ratio)**-1
        rdiagbuf(13,ii) = errinv_adjst         ! read_prepbufr inverse obs error (ratio)**-1
        rdiagbuf(14,ii) = errinv_final         ! final inverse observation error (ratio)**-1
        rdiagbuf(15,ii) = data(ifuvob,i)       ! fuv ratio observation 
        rdiagbuf(16,ii) = ddiff                ! obs-ges 
        rdiagbuf(17,ii) = data(ifuvob,i)-rfuv  ! obs-ges w/o bias correction (future slot)
        rdiagbuf(18,ii)=data(ioza,i)*rad2deg   ! observing zenith angle
        rdiagbuf(19,ii)=data(isza,i)*rad2deg  ! solar zenith angle

        rdiagbuf(20,ii) = 1.e+10_r_single    ! ges ensemble spread (filled in EnKF)
        rdiagbuf(21,ii) = 1.e+10_r_single    ! ges ensemble spread (filled in EnKF)

        if (lobsdiagsave) then
            write(6,*)'wrong here, stop in setupfuv.f90 '
            stop
           ioff=nreal
           do jj=1,miter
              ioff=ioff+1
              if (odiag%muse(jj)) then
                 rdiagbuf(ioff,ii) = one
              else
                 rdiagbuf(ioff,ii) = -one
              endif
           enddo
           do jj=1,miter+1
              ioff=ioff+1
              rdiagbuf(ioff,ii) = odiag%nldepart(jj)
           enddo
           do jj=1,miter
              ioff=ioff+1
              rdiagbuf(ioff,ii) = odiag%tldepart(jj)
           enddo
           do jj=1,miter
              ioff=ioff+1
              rdiagbuf(ioff,ii) = odiag%obssen(jj)
           enddo
        endif

  end subroutine contents_binary_diag_
  subroutine contents_netcdf_diag_(odiag)
  type(obs_diag),pointer,intent(in):: odiag
! Observation class
  character(7),parameter     :: obsclass = '    gold'
  real(r_kind),dimension(miter) :: obsdiag_iuse
           call nc_diag_metadata("Satellite_ID",            satellite_id        )
           call nc_diag_metadata("Observation_Type",        itype_gold(ikx)     )
           call nc_diag_metadata("Observation_Subtype",     isubtype_gold(ikx)  )
           call nc_diag_metadata("Latitude",                sngl(data(ilate,i)) )
           call nc_diag_metadata("Longitude",               sngl(data(ilone,i)) )
           call nc_diag_metadata("Pressure",                sngl(presw)         )
           call nc_diag_metadata("Height",                  sngl(data(ihgt,i))  )
           call nc_diag_metadata("Time",                    sngl(dtime-time_offset))
           call nc_diag_metadata("Prep_QC_Mark",            sngl(zero)          )
           call nc_diag_metadata("Prep_Use_Flag",           sngl(data(iuse,i))  )
           if(muse(i)) then
              call nc_diag_metadata("Analysis_Use_Flag",    sngl(one)           )
           else
              call nc_diag_metadata("Analysis_Use_Flag",    sngl(-one)          )
           endif

           call nc_diag_metadata("Errinv_Input",            sngl(errinv_input)  )
           call nc_diag_metadata("Errinv_Adjust",           sngl(errinv_adjst)  )
           call nc_diag_metadata("Errinv_Final",            sngl(errinv_final)  )

           call nc_diag_metadata("Observation",             sngl(data(ifuvob,i)) )
           call nc_diag_metadata("Obs_Minus_Forecast_adjusted",   sngl(ddiff)   )
           call nc_diag_metadata("Obs_Minus_Forecast_unadjusted", sngl(data(ifuvob,i)-rfuv) )

           if (lobsdiagsave) then
              do jj=1,miter
                 if (odiag%muse(jj)) then
                       obsdiag_iuse(jj) =  one
                 else
                       obsdiag_iuse(jj) = -one
                 endif
              enddo

              call nc_diag_data2d("ObsDiagSave_iuse",     obsdiag_iuse              )
              call nc_diag_data2d("ObsDiagSave_nldepart", odiag%nldepart )
              call nc_diag_data2d("ObsDiagSave_tldepart", odiag%tldepart )
              call nc_diag_data2d("ObsDiagSave_obssen"  , odiag%obssen   )
           endif

  end subroutine contents_netcdf_diag_

  subroutine final_vars_
    if(allocated(ges_z )) deallocate(ges_z )
    if(allocated(ges_ps)) deallocate(ges_ps)
    if(allocated(ges_tv)) deallocate(ges_tv)
    if(allocated(ges_o)) deallocate(ges_o)
    if(allocated(ges_o2)) deallocate(ges_o2)
    if(allocated(der_n2)) deallocate(der_n2)
  end subroutine final_vars_

end subroutine setupfuv
end module fuv_setup