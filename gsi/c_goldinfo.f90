module goldinfo
!$$$ module documentation block
!
! module:  goldinfo
!    prgmmr: cantrall                       date: 2020-07-24
! 
! abstract: This module contains variables and routines realted to the
!           assimilation of gold observations
!
! Subroutines included:
!   sub init_gold           - set gold related variables to defaults
!   sub goldinfo_read       - read in gold info
!   sub goldinfo_cleanup    - allocate and dellaocate public variables (copied from Chih-Ting Hsu EnKF code for sTEC)
!    
! Variable Definitions: (Some may not be neccessary)
!   def ctype_gold          - observation type - string
!   def ntype_gold          - number of input gold types
!   def itype_gold          - observation type - code
!   def isubtype_gold       - observation subtype - code
!   def typeuse_gold        - use flag
!   def grosserr_gold       - gross error parameter - gross error
!
!$$$ end documentation block

 use kinds, only: r_kind,i_kind
 implicit none

! set default to private
private
! set subroutines to public
public :: init_gold
public :: goldinfo_read
! set passed variables to public
public :: ctype_gold, ntype_gold, itype_gold, isubtype_gold, typeuse_gold, grosserr_gold
public :: mype_gold
public :: diag_gold

character(len=80)  :: fname = 'goldinfo'
logical diag_gold
integer(i_kind) :: mype_gold, ntype_gold
real(r_kind),allocatable,dimension(:) :: grosserr_gold
integer(i_kind),allocatable,dimension(:) :: typeuse_gold, itype_gold, isubtype_gold
character(len=20),allocatable,dimension(:) :: ctype_gold

contains

subroutine init_gold
!$$$ subprogram documentation block
!
! subprogram:   init_gold   initialize parameters for gold data
!   prgmmer:    cantrall                        date: 2020-07-27
!
! abstract: This routine sets default values for variables used in the gold
!           processing routines
!
!   input argument list:
!
!   output argument list:
!
!$$$
    use mpimod, only: npe       ! contains the number of mpi tasks, variable "npe" 
    implicit none

    ntype_gold = 0               ! number of entries read from goldinfo
    diag_gold = .true.          ! true = generate diag file
    mype_gold = 0   ! mpi task to write gold summary report
    
end subroutine init_gold

subroutine goldinfo_cleanup
    if (allocated(ctype_gold   )) deallocate(ctype_gold)
    if (allocated(grosserr_gold    )) deallocate(grosserr_gold)
    if (allocated(typeuse_gold )) deallocate(typeuse_gold)
    if (allocated(itype_gold   )) deallocate(itype_gold)
    if (allocated(isubtype_gold)) deallocate(isubtype_gold)
end subroutine goldinfo_cleanup

subroutine goldinfo_read
!$$$ subprogram documentation block
!
! subprogram:   goldinfo_read       read gold information file
!    !prgmmr:   cantrall            date: 2020-07-27
! 
! abstract: This routine reads the gold information file, global_goldinfo.txt
!
! input argument list:
!   mype = mpi task id
!
! output argument list:
!
!$$$
    use mpimod, only: mype 
    use obsmod, only: iout_gold 
    implicit none
    
    logical lexist
    character(len=1):: cflg
    character(len=120):: crecord
    integer(i_kind) lunin,j,k,istat,nlines
    data lunin / 47 /
 

!   Check the status of input file
    inquire(file=trim(fname),exist=lexist)
    !inquire(file=trim('global_goldinfo.txt'),exist=lexist) !NTC - cantrall (comment this line and uncomment above line)
    if ( lexist ) then

!       Determine number of entries in gold information file
        open(lunin,file=fname,form='formatted')
        !open(lunin,file='global_goldinfo.txt',form='formatted') !NTC-cantrall (comment this line and uncomment above line)
        j=0
        nlines=0
        readl: do
            read(lunin,100,iostat=istat,end=123) cflg,crecord
            if (istat /= 0) exit
            nlines = nlines+1
            if (cflg == '!') cycle
            j=j+1
        end do readl
123     continue
        if (istat>0) then
            write(6,*)'GOLDINFO_READ: ***ERROR*** reading goldinfo, istat=',istat
            close(lunin)
            write(6,*)'GOLDINFO_READ: stop program execution'
            call stop2(79)
        endif
        ntype_gold=j


!       Allocate arrays to hold gold information
        allocate(ctype_gold     (ntype_gold),&
                 itype_gold     (ntype_gold),&
                 isubtype_gold  (ntype_gold),&
                 typeuse_gold   (ntype_gold),&
                 grosserr_gold  (ntype_gold))

!       All mpi tasks open and read gold information file. 
!       Task mype_gold writes  information to gold runtime file
        if (mype==mype_gold) then
            open(iout_gold)
            !open(iout_gold,file='test_writeout.txt') !NTC - cantrall (remove line)
            write(iout_gold,110) ntype_gold
110         format('GOLDINFO_READ:  ntype_gold=',1x,i6)
        endif
        rewind(lunin)
        
!       Read info file
        j=0
        do k=1,nlines
            read(lunin,100) cflg, crecord
            if (cflg == '!') cycle
            j=j+1
            read(crecord,*) ctype_gold(j), itype_gold(j), isubtype_gold(j),&
                            typeuse_gold(j), grosserr_gold(j)
            if (mype==mype_gold) write(iout_gold,130) ctype_gold(j),&
                            itype_gold(j), isubtype_gold(j),&
                            typeuse_gold(j), grosserr_gold(j)
        enddo
        close(lunin)
        
        if (mype==mype_gold) close(iout_gold)
    
100     format(a1,a120)
130     format(a20,' itype = ',i2,' isubtype = ',i2, &
               ' typeuse = ',i2,' error = ',f7.3)
    end if
!   Successful read, return to calling routine
    return
    
end subroutine goldinfo_read

end module goldinfo



