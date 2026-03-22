module params
real(8),parameter:: timemax = 3.0d0 ! simulation end time (Mattia & Mignone 2022 Fig.12)

! option
integer, parameter :: flag_flux  = 2 ! 1 (HLL), 2 (HLLC), 3 (HLLD, not yet implemented)
integer, parameter :: flag_recon = 2 ! 1 (flat/1st-order), 2 (PLM/2nd-order van Leer)

! coordinate
integer,parameter::nx = 32   ! the number of grids in the simulation box
integer,parameter::ny = 2*nx ! the number of grids in the simulation box
integer,parameter::ngh = 2         ! the number of ghost cells
integer,parameter::nxtot = nx+2*ngh+1 ! the total number of grids including ghost cells
integer,parameter::nytot = ny+2*ngh+1 ! the total number of grids including ghost cells
integer,parameter::is = ngh+1         ! the index of the leftmost grid
integer,parameter::js = ngh+1         ! the index of the leftmost grid
integer,parameter::ie = nx+ngh     ! the index of the rightmost grid
integer,parameter::je = ny+ngh     ! the index of the rightmost grid
real(8),parameter::xmin = -0.5d0,xmax = 0.5d0
real(8),parameter::ymin = -1.0d0,ymax = 1.0d0

real(8),parameter::cfl_number = 0.4d0

! indices of the conservative variables
integer, parameter :: IDN = 1
integer, parameter :: IMX = 2
integer, parameter :: IMY = 3
integer, parameter :: IMZ = 4
integer, parameter :: IPR = 5
integer, parameter :: ISC = 6
integer, parameter :: IBX = 7
integer, parameter :: IBY = 8
integer, parameter :: IBZ = 9
integer, parameter :: NVAR = 6
integer, parameter :: NFLX = 9

! indices of the primitive variables
integer, parameter :: IVX = 2
integer, parameter :: IVY = 3
integer, parameter :: IVZ = 4
integer, parameter :: IEN = 5

real(8),parameter::gam=4.0d0/3.0d0 !! adiabatic index (Mattia & Mignone 2022, KHI test)

! output
character(20),parameter::dirname="ct" ! directory name
logical, parameter :: flag_binary = .true.

! snapshot
integer, parameter :: unitsnap = 17
real(8), parameter:: dtsnap=2.0d-1

! realtime analysis
integer, parameter :: nevo = 3
integer, parameter :: unitevo =11
integer, parameter :: unitbin =13
end module params
