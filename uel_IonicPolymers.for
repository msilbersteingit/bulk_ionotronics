************************************************************************
!     
!     User element subroutine for coupled electrochemical and mechanical
!     behavior of soft ionic polymers
!
!     Nikola Bosnjak, March 2022
!
!     This finite element code is used to solve coupled multiphysics problems
!     in Bosnjak, N., Tepermeister, M., & Silberstein, M.N (2022). Modeling  
!     coupled electrochemical and mechanical behavior of soft ionic materials
!     and ionotronic devices. Journal of the Mechanics and Physics of Solids,
!     168, 105014.
!
!     
!     The numerical approach follows the one laid out in Chester, S. A.,
!     Di Leo, C. V., & Anand, L. (2015). A finite element implementation
!     of a coupled diffusion-deformation theory for elastomeric gels.
!     International Journal of Solids and Structures, 52, 1-18.
!      
!     Nodal degrees of freedom are the chemical potential of two mobile species,
!     the electric potential and displacements.
! 
!     This subroutine is for the following element types:
!     > three-dimensional 8 node isoparametric element as shown below
!     with 1pt (reduced) or 8pt (full) gauss integration.
!
!     In order to avoid locking for the fully-integrated element, we
!     use the F-bar method of de Souza Neto (1996).
!
!     Mechanical, traction- and pressure-type boundary conditions 
!     may be applied to the dummy mesh using the Abaqus built-in 
!     commands *Dload or *Dsload.
!
!     Surface flux boundary conditions are supported for the following
!     elements.  Based on our convention, the face on which the fliud
!     flux is applied is the "label", i.e.
!     - U1,U2,U3,U4,... refer to fluid fluxes applied to faces 
!     1,2,3,4,... respectively,
!
!
!
!  8-node     8-----------7
!  brick     /|          /|       zeta
!           / |         / |       
!          5-----------6  |       |     eta
!          |  |        |  |       |   /
!          |  |        |  |       |  /
!          |  4--------|--3       | /
!          | /         | /        |/
!          |/          |/         O--------- xi
!          1-----------2        origin at cube center
!
!     Face numbering follows:
!       Face 1 = nodes 1,2,3,4
!       Face 2 = nodes 5,8,7,6
!       Face 3 = nodes 1,5,6,2
!       Face 4 = nodes 2,6,7,3
!       Face 5 = nodes 3,7,8,4
!       Face 6 = nodes 4,8,5,1
!
!
!     User element statement in the input file (set ? values as needed):
!
!     3D elements
!     *User Element,Nodes=8,Type=U3,Iproperties=2,Properties=14,Coordinates=3,
!     Variables=?,Unsymm
!     1,2,3,11,12
!
!     State Variables
!     --------------------------------------------------------------
!     Global SDV's (used for visualization)
!     1) rrho
!     2) ccRp_tau
!     3) ccRn_tau
!     4) ccRfix
!
!     Local SDV's (used for the solution procedure)
!       j = 0
!       do k = 1,nIntPt
!         svars(1+j) = Efield_tau(1,1)
!         svars(2+j) = Efield_tau(2,1(      
!         svars(3+j) = Efield_tau(3,1) 
!         svars(4+j) = ccRp_tau 
!         svars(5+j) = ccRn_tau       
!         j = j + nlSdv
!       end loop over k
!
!     In the input file, set 'User output variables'= number of global SDV's
!
!     In the input file, set 'ngSdv'= number of global SDV's
!
!     In the input file, set 'nlSdv'= number of local SDV's
!
!     In the input file, set 'varibles'=(nlSdv*nIntPt)
!
!
!     Material properties 
!     --------------------------------------------------------------
!     eNa     = props(1)  ! elementary charge per mol (C/mol)
!     RT      = props(2)  ! universal gas constant times temperature (J/mol) 
!     epsilon = props(3)  ! dielectric permittivity (F/m)
!     zP      = props(4)  ! charge number of positive ions
!     OmegaP  = props(5)  ! molar volume of positive ions (m3/mol)
!     ccRp0   = props(6)  ! initial concentration of positive ions (mol/m3)
!     ccRpBar = props(7)  ! maximum concentration of positive ions (mol/m3)
!     zN      = props(8)  ! charge number of negative ions
!     OmegaN  = props(9)  ! molar volume of negative ions (m3/mol)
!     ccRn0   = props(10) ! initial concentration of negative ions (mol/m3)
!     ccRnBar = props(11) ! maximum concentration of negative ions (mol/m3)
!     zFix    = props(12) ! charge number of fixed charge groups 
!     ccRfix  = props(13) ! concentration of fixed charge groups (mol/m3)
!     lambdaD = props(14) ! Debye length (m)
!     Gshear  = props(15) ! shear modulus (Pa)
!     Kbulk   = props(16) ! bulk modulus (Pa)
!     Length  = props(17) ! length of the body under consideration (m)
!     cE      = props(18) ! electrolyte strength (mol/m3)
!     
!     Element properties
!     --------------------------------------------------------------     
!     nlSdv   = jprops(1) ! Number of local sdv's per integ pt
!     ngSdv   = jprops(2) ! Number of global sdv's per integ pt
!     nInt    = jprops(3) ! Number of volume integ. points
!     nIntS   = jprops(4) ! Number of surface integ. points 
!     
!**********************************************************************

      module global

      ! This module is used to transfer SDV's from the UEL
      !  to the UVARM so that SDV's can be visualized on a
      !  dummy mesh
      !
      !  globalSdv(X,Y,Z)
      !   X - element pointer
      !   Y - integration point pointer
      !   Z - SDV pointer
      !
      !  numElem
      !   Total number of elements in the real mesh, the dummy
      !   mesh needs to have the same number of elements, and 
      !   the dummy mesh needs to have the same number of integ
      !   points.  You must set that parameter value here.
      !
      !  ElemOffset
      !   Offset between element numbers on the real mesh and
      !    dummy mesh.  That is set in the input file, and 
      !    that value must be set here the same.

      integer numElem,ElemOffset,err

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Set the number of UEL elements used here
      parameter(numElem=4100)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Set the offset here for UVARM plotting, must match input file!
      parameter(ElemOffset=10000)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real*8, allocatable :: globalSdv(:,:,:)

      integer test

      end module global

***********************************************************************

      SUBROUTINE SDVINI(STATEV,COORDS,NSTATV,NCRDS,NOEL,NPT,
     1     LAYER,KSPT)

      !DEC$ ATTRIBUTES ALIAS:"sdvini"::SDVINI
      
      use global

          
      INCLUDE 'ABA_PARAM.INC'

      DIMENSION STATEV(NSTATV),COORDS(NCRDS)


      statev = 0.999
      test = max(test,noel)


      RETURN
      END


      
***********************************************************************

      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     1 NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     2 JMAC,JMATYP,MATLAYO,LACCFLA)

      ! This subroutine is used to transfer SDV's from the UEL
      !  onto the dummy mesh for viewing.  Note that an offset of
      !  ElemOffset is used between the real mesh and the dummy mesh.
      !  If your model has more than ElemOffset UEL elements, then
      !  this will need to be modified.

      !DEC$ ATTRIBUTES ALIAS:"uvarm"::UVARM

      
      use global
     
      include 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3 FLGRAY(15)
      DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)

C     The dimensions of the variables FLGRAY, ARRAY and JARRAY
C     must be set equal to or greater than 15.

      uvar(1) = globalSdv(noel-ElemOffset,npt,1)
      uvar(2) = globalSdv(noel-ElemOffset,npt,2)
      uvar(3) = globalSdv(noel-ElemOffset,npt,3)
      uvar(4) = globalSdv(noel-ElemOffset,npt,4)
      
c      for example
c      uvar(2) = globalSdv(noel-ElemOffset,npt,2)
c      uvar(3) = globalSdv(noel-ElemOffset,npt,3)
c      uvar(4) = globalSdv(noel-ElemOffset,npt,4)

      return
      end subroutine uvarm

****************************************************************************

      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +     PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
     +     TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +     PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +     NJPROP,PERIOD)

            !DEC$ ATTRIBUTES ALIAS:"uel"::UEL


      use global
*
      IMPLICIT NONE
*
*     VARIABLES DEFINED IN UEL, PASSED BACK TO ABAQUS
*
      REAL(8) :: RHS,AMATRX,SVARS,ENERGY
*
*     VARIABLES PASSED INTO UEL 
*
      REAL(8) :: PROPS,coords,Uall,DUall,Vel,Accn,TIME,
     1  DTIME,PARAMS,ADLMAG,PREDEF,DDLMAG,PNEWDT,PERIOD
      INTEGER :: NDOFEL,NRHS,NSVARS,NPROPS,MCRD,NNODE,JTYPE,KSTEP,KINC,
     1  JELEM,NDLOAD,JDLTYP,NPREDF,LFLAGS,MLVARX,MDLOAD,JPROPS,NJPROP
*
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(*),ENERGY(8),coords(MCRD,NNODE),Uall(NDOFEL),
     2 DUall(MLVARX,*),Vel(NDOFEL),Accn(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

      integer lenJobName,lenOutDir,nDim,nInt,nIntS
      real*8 transient
      character*256 jobName,outDir,fileName


      

      !----------------------------------------------------------------
      ! 
      ! Perform initial checks
      !
      !
      ! Open the debug/error message file
      !
      !DEC$ ATTRIBUTES ALIAS:"getjobname"::GETJOBNAME
      !DEC$ ATTRIBUTES ALIAS:"getoutdir"::GETOUTDIR      
      !DEC$ ATTRIBUTES ALIAS:"xit"::XIT
      !
      call getJobName(jobName,lenJobName)
      call getOutDir(outDir,lenOutDir)
      fileName = outDir(1:lenOutDir)//'\aaMSGS_'//
     +     jobName(1:lenJobName)//'.dat'
      open(unit=80,file=fileName,status='unknown')


      ! Check the procedure type, this should be a coupled
      !  temperature displacement or pore pressure displacement
      !
      if((lflags(1).eq.64).or.(lflags(1).eq.65).or.
     +     (lflags(1).eq.62).or.(lflags(1).eq.63).or.
     +     (lflags(1).eq.71).or.
     +     (lflags(1).eq.72).or.(lflags(1).eq.73)) then
         !
         ! all is good
         !
      else
         write(*,*) 'Abaqus does not have the right procedure'
         write(*,*) 'go back and chekc the procedure type'
         write(*,*) 'lflags(1)=',lflags(1)
         write(80,*) 'Abaqus does not have the right procedure'
         write(80,*) 'go back and chekc the procedure type'
         write(80,*) 'lflags(1)=',lflags(1)
         call xit
      endif


      ! Make sure Abaqus knows you are doing a large
      !  deformation problem, I think this only matters
      !  when it comes to output in viewer
      !
      if(lflags(2).eq.0) then
         !
         ! lflags(2)=0 -> small disp.
         ! lflags(2)=1 -> large disp.
         !
         write(*,*) 'Abaqus thinks you are doing'
         write(*,*) 'a small displacement analysis'
         write(*,*) 'go in and set nlgeom=yes'
         write(80,*) 'Abaqus thinks you are doing'
         write(80,*) 'a small displacement analysis'
         write(80,*) 'go in and set nlgeom=yes'
         call xit
      endif


      ! Check to see if you are doing a general
      !  step or a linear purturbation step
      !
      if(lflags(4).eq.1) then
         !
         ! lflags(4)=0 -> general step
         ! lflags(4)=1 -> linear purturbation step
         !
         write(*,*) 'Abaqus thinks you are doing'
         write(*,*) 'a linear purturbation step'
         write(80,*) 'Abaqus thinks you are doing'
         write(80,*) 'a linear purturbation step'
         call xit         
      endif


      ! Do nothing if a ``dummy'' step
      !
      if(dtime.eq.0.0) return
      !
      ! Done with initial checks
      !
      !----------------------------------------------------------------


      !----------------------------------------------------------------
      ! 
      ! Obtain the number of integration points
      !
      nInt = jprops(3)
      nIntS = jprops(4)
      !
      ! Done with integ points
      !
      !----------------------------------------------------------------


      !----------------------------------------------------------------
      ! 
      ! Check for steady state or transient analysis
      !
      if((lflags(1).eq.62).or.(lflags(1).eq.63).or.
     +     (lflags(1).eq.71)) then
         !
         ! This is steady state
         !
         transient = 0.d0
      else
         !
         ! This is transient
         !
         transient = 1.d0
      endif
      !
      ! Done with check for steady state or transient
      !
      !----------------------------------------------------------------


      !----------------------------------------------------------------
      ! 
      ! Call the paricular element to perform the analysis
      !
      if(jtype.eq.1) then
         !
         ! This is a plane strain analysis
         !
         !
      elseif(jtype.eq.3) then
         !
         ! This is a 3D analysis
         !
         nDim = 3
         call U3D8(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +        PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
     +        TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +        PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +        NJPROP,PERIOD,
     +        nDim,nInt,nIntS,transient)
         !
         !
      else
         !
         ! We have a problem...
         !
         write(*,*) 'Element type not supported, jtype=',jtype
         write(80,*) 'Element type not supported, jtype=',jtype
         call xit
         !
      endif
      !
      ! Done with this element, RHS and AMATRX already returned
      !  as output from the specific element routine called
      !
      !----------------------------------------------------------------


      return
      end subroutine uel

************************************************************************

      subroutine U3D8(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +     PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
     +     TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +     PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +     NJPROP,PERIOD,
     +     nDim,nInt,nIntS,transient)

      !DEC$ ATTRIBUTES ALIAS:"xit"::XIT
      
      use global

      IMPLICIT NONE
*
*     VARIABLES DEFINED IN UEL, PASSED BACK TO ABAQUS
*
      REAL(8) :: RHS,AMATRX,SVARS,ENERGY
*
*     VARIABLES PASSED INTO UEL 
*
      REAL(8) :: PROPS,coords,Uall,DUall,Vel,Accn,TIME,
     1  DTIME,PARAMS,ADLMAG,PREDEF,DDLMAG,PNEWDT,PERIOD
      INTEGER :: NDOFEL,NRHS,NSVARS,NPROPS,MCRD,NNODE,JTYPE,KSTEP,KINC,
     1  JELEM,NDLOAD,JDLTYP,NPREDF,LFLAGS,MLVARX,MDLOAD,JPROPS,NJPROP
*
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(*),ENERGY(8),coords(MCRD,NNODE),Uall(NDOFEL),
     2 DUall(MLVARX,*),Vel(NDOFEL),Accn(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)


      real*8 u(nNode,3),du(nNode,ndofel),uNew(nNode,ndofel)
      real*8 uOld(nNode,ndofel),u_t(nNode,ndofel),v(nNode,3)
      real*8 coordsC(mcrd,nNode),PhiNew(nNode),PhiOld(nNode)
      real*8 dPhi(nNode),MuPNew(nNode),MuPOld(nNode),dMuP(nNode)
      real*8 MuNNew(nNode),MuNOld(nNode),dMuN(nNode)

      
      integer i,j,k,l,m,n,nIntPt,nDim,intpt,pOrder,a1,b1,a11,b11,face
      integer nInt,ii,jj,a12,b12,pe,nSdv,stat,q,nIntV,nIntPtV,p,ngSdv
      integer nlSdv,kk,lenJobName,lenOutDir,nIntS,faceFlag,nFibFam

      real*8 Iden(3,3),Le,theta0,phi0,Ru(3*nNode,1),Re(nNode,1),body(3)
      real*8 Kuu(3*nNode,3*nNode),Kee(nNode,nNode),sh0(nNode),detMapJ0
      real*8 dshxi(nNode,3),dsh0(nNode,3),dshC0(nNode,3),detMapJ0C
      real*8 Fc_tau(3,3),Fc_t(3,3),detFc_tau,detFc_t,w(nInt)
      real*8 sh(nNode),detMapJ,dsh(nNode,3),detMapJC,phiLmt,umeror
      real*8 dshC(nNode,3),F_tau(3,3)
      real*8 F_t(3,3),detF_tau,xi(nInt,3),detF,TR_tau(3,3),T_tau(3,3)
      real*8 SpTanMod(3,3,3,3)
      real*8 Smat(6,1),Bmat(6,3*nNode),BodyForceRes(3*nNode,1)
      real*8 Gmat(9,3*nNode),G0mat(9,3*nNode),Amat(9,9),Qmat(9,9),dA
      real*8 xLocal(nIntS),yLocal(nIntS),zLocal(nIntS),wS(nIntS),detF_t
      real*8 Kue(3*nNode,nNode),Keu(nNode,3*nNode),Nvec(1,nNode),ResFac
      real*8 AmatUC(6,1),TanFac
      real*8 Rp(nNode,1),Kpu(nNode,nDim*nNode)
      real*8 Kup(nDim*nNode,nNode),Kpp(nNode,nNode),Kpe(nNode,nNode)
      real*8 Kep(nNode,nNode),transient
      real*8 beta0,beta_t,beta_tau,Knn(nNode,nNode),Rn(nNode,1)
      real*8 Kun(nDim*nNode,nNode),Knu(nNode,nDim*nNode)
      real*8 Ken(nNode,nNode),Kne(nNode,nNode),Kpn(nNode,nNode)
      real*8 Knp(nNode,nNode),Efield_t(nDim,1),Efield0(nDim,1)
      real*8 dEdt_t(nDim,1),dEdt0(nDim,1)
      real*8 chargePC,Phi_tau,Phi_t,Cp_tau,Cp_t,Cn_tau,Cn_t
      real*8 ccRp0,ccRn0
      real*8 CpLmt,dPhidt,dPhidX(nDim,1),dCpdt,dCpdX(nDim,1)
      real*8 dCndt,dCndX(nDim,1),lambdaD,Length,flux(nDim,1)
      real*8 normal(nDim,1),Efield_tau(nDim,1),SpUPhiMod(3,3,3),volTerm
      real*8 SpPhiUMod(3,3,3),CpFlux(nDim,1),dEdt_1,dEdt_2,dEdt_3,rho
      real*8 RT,ccRp_t,ccRn_t,ccRp_tau,ccRn_tau,MuP_tau,MuN_tau,MuP_t
      real*8 MuN_t,dMuPdX(nDim,1),dMuNdX(nDim,1),dMuPdt,dMuNdt,dPhidMuP
      real*8 dPhidMuN,CRp_tau,CRn_tau,fluxRp(nDim,1),fluxRn(nDim,1)
      real*8 JvecP(nDim,1),JvecN(nDim,1),dCRpDotdPhi,dCRnDotdPhi
      real*8 dCRpDotdMu,dCRnDotdMu,dRhodPhi,dCRpdMu,dCRndMu
      real*8 dCRpdPhi,dCRndPhi,dMuPdPhi,dMuNdPhi,dCRpdt,dCRndt
      real*8 tempPP(8,3),tempPP2(8,3),tempNN(8,3),tempNN2(8,3)      
      real*8 dCpdMu,dCndMu,zP,zN,dCpdPhi,dCndPhi,ccRfix,rrho
      real*8 detFe,detFs,pressureP,pressureN,CRfix,zFix,eNa

      real*8 zero,one,two,half,Pi,three,third
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,Pi=3.141592653d0,
     +     three=3.d0,third=0.3333333)

      character*256 jobName,outDir,fileName

      ! Get element parameters
      !
      nlSdv  = jprops(1) 	! number of local sdv's per integ point
      ngSdv  = jprops(2) 	! number of global sdv's per integ point
     
      

      ! Allocate memory for the globalSdv's used for viewing
      !  results on the dummy mesh
      !
      if(.not.allocated(globalSdv)) then
         !
         ! allocate memory for the globalSdv's
         !
         ! numElem needs to be set in the MODULE
         ! nInt needs to be set in the UEL
         !
         stat=0
c         allocate(globalSdv(numElem,nInt,ngSdv))
c         deallocate(globalSdv)
         allocate(globalSdv(numElem,nInt,ngSdv),stat=err)
         if(stat.ne.0) then
            write(*,*) '//////////////////////////////////////////////'
            write(*,*) 'error when allocating globalSdv'
            write(*,*) '//////////////////////////////////////////////'
            write(*,*) '   stat=',stat
            write(*,*) '  ngSdv=',ngSdv
            write(*,*) '   nInt=',nInt
            write(*,*) 'numElem=',numElem
            write(*,*) '  nNode=',nNode
            write(*,*) 'lbound(globalSdv)',lbound(globalSdv)
            write(*,*) 'ubound(globalSdv)',ubound(globalSdv)
            write(*,*) '//////////////////////////////////////////////'
            write(80,*) '//////////////////////////////////////////////'
            write(80,*) 'error when allocating globalSdv'
            write(80,*) '//////////////////////////////////////////////'
            write(80,*) '   stat=',stat
            write(80,*) '  ngSdv=',ngSdv
            write(80,*) '   nInt=',nInt
            write(80,*) 'numElem=',numElem
            write(80,*) '  nNode=',nNode
            write(80,*) 'lbound(globalSdv)=',lbound(globalSdv)
            write(80,*) 'ubound(globalSdv)=',ubound(globalSdv)
            write(80,*) '//////////////////////////////////////////////'
            call xit
         endif
         write(*,*) '-------------------------------------------------'
         write(*,*) '----------- globalSDV ALLOCATED -----------------'
         write(*,*) '-------------------------------------------------'
         write(*,*) '---------- YOU PUT NUMBER OF ELEMENTS -----------'
         write(*,*) '---------- numElem=',numElem
         write(*,*) '---------- U3D8 ELEMENTS ------------------------'
         write(*,*) '-------------------------------------------------'
         write(*,*) '---------- YOU PUT NUMBER OF POINTS -------------'
         write(*,*) '---------- nInt =',nInt
         write(*,*) '---------- nIntS=',nIntS
         write(*,*) '-------------------------------------------------'
         write(*,*) '---------- YOU PUT NUMBER OF SDVs ---------------'
         write(*,*) '---------- ngSdv=',ngSdv
         write(*,*) '-------------------------------------------------'
      endif


      ! Identity tensor
      !
      call onem(Iden)


      ! Obtain initial conditions
      !
      Efield0 = zero
      dEdt0   = zero
      
      RT      = props(2)
      zP      = props(4)
      ccRp0   = props(6)
      zN      = props(8)
      zFix    = props(12)
      ccRn0   = props(10)
      lambdaD = props(16)
      Length  = props(17)
      

      ! Initialize the residual and tangent matrices to zero.
      !
      Ru  = zero
      Re  = zero
      Rp  = zero
      Rn  = zero
      Kuu = zero
      Kee = zero
      Kpp = zero
      Knn = zero
      Kue = zero
      Kup = zero
      Kun = zero
      Keu = zero
      Kep = zero
      Ken = zero
      Kpu = zero
      Kpe = zero
      Kpn = zero
      Knu = zero
      Kne = zero
      Knp = zero
      Energy = zero


      ! Body forces
      !
      body(1:3) = zero


      ! Obtain nodal displacements,electric potentials and chemical potentials
      !
      k = 0
      do i=1,nNode
         do j=1,nDim
            k = k + 1
            u(i,j) = Uall(k)
            du(i,j) = DUall(k,1)
            uOld(i,j) = u(i,j) - du(i,j)           
         enddo
         k = k + 1
         PhiNew(i) = Uall(k)
         dPhi(i) = DUall(k,1)
         PhiOld(i) = PhiNew(i) - dPhi(i)
         k = k + 1
         MuPNew(i) = Uall(k)
         dMuP(i) = DUall(k,1)
         MuPOld(i) = MuPNew(i) - dMuP(i)
         k=k+1
         MuNNew(i) = Uall(k)
         dMuN(i) = DUall(k,1)
         MuNOld(i) = MuNNew(i) - dMuN(i)
      enddo

      
      
      ! Obtain current nodal coordinates
      !
      do i=1,nNode
         do j=1,nDim
            coordsC(j,i) = coords(j,i) + u(i,j)
         enddo
      enddo

      
      !
      ! displacement increment, based on element diagonal
      !
      Le = dsqrt(((coordsC(1,1)-coordsC(1,7))**two) + 
     +     ((coordsC(2,1)-coordsC(2,7))**two) +
     +     ((coordsC(3,1)-coordsC(3,7))**two))
      !
      do i=1,nNode
         do j=1,nDim
            if(dabs(du(i,j)).gt.10.d0*Le) then
               pnewdt = 0.5
               return
            endif
         enddo
      enddo



      !----------------------------------------------------------------
      ! 
      ! Take this opportunity to perform calculations at the element
      !  centroid.  Here, check for hourglass stabilization and get
      !  the deformation gradient for use in the `F-bar' method.
      !
      ! Reference for the F-bar method:
      !  de Souza Neto, E.A., Peric, D., Dutko, M., Owen, D.R.J., 1996.
      !  Design of simple low order finite elements for large strain
      !  analysis of nearly incompressible solids. International Journal
      !  of Solids and Structures, 33, 3277-3296.
      !
      !
      ! Obtain shape functions and their local gradients at the element
      !  centriod, that means xi=eta=zeta=0.0, and nIntPt=1
      !
      if(nNode.eq.8) then
         call calcShape3DLinear(1,zero,1,sh0,dshxi)
      else
         write(*,*) 'Incorrect number of nodes: nNode.ne.8'
         write(80,*) 'Incorrect number of nodes: nNode.ne.8'
         call xit
      endif


      ! Map shape functions from local to global reference coordinate system
      !
      call mapShape3D(nNode,dshxi,coords,dsh0,detMapJ0,stat)
      if(stat.eq.0) then
         pnewdt = 0.5
         return
      endif



      ! Map shape functions from local to global current coordinate system
      !
      call mapShape3D(nNode,dshxi,coordsC,dshC0,detMapJ0C,stat)
      if(stat.eq.0) then
         pnewdt = 0.5
         return
      endif


      ! Calculate the deformation gradient at the element centriod
      !  at the the begining and end of the increment for use in 
      !  the `F-bar' method
      !
      Fc_tau = Iden
      Fc_t = Iden
      do i=1,nDim
         do j=1,nDim
            do k=1,nNode
               Fc_tau(i,j) = Fc_tau(i,j) + dsh0(k,j)*u(k,i)
               Fc_t(i,j) = Fc_t(i,j) + dsh0(k,j)*uOld(k,i)
            enddo
         enddo
      enddo
      call mdet(Fc_tau,detFc_tau)
      call mdet(Fc_t,detFc_t)
      !
      ! With the deformation gradient known at the element centriod
      !  we are now able to implement the `F-bar' method later
      !
      !----------------------------------------------------------------




      !----------------------------------------------------------------
      ! Begin the loop over integration points
      !
      ! Obtain integration point local coordinates and weights
      !
      if(nInt.eq.1) then
         call xint3D1pt(xi,w,nIntPt) ! 1-pt integration, nInt=1 above
      elseif(nInt.eq.8) then
         call xint3D8pt(xi,w,nIntPt) ! 8-pt integration, nInt=8 above
      else
         write(*,*) 'Invalid number of int points, nInt=',nInt
         write(80,*) 'Invalid number of int points, nInt=',nInt
         call xit
      endif


      ! Loop over integration points
      !
      jj = 0 ! jj is used for tracking the state variables
      do intpt=1,nIntPt


         ! Obtain state variables from previous increment
         !
         if((kinc.le.1).and.(kstep.eq.1)) then
            !
            ! this is the first increment, of the first step
            ! give initial conditions (or just anything)
            !
            Efield_t = Efield0
            ccRp_t   = ccRp0
            ccRn_t   = ccRn0
            !
         else
            !
            ! this is not the first increment, read old values
            !
            Efield_t(1,1) = svars(1+jj)
            Efield_t(2,1) = svars(2+jj)
            Efield_t(3,1) = svars(3+jj)
            ccRp_t        = svars(4+jj)
            ccRn_t        = svars(5+jj)
            !
         endif

         

         ! Obtain shape functions and their local gradients
         !
         if(nNode.eq.8) then
            call calcShape3DLinear(nIntPt,xi,intpt,sh,dshxi)
         else
            write(*,*) 'Incorrect number of nodes: nNode.ne.8'
            write(80,*) 'Incorrect number of nodes: nNode.ne.8'
            call xit
         endif


         ! Map shape functions from local to global reference coordinate system
         !
         call mapShape3D(nNode,dshxi,coords,dsh,detMapJ,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif


         ! Map shape functions from local to global current coordinate system
         !
         call mapShape3D(nNode,dshxi,coordsC,dshC,detMapJC,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif


         ! Obtain the electric potential and its derivatives at 
         ! this intPt at the begining and end of the increment
         !
         Phi_tau = zero
         Phi_t = zero
         dPhidt = zero
         dPhidX = zero
         do k=1,nNode
            Phi_tau = Phi_tau + PhiNew(k)*sh(k)
            Phi_t   = Phi_t   + PhiOld(k)*sh(k)
            do i=1,nDim
               dPhidX(i,1) = dPhidX(i,1) + PhiNew(k)*dshC(k,i)
            enddo
         enddo
         dPhidt = (Phi_tau - Phi_t)/dtime



         ! Obtain the chemical potential of positive ions and its derivatives
         ! at this intPt at the begining and end of the increment
         !
         MuP_tau = zero
         MuP_t = zero
         dMuPdt = zero
         dMuPdX = zero
         do k=1,nNode
            MuP_tau = MuP_tau + MuPNew(k)*sh(k)        
            MuP_t   = MuP_t   + MuPOld(k)*sh(k)
            do i=1,nDim
               dMuPdX(i,1) = dMuPdX(i,1) + MuPNew(k)*dshC(k,i)                
            enddo
         enddo
         dMuPdt = (MuP_tau - MuP_t)/dtime

         
         ! Obtain the chemical potential of negative ions and its derivatives
         ! at this intPt at the begining and end of the increment
         !
         MuN_tau = zero
         MuN_t   = zero
         dMuNdt  = zero
         dmuNdX  = zero
         do k=1,nNode
            MuN_tau = MuN_tau + MuNNew(k)*sh(k)
            MuN_t   = MuN_t   + MuNOld(k)*sh(k)
            do i=1,nDim
               dMuNdX(i,1) = dMuNdX(i,1) + MuNNew(k)*dshC(k,i)
            enddo
         enddo
         dMuNdt = (MuN_tau - MuN_t)/dtime
         
         
         
         ! Obtain, and modify the deformation gradient at this integration
         !  point.  Modify the deformation gradienet for use in the `F-bar'
         !  method.
         !
         F_tau = Iden
         F_t = Iden
         do i=1,nDim
            do j=1,nDim
               do k=1,nNode
                  F_tau(i,j) = F_tau(i,j) + dsh(k,j)*u(k,i)
                  F_t(i,j) = F_t(i,j) + dsh(k,j)*uOld(k,i)
               enddo
            enddo
         enddo
         !
         ! Modify the deformation gradient for the `F-bar' method
         !  only when using the 8 node fully integrated linear
         !  element, do not use the `F-bar' method for any other element
         !
         if((nNode.eq.8).and.(nInt.eq.8)) then
            call mdet(F_tau,detF_tau)
            call mdet(F_t,detF_t)
            F_tau = ((detFc_tau/detF_tau)**third)*F_tau
            F_t = ((detFc_tau/detF_tau)**third)*F_t
         endif
         call mdet(F_tau,detF)

         
         
         ! Obtain the electric field
         !
         do i=1,nDim
            Efield_tau(i,1) = -dPhidX(i,1)
         enddo

         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         !
         ! Perform the time integration at this integ. point to compute
         !  all the specific forms and parameters needed for the solution
!
         call  integ(
              ! inputs:
     +        props,nprops,dtime,F_tau,
     +        Phi_tau,Efield_tau,
     +        MuP_tau,ccRp_t,dMuPdX,
     +        MuN_tau,ccRn_t,dMuNdX,
              ! outputs:
     +        dPhidMuP,dPhidMuN,
     +        CRp_tau,Cp_tau,ccRp_tau,fluxRp,dCRpdPhi,dCRpdMu,
     +        dCRpDotdPhi,dCRpDotdMu,dmuPdPhi,dCRpdt,
     +        CRn_tau,Cn_tau,ccRn_tau,fluxRn,dCRndPhi,dCRndMu,
     +        dCRnDotdPhi,dCRnDotdMu,dMuNdPhi,dCRndt,
     +        rrho,Rho,dRhodPhi,CRfix,
     +        detFe,detFs,pressureP,pressureN,
     +        T_tau,SpTanMod)
         !
         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


         ! Save the state variables at this integ point
         ! at the end of the increment
         !
         svars(1+jj) = Efield_tau(1,1)
         svars(2+jj) = Efield_tau(2,1)
         svars(3+jj) = Efield_tau(3,1)         
         svars(4+jj) = ccRp_tau
         svars(5+jj) = ccRn_tau

         
         jj = jj + nlSdv ! setup for the next intPt


         ! Save the state variables at this integ point in the
         ! global array used for plotting field output
         !        
         globalSdv(jelem,intPt,1) = rrho
         globalSdv(jelem,intPt,2) = ccRp_tau/detF
         globalSdv(jelem,intPt,3) = ccRn_tau/detF
         globalSdv(jelem,intPt,4) = ccRfix/detF
       
        

         ! Time stepping algorithim based on the constitutive response
         !
         !CpLmt = 0.005d0
         !umeror = dabs((Cp_tau - Cp_t)/CpLmt)
         !if(umeror.le.half) then
         !   pnewdt = 1.5d0
         !elseif((umeror.gt.half).and.(umeror.le.0.8d0)) then
         !   pnewdt = 1.25d0
         !elseif((umeror.gt.0.8d0).and.(umeror.le.1.25d0)) then
         !   pnewdt = 0.75d0
         !else
         !   pnewdt = half
         !endif


         ! Compute/update the displacement residual vector
         !
         Smat(1,1) = T_tau(1,1)
         Smat(2,1) = T_tau(2,2)
         Smat(3,1) = T_tau(3,3)
         Smat(4,1) = T_tau(1,2)
         Smat(5,1) = T_tau(2,3)
         Smat(6,1) = T_tau(1,3)
         !
         Bmat = zero
         do kk=1,nNode
            Bmat(1,1+nDim*(kk-1)) = dshC(kk,1)
            Bmat(2,2+nDim*(kk-1)) = dshC(kk,2)
            Bmat(3,3+nDim*(kk-1)) = dshC(kk,3)
            Bmat(4,1+nDim*(kk-1)) = dshC(kk,2)
            Bmat(4,2+nDim*(kk-1)) = dshC(kk,1)
            Bmat(5,2+nDim*(kk-1)) = dshC(kk,3)
            Bmat(5,3+nDim*(kk-1)) = dshC(kk,2)
            Bmat(6,1+nDim*(kk-1)) = dshC(kk,3)
            Bmat(6,3+nDim*(kk-1)) = dshC(kk,1)
         enddo
         !
         BodyForceRes = zero
         do kk=1,nNode
            BodyForceRes(1+nDim*(kk-1),1) = sh(kk)*body(1)
            BodyForceRes(2+nDim*(kk-1),1) = sh(kk)*body(2)
            BodyForceRes(3+nDim*(kk-1),1) = sh(kk)*body(3)
         enddo
         !

                
         do kk=1,nNode
            Nvec(1,kk) = sh(kk)
         enddo
         
         Ru = Ru + detmapJC*w(intpt)*
     +        (
     +        -matmul(transpose(Bmat),Smat)
     +        + BodyForceRes
     +        )
         

         

 
         
         
         ! Compute/update the displacement tangent matrix
         !
         Gmat = zero
         do kk=1,nNode
            Gmat(1,1+nDim*(kk-1)) = dshC(kk,1)
            Gmat(2,2+nDim*(kk-1)) = dshC(kk,1)
            Gmat(3,3+nDim*(kk-1)) = dshC(kk,1)
            Gmat(4,1+nDim*(kk-1)) = dshC(kk,2)
            Gmat(5,2+nDim*(kk-1)) = dshC(kk,2)
            Gmat(6,3+nDim*(kk-1)) = dshC(kk,2)
            Gmat(7,1+nDim*(kk-1)) = dshC(kk,3)
            Gmat(8,2+nDim*(kk-1)) = dshC(kk,3)
            Gmat(9,3+nDim*(kk-1)) = dshC(kk,3)
         enddo

         G0mat = zero
         do kk=1,nNode
            G0mat(1,1+nDim*(kk-1)) = dshC0(kk,1)
            G0mat(2,2+nDim*(kk-1)) = dshC0(kk,1)
            G0mat(3,3+nDim*(kk-1)) = dshC0(kk,1)
            G0mat(4,1+nDim*(kk-1)) = dshC0(kk,2)
            G0mat(5,2+nDim*(kk-1)) = dshC0(kk,2)
            G0mat(6,3+nDim*(kk-1)) = dshC0(kk,2)
            G0mat(7,1+nDim*(kk-1)) = dshC0(kk,3)
            G0mat(8,2+nDim*(kk-1)) = dshC0(kk,3)
            G0mat(9,3+nDim*(kk-1)) = dshC0(kk,3)
         enddo

         Amat = zero
         Amat(1,1) = SpTanMod(1,1,1,1)
         Amat(1,2) = SpTanMod(1,1,2,1)
         Amat(1,3) = SpTanMod(1,1,3,1)
         Amat(1,4) = SpTanMod(1,1,1,2)
         Amat(1,5) = SpTanMod(1,1,2,2)
         Amat(1,6) = SpTanMod(1,1,3,2)
         Amat(1,7) = SpTanMod(1,1,1,3)
         Amat(1,8) = SpTanMod(1,1,2,3)
         Amat(1,9) = SpTanMod(1,1,3,3)
         Amat(2,1) = SpTanMod(2,1,1,1)
         Amat(2,2) = SpTanMod(2,1,2,1)
         Amat(2,3) = SpTanMod(2,1,3,1)
         Amat(2,4) = SpTanMod(2,1,1,2)
         Amat(2,5) = SpTanMod(2,1,2,2)
         Amat(2,6) = SpTanMod(2,1,3,2)
         Amat(2,7) = SpTanMod(2,1,1,3)
         Amat(2,8) = SpTanMod(2,1,2,3)
         Amat(2,9) = SpTanMod(2,1,3,3)
         Amat(3,1) = SpTanMod(3,1,1,1)
         Amat(3,2) = SpTanMod(3,1,2,1)
         Amat(3,3) = SpTanMod(3,1,3,1)
         Amat(3,4) = SpTanMod(3,1,1,2)
         Amat(3,5) = SpTanMod(3,1,2,2)
         Amat(3,6) = SpTanMod(3,1,3,2)
         Amat(3,7) = SpTanMod(3,1,1,3)
         Amat(3,8) = SpTanMod(3,1,2,3)
         Amat(3,9) = SpTanMod(3,1,3,3)
         Amat(4,1) = SpTanMod(1,2,1,1)
         Amat(4,2) = SpTanMod(1,2,2,1)
         Amat(4,3) = SpTanMod(1,2,3,1)
         Amat(4,4) = SpTanMod(1,2,1,2)
         Amat(4,5) = SpTanMod(1,2,2,2)
         Amat(4,6) = SpTanMod(1,2,3,2)
         Amat(4,7) = SpTanMod(1,2,1,3)
         Amat(4,8) = SpTanMod(1,2,2,3)
         Amat(4,9) = SpTanMod(1,2,3,3)
         Amat(5,1) = SpTanMod(2,2,1,1)
         Amat(5,2) = SpTanMod(2,2,2,1)
         Amat(5,3) = SpTanMod(2,2,3,1)
         Amat(5,4) = SpTanMod(2,2,1,2)
         Amat(5,5) = SpTanMod(2,2,2,2)
         Amat(5,6) = SpTanMod(2,2,3,2)
         Amat(5,7) = SpTanMod(2,2,1,3)
         Amat(5,8) = SpTanMod(2,2,2,3)
         Amat(5,9) = SpTanMod(2,2,3,3)
         Amat(6,1) = SpTanMod(3,2,1,1)
         Amat(6,2) = SpTanMod(3,2,2,1)
         Amat(6,3) = SpTanMod(3,2,3,1)
         Amat(6,4) = SpTanMod(3,2,1,2)
         Amat(6,5) = SpTanMod(3,2,2,2)
         Amat(6,6) = SpTanMod(3,2,3,2)
         Amat(6,7) = SpTanMod(3,2,1,3)
         Amat(6,8) = SpTanMod(3,2,2,3)
         Amat(6,9) = SpTanMod(3,2,3,3)
         Amat(7,1) = SpTanMod(1,3,1,1)
         Amat(7,2) = SpTanMod(1,3,2,1)
         Amat(7,3) = SpTanMod(1,3,3,1)
         Amat(7,4) = SpTanMod(1,3,1,2)
         Amat(7,5) = SpTanMod(1,3,2,2)
         Amat(7,6) = SpTanMod(1,3,3,2)
         Amat(7,7) = SpTanMod(1,3,1,3)
         Amat(7,8) = SpTanMod(1,3,2,3)
         Amat(7,9) = SpTanMod(1,3,3,3)
         Amat(8,1) = SpTanMod(2,3,1,1)
         Amat(8,2) = SpTanMod(2,3,2,1)
         Amat(8,3) = SpTanMod(2,3,3,1)
         Amat(8,4) = SpTanMod(2,3,1,2)
         Amat(8,5) = SpTanMod(2,3,2,2)
         Amat(8,6) = SpTanMod(2,3,3,2)
         Amat(8,7) = SpTanMod(2,3,1,3)
         Amat(8,8) = SpTanMod(2,3,2,3)
         Amat(8,9) = SpTanMod(2,3,3,3)
         Amat(9,1) = SpTanMod(3,3,1,1)
         Amat(9,2) = SpTanMod(3,3,2,1)
         Amat(9,3) = SpTanMod(3,3,3,1)
         Amat(9,4) = SpTanMod(3,3,1,2)
         Amat(9,5) = SpTanMod(3,3,2,2)
         Amat(9,6) = SpTanMod(3,3,3,2)
         Amat(9,7) = SpTanMod(3,3,1,3)
         Amat(9,8) = SpTanMod(3,3,2,3)
         Amat(9,9) = SpTanMod(3,3,3,3)


         Qmat = zero
         Qmat(1,1) = third*(Amat(1,1)+Amat(1,5)+Amat(1,9)) 
     +        - (two/three)*T_tau(1,1)
         Qmat(2,1) = third*(Amat(2,1)+Amat(2,5)+Amat(2,9))
     +        - (two/three)*T_tau(2,1)
         Qmat(3,1) = third*(Amat(3,1)+Amat(3,5)+Amat(3,9))
     +        - (two/three)*T_tau(3,1)
         Qmat(4,1) = third*(Amat(4,1)+Amat(4,5)+Amat(4,9))
     +        - (two/three)*T_tau(1,2)
         Qmat(5,1) = third*(Amat(5,1)+Amat(5,5)+Amat(5,9))
     +        - (two/three)*T_tau(2,2)
         Qmat(6,1) = third*(Amat(6,1)+Amat(6,5)+Amat(6,9))
     +        - (two/three)*T_tau(3,2)
         Qmat(7,1) = third*(Amat(7,1)+Amat(7,5)+Amat(7,9))
     +        - (two/three)*T_tau(1,3)
         Qmat(8,1) = third*(Amat(8,1)+Amat(8,5)+Amat(8,9))
     +        - (two/three)*T_tau(2,3)
         Qmat(9,1) = third*(Amat(9,1)+Amat(9,5)+Amat(9,9))
     +        - (two/three)*T_tau(3,3)
         Qmat(1,5) = Qmat(1,1)
         Qmat(2,5) = Qmat(2,1)
         Qmat(3,5) = Qmat(3,1)
         Qmat(4,5) = Qmat(4,1)
         Qmat(5,5) = Qmat(5,1)
         Qmat(6,5) = Qmat(6,1)
         Qmat(7,5) = Qmat(7,1)
         Qmat(8,5) = Qmat(8,1)
         Qmat(9,5) = Qmat(9,1)
         Qmat(1,9) = Qmat(1,1)
         Qmat(2,9) = Qmat(2,1)
         Qmat(3,9) = Qmat(3,1)
         Qmat(4,9) = Qmat(4,1)
         Qmat(5,9) = Qmat(5,1)
         Qmat(6,9) = Qmat(6,1)
         Qmat(7,9) = Qmat(7,1)
         Qmat(8,9) = Qmat(8,1)
         Qmat(9,9) = Qmat(9,1)
         

         if((nNode.eq.8).and.(nInt.eq.8)) then
            !
            ! This is the tangent using the F-bar method with the
            !  8 node fully integrated linear element
            !
            Kuu = Kuu + detMapJC*w(intpt)*
     +           (
     +           matmul(matmul(transpose(Gmat),Amat),Gmat)
     +           + matmul(transpose(Gmat),matmul(Qmat,(G0mat-Gmat)))
     +           )
         else
            !
            ! This is the tangent NOT using the F-bar method with all
            !  other elements
!
            Kuu = Kuu + detMapJC*w(intpt)*
     +           (
     +           matmul(matmul(transpose(Gmat),Amat),Gmat)
     +           )
         endif


         ! Compute the Debye length based on the eletrolyte strength
         ! per unit spat. vol. in (m)
         lambdaD = lambdaD*detF**half
         
         
         ! Compute/update the non-dim. electric potential residual vector
         !
         Re = Re + detmapJC*w(intpt)*
     +        (
     +        transpose(Nvec)*Rho -
     +        two*lambdaD**two/Length**two*
     +        matmul(dshC,dPhidX)              
     +        )
         
         
         
         ! Compute/update the positive ion non-dim. chem. pot. residual vector
         !
         Rp = Rp + detmapJC*w(intpt)*
     +        (
     +        transpose(Nvec)*dCRpdt
     +        + Cp_tau*lambdaD/Length*detF*
     +        (
     +        matmul(dshC,dmuPdX)
     +        + zP*matmul(dshC,dPhidX)
     +        )
     +        )

         
         ! Compute/update the negative ion non-dim. chem. pot. residual vector
         !
         Rn = Rn + detmapJC*w(intpt)*
     +        (
     +        transpose(Nvec)*dCRndt
     +        + Cn_tau*lambdaD/Length*detF*
     +        (
     +        matmul(dshC,dmuNdX)
     +        + zN*matmul(dshC,dPhidX)
     +        )
     +        )

         
         ! Compute/update non-dim. electric potential tangent matrix
         !
         Kee = Kee + detmapJC*w(intPt)*
     +        (
     +        -matmul(transpose(Nvec),Nvec)*(dRhodPhi)/detF
     +        +two*LambdaD**two/Length**two*matmul(dshC,transpose(dshC))
     +        )
        
         
         ! Compute/update the positive ion non-dim. chem. pot. tangent matrix
         !
         tempPP  = matmul(transpose(Nvec),transpose(dmuPdX))
         tempPP2 = matmul(transpose(Nvec),transpose(dPhidX))
         dCpdmu  = dCRpdmu/detF
         
         Kpp = Kpp + detmapJC*w(intPt)*
     +        (
     +        -matmul(transpose(Nvec),Nvec)*dCRpdotdmu
     +        -lambdaD/Length*
     +        (
     +        (
     +        dCpdmu*matmul(dshC,transpose(tempPP))
     +        + Cp_tau*matmul(dshC,transpose(dshC))
     +        )
     +        + zP*
     +        (
     +        dCpdmu*matmul(dshC,transpose(tempPP2))
     +        +Cp_tau*dPhidmuP*matmul(dshC,transpose(dshC))
     +        )
     +        )
     +        )

         
          
         ! Compute/update the negative ion non-dim. chem. pot. tangent matrix
         !
         tempNN  = matmul(transpose(Nvec),transpose(dmuNdX))
         tempNN2 = matmul(transpose(Nvec),transpose(dPhidX))
         dCndmu  = dCRndmu/detF
         
         Knn = Knn + detmapJC*w(intPt)*
     +        (
     +        -matmul(transpose(Nvec),Nvec)*dCRndotdmu
     +        -lambdaD/Length*
     +        (
     +        (
     +        dCndmu*matmul(dshC,transpose(tempNN))
     +        + Cn_tau*matmul(dshC,transpose(dshC))
     +        )
     +        + zN*
     +        (
     +        dCndmu*matmul(dshC,transpose(tempNN2))
     +        +Cn_tau*dPhidmuN*matmul(dshC,transpose(dshC))
     +        )
     +        )
     +        )

         
         ! Compute/update the non-dim. el. pot. - positive ion chem. pot.
         ! tangent matrix
         !              
         Kep = Kpn + detmapJC*w(intPt)*
     +        (
     +        - matmul(transpose(Nvec),Nvec)*zP*dCpdmu/detF
     +        + two*lambdaD**two/Length**two*
     +        dPhidmuP*matmul(dshC,transpose(dshC))
     +        )

         
         ! Compute/update the non-dim. el. pot. - negative ion chem. pot.
         ! tangent matrix
         !                  
         Ken = Ken + detmapJC*w(intPt)*
     +        (
     +        - matmul(transpose(Nvec),Nvec)*zN*dCndmu/detF
     +        + two*lambdaD**two/Length**two*
     +        dPhidmuN*matmul(dshC,transpose(dshC))
     +        )
         
          
         ! Compute/update the posititve ion non-dim. chem pot. - el. pot. 
         ! tangent matrix
         !    
         dCpdPhi = dCRpdPhi/detF
         
         Kpe = Kpe + detmapJC*w(intPt)*
     +        (
     +        - matmul(transpose(Nvec),Nvec)*dCRpdotdPhi
     +        + lambdaD/Length*
     +        ( 
     +        (
     +        dCpdPhi*matmul(dshC,transpose(tempPP))
     +        + Cp_tau*dmuPdPhi*matmul(dshC,transpose(dshC))
     +        )
     +        +zP*
     +        (
     +        dCpdPhi*matmul(dshC,transpose(tempPP2))
     +        + Cp_tau*matmul(dshC,transpose(dshC))
     +        )
     +        )
     +        )
     
       
          
         ! Compute/update the negative ion non-dim. chem pot. - el. pot. 
         ! tangent matrix
         ! 
         dCndPhi = dCRndPhi/detF
         
         Kne = Kne + detmapJC*w(intPt)*
     +        (
     +        - matmul(transpose(Nvec),Nvec)*dCRndotdPhi
     +        + lambdaD/Length*
     +        ( 
     +        (
     +        dCndPhi*matmul(dshC,transpose(tempNN))
     +        + Cn_tau*dmuNdPhi*matmul(dshC,transpose(dshC))
     +        )
     +        +zN*
     +        (
     +        dCndPhi*matmul(dshC,transpose(tempNN2))
     +        + Cn_tau*matmul(dshC,transpose(dshC))
     +        )
     +        )
     +        )
         

         
      enddo
      !
      ! End the loop over integration points
      !----------------------------------------------------------------




      !----------------------------------------------------------------
      ! Return Abaqus the RHS vector and the Stiffness matrix.
      !
      call AssembleElement(nDim,nNode,ndofel,
     +     Ru,Re,Rp,Rn,
     +     Kuu,Kue,Kup,Kun,
     +     Keu,Kee,Kep,Ken,
     +     Kpu,Kpe,Kpp,Kpn,
     +     Knu,Kne,Knp,Knn,
     +     rhs,amatrx)
      !
      ! End return of RHS and AMATRX
      !----------------------------------------------------------------

      
      return
      end subroutine U3D8

************************************************************************

      subroutine integ(
           ! inputs:
     +     props,nprops,dtime,F_tau,
     +     Phi_tau,Efield_tau,
     +     MuP_tau,ccRp_t,dMuPdX,
     +     MuN_tau,ccRn_t,dMuNdX,
           ! outputs:
     +     dPhidMuP,dPhidMuN,
     +     CRp_tau,Cp_tau,ccRp_tau,fluxP,dCRpdPhi,dCRpdMu,
     +     dCRpDotdPhi,dCRpDotdMu,dMuPdPhi,dCRpdt,
     +     CRn_tau,Cn_tau,ccRn_tau,fluxN,dCRndPhi,dCRndMu,
     +     dCRnDotdPhi,dCRnDotdMu,dMuNdPhi,dCRndt,
     +     rrho,Rho,dRhodPhi,CRfix,
     +     detFe,detFs,pressureP,pressureN,
     +     T_tau,SpTanMod)

      !DEC$ ATTRIBUTES ALIAS:"xit"::XIT

      ! This subroutine computes everything required for the time integration
      ! of the problem.
      !
      ! Inputs:
      !  1)  material parameters, props(nprops)
      !  2)  number of material properties, nprops
      !  3)  non-dim. time increment, dtime
      !  4)  deformation gradient, F_tau(3,3)
      !  5)  non-dim. electric potential, Phi_tau 
      !  6)  non-dim. electric field, Efield_tau
      !  7)  non-dim. chemical potential of positive ions, MuP_tau
      !  8)  old positive ion conc. in ref. configuration, ccRp_t (mol/m3)
      !  9)  non-dim. positive ion chemical potential gradient, dMuPdX
      !  10)  non-dim chemical potential of negative ions, MuN_tau
      !  11)  old negative ion conc. in ref. configuration, ccRn_t (mol/m3)
      !  12) negative ions non-dim. chemical potential gradient, dMuNdX


      
      ! Outputs:
          
      ! Electric potential:
      !  1)  der. of non-dim. el. pot. wrt chem pot. of pos. i., dPhidMuP 
      !  2)  der. of non-dim. el. pot. wrt chem pot. of neg. i., dPhidMuN 
      
      ! Positive ions:
      !  3)  non-dim. concentration per unit ref. vol., CRp_tau
      !  4)  non-dim. concentration per unit spat. vol., Cp_tau    
      !  4)  concentration per unit ref. config., ccRp_tau (mol/m3)
      !  5)  non-dim. flux of positive ions in spat. config., fluxP    
      !  6)  der. of non-dim. conc. wrt norm. electric pot., dCRpdPhi
      !  7)  der. of non-dim. conc. wrt chem pot., dCRpdMu
      !  8)  der. of non-dim. conc. t. der. wrt norm. el. pot., dCRpDotdPhi
      !  9)  der. of norm. conc. t. der. wrt echem pot., dCRpDotdmuP (mol/J)
      !  10) der. of echem pot. wrt norm. el. pot., dMuPdPhi (mol/J)
      !  11) non-dim. time derivative of ion conc. per unit ref. vol., dCRpdt
      
      ! Negative ions:
      !  3)  non-dim. concentration per unit ref. vol., CRn_tau
      !  4)  non-dim. concentration per unit spat. vol., Cn_tau    
      !  4)  concentration per unit ref. config., ccRn_tau (mol/m3)
      !  5)  non-dim. flux of positive ions in spat. config., fluxN    
      !  6)  der. of non-dim. conc. wrt norm. electric pot., dCRndPhi
      !  7)  der. of non-dim. conc. wrt chem pot., dCRndMu
      !  8)  der. of non-dim. conc. t. der. wrt norm. el. pot., dCRnDotdPhi
      !  9)  der. of norm. conc. t. der. wrt echem pot., dCRnDotdmuN (mol/J)
      !  10) der. of echem pot. wrt norm. el. pot., dMuNdPhi (mol/J)
      !  11) non-dim. time derivative of ion conc. per unit ref. vol., dCRndt

      ! Fixed and total charge concentration:
      !  3)  charge concentration per unit spat. vol., rrho (C/m3)
      !  4)  charge concentration per unit spat. vol., Rho 
      !  4)  der. of non-dim. charge conc. wrt el. pot, dRhodPhi
      !  5)  non-dim. conc. of fixed charge per unit ref. vol., CRfix

      
      ! Mechanical:
      !  12) volumetric deformation due to elasticity, detFe
      !  15) volumetric deformation due to swelling, detFs
      !  15) pressure term in positive ion chem. pot., pressureP
      !  15) pressure term in positive ion chem. pot., pressureN
      !  19)  Cauchy stress, T_tau(3,3)
      !  20)  spatial tangent modulus, SpTanMod(3,3,3,3)
      

      implicit none

      integer i,j,k,l,m,n,nprops,nargs,stat,nFibFam,kstep,jelem
      parameter(nargs=8)

      real*8 Iden(3,3),props(nprops),F_tau(3,3)
      real*8 theta,TR_tau(3,3),T_tau(3,3),dTRdF(3,3,3,3),Gshear,Kbulk
      real*8 spTanMod(3,3,3,3),chi,mu0,Vmol,Rgas,detF,FinvT(3,3),dPdt
      real*8 B_tau(3,3),trB_tau,C_tau(3,3),trC_tau,args(nargs),detFe,Nr
      real*8 detFs,dtime,Finv(3,3),G0,lamL,dGdF(3,3),lamBar,auxL
      real*8 T_tau_el(3,3),epsilon,Efield(3,1),SpUPhiMod(3,3,3)
      real*8 SpPhiUMod(3,3,3),dPhidX(3,1),Efield_tau(3,1)
      real*8 Phi_tau,MuP_tau,ccRp_t,dMuPdX(3,1),MuN_tau,ccRn_t
      real*8 dMuNdX(3,1),dPhidMuP,dPhidMuN
      real*8 CRp_tau,ccRp_tau,fluxP(3,1),dCRpdPhi,dCRpdMu,dCRpDotdPhi
      real*8 dCRpDotdMu,dMuPdPhi,CRn_tau,ccRn_tau,fluxN(3,1),dCRndPhi
      real*8 dCRndMu,dCRnDotdPhi,dCRnDotdMu,dMuNdPhi,RT,eNa,OmegaP
      real*8 OmegaN,zP,zN,ccRp0,ccRn0,ccRpBar,ccRnBar,CRfix,zFix
      real*8 pphi_tau,pressure,CRp_t,CRn_t,deltaMU
      real*8 deltaPHI,MuP_p,MuP_m,ccRp_p,ccRp_m,pphi_p,pphi_m
      real*8 dCRpdt,dCRndt,pressureP,pressureN,MuN_p,MuN_m,ccRn_p,ccRn_m
      real*8 Rho,ccRpdot_p,ccRpdot_m,ccRndot_p,ccRndot_m,dRhodPhi,cE
      real*8 factor,ccRfix,ccP_tau,ccN_tau,rrho,Cp_tau,Cn_tau
      
      
      
      real*8 zero,one,two,three,third,half,nine,twoThirds
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,third=0.3333333,
     +     half=1.d0/2.d0,nine=9.d0,twoThirds=0.6666666667)

      ! Identity tensor
      !
      call onem(Iden)

      
      ! Obtain material properties
    
      ! Physical constants
      eNa = props(1)
      RT  = props (2)
      
      ! Electrical
      epsilon = props(3)
      
      ! Positive ions
      zP      = props(4)
      OmegaP  = props(5)
      ccRp0   = props(6)
      ccRpBar = props(7)
      
      ! Negative ions
      zN      = props(8)
      omegaN  = props(9)
      ccRn0   = props(10)
      ccRnBar = props(11)
      
      ! Fixed charge groups
      zFix   = props(12) 
      ccRfix = props(13)
      
      ! Electrolyte strength
      cE     = props(18)
      
      ! Mechanical
      G0      = props(14)
      Kbulk   = props(15)


      
      ! Compute the inverse of F, its determinant, and its transpose
      !
      call matInv3D(F_tau,Finv,detF,stat)
      if(stat.eq.0) then
         write(*,*) 'Problem: detF.lt.zero'
         call xit
      endif
      FinvT = transpose(Finv)

      
      
      
      ! Compute the left Cauchy-Green tensor and its trace
      !
      B_tau = matmul(F_tau,transpose(F_tau))
      trB_tau = B_tau(1,1) + B_tau(2,2) + B_tau(3,3)


      ! Compute the right Cauchy-Green tensor and its trace
      !
      C_tau = matmul(transpose(F_tau),F_tau)
      trC_tau = C_tau(1,1) + C_tau(2,2) + C_tau(3,3)
	

      ! Compute the swelling volume ratio, detFs
      !
      detFs=one+omegaP*(ccRp_t - ccRp0)+omegaN*(ccRn_t - ccRn0)
 
     
      ! Compute the elastic volume ratio, detFe
      !
      detFe = detF/detFs
      
      
      ! Compute the concentration and flux of ions based on the
      ! electrochemical potential, electric potential and mechanical pressure.
      !
      ! Electrochemical potential is in J/mol, so we need to convert
      ! non-dimensional electric potential Phi to dimensional phi (V)
      ! to compute concentration from electrochemical potential
      !
      pphi_tau = Phi_tau*RT/eNa

      ! Ionic strength of the electrolyte per unit spat. volume (mol/m3)
      !
      cE = cE/detF
     

      ! Pressure terms in the chemical potential relations
      !
      pressureP = ( - OmegaP*Kbulk*dlog(detFe) + 
     +     half*Kbulk*OmegaP*dlog(detFe)**two)

      pressureN = (- OmegaN*Kbulk*dlog(detFe) + 
     +     half*Kbulk*OmegaN*dlog(detFe)**two) 
      

      
      ! Positive ion-related quantities

      ! Concentration of positive ions per unit spat. vol. (mol/m3)
      !
      ccP_tau = exp(MuP_tau - pressureP/RT)*ccRpBar

      ! Concentration of positive ions per unit ref. vol. (mol/m3)
      !
      ccRp_tau = ccP_tau*detF

      ! Non-dim. concentration of positive ions per unit spat. vol.
      !
      Cp_tau  = ccP_tau/cE
      
      ! Non-dim. concentration of positive ions per unit ref. vol.
      !
      CRp_tau  = Cp_tau*detF
      
      ! Non-dim. concentration of positive ions at the previous time increment
      !
      CRp_t  = ccRp_t/cE

      ! Non-dim. time derivative of positive ions
      !
      dCRpdt = (CRp_tau - CRp_t)/dtime
      
     
      ! Non-dim. flux of positive ions
      !
      fluxP = -Cp_tau*dMuPdX
      
      ! Compute the derivative of CRp wrt non-dim chem. pot.
      !
      dCRpdMu = ccRp_tau/cE
    
      ! Compute the derivative of CRp wrt non-dim. non-dim. el. pot. 
      !
      dCRpdPhi = zero
      
      
      ! Compute the derivative of dCRpdt wrt non-dim. chem. pot.
      ! obtained using implicit differentiation on the non-dim. chem. pot.

      ! Compute the perturbation on the chemical potential
      !
      if(dabs(MuP_tau).gt.one) then
         deltaMU = dabs(MuP_tau)*1.d-8
      else
         deltaMU = 1.d-8
      endif
      !
      MuP_p  = MuP_tau + deltaMU
      ccRp_p = exp(MuP_p  - pressureP/RT)*ccRpBar 
      ccRpDot_p = (ccRp_p-ccRp_t)/dtime 
      !
      MuP_m  = MuP_tau - deltaMU
      ccRp_m = exp(MuP_m - pressureP/RT)*ccRpBar 
      ccRpDot_m = (ccRp_m-ccRp_t)/dtime 
      !
      dCRpDotdMu = (ccRpDot_p-ccRpDot_m)/(two*deltaMU*cE)
      

      ! Compute the derivative of dCRpdt wrt non-dim. el. pot.
      !
      dCRpDotdPhi = zero
      
      ! Compute the derivative of non-dim. chem. pot. wrt non-dim. el. pot.
      !
      dMuPdPhi = zero
      
      ! Compute the derivative of el. pot. wrt non-dim. chem. pot.
      !
      dPhidMuP = zero
      
      
      
      
      ! Negative ion-related quantities
     
      ! Concentration of negative ions per unit spat. vol. (mol/m3)
      !
      ccN_tau = exp(MuN_tau - pressureN/RT)*ccRnBar

      ! Concentration of negative ions per unit ref. vol. (mol/m3)
      !
      ccRn_tau = ccN_tau*detF

      ! Non-dim. concentration of negative ions per unit spat. vol.
      !
      Cn_tau  = ccN_tau/cE
      
      ! Non-dim. concentration of negative ions per unit ref. vol.
      !
      CRn_tau  = Cn_tau*detF
      
      ! Non-dim. concentration of negative ions at the previous increment
      !
      CRn_t  = ccRn_t/cE
      
      ! Non-dim. time derivative of negative ions
      !
      dCRndt = (CRn_tau - CRn_t)/dtime
      
      ! Normalized flux of negative ions
      !
      fluxN = -Cn_tau*dMuNdX
      
      ! Compute the derivative of CRn wrt non-dim. chem. pot.
      !
      dCRndMu = ccRn_tau/cE
      
      ! Compute the derivative of CRp wrt non-dim. el. pot. 
      !
      dCRndPhi = zero
      
      ! Compute the derivative of dCRndt wrt non-dim. chem. pot.
      ! obtained using implicit differentiation on the non-dim. chem. pot.
      !
      MuN_p     = MuN_tau + deltaMU
      ccRn_p    = exp(MuN_p - pressureN/RT)*ccRnBar
      ccRnDot_p = (ccRn_p-ccRn_t)/dtime 
      !
      MuN_m     = MuN_tau - deltaMU
      ccRn_m    = exp(MuN_m - pressureN/RT)*ccRnBar 
      ccRnDot_m = (ccRn_m-ccRn_t)/dtime 
      !
      dCRnDotdMu = (ccRnDot_p-ccRnDot_m)/(two*deltaMU*cE)
      
      
      ! Compute the derivative of dCRndt wrt non-dim. el. pot.
      !
      dCRnDotdPhi = zero
      
      ! Compute the derivative of non-dim chem. pot. wrt non-dim. el. pot.
      !
      dMuNdPhi = zero
      
      ! Compute the derivative of non-dim el. pot. wrt non-dim. chem. pot.
      !
      dPhidMuN = zero
      
      
      
      ! Compute the charge density per unit spat. volume (mol/m3)
      !
      rrho =  one/detF*(zP*ccRp_tau+zN*ccRn_tau+zFix*ccRfix)

      ! Non-dim. conc. of fixed charges per unit ref. volume 
      !
      CRfix = ccRfix/cE

      ! Compute the non-dim. charge density per unit spat. volume
      !
      Rho = one/detF*(zP*CRp_tau+zN*CRn_tau+zFix*CRfix)

      ! Compute the derivative of non-dim charge density wrt non-dim el. pot.
      !
      dRhodPhi =   zero
	        
      if(lamL.le.zero) then
         Gshear = G0
      else
         lamBar = dsqrt(third*trC_tau)
         auxL = (three-(lamBar/lamL)**two)/(one - (lamBar/lamL)**two)
         Gshear = third*G0*auxL
      endif
      
 
      ! Compute the 1st Piola stress
      !
      TR_tau = 
     +     (
     +     Gshear*(F_tau - FinvT)
     +     + detFs*Kbulk*dlog(detFe)*FinvT
     +     )

      
      ! Compute the Cauchy stress
      !
      T_tau = (Gshear*(B_tau-Iden)
     +     +detFs*Kbulk*dlog(detFe)*Iden)/detF


      ! Compute the electrostatic contribution to the Cauchy
      ! stress (Maxwell stress)
      !
      T_tau_el = epsilon*RT*two/eNa**two*
     +     (
     +     matmul(Efield_tau,transpose(Efield_tau))
     +     - half*Iden*
     +     (
     +       Efield_tau(1,1)*Efield_tau(1,1)
     +     + Efield_tau(2,1)*Efield_tau(2,1)
     +     + Efield_tau(3,1)*Efield_tau(3,1)
     +     )
     +     )



      
      ! Compute dTRdF, the so-called material tangent modulus
      !
      dTRdF = zero
      do i=1,3
         do j = 1,3
            do k = 1,3
               do l = 1,3          
                  dTRdF(i,j,k,l) = dTRdF(i,j,k,l)
     +                 + (Gshear*
     +                 (
     +                 Iden(i,k)*Iden(j,l)
     +                 + Finv(l,i)*Finv(j,k)
     +                 )
     +                 + dGdF(k,l)*
     +                 (
     +                 (detF**(-twoThirds))*
     +                 (F_tau(i,j) - Finv(j,i))
     +                 )
     +                 + detFs*Kbulk*
     +                 (
     +                 Finv(j,i)*Finv(l,k)
     +                 - dlog(detFe)*Finv(l,i)*Finv(j,k)
     +                 )
     +                 )
               enddo
            enddo
         enddo
      enddo
      !
      ! Calculate the so-called spatial tangent modulus, based
      !  on the push forward of the material tangent modulus
      !
      SpTanMod = zero
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  do m=1,3
                     do n=1,3
                        SpTanMod(i,j,k,l) = SpTanMod(i,j,k,l) + 
     +                       (dTRdF(i,m,k,n)*F_tau(j,m)*F_tau(l,n))/detF
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      ! Calculate the spatial stress/non-dim. el. pot. modulus
      !
      SpUPhiMod = zero


      ! Calculate the non-dim. el. pot./displacement modulus
      !
      SpPhiUMod = zero


      return
      end subroutine integ

******************************************************************************
      
      subroutine AssembleElement(nDim,nNode,ndofel,
     +     Ru,Re,Rp,Rn,
     +     Kuu,Kue,Kup,Kun,
     +     Keu,Kee,Kep,Ken,
     +     Kpu,Kpe,Kpp,Kpn,
     +     Knu,Kne,Knp,Knn,
     +     rhs,amatrx)

      !DEC$ ATTRIBUTES ALIAS:"xit"::XIT
      
      !
      ! Subroutine to assemble the local elements residual and tangent
      !

      implicit none

      integer i,j,k,l,m,n,A11,A12,B11,B12,nDim,nNode,nDofEl,nDofN

      real*8 Ru(nDim*nNode,1),Re(nNode,1),Rp(nNode,1),Rn(nNode,1)
      real*8 Kuu(nDim*nNode,nDim*nNode),Kee(nNode,nNode)
      real*8 Kpp(nNode,nNode),Knn(nNode,nNode)
      real*8 Kue(nDim*nNode,nNode),Kup(nDim*nNode,nNode)
      real*8 Kun(nDim*nNode,nNode)
      real*8 Keu(nNode,nDim*nNode),Kep(nNode,nNode),Ken(nNode,nNode)
      real*8 Kpu(nNode,nDim*nNode),Kpe(nNode,nNode),Kpn(nNode,nNode)
      real*8 Knu(nNode,nDim*nNode),Kne(nNode,nNode),Knp(nNode,nNode)
      real*8 amatrx(ndofel,ndofel),rhs(ndofel,1)
      

      
      ! Total number of degrees of freedom per node
      !
      nDofN = nDofEl/nNode


      ! init
      !
      rhs(:,1) = 0.d0
      amatrx = 0.d0

      if(nDim.eq.2) then
         !
         ! Assemble the element level residual
         !
         do i=1,nNode
            A11 = nDofN*(i-1)+1
            A12 = nDim*(i-1)+1
            !
            ! displacement
            !
            rhs(A11,1)   = Ru(A12,1)
            rhs(A11+1,1) = Ru(A12+1,1)
            !
            ! electric potential
            !
            rhs(A11+2,1) = Re(i,1)
            !
            ! concentration of + ions
            !
            rhs(A11+3,1) = Rp(i,1)
            !
            ! concentration of - ions
            !
            rhs(A11+4,1) = Rn(i,1)
         enddo
         !
         ! Assemble the element level tangent matrix
         !
         do i=1,nNode
            do j=1,nNode
               A11 = nDofN*(i-1)+1
               A12 = nDim*(i-1)+1
               B11 = nDofN*(j-1)+1
               B12 = nDim*(j-1)+1
               !
               ! displacement
               !
               amatrx(A11,B11) = Kuu(A12,B12)
               amatrx(A11,B11+1) = Kuu(A12,B12+1)
               amatrx(A11+1,B11) = Kuu(A12+1,B12)
               amatrx(A11+1,B11+1) = Kuu(A12+1,B12+1)
               !
               ! electric potential
               !
               amatrx(A11+2,B11+2) = Kee(i,j)
               !
               ! concentration of + ions
               !
               amatrx(A11+3,B11+3) = Kpp(i,j)
               !
               ! concentration of - ions
               !
               amatrx(A11+4,B11+4) = Knn(i,j)
               !
               ! displacement - electric potential
               !
               amatrx(A11,B11+2) = Kue(A12,j)
               amatrx(A11+1,B11+2) = Kue(A12+1,j)
               !
               ! displacement - concentration of + ions
               !
               amatrx(A11,B11+3) = Kup(A12,j)
               amatrx(A11+1,B11+3) = Kup(A12+1,j)
               !
               ! displacement - concentration of - ions
               !
               amatrx(A11,B11+4)   = Kun(A12,j)
               amatrx(A11+1,B11+4) = Kun(A12+1,j)
               !
               ! electric potential - displacement
               !
               amatrx(A11+2,B11) = Keu(i,B12)
               amatrx(A11+2,B11+1) = Keu(i,B12+1)
               !
               ! electric potential - concentration of + ions
               !
               amatrx(A11+2,B11+3) = Kep(i,j)
               !
               ! electric potential - concentration of - ions
               !
               amatrx(A11+2,B11+4) = Ken(i,j)
               !
               ! concentration of + ions - displacement
               !
               amatrx(A11+3,B11) = Kpu(i,B12)
               amatrx(A11+3,B11+1) = Kpu(i,B12+1)
               !
               ! concentration of + ions - electric potential
               !
               amatrx(A11+3,B11+2) = Kpe(i,j)
               !
               ! concentration of + ions - concentration of - ions
               !
               amatrx(A11+3,B11+4) = Kpn(i,j)
               !
               ! concentration of - ions - displacement
               !
               amatrx(A11+4,B11) = Knu(i,B12)
               amatrx(A11+4,B11+1) = Knu(i,B12+1)
               !
               ! concentration of - ions - electric potential
               !
               amatrx(A11+4,B11+2) = Kne(i,j)
               !
               ! concentration of - ions - concentration of + ions
               !
               amatrx(A11+4,B11+3) = Knp(i,j)            
            enddo
         enddo
         !
      elseif(nDim.eq.3) then
         !
         ! Assemble the element level residual
         !
         do i=1,nNode
            A11 = nDofN*(i-1)+1
            A12 = nDim*(i-1)+1
            !
            ! displacement
            !
            rhs(A11,1)   = Ru(A12,1)
            rhs(A11+1,1) = Ru(A12+1,1)
            rhs(A11+2,1) = Ru(A12+2,1)
            !
            ! electric potential
            !
            rhs(A11+3,1) = Re(i,1)
            !
            ! concentration of + ions
            !
            rhs(A11+4,1) = Rp(i,1)
            !
            ! concentration of - ions
            !
            rhs(A11+5,1) = Rn(i,1)
            !
         enddo
         !
         ! Assembly the element level tangent matrix
         !
         do i=1,nNode
            do j=1,nNode
               A11 = nDofN*(i-1)+1
               A12 = nDim*(i-1)+1
               B11 = nDofN*(j-1)+1
               B12 = nDim*(j-1)+1
               !
               ! displacement
               !
               amatrx(A11,B11)     = Kuu(A12,B12)
               amatrx(A11,B11+1)   = Kuu(A12,B12+1)
               amatrx(A11,B11+2)   = Kuu(A12,B12+2)
               amatrx(A11+1,B11)   = Kuu(A12+1,B12)
               amatrx(A11+1,B11+1) = Kuu(A12+1,B12+1)
               amatrx(A11+1,B11+2) = Kuu(A12+1,B12+2)
               amatrx(A11+2,B11)   = Kuu(A12+2,B12)
               amatrx(A11+2,B11+1) = Kuu(A12+2,B12+1)
               amatrx(A11+2,B11+2) = Kuu(A12+2,B12+2)
               !
               ! electric potential
               !
               amatrx(A11+3,B11+3) = Kee(i,j)
               !
               ! concentration of + ions
               !
               amatrx(A11+4,B11+4) = Kpp(i,j)
               !
               ! concentration of - ions
               !
               amatrx(A11+5,B11+5) = Knn(i,j)
               !
               ! displacement - electric potential
               !
               amatrx(A11,B11+3)   = Kue(A12,j)
               amatrx(A11+1,B11+3) = Kue(A12+1,j)
               amatrx(A11+2,B11+3) = Kue(A12+2,j)
               !
               ! displacement - concentration of + ions
               !
               amatrx(A11,B11+4)   = Kup(A12,j)
               amatrx(A11+1,B11+4) = Kup(A12+1,j)
               amatrx(A11+2,B11+4) = Kup(A12+2,j)
               !
               ! displacement - concentration of - ions
               !
               amatrx(A11,B11+5)   = Kun(A12,j)
               amatrx(A11+1,B11+5) = Kun(A12+1,j)
               amatrx(A11+2,B11+5) = Kun(A12+2,j)
               !
               ! electric potential - displacement
               !
               amatrx(A11+3,B11)   = Keu(i,B12)
               amatrx(A11+3,B11+1) = Keu(i,B12+1)
               amatrx(A11+3,B11+2) = Keu(i,B12+2)
               !
               ! electric potential - concentration of + ions
               !
               amatrx(A11+3,B11+4) = Kep(i,j)
               !
               ! electric potential - concentration of - ions
               !
               amatrx(A11+3,B11+5) = Ken(i,j)
               !
               ! concentration of + ions - displacement
               !
               amatrx(A11+4,B11)   = Kpu(i,B12)
               amatrx(A11+4,B11+1) = Kpu(i,B12+1)
               amatrx(A11+4,B11+2) = Kpu(i,B12+2)
               !
               ! concentration of + ions - electric potential
               !
               amatrx(A11+4,B11+3) = Kpe(i,j)
               !
               ! concentration of + ions - concentration of - ions
               !
               amatrx(A11+4,B11+5) = Kpn(i,j)
               !
               ! concentration of - ions - displacement
               !
               amatrx(A11+5,B11)   = Knu(i,B12)
               amatrx(A11+5,B11+1) = Knu(i,B12+1)
               amatrx(A11+5,B11+2) = Knu(i,B12+2)
               !
               ! concentration of - ions - electric potential
               !
               amatrx(A11+5,B11+3) = Kne(i,j)
               !
               ! concentration of - ions - concentration of + ions
               !
               amatrx(A11+5,B11+4) = Knp(i,j)
               !
            enddo
         enddo
         !
      else
         write(*,*) 'How did you get nDim=',nDim
         call xit
      endif

      return
      end subroutine AssembleElement

************************************************************************

      subroutine computeSUPG(nDim,nNode,coordsC,direction,epsilon,SUPG)

      !DEC$ ATTRIBUTES ALIAS:"xit"::XIT
      
      ! This subroutine computes the SUPG stabilization parameter often
      !  called ``tau'' in the literature.  Since we have coth(a)-1/a as
      !  part of this calculation, we use the taylor expansion near zero
      !  to avoid the divide by zero.

      implicit none
      
      integer nDim,nNode

      real*8 hXi,hEta,hZeta,eXi(nDim,1),eEta(nDim,1),eZeta(nDim,1),dXi
      real*8 dEta,dZeta,alphaXi,alphaEta,alphaZeta,xiBar,etaBar,zetaBar
      real*8 coordsC(nDim,nNode),direction(nDim,1),epsilon,tmpX,tmpY
      real*8 tmpZ,SUPG

      if(nDim.eq.2) then
         !
         ! This is a 2D problem
         !
         hXi = 0.5d0*dsqrt(
     +       (coordsC(1,3)+coordsC(1,2)-coordsC(1,1)-coordsC(1,4))**2.d0
     +      +(coordsC(2,3)+coordsC(2,2)-coordsC(2,1)-coordsC(2,4))**2.d0
     +        )
         hEta = 0.5d0*dsqrt(
     +       (coordsC(1,3)+coordsC(1,4)-coordsC(1,1)-coordsC(1,2))**2.d0
     +      +(coordsC(2,3)+coordsC(2,4)-coordsC(2,1)-coordsC(2,2))**2.d0
     +        )

      
         eXi(1,1) = 
     +        (coordsC(1,3)+coordsC(1,2)-coordsC(1,1)-coordsC(1,4))
     +        /(2.d0*hXi)
         eXi(2,1) = 
     +        (coordsC(2,3)+coordsC(2,2)-coordsC(2,1)-coordsC(2,4))
     +        /(2.d0*hXi)
         eEta(1,1) = 
     +        (coordsC(1,3)+coordsC(1,4)-coordsC(1,1)-coordsC(1,2))
     +        /(2.d0*hEta)
         eEta(2,1) = 
     +        (coordsC(2,3)+coordsC(2,4)-coordsC(2,1)-coordsC(2,2))
     +        /(2.d0*hEta)
         

         dXi = direction(1,1)*eXi(1,1) + direction(2,1)*eXi(2,1)
         dEta = direction(1,1)*eEta(1,1) + direction(2,1)*eEta(2,1)
      

         if(epsilon.lt.1.d-10) then
            !
            ! we have zero diffusion in the PDE, but we assume
            !  a very small value only to compute the SUPG parameter
            !
            epsilon = 1.d-9
            !
         endif

         alphaXi = (dXi*hXi)/(2.d0*epsilon)
         alphaEta = (dEta*hEta)/(2.d0*epsilon)
         !
         if(abs(alphaXi).le.0.1d0) then
            xiBar = (1.d0/3.d0)*alphaXi 
     +           - (1.d0/45.d0)*alphaXi**3.d0 
     +           + (2.d0/945.d0)*alphaXi**5.d0
     +           - (1.d0/4725.d0)*alphaXi**7.d0
         else
            xiBar = 1.d0/dtanh(alphaXi) - 1.d0/alphaXi   
         endif
         !
         if(abs(alphaEta).le.0.1d0) then
            etaBar = (1.d0/3.d0)*alphaEta 
     +           - (1.d0/45.d0)*alphaEta**3.d0 
     +           + (2.d0/945.d0)*alphaEta**5.d0
     +           - (1.d0/4725.d0)*alphaEta**7.d0
         else
            etaBar = 1.d0/dtanh(alphaEta) - 1.d0/alphaEta
         endif
         !
         SUPG = (xiBar*dXi*hXi + etaBar*dEta*hEta)/2.d0
         !
      elseif(nDim.eq.3) then
         !
         ! This is a 3D problem
         !
         tmpX = coordsC(1,2)+coordsC(1,3)+coordsC(1,6)+coordsC(1,7)
     +        - coordsC(1,1)-coordsC(1,4)-coordsC(1,5)-coordsC(1,8)
         tmpY = coordsC(2,2)+coordsC(2,3)+coordsC(2,6)+coordsC(2,7)
     +        - coordsC(2,1)-coordsC(2,4)-coordsC(2,5)-coordsC(2,8)
         tmpZ = coordsC(3,2)+coordsC(3,3)+coordsC(3,6)+coordsC(3,7)
     +        - coordsC(3,1)-coordsC(3,4)-coordsC(3,5)-coordsC(3,8)
         hXi = 0.25d0*dsqrt(tmpX**2.d0 + tmpY**2.d0 + tmpZ**2.d0)
         eXi(1,1) = 0.25d0*tmpX/hXi
         eXi(2,1) = 0.25d0*tmpY/hXi
         eXi(3,1) = 0.25d0*tmpZ/hXi


         tmpX = coordsC(1,3)+coordsC(1,4)+coordsC(1,7)+coordsC(1,8)
     +        - coordsC(1,1)-coordsC(1,2)-coordsC(1,5)-coordsC(1,6)
         tmpY = coordsC(2,3)+coordsC(2,4)+coordsC(2,7)+coordsC(2,8)
     +        - coordsC(2,1)-coordsC(2,2)-coordsC(2,5)-coordsC(2,6)
         tmpZ = coordsC(3,3)+coordsC(3,4)+coordsC(3,7)+coordsC(3,8)
     +        - coordsC(3,1)-coordsC(3,2)-coordsC(3,5)-coordsC(3,6)
         hEta = 0.25d0*dsqrt(tmpX**2.d0 + tmpY**2.d0 + tmpZ**2.d0)
         eEta(1,1) = 0.25d0*tmpX/hEta
         eEta(2,1) = 0.25d0*tmpY/hEta
         eEta(3,1) = 0.25d0*tmpZ/hEta

         tmpX = coordsC(1,5)+coordsC(1,6)+coordsC(1,7)+coordsC(1,8)
     +        - coordsC(1,1)-coordsC(1,2)-coordsC(1,3)-coordsC(1,4)
         tmpY = coordsC(2,5)+coordsC(2,6)+coordsC(2,7)+coordsC(2,8)
     +        - coordsC(2,1)-coordsC(2,2)-coordsC(2,3)-coordsC(2,4)
         tmpZ = coordsC(3,5)+coordsC(3,6)+coordsC(3,7)+coordsC(3,8)
     +        - coordsC(3,1)-coordsC(3,2)-coordsC(3,3)-coordsC(3,4)
         hZeta = 0.25d0*dsqrt(tmpX**2.d0 + tmpY**2.d0 + tmpZ**2.d0)
         eZeta(1,1) = 0.25d0*tmpX/hZeta
         eZeta(2,1) = 0.25d0*tmpY/hZeta
         eZeta(3,1) = 0.25d0*tmpZ/hZeta


         dXi = direction(1,1)*eXi(1,1) + direction(2,1)*eXi(2,1) 
     +        + direction(3,1)*eXi(3,1) 
         dEta = direction(1,1)*eEta(1,1) + direction(2,1)*eEta(2,1)
     +        + direction(3,1)*eEta(3,1)
         dZeta = direction(1,1)*eZeta(1,1) + direction(2,1)*eZeta(2,1)
     +        + direction(3,1)*eZeta(3,1)



         if(epsilon.lt.1.d-10) then
            !
            ! we have zero diffusion in the PDE, but we assume
            !  a very small value only to compute the SUPG parameter
            !
            epsilon = 1.d-9
            !
         endif

         alphaXi = (dXi*hXi)/(2.d0*epsilon)
         alphaEta = (dEta*hEta)/(2.d0*epsilon)
         alphaZeta = (dZeta*hZeta)/(2.d0*epsilon)
            
         if(abs(alphaXi).le.0.1d0) then
            xiBar = (1.d0/3.d0)*alphaXi 
     +           - (1.d0/45.d0)*alphaXi**3.d0 
     +           + (2.d0/945.d0)*alphaXi**5.d0
     +           - (1.d0/4725.d0)*alphaXi**7.d0
         else
            xiBar = 1.d0/dtanh(alphaXi) - 1.d0/alphaXi   
         endif
         
         if(abs(alphaEta).le.0.1d0) then
            etaBar = (1.d0/3.d0)*alphaEta 
     +           - (1.d0/45.d0)*alphaEta**3.d0 
     +           + (2.d0/945.d0)*alphaEta**5.d0
     +           - (1.d0/4725.d0)*alphaEta**7.d0
         else
            etaBar = 1.d0/dtanh(alphaEta) - 1.d0/alphaEta
         endif
         
         if(abs(alphaZeta).le.0.1d0) then
            zetaBar = (1.d0/3.d0)*alphaZeta 
     +           - (1.d0/45.d0)*alphaZeta**3.d0 
     +           + (2.d0/945.d0)*alphaZeta**5.d0
     +           - (1.d0/4725.d0)*alphaZeta**7.d0
         else
            zetaBar = 1.d0/dtanh(alphaZeta) - 1.d0/alphaZeta
         endif
         !
         SUPG = (xiBar*dXi*hXi + etaBar*dEta*hEta + zetaBar*dZeta*hZeta)
     +        /2.d0
         !
      else
         write(*,*) 'SUPG, bad nDim=',nDim
         write(80,*) 'SUPG: bad nDim=',nDim
         call xit
      endif


      return
      end subroutine computeSUPG
        
!****************************************************************************
!     Element subroutines
!****************************************************************************

      subroutine xint2D1pt(xi,w,nIntPt)

       
      !
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 1 gauss point for integration
      !
      !  xi(nIntPt,2): xi,eta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      !
      implicit none
      !
      integer nIntPt,nDim
      !
      real*8 xi(1,2), w(1)


      ! Initialize
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 1


      ! Gauss weights
      !
      w = 4.d0
      

      ! Gauss pt location in master element
      !
      xi(1,1) = 0.d0
      xi(1,2) = 0.d0


      return
      end subroutine xint2D1pt
      
!************************************************************************

      subroutine xint2D4pt(xi,w,nIntPt)

      
      !
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 4 gauss points for integration
      !
      !  xi(nIntPt,2): xi,eta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      !
      implicit none
      !
      integer nIntPt,nDim
      !
      real*8 xi(4,2), w(4)


      ! Initialize
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 4


      ! Gauss weights
      !
      w(1) = 1.d0
      w(2) = 1.d0
      w(3) = 1.d0
      w(4) = 1.d0
      

      ! Gauss pt locations in master element
      !
      xi(1,1) = -dsqrt(1.d0/3.d0)
      xi(1,2) = -dsqrt(1.d0/3.d0)
      xi(2,1) = dsqrt(1.d0/3.d0)
      xi(2,2) = -dsqrt(1.d0/3.d0)
      xi(3,1) = -dsqrt(1.d0/3.d0)
      xi(3,2) = dsqrt(1.d0/3.d0)
      xi(4,1) = dsqrt(1.d0/3.d0)
      xi(4,2) = dsqrt(1.d0/3.d0)


      return
      end subroutine xint2D4pt

************************************************************************

      subroutine xint3D1pt(xi,w,nIntPt)
  
      
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 3D elements
      !  using a 2 gauss points for integration
      !
      !  xi(nIntPt,3): xi,eta,zeta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      
      implicit none

      integer nIntPt,nDim

      real*8 xi(1,3),w(1)


      ! Init
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 1


      ! Gauss weights
      !
      w(1) = 8.d0
      

      ! Gauss pt locations in master element
      !
      xi(1,1) = 0.d0
      xi(1,2) = 0.d0
      xi(1,3) = 0.d0

      return
      end subroutine xint3D1pt
     
************************************************************************

      subroutine xint3D8pt(xi,w,nIntPt)
     
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 3D elements
      !  using 8 gauss points for integration
      !
      !  xi(nIntPt,3): xi,eta,zeta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      
      implicit none

      integer nIntPt,nDim

      real*8 xi(8,3),w(8)


      ! Init
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 8


      ! Gauss weights
      !
      w(1) = 1.d0
      w(2) = 1.d0
      w(3) = 1.d0
      w(4) = 1.d0
      w(5) = 1.d0
      w(6) = 1.d0
      w(7) = 1.d0
      w(8) = 1.d0
      

      ! Gauss pt locations in master element
      !
      xi(1,1) = -dsqrt(1.d0/3.d0)
      xi(1,2) = -dsqrt(1.d0/3.d0)
      xi(1,3) = -dsqrt(1.d0/3.d0)
      xi(2,1) = dsqrt(1.d0/3.d0)
      xi(2,2) = -dsqrt(1.d0/3.d0)
      xi(2,3) = -dsqrt(1.d0/3.d0)
      xi(3,1) = -dsqrt(1.d0/3.d0)
      xi(3,2) = dsqrt(1.d0/3.d0)
      xi(3,3) = -dsqrt(1.d0/3.d0)
      xi(4,1) = dsqrt(1.d0/3.d0)
      xi(4,2) = dsqrt(1.d0/3.d0)
      xi(4,3) = -dsqrt(1.d0/3.d0)
      xi(5,1) = -dsqrt(1.d0/3.d0)
      xi(5,2) = -dsqrt(1.d0/3.d0)
      xi(5,3) = dsqrt(1.d0/3.d0)
      xi(6,1) = dsqrt(1.d0/3.d0)
      xi(6,2) = -dsqrt(1.d0/3.d0)
      xi(6,3) = dsqrt(1.d0/3.d0)
      xi(7,1) = -dsqrt(1.d0/3.d0)
      xi(7,2) = dsqrt(1.d0/3.d0)
      xi(7,3) = dsqrt(1.d0/3.d0)
      xi(8,1) = dsqrt(1.d0/3.d0)
      xi(8,2) = dsqrt(1.d0/3.d0)
      xi(8,3) = dsqrt(1.d0/3.d0)


      return
      end subroutine xint3D8pt

************************************************************************

      subroutine xintSurf2D1pt(face,xLocal,yLocal,w)

      !DEC$ ATTRIBUTES ALIAS:"xit"::XIT
      
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 2 gauss points for surface integration
      !
      !  xLocal(nIntPt): x coordinates for the integration pts
      !  yLocal(nIntPt): y coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights

      implicit none

      integer face

      real*8 xLocal(1),yLocal(1),w(1),zero,one,two
      parameter(zero=0.d0,one=1.d0,two=2.d0)


      ! Gauss weights
      !
      w(1) = two
      

      ! Gauss pt locations in master element
      !
      if(face.eq.1) then
         xLocal(1) = zero
         yLocal(1) = -one
      elseif(face.eq.2) then
         xLocal(1) = one
         yLocal(1) = zero
      elseif(face.eq.3) then
         xLocal(1) = zero
         yLocal(1) = one
      elseif(face.eq.4) then
         xLocal(1) = -one
         yLocal(1) = zero
      else
         write(*,*) 'face.ne.1,2,3,4'
         write(80,*) 'face.ne.1,2,3,4'
         call xit
      endif

      end subroutine xintSurf2D1pt

************************************************************************

      subroutine xintSurf2D2pt(face,xLocal,yLocal,w)

      !DEC$ ATTRIBUTES ALIAS:"xit"::XIT
      
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 2 gauss points for surface integration
      !
      !  xLocal(nIntPt): x coordinates for the integration pts
      !  yLocal(nIntPt): y coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights

      implicit none

      integer face

      real*8 xLocal(2),yLocal(2),w(2),one,three
      parameter(one=1.d0,three=3.d0)


      ! Gauss weights
      !
      w(1) = one
      w(2) = one
      

      ! Gauss pt locations in master element
      !
      if(face.eq.1) then
         xLocal(1) = -dsqrt(one/three)
         yLocal(1) = -one
         xLocal(2) = dsqrt(one/three)
         yLocal(2) = -one
      elseif(face.eq.2) then
         xLocal(1) = one
         yLocal(1) = -dsqrt(one/three)
         xLocal(2) = one
         yLocal(2) = dsqrt(one/three)
      elseif(face.eq.3) then
         xLocal(1) = -dsqrt(one/three)
         yLocal(1) = one
         xLocal(2) = dsqrt(one/three)
         yLocal(2) = one
      elseif(face.eq.4) then
         xLocal(1) = -one
         yLocal(1) = dsqrt(one/three)
         xLocal(2) = -one
         yLocal(2) = -dsqrt(one/three)
      else
         write(*,*) 'face.ne.1,2,3,4'
         write(80,*) 'face.ne.1,2,3,4'
         call xit
      endif

      end subroutine xintSurf2D2pt

************************************************************************

      subroutine xintSurf2D3pt(face,xLocal,yLocal,w)

      !DEC$ ATTRIBUTES ALIAS:"xit"::XIT
      
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 2 gauss points for surface integration
      !
      !  xLocal(nIntPt): x coordinates for the integration pts
      !  yLocal(nIntPt): y coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights

      implicit none

      integer face

      real*8 xLocal(3),yLocal(3),w(3),zero,one,two,three,five,eight,nine
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,five=5.d0,
     +     eight=8.d0,nine=9.d0)


      ! Gauss weights
      !
      w(1) = five/nine
      w(2) = eight/nine
      w(3) = five/nine
      

      ! Gauss pt locations in master element
      !
      if(face.eq.1) then
         xLocal(1) = -dsqrt(three/five)
         yLocal(1) = -one
         xLocal(2) = zero
         yLocal(2) = -one
         xLocal(2) = dsqrt(three/five)
         yLocal(2) = -one
      elseif(face.eq.2) then
         xLocal(1) = one
         yLocal(1) = -dsqrt(three/five)
         xLocal(2) = one
         yLocal(2) = zero
         xLocal(3) = one
         yLocal(3) = dsqrt(three/five)
      elseif(face.eq.3) then
         xLocal(1) = -dsqrt(three/five)
         yLocal(1) = one
         xLocal(2) = zero
         yLocal(2) = one
         xLocal(3) = dsqrt(three/five)
         yLocal(3) = one
      elseif(face.eq.4) then
         xLocal(1) = -one
         yLocal(1) = dsqrt(three/five)
         xLocal(2) = -one
         yLocal(2) = zero
         xLocal(3) = -one
         yLocal(3) = -dsqrt(three/five)
      else
         write(*,*) 'face.ne.1,2,3,4'
         write(80,*) 'face.ne.1,2,3,4'
         call xit
      endif

      end subroutine xintSurf2D3pt

************************************************************************

      subroutine xintSurf3D1pt(face,xLocal,yLocal,zLocal,w)

      !DEC$ ATTRIBUTES ALIAS:"xit"::XIT
      
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 3D elements
      !  using 1 gauss point for surface integration
      !
      !  xLocal(nIntPt): x coordinates for the integration pts
      !  yLocal(nIntPt): y coordinates for the integration pts
      !  zLocal(nIntPt): z coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights

      implicit none

      integer face

      real*8 xLocal(1),yLocal(1),zLocal(1),w(1),zero,one,four
      parameter(zero=0.d0,one=1.d0,four=4.d0)


      ! Gauss weights
      !
      w(1) = four
      

      ! Gauss pt locations in master element
      !
      if(face.eq.1) then
         xLocal(1) = zero
         yLocal(1) = zero
         zLocal(1) = -one
      elseif(face.eq.2) then
         xLocal(1) = zero
         yLocal(1) = zero
         zLocal(1) = one
      elseif(face.eq.3) then
         xLocal(1) = zero
         yLocal(1) = -one
         zLocal(1) = zero
      elseif(face.eq.4) then
         xLocal(1) = one
         yLocal(1) = zero
         zLocal(1) = zero
      elseif(face.eq.5) then
         xLocal(1) = zero
         yLocal(1) = one
         zLocal(1) = zero
      elseif(face.eq.6) then
         xLocal(1) = -one
         yLocal(1) = zero
         zLocal(1) = zero
      else
         write(*,*) 'face.ne.1,2,3,4,5,6'
         write(80,*) 'face.ne.1,2,3,4,5,6'
         call xit
      endif

      end subroutine xintSurf3D1pt

************************************************************************

      subroutine xintSurf3D4pt(face,xLocal,yLocal,zLocal,w)

      !DEC$ ATTRIBUTES ALIAS:"xit"::XIT
      
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 3D elements
      !  using 4 gauss points for surface integration
      !
      !  xLocal(nIntPt): x coordinates for the integration pts
      !  yLocal(nIntPt): y coordinates for the integration pts
      !  yLocal(nIntPt): z coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights

      implicit none

      integer face

      real*8 xLocal(4),yLocal(4),zLocal(4),w(4),one,three
      parameter(one=1.d0,three=3.d0)


      ! Gauss weights
      !
      w(1) = one
      w(2) = one
      w(3) = one
      w(4) = one
      

      ! Gauss pt locations in master element
      !
      if(face.eq.1) then
         xLocal(1) = -dsqrt(one/three)
         yLocal(1) = -dsqrt(one/three)
         zLocal(1) = -one
         xLocal(2) = dsqrt(one/three)
         yLocal(2) = -dsqrt(one/three)
         zLocal(2) = -one
         xLocal(3) = dsqrt(one/three)
         yLocal(3) = dsqrt(one/three)
         zLocal(3) = -one
         xLocal(4) = -dsqrt(one/three)
         yLocal(4) = dsqrt(one/three)
         zLocal(4) = -one
      elseif(face.eq.2) then
         xLocal(1) = -dsqrt(one/three)
         yLocal(1) = -dsqrt(one/three)
         zLocal(1) = one
         xLocal(2) = dsqrt(one/three)
         yLocal(2) = -dsqrt(one/three)
         zLocal(2) = one
         xLocal(3) = dsqrt(one/three)
         yLocal(3) = dsqrt(one/three)
         zLocal(3) = one
         xLocal(4) = -dsqrt(one/three)
         yLocal(4) = dsqrt(one/three)
         zLocal(4) = one
      elseif(face.eq.3) then
         xLocal(1) = -dsqrt(one/three)
         yLocal(1) = -one
         zLocal(1) = -dsqrt(one/three)
         xLocal(2) = dsqrt(one/three)
         yLocal(2) = -one
         zLocal(2) = -dsqrt(one/three)
         xLocal(3) = dsqrt(one/three)
         yLocal(3) = -one
         zLocal(3) = dsqrt(one/three)
         xLocal(4) = -dsqrt(one/three)
         yLocal(4) = -one
         zLocal(4) = dsqrt(one/three)
      elseif(face.eq.4) then
         xLocal(1) = one
         yLocal(1) = -dsqrt(one/three)
         zLocal(1) = -dsqrt(one/three)
         xLocal(2) = one
         yLocal(2) = dsqrt(one/three)
         zLocal(2) = -dsqrt(one/three)
         xLocal(3) = one
         yLocal(3) = dsqrt(one/three)
         zLocal(3) = dsqrt(one/three)
         xLocal(4) = one
         yLocal(4) = -dsqrt(one/three)
         zLocal(4) = dsqrt(one/three)
      elseif(face.eq.5) then
         xLocal(1) = -dsqrt(one/three)
         yLocal(1) = one
         zLocal(1) = -dsqrt(one/three)
         xLocal(2) = dsqrt(one/three)
         yLocal(2) = one
         zLocal(2) = -dsqrt(one/three)
         xLocal(3) = dsqrt(one/three)
         yLocal(3) = one
         zLocal(3) = dsqrt(one/three)
         xLocal(4) = -dsqrt(one/three)
         yLocal(4) = one
         zLocal(4) = dsqrt(one/three)
      elseif(face.eq.6) then
         xLocal(1) = -one
         yLocal(1) = -dsqrt(one/three)
         zLocal(1) = -dsqrt(one/three)
         xLocal(2) = -one
         yLocal(2) = dsqrt(one/three)
         zLocal(2) = -dsqrt(one/three)
         xLocal(3) = -one
         yLocal(3) = dsqrt(one/three)
         zLocal(3) = dsqrt(one/three)
         xLocal(4) = -one
         yLocal(4) = -dsqrt(one/three)
         zLocal(4) = dsqrt(one/three)
      else
         write(*,*) 'face.ne.1,2,3,4,5,6'
         write(80,*) 'face.ne.1,2,3,4,5,6'
         call xit
      endif

      end subroutine xintSurf3D4pt
     
!************************************************************************

      subroutine calcShape2DLinear(nIntPt,xi_int,intpt,sh,dshxi)
      !
      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element


      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element
      !
      !                          eta
      !   4-----------3          |
      !   |           |          |
      !   |           |          |
      !   |           |          |
      !   |           |          |
      !   |           |          O--------- xi
      !   1-----------2        origin at center
      !
      !
      ! sh(i) = shape function of node i at the intpt.
      ! dshxi(i,j) = derivative wrt j direction of shape fn of node i
      !
      implicit none
      !
      integer intpt,nDim,nIntPt
      !
      real*8 xi_int(nIntPt,2),sh(4),dshxi(4,2),xi,eta
      !
      real*8 zero,one,fourth
      parameter(zero=0.d0,one=1.d0,fourth=1.d0/4.d0)
      

      ! Location in the master element
      !
      xi = xi_int(intpt,1)
      eta = xi_int(intpt,2)
      
      
      ! The shape functions
      !
      sh(1) = fourth*(one - xi)*(one - eta)
      sh(2) = fourth*(one + xi)*(one - eta)
      sh(3) = fourth*(one + xi)*(one + eta)
      sh(4) = fourth*(one - xi)*(one + eta)
      
      
      ! The first derivatives
      !
      dshxi(1,1) = -fourth*(one - eta)
      dshxi(1,2) = -fourth*(one - xi)
      dshxi(2,1) = fourth*(one - eta)
      dshxi(2,2) = -fourth*(one + xi)
      dshxi(3,1) = fourth*(one + eta)
      dshxi(3,2) = fourth*(one + xi)
      dshxi(4,1) = -fourth*(one + eta)
      dshxi(4,2) = fourth*(one - xi)
      

      return
      end subroutine calcShape2DLinear

************************************************************************

      subroutine calcShape3DLinear(nIntPt,xi_int,intpt,sh,dshxi)
      !
      !
      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element
      !
      ! This subroutine uses a 8-node linear 3D element as shown
      !
      !      8-----------7
      !     /|          /|       zeta
      !    / |         / |       
      !   5-----------6  |       |     eta
      !   |  |        |  |       |   /
      !   |  |        |  |       |  /
      !   |  4--------|--3       | /
      !   | /         | /        |/
      !   |/          |/         O--------- xi
      !   1-----------2        origin at cube center
      !
      !
      ! sh(i) = shape function of node i at the intpt.
      ! dshxi(i,j) = derivative wrt j direction of shape fn of node i
      ! d2shxi(i,j,k) = derivatives wrt j and k of shape fn of node i

      implicit none

      integer intpt,nDim,nIntPt,i,j

      real*8 xi_int(nIntPt,3),sh(8),dshxi(8,3)
      real*8 d2shxi(8,3,3),xi,eta,zeta

      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)
      

      ! Location in the master element
      !
      xi = xi_int(intpt,1)
      eta = xi_int(intpt,2)
      zeta = xi_int(intpt,3)
      !
      ! The shape functions
      !
      sh(1) = eighth*(one - xi)*(one - eta)*(one - zeta)
      sh(2) = eighth*(one + xi)*(one - eta)*(one - zeta)
      sh(3) = eighth*(one + xi)*(one + eta)*(one - zeta)
      sh(4) = eighth*(one - xi)*(one + eta)*(one - zeta)
      sh(5) = eighth*(one - xi)*(one - eta)*(one + zeta)
      sh(6) = eighth*(one + xi)*(one - eta)*(one + zeta)
      sh(7) = eighth*(one + xi)*(one + eta)*(one + zeta)
      sh(8) = eighth*(one - xi)*(one + eta)*(one + zeta)
      !
      ! The first derivatives
      !


      dshxi(1,1) = -eighth*(one - eta)*(one - zeta)
      dshxi(1,2) = -eighth*(one - xi)*(one - zeta)
      dshxi(1,3) = -eighth*(one - xi)*(one - eta)
      dshxi(2,1) = eighth*(one - eta)*(one - zeta)
      dshxi(2,2) = -eighth*(one + xi)*(one - zeta)
      dshxi(2,3) = -eighth*(one + xi)*(one - eta)
      dshxi(3,1) = eighth*(one + eta)*(one - zeta)
      dshxi(3,2) = eighth*(one + xi)*(one - zeta)
      dshxi(3,3) = -eighth*(one + xi)*(one + eta)
      dshxi(4,1) = -eighth*(one + eta)*(one - zeta)
      dshxi(4,2) = eighth*(one - xi)*(one - zeta)
      dshxi(4,3) = -eighth*(one - xi)*(one + eta)
      dshxi(5,1) = -eighth*(one - eta)*(one + zeta)
      dshxi(5,2) = -eighth*(one - xi)*(one + zeta)
      dshxi(5,3) = eighth*(one - xi)*(one - eta)
      dshxi(6,1) = eighth*(one - eta)*(one + zeta)
      dshxi(6,2) = -eighth*(one + xi)*(one + zeta)
      dshxi(6,3) = eighth*(one + xi)*(one - eta)
      dshxi(7,1) = eighth*(one + eta)*(one + zeta)
      dshxi(7,2) = eighth*(one + xi)*(one + zeta)
      dshxi(7,3) = eighth*(one + xi)*(one + eta)
      dshxi(8,1) = -eighth*(one + eta)*(one + zeta)
      dshxi(8,2) = eighth*(one - xi)*(one + zeta)
      dshxi(8,3) = eighth*(one - xi)*(one + eta)
      !
      ! The second derivatives
      !
      d2shxi = zero
      d2shxi(1,1,2) = eighth*(one - zeta)
      d2shxi(1,2,1) = d2shxi(1,1,2)
      d2shxi(1,1,3) = eighth*(one - eta)
      d2shxi(1,3,1) = d2shxi(1,1,3)
      d2shxi(1,2,3) = eighth*(one - xi)
      d2shxi(1,3,2) = d2shxi(1,2,3)
      d2shxi(2,1,2) = -eighth*(one - zeta)
      d2shxi(2,2,1) = d2shxi(2,1,2)
      d2shxi(2,1,3) = -eighth*(one - eta)
      d2shxi(2,3,1) = d2shxi(2,1,3)
      d2shxi(2,2,3) = eighth*(one + xi)
      d2shxi(2,3,2) = d2shxi(2,2,3)
      d2shxi(3,1,2) = eighth*(one - zeta)
      d2shxi(3,2,1) = d2shxi(2,1,2)
      d2shxi(3,1,3) = -eighth*(one + eta)
      d2shxi(3,3,1) = d2shxi(2,1,3)
      d2shxi(3,2,3) = -eighth*(one + xi)
      d2shxi(3,3,2) = d2shxi(2,2,3)
      d2shxi(4,1,2) = -eighth*(one - zeta)
      d2shxi(4,2,1) = d2shxi(2,1,2)
      d2shxi(4,1,3) = eighth*(one + eta)
      d2shxi(4,3,1) = d2shxi(2,1,3)
      d2shxi(4,2,3) = -eighth*(one - xi)
      d2shxi(4,3,2) = d2shxi(2,2,3)
      d2shxi(5,1,2) = eighth*(one + zeta)
      d2shxi(5,2,1) = d2shxi(2,1,2)
      d2shxi(5,1,3) = -eighth*(one - eta)
      d2shxi(5,3,1) = d2shxi(2,1,3)
      d2shxi(5,2,3) = -eighth*(one - xi)
      d2shxi(5,3,2) = d2shxi(2,2,3)
      d2shxi(6,1,2) = eighth*(one + zeta)
      d2shxi(6,2,1) = d2shxi(2,1,2)
      d2shxi(6,1,3) = eighth*(one - eta)
      d2shxi(6,3,1) = d2shxi(2,1,3)
      d2shxi(6,2,3) = -eighth*(one + xi)
      d2shxi(6,3,2) = d2shxi(2,2,3)
      d2shxi(7,1,2) = eighth*(one + zeta)
      d2shxi(7,2,1) = d2shxi(2,1,2)
      d2shxi(7,1,3) = eighth*(one + eta)
      d2shxi(7,3,1) = d2shxi(2,1,3)
      d2shxi(7,2,3) = eighth*(one + xi)
      d2shxi(7,3,2) = d2shxi(2,2,3)
      d2shxi(8,1,2) = -eighth*(one + zeta)
      d2shxi(8,2,1) = d2shxi(2,1,2)
      d2shxi(8,1,3) = -eighth*(one + eta)
      d2shxi(8,3,1) = d2shxi(2,1,3)
      d2shxi(8,2,3) = eighth*(one - xi)
      d2shxi(8,3,2) = d2shxi(2,2,3)
      
      return
      end subroutine calcShape3DLinear

!************************************************************************


      subroutine computeSurf(xLocal,yLocal,face,coords,sh,ds,normal)
      
      ! This subroutine computes the shape functions, derivatives
      !  of shape functions, and the length ds, so that one can
      !  do the numerical integration on the boundary for fluxes 
      !  on the 4-node quadrilateral elements


      !DEC$ ATTRIBUTES ALIAS:"xit"::XIT
      
      implicit none

      integer face

      real*8 xLocal,yLocal,ds,dshxi(4,2),sh(4),dXdXi,dXdEta,dYdXi
      real*8 dYdEta,one,coords(2,4),fourth,shape,normal(2,1)
      parameter(one=1.d0,fourth=1.d0/4.d0)

      sh(1) = fourth*(one - xLocal)*(one - yLocal)
      sh(2) = fourth*(one + xLocal)*(one - yLocal)
      sh(3) = fourth*(one + xLocal)*(one + yLocal)
      sh(4) = fourth*(one - xLocal)*(one + yLocal)
      
      dshxi(1,1) = -fourth*(one - yLocal)
      dshxi(1,2) = -fourth*(one - xLocal)
      dshxi(2,1) = fourth*(one - yLocal)
      dshxi(2,2) = -fourth*(one + xLocal)
      dshxi(3,1) = fourth*(one + yLocal)
      dshxi(3,2) = fourth*(one + xLocal)
      dshxi(4,1) = -fourth*(one + yLocal)
      dshxi(4,2) = fourth*(one - xLocal)

      dXdXi = dshxi(1,1)*coords(1,1)+dshxi(2,1)*coords(1,2)
     +     + dshxi(3,1)*coords(1,3)+dshxi(4,1)*coords(1,4)
      dXdEta = dshxi(1,2)*coords(1,1)+dshxi(2,2)*coords(1,2)
     +     + dshxi(3,2)*coords(1,3)+dshxi(4,2)*coords(1,4)
      dYdXi = dshxi(1,1)*coords(2,1)+dshxi(2,1)*coords(2,2)
     +     + dshxi(3,1)*coords(2,3)+dshxi(4,1)*coords(2,4)
      dYdEta = dshxi(1,2)*coords(2,1)+dshxi(2,2)*coords(2,2)
     +     + dshxi(3,2)*coords(2,3)+dshxi(4,2)*coords(2,4)


      ! Jacobian of the mapping
      !
      if((face.eq.2).or.(face.eq.4)) then
         ds = dsqrt(dXdEta*dXdEta + dYdEta*dYdEta)
      elseif((face.eq.1).or.(face.eq.3)) then
         ds = dsqrt(dXdXi*dXdXi + dYdXi*dYdXi)
      else
         write(*,*) 'never should get here'
         call xit
      endif


      ! Surface normal, outward pointing in this case. Useful for
      !  ``follower'' type loads. The normal is referential or spatial
      !  depending on which coords were supplied to this subroutine
      !  (NOT fully tested)
      !
      if((face.eq.2).or.(face.eq.4)) then
         normal(1,1) = dYdEta/dsqrt(dXdEta*dXdEta + dYdEta*dYdEta)
         normal(2,1) = -dXdEta/dsqrt(dXdEta*dXdEta + dYdEta*dYdEta)
         if(face.eq.4) normal = -normal
      elseif((face.eq.1).or.(face.eq.3)) then
         normal(1,1) = dYdXi/dsqrt(dXdXi*dXdXi + dYdXi*dYdXi)
         normal(2,1) = -dXdXi/dsqrt(dXdXi*dXdXi + dYdXi*dYdXi)
         if(face.eq.3) normal = -normal
      else
         write(*,*) 'never should get here'
         call xit
      endif

      return
      end subroutine computeSurf

************************************************************************

      subroutine computeSurf3D(xLocal,yLocal,zLocal,face,coords,sh,
     +     dshxi,dA,normal)

      ! This subroutine computes the shape functions, derivatives
      !  of shape functions, and the area dA, so that one can
      !  do the numerical integration on the boundary for fluxes 
      !  on the 8-node brick elements


      !DEC$ ATTRIBUTES ALIAS:"xit"::XIT
      
      implicit none

      integer face,stat,i,j,k

      real*8 xLocal,yLocal,zLocal,dA,dshxi(8,3),sh(8),zero,dsh(8,3),one
      real*8 coords(3,8),two,eighth,mapJ(3,3),mag,normal(3,1)

      real*8 dXdXi,dXdEta,dXdZeta,dYdXi,dYdEta,dYdZeta,dZdXi,dZdEta
      real*8 dZdZeta

      parameter(one=1.d0,two=2.d0,eighth=1.d0/8.d0,zero=0.d0)

      ! The shape functions
      !
      sh(1) = eighth*(one - xLocal)*(one - yLocal)*(one - zLocal)
      sh(2) = eighth*(one + xLocal)*(one - yLocal)*(one - zLocal)
      sh(3) = eighth*(one + xLocal)*(one + yLocal)*(one - zLocal)
      sh(4) = eighth*(one - xLocal)*(one + yLocal)*(one - zLocal)
      sh(5) = eighth*(one - xLocal)*(one - yLocal)*(one + zLocal)
      sh(6) = eighth*(one + xLocal)*(one - yLocal)*(one + zLocal)
      sh(7) = eighth*(one + xLocal)*(one + yLocal)*(one + zLocal)
      sh(8) = eighth*(one - xLocal)*(one + yLocal)*(one + zLocal)


      ! Shape function derivatives
      !
      dshxi(1,1) = -eighth*(one - yLocal)*(one - zLocal)
      dshxi(1,2) = -eighth*(one - xLocal)*(one - zLocal)
      dshxi(1,3) = -eighth*(one - xLocal)*(one - yLocal)
      dshxi(2,1) = eighth*(one - yLocal)*(one - zLocal)
      dshxi(2,2) = -eighth*(one + xLocal)*(one - zLocal)
      dshxi(2,3) = -eighth*(one + xLocal)*(one - yLocal)
      dshxi(3,1) = eighth*(one + yLocal)*(one - zLocal)
      dshxi(3,2) = eighth*(one + xLocal)*(one - zLocal)
      dshxi(3,3) = -eighth*(one + xLocal)*(one + yLocal)
      dshxi(4,1) = -eighth*(one + yLocal)*(one - zLocal)
      dshxi(4,2) = eighth*(one - xLocal)*(one - zLocal)
      dshxi(4,3) = -eighth*(one - xLocal)*(one + yLocal)
      dshxi(5,1) = -eighth*(one - yLocal)*(one + zLocal)
      dshxi(5,2) = -eighth*(one - xLocal)*(one + zLocal)
      dshxi(5,3) = eighth*(one - xLocal)*(one - yLocal)
      dshxi(6,1) = eighth*(one - yLocal)*(one + zLocal)
      dshxi(6,2) = -eighth*(one + xLocal)*(one + zLocal)
      dshxi(6,3) = eighth*(one + xLocal)*(one - yLocal)
      dshxi(7,1) = eighth*(one + yLocal)*(one + zLocal)
      dshxi(7,2) = eighth*(one + xLocal)*(one + zLocal)
      dshxi(7,3) = eighth*(one + xLocal)*(one + yLocal)
      dshxi(8,1) = -eighth*(one + yLocal)*(one + zLocal)
      dshxi(8,2) = eighth*(one - xLocal)*(one + zLocal)
      dshxi(8,3) = eighth*(one - xLocal)*(one + yLocal)


      dXdXi = zero
      dXdEta = zero
      dXdZeta = zero
      dYdXi = zero
      dYdEta = zero
      dYdZeta = zero
      dZdXi = zero
      dZdEta = zero
      dZdZeta = zero
      do k=1,8
         dXdXi = dXdXi + dshxi(k,1)*coords(1,k)
         dXdEta = dXdEta + dshxi(k,2)*coords(1,k)
         dXdZeta = dXdZeta + dshxi(k,3)*coords(1,k)
         dYdXi = dYdXi + dshxi(k,1)*coords(2,k)
         dYdEta = dYdEta + dshxi(k,2)*coords(2,k)
         dYdZeta = dYdZeta + dshxi(k,3)*coords(2,k)
         dZdXi = dZdXi + dshxi(k,1)*coords(3,k)
         dZdEta = dZdEta + dshxi(k,2)*coords(3,k)
         dZdZeta = dZdZeta + dshxi(k,3)*coords(3,k)
      enddo


      ! Jacobian of the mapping
      !
      if((face.eq.1).or.(face.eq.2)) then
         ! zeta = constant on this face
         dA = dsqrt(
     +          (dYdXi*dZdEta - dYdEta*dZdXi)**two
     +        + (dXdXi*dZdEta - dXdEta*dZdXi)**two
     +        + (dXdXi*dYdEta - dXdEta*dYdXi)**two
     +        )
      elseif((face.eq.3).or.(face.eq.5)) then
         ! eta = constant on this face
         dA = dsqrt(
     +          (dYdXi*dZdZeta - dYdZeta*dZdXi)**two
     +        + (dXdXi*dZdZeta - dXdZeta*dZdXi)**two
     +        + (dXdXi*dYdZeta - dXdZeta*dYdXi)**two
     +        )
      elseif((face.eq.4).or.(face.eq.6)) then
         ! xi = constant on this face
         dA = dsqrt(
     +          (dYdEta*dZdZeta - dYdZeta*dZdEta)**two
     +        + (dXdEta*dZdZeta - dXdZeta*dZdEta)**two
     +        + (dXdEta*dYdZeta - dXdZeta*dYdEta)**two
     +        )
         else
            write(*,*) 'never should get here'
            call xit
      endif


      ! Surface normal, outward pointing in this case. Useful for
      !  ``follower'' type loads. The normal is referential or spatial
      !  depending on which coords were supplied to this subroutine
      !  (NOT fully tested)
      !
      if((face.eq.1).or.(face.eq.2)) then
         ! zeta = constant on this face
         normal(1,1) = dYdXi*dZdEta - dYdEta*dZdXi
         normal(2,1) = dXdXi*dZdEta - dXdEta*dZdXi
         normal(3,1) = dXdXi*dYdEta - dXdEta*dYdXi
         if(face.eq.1) normal = -normal
      elseif((face.eq.3).or.(face.eq.5)) then
         ! eta = constant on this face
         normal(1,1) = dYdXi*dZdZeta - dYdZeta*dZdXi
         normal(2,1) = dXdXi*dZdZeta - dXdZeta*dZdXi
         normal(3,1) = dXdXi*dYdZeta - dXdZeta*dYdXi
         if(face.eq.5) normal = -normal
      elseif((face.eq.4).or.(face.eq.6)) then
         ! xi = constant on this face
         normal(1,1) = dYdEta*dZdZeta - dYdZeta*dZdEta
         normal(2,1) = dXdEta*dZdZeta - dXdZeta*dZdEta
         normal(3,1) = dXdEta*dYdZeta - dXdZeta*dYdEta
         if(face.eq.4) normal = -normal
      else
         write(*,*) 'never should get here'
         call xit
      endif
      mag = dsqrt(normal(1,1)**two+normal(2,1)**two+normal(3,1)**two)
      normal(1,1) = normal(1,1)/mag
      normal(2,1) = normal(2,1)/mag
      normal(3,1) = normal(3,1)/mag

      end subroutine computeSurf3D

************************************************************************

      subroutine mapShape2D(nNode,dshxi,coords,dsh,detMapJ,stat)
      !
      ! Map derivatives of shape fns from xi-eta-zeta domain
      !  to x-y-z domain.
      !
      implicit none
      !
      integer i,j,k,nNode,ieror,stat
      !
      real*8 dshxi(nNode,2),dsh(nNode,2),coords(3,nNode),mapJ(2,2),
     +  mapJ_inv(2,2),detmapJ
      !
      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)


      ! Calculate the mapping Jacobian matrix:
      !
      mapJ = zero
      do i=1,2
        do j=1,2
          do k=1,nNode
              mapJ(i,j) = mapJ(i,j) + dshxi(k,i)*coords(j,k)
          end do
        end do
      end do


      ! Calculate the inverse and the determinant of Jacobian
      !
      call matInv2D(mapJ,mapJ_inv,detMapJ,stat)


      ! Calculate first derivatives wrt x, y, z
      !
      dsh = transpose(matmul(mapJ_inv,transpose(dshxi)))
      

      return
      end subroutine mapShape2D

!*************************************************************************

      subroutine mapShape2Da(nNode,dshxi,coords,dsh,detMapJ,stat)
      !
      ! Map derivatives of shape fns from xi-eta-zeta domain
      !  to x-y-z domain.
      !
      ! This subroutine is exactly the same as the regular mapShape2D
      !  with the exception that coords(2,nNode) here and coords(3,nNode)
      !  in the regular.  I have noticed that a "heat transfer" and 
      !  "static" step uses MCRD=2, but for "coupled-temperature-displacement"
      !  you will get MCRD=3, even for a plane analysis.
      !
      implicit none
      !
      integer i,j,k,nNode,ieror,stat
      !
      real*8 dshxi(nNode,2),dsh(nNode,2),coords(2,nNode),mapJ(2,2),
     +  mapJ_inv(2,2),detmapJ
      !
      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)


      ! Calculate the mapping Jacobian matrix:
      !
      mapJ = zero
      do i=1,2
        do j=1,2
          do k=1,nNode
              mapJ(i,j) = mapJ(i,j) + dshxi(k,i)*coords(j,k)
          end do
        end do
      end do


      ! Calculate the inverse and the determinant of Jacobian
      !
      call matInv2D(mapJ,mapJ_inv,detMapJ,stat)


      ! Calculate first derivatives wrt x, y, z
      !
      dsh = transpose(matmul(mapJ_inv,transpose(dshxi)))
      

      return
      end subroutine mapShape2Da

************************************************************************

      subroutine mapShape3D(nNode,dshxi,coords,dsh,detMapJ,stat)
      !
      ! Map derivatives of shape fns from xi-eta-zeta domain
      !  to x-y-z domain.  This subroutine works for both 8-node
      !  linear and 20-node quadratic 3D elements.
      !
      implicit none

      integer i,j,k,nNode,ieror,stat

      real*8 dshxi(nNode,3),dsh(nNode,3),coords(3,nNode)
      real*8 mapJ(3,3),mapJ_inv(3,3),detmapJ

      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)
      

      ! Calculate the mapping Jacobian matrix:
      !
      mapJ = zero
      do i=1,3
        do j=1,3
          do k=1,nNode
              mapJ(i,j) = mapJ(i,j) + dshxi(k,i)*coords(j,k)
          end do
        end do
      end do


      ! Calculate the inverse and the determinant of Jacobian
      !
      call matInv3D(mapJ,mapJ_inv,detMapJ,stat)


      ! Calculate first derivatives wrt x, y, z
      !
      dsh = transpose(matmul(mapJ_inv,transpose(dshxi)))


      ! The second derivatives may be calculated.
      !

      return
      end subroutine mapShape3D

!****************************************************************************
!     Utility subroutines
!****************************************************************************

      subroutine matInv3D(A,A_inv,det_A,istat)
      !
      ! Returns A_inv, the inverse and det_A, the determinant
      ! Note that the det is of the original matrix, not the
      ! inverse
      !
      implicit none
      !
      integer istat
      !
      real*8 A(3,3),A_inv(3,3),det_A,det_A_inv


      istat = 1
      
      det_A = A(1,1)*(A(2,2)*A(3,3) - A(3,2)*A(2,3)) -
     +        A(2,1)*(A(1,2)*A(3,3) - A(3,2)*A(1,3)) +
     +        A(3,1)*(A(1,2)*A(2,3) - A(2,2)*A(1,3))
      
      if (det_A .le. 0.d0) then
        write(*,*) 'WARNING: subroutine matInv3D:'
        write(*,*) 'WARNING: det of mat=',det_A
        istat = 0
        return
      end if
          
      det_A_inv = 1.d0/det_A
        
      A_inv(1,1) = det_A_inv*(A(2,2)*A(3,3)-A(3,2)*A(2,3))
      A_inv(1,2) = det_A_inv*(A(3,2)*A(1,3)-A(1,2)*A(3,3))
      A_inv(1,3) = det_A_inv*(A(1,2)*A(2,3)-A(2,2)*A(1,3))
      A_inv(2,1) = det_A_inv*(A(3,1)*A(2,3)-A(2,1)*A(3,3))
      A_inv(2,2) = det_A_inv*(A(1,1)*A(3,3)-A(3,1)*A(1,3))
      A_inv(2,3) = det_A_inv*(A(2,1)*A(1,3)-A(1,1)*A(2,3))
      A_inv(3,1) = det_A_inv*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
      A_inv(3,2) = det_A_inv*(A(3,1)*A(1,2)-A(1,1)*A(3,2))
      A_inv(3,3) = det_A_inv*(A(1,1)*A(2,2)-A(2,1)*A(1,2))
      

      return
      end subroutine matInv3D

!****************************************************************************

      subroutine matInv2D(A,A_inv,det_A,istat)
      !
      ! Returns A_inv, the inverse, and det_A, the determinant
      ! Note that the det is of the original matrix, not the
      ! inverse
      !
      implicit none
      !
      integer istat
      !
      real*8 A(2,2),A_inv(2,2),det_A,det_A_inv

      
      istat = 1
      
      det_A = A(1,1)*A(2,2) - A(1,2)*A(2,1)
        
      if (det_A .le. 0.d0) then
        write(*,*) 'WARNING: subroutine matInv2D:'
        write(*,*) 'WARNING: det of mat=',det_A
        istat = 0
        return
      end if
            
      det_A_inv = 1.d0/det_A
          
      A_inv(1,1) =  det_A_inv*A(2,2)
      A_inv(1,2) = -det_A_inv*A(1,2)
      A_inv(2,1) = -det_A_inv*A(2,1)
      A_inv(2,2) =  det_A_inv*A(1,1)


      return
      end subroutine matInv2D

!****************************************************************************

      subroutine mdet(A,det)
      !
      ! This subroutine calculates the determinant
      ! of a 3 by 3 matrix [A]
      !
      implicit none
      !
      real*8  A(3,3),det


      det = A(1,1)*A(2,2)*A(3,3) 
     +	  + A(1,2)*A(2,3)*A(3,1)
     +	  + A(1,3)*A(2,1)*A(3,2)
     +	  - A(3,1)*A(2,2)*A(1,3)
     +	  - A(3,2)*A(2,3)*A(1,1)
     +	  - A(3,3)*A(2,1)*A(1,2)


      return
      end subroutine mdet
	
!****************************************************************************

      subroutine onem(A)
      !
      ! This subroutine stores the identity matrix in the
      ! 3 by 3 matrix [A]
      !
      implicit none
      !
      integer i,j
      !
      real*8 A(3,3)


      do i=1,3
         do J=1,3
	    if (i .eq. j) then
              A(i,j) = 1.0
            else
              A(i,j) = 0.0
            end if
         end do
      end do


      return
      end subroutine onem

****************************************************************************

    
