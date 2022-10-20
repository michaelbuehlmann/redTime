!     Code for Anisotropies in the Microwave Background
!     by Antony Lewis (http://cosmologist.info/) and Anthony Challinor
!     See readme.html for documentation. This is a sample driver routine that reads
!     in one set of parameters and produdes the corresponding output. 

    program driver
        use IniFile
        use CAMB
        use LambdaGeneral
        use Lensing
        use AMLUtils
        use Transfer
        use constants
        use Bispectrum

        implicit none

        Type(CAMBparams) P
        
        character(LEN=Ini_max_string_len) numstr, VectorFileName, &
            InputFile, ScalarFileName, TensorFileName, TotalFileName, LensedFileName,&
            LensedTotFileName, LensPotentialFileName, PowerFileName
        integer, parameter :: nmodels = 150; 
        integer :: i,npoints,newpts,modelno
        character(LEN=Ini_max_string_len) TransferFileNames(max_transfer_redshifts), &
               MatterPowerFileNames(max_transfer_redshifts), outroot, version_check
        logical bad

        real :: minkh, dlnkh, maxkh
	real, dimension(1000) :: outpower
        real(dl),dimension(1000) :: k,outpower_double, outpower_new, d2outpower

        real(dl) output_factor, Age, nmassive, omh2
	real(dl), dimension(nmodels) :: ombh2, omch2, omnuh2, omk, w0, wa, running, n_s, hubble,sigma8
        real(dl), parameter :: RECFAST_fudge_default = 1.14_dl
        real(dl), parameter :: RECFAST_fudge_default2 = 1.105d0 
        real(dl) kh, h, x, win, delta, k_new 
        real(dl) lnk, dlnk, lnko
        real(dl) dsig8, dsig8o, sig8, sig8o, sigr8
	real(dl),dimension(5) :: any
        real(dl), dimension (8) :: maxparams=(/0.155, 0.0235, 0.9, 0.85, 1.05, -0.7, 0.6, 0.01/)
        real(dl), dimension (8) :: minparams=(/0.12, 0.0215, 0.7, 0.55, 0.85,  -1.3, -1.0, 0.0/)

        open(unit=11, file="SuperCoyote/s-lhs.150.13_1", action="read")

	do i = 1, nmodels
           read(11,*)  omh2, ombh2(i), sigma8(i), hubble(i), n_s(i), w0(i), wa(i), omnuh2(i), any(1:5)
           omh2 = omh2*(maxparams(1)-minparams(1))+minparams(1)
	   ombh2(i) = ombh2(i)*(maxparams(2)-minparams(2))+minparams(2)
	   omnuh2(i) = omnuh2(i)*(maxparams(8)-minparams(8))+minparams(8)
	   sigma8(i) = sigma8(i)*(maxparams(3)-minparams(3))+minparams(3)
	   hubble(i) = hubble(i)*(maxparams(4)-minparams(4))+minparams(4)
	   n_s(i) = n_s(i)*(maxparams(5)-minparams(5))+minparams(5)
	   w0(i) = w0(i)*(maxparams(6)-minparams(6))+minparams(6)
	   wa(i) = wa(i)*(maxparams(7)-minparams(7))+minparams(7)
	   omch2(i) = omh2-ombh2(i)-omnuh2(i)
        enddo
        
        close(11)

!!$        omh2 = 0.5*(maxparams(1)-minparams(1))+minparams(1)
!!$	ombh2(1) = 0.5*(maxparams(2)-minparams(2))+minparams(2)
!!$	omnuh2(1) = 0.5*(maxparams(8)-minparams(8))+minparams(8)
!!$        omch2(1) = omh2-ombh2(1)-omnuh2(1)
!!$	w0(1) = 0.5*(maxparams(6)-minparams(6))+minparams(6)
!!$	wa(1) = 0.5*(maxparams(7)-minparams(7))+minparams(7)
!!$	n_s(1) = 0.5*(maxparams(5)-minparams(5))+minparams(5)
!!$	hubble(1) = 0.5*(maxparams(4)-minparams(4))+minparams(4)
!!$	sigma8(1) = 0.5*(maxparams(3)-minparams(3))+minparams(3)

       do modelno = 1, nmodels
          P%thismodel = modelno
          
          
          outroot = 'testdriver'
          if (outroot /= '') outroot = trim(outroot) // '_'

          call CAMB_SetDefParams(P)

          P%WantScalars = .false.
          P%WantVectors = .false.
          P%WantTensors = .false.

          P%OutputNormalization=outNone

          output_factor = 7.4311e12

          P%tcmb   = 2.726_dl
          P%yhe    = 0.24_dl

          P%WantCls= P%WantScalars .or. P%WantTensors .or. P%WantVectors

          P%WantTransfer=.true.

          P%NonLinear = NonLinear_none

          P%DoLensing = .false.
          if (P%WantCls) then
             if (P%WantScalars  .or. P%WantVectors) then
                P%Max_l = 2200
                P%Max_eta_k = P%Max_l*2._dl
                if (P%WantScalars) then
                   P%DoLensing = .true.
                   if (P%DoLensing) lensing_method = 1
                end if
             end if
          endif
          vec_sig0 = 1
          Magnetic = 0


          if (P%WantTensors) then
             P%Max_l_tensor = 1500
             P%Max_eta_k_tensor =  P%Max_l_tensor*2._dl
          endif






!!!! Transfer parameters

          P%Transfer%high_precision=  .true. 
          P%transfer%kmax          =  10
          P%transfer%k_per_logint  =  0 !!! this needs to be changed for high precision
          P%transfer%num_redshifts =  1

          transfer_interp_matterpower = .true.
!!!        transfer_power_var = Ini_read_int('transfer_power_var',transfer_power_var)       ???     what is this ???
          P%transfer%num_redshifts = 1
          P%transfer%redshifts(1)  = 0._dl
          transferFileNames(1)     = trim(outroot)//"transfer_out.dat"
          MatterPowerFilenames(1)  = trim(outroot)//"matterpower_out.dat"


          P%transfer%kmax=P%transfer%kmax*(P%h0/100._dl)

          Ini_fail_on_not_found = .false. 

!!! Re-ionization parameters
          P%Reion%Reionization = .true.
          P%Reion%use_optical_depth = .true.
          P%Reion%optical_depth = 0.09_dl
          P%Reion%redshift = 11._dl
          P%Reion%delta_redshift = 1.5_dl
          P%Reion%fraction = -1._dl

!!!!! Initial power spectrum parameters 

          P%InitPower%nn = 1 !!! I only want one power spectrum 
          P%InitPower%k_0_scalar = 0.05_dl
          P%InitPower%k_0_tensor = 0.05_dl
          P%InitPower%ScalarPowerAmp(1) = 2.5e-9





!!! Recombination parameters
          P%Recomb%RECFAST_fudge_He = 0.86_dl
          P%Recomb%RECFAST_Heswitch = 6._dl
          P%Recomb%RECFAST_Hswitch = .true.
          P%Recomb%RECFAST_fudge = 1.14_dl
          P%Recomb%RECFAST_fudge = P%Recomb%RECFAST_fudge - (RECFAST_fudge_default - RECFAST_fudge_default2)



!!! Bi-spectrum parameters - I think this is enough 
          !!      B%do_lensing_bispectrum = .false.
          !!      B%do_primordial_bispectrum = .false.

          if (P%WantScalars .or. P%WantTransfer) then
             P%Scalar_initial_condition = 1
             if (P%Scalar_initial_condition == initial_vector) then
                P%InitialConditionVector=0
                P%InitialConditionVector(1)  = -1
             end if

          end if


          if (P%WantScalars) then
             ScalarFileName = trim(outroot)//"scalCls.dat"
             LensedFileName =  trim(outroot) //"lensedCls.dat"
             LensPotentialFileName =  "lenspotentialCls.dat"
             if (LensPotentialFileName/='') LensPotentialFileName = concat(outroot,LensPotentialFileName)
          end if

          Ini_fail_on_not_found = .false. 

          !optional parameters controlling the computation
          !Accuracy parameters - these are important!!!

          P%AccuratePolarization = .true.
          P%AccurateReionization = .false.
          P%AccurateBB = .false.


          !Mess here to fix typo with backwards compatibility
          DoLateRadTruncation = .true.
          DoTensorNeutrinos = .true.
          FeedbackLevel = 1

          P%MassiveNuMethod  = Nu_best

          ThreadNum      = 1
          AccuracyBoost  = 1
          lAccuracyBoost = 1
          HighAccuracyDefault = .false.
          !!       use_spline_template = Ini_Read_Logical('use_spline_template',use_spline_template)
          if (HighAccuracyDefault) then
             P%Max_eta_k=max(min(P%max_l,3000)*2.5_dl,P%Max_eta_k)
          end if


!  Read initial parameters.
!  Set cosmo here

       
          w_lam = w0(modelno)
          w_lamprime = -1._dl*wa(modelno) !!! this is because w=w0+(a-1)*wa in equations.f90 are reversed


          cs2_lam = 1.d0
          P%H0     = hubble(modelno)*100._dl
          write(*,*) 'w0 = ', w_lam, 'wa = ', w_lamprime, 'h0 = ', hubble(modelno), &
               &'ns = ', n_s(modelno), 'om = ', omch2(modelno), 'omb= ', ombh2(modelno) 
          P%omegab = ombh2(modelno)/(P%H0/100)**2
          P%omegac = omch2(modelno)/(P%H0/100)**2
          P%omegan = omnuh2(modelno)/(P%H0/100)**2
!!$          P%omegan = 0._dl
          P%omegav = 1- P%omegab-P%omegac - P%omegan
          

          P%InitPower%an(1) = n_s(modelno)
          P%InitPower%n_run(1) = 0; 

          if (P%WantTensors) then
             P%InitPower%ant(1) = 0
             P%InitPower%rat(1) = 1
          end if
          if (P%omegan > 0) then
             P%Num_Nu_massless  = 1.046_dl!!! I don't know about this one...
             nmassive = 2._dl 
             !Store fractional numbers in the massless total
             P%Num_Nu_massive   = int(nmassive+1e-6)
             P%Num_Nu_massless  = P%Num_Nu_massless + nmassive-P%Num_Nu_massive
             P%nu_mass_splittings = .true.
             P%Nu_mass_eigenstates = 1!P%Num_Nu_massive
             if (P%Nu_mass_eigenstates > max_nu) stop 'too many mass eigenstates'
             P%Nu_mass_degeneracies(1)=2
             P%Nu_mass_degeneracies(2)=1
             P%Nu_mass_fractions(:)=1
!!$             P%Nu_mass_fractions(1) =0.03
!!$             P%Nu_mass_fractions(2) =0.07
!!$             P%Nu_mass_fractions(3) =0.9
!!$
!!$             P%Nu_mass_degeneracies(:)= 1
!!$             P%Nu_mass_fractions(:) =0.03
!!$             P%Nu_mass_fractions(2) =0.07
!!$             P%Nu_mass_fractions(3) =0.9
          else
             P%Num_Nu_massless  = 3.046  !!! I don't know about this one...
             nmassive = 0.0_dl
             !Store fractional numbers in the massless total
             P%Num_Nu_massive   = int(nmassive+1e-6)
             P%Num_Nu_massless  = P%Num_Nu_massless + nmassive-P%Num_Nu_massive
             P%nu_mass_splittings = .true.
             P%Nu_mass_eigenstates = 1
             P%Nu_mass_degeneracies(1)= P%Num_nu_massive !! ok for now
             P%Nu_mass_fractions(1)=1  !!! again this is ok for massless neutrinos
          endif


          if (.not. CAMB_ValidateParams(P)) stop 'Stopped due to parameter error'

          if (FeedbackLevel > 0) then
             Age = CAMB_GetAge(P) 
             write (*,'("Age of universe/GYr  = ",f7.3)') Age  
          end if

          if (global_error_flag==0) call CAMB_GetResults(P)
          if (global_error_flag/=0) then
             write (*,*) 'Error result '//trim(global_error_message)
             stop
          endif

!!$          call Transfer_SaveMatterPower(MT,MatterPowerFileNames)
!!$          call Transfer_output_sig8(MT)

          minkh = 1e-4
          dlnkh = 0.02
          npoints = log(MT%TransferData(1,MT%num_q_trans,1)/minkh)/dlnkh+1	

          call Transfer_GetMatterPower(MT,outpower, 1, 1, minkh, dlnkh, npoints)
          outpower(:) = outpower(:)*sigma8(modelno)**2.d0/MT%sigma_8(1,1)**2.d0

          lnko=0
          dsig8o=0
          sig8=0
          sig8o=0
          sigr8=8._dl
          if (modelno < 10) then
             write (PowerFileName, "(A22,I1,A4)") "SuperCoyote/pk_lin_M00", modelno, ".dat"
          else if (modelno < 100) then
             write (PowerFileName, "(A21,I2,A4)") "SuperCoyote/pk_lin_M0", modelno, ".dat"
          else 
             write (PowerFileName, "(A20,I3,A4)") "SuperCoyote/pk_lin_M", modelno, ".dat"
          endif
          open(unit=10, file=PowerFileName, action='write')

!!$          open(unit=10, file='SuperCoyote/pk_test.dat', action='write')
          do i=1, npoints
             k(i) = exp(log(minkh) + dlnkh*(i-1))
             x= k(i) *sigr8
             win =3._dl*(sin(x)-x*cos(x))/x**3._dl
             lnk=k(i)!log(k)
             if (i==1) then
                dlnk=0.5_dl 
             else
                dlnk=lnk-lnko
             end if
             dsig8=(win*lnk)**2*outpower(i)
             sig8=sig8+(dsig8+dsig8o)*dlnk/2
             dsig8o=dsig8
             lnko=lnk
          end do
          sig8 = sig8/acos(-1.)**2._dl/2._dl
          write(*,*) "New sigma8 = ", sqrt(sig8), "sigma8 wanted = ", sigma8(modelno)
          MT%sigma_8(1,1) = sqrt(sig8)

          minkh = 5e-3
          maxkh = 1
          newpts = 100
          dlnkh = (log(maxkh)-log(minkh))/(newpts-1)
!!$          dlnkh = 0.02

          k(:) = k(:)*hubble(modelno) !!! this is the k spacing that we have in 1/Mpc
          outpower_double(:) = outpower(:)/hubble(modelno)**3. !! this is the P(k) in units of Mpc^3

          call spline(k,outpower_double,npoints,1d30,1d30,d2outpower)

          do i=1, newpts
             k_new = exp(log(minkh) + dlnkh*(i-1)) !!! this is the k spacing that we want in 1/Mpc
             call my_splint(k,outpower_double,d2outpower,npoints,k_new, outpower_new(i))
             write(10,*) k_new, outpower_new(i)
          end do

   end do
!!$        if (P%WantCls) then
!!$  
!!$         if (P%OutputNormalization == outCOBE) then
!!$
!!$            if (P%WantTransfer) call Transfer_output_Sig8AndNorm(MT)
!!$           
!!$          end if
!!$
!!$         call output_cl_files(ScalarFileName, TensorFileName, TotalFileName, &
!!$              LensedFileName, LensedTotFilename, output_factor)
!!$              
!!$         call output_lens_pot_files(LensPotentialFileName, output_factor)
!!$
!!$         if (P%WantVectors) then
!!$           call output_veccl_files(VectorFileName, output_factor)
!!$         end if
!!$
!!$        end if

        call CAMB_cleanup             
        end program driver

        SUBROUTINE my_splint(xa,ya,y2a,n,x,y)
          INTEGER n
          DOUBLE PRECISION x,y,xa(n),y2a(n),ya(n)
          INTEGER k,khi,klo
          DOUBLE PRECISION a,b,h
          klo=1
          khi=n
1         if (khi-klo.gt.1) then
             k=(khi+klo)/2
             if(xa(k).gt.x)then
                khi=k
             else
                klo=k
             endif
             goto 1
          endif
          h=xa(khi)-xa(klo)
          if (h.eq.0.) pause 'bad xa input in splint'
          a=(xa(khi)-x)/h
          b=(x-xa(klo))/h
          y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
          return
        END SUBROUTINE my_splint

