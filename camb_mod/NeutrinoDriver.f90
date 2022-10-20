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
        integer :: i,npoints,nmodels, modelno
        character(LEN=Ini_max_string_len) TransferFileNames(max_transfer_redshifts), &
               MatterPowerFileNames(max_transfer_redshifts), outroot, version_check
        logical bad

        real :: minkh, dlnkh, scalef
        integer, parameter ::  itf = 6
        real :: outpower(itf,1000), growthrate(1000)
        real(dl) output_factor, Age, nmassive, nu_total_mass
	real(dl), dimension(100) :: ombh2, omch2, w0, wa, n_s, hubble,sigma8, nu1, nu2, nu3
        real(dl), parameter :: RECFAST_fudge_default = 1.14_dl
        real(dl), parameter :: RECFAST_fudge_default2 = 1.105d0 
        real(dl) kh, k, h, x, win, delta
        real(dl) lnk, dlnk, lnko
        real(dl) dsig8, dsig8o, sig8, sig8o, sigr8
	real(dl),dimension(5) :: any
        real(dl), dimension (9) :: minparams=(/0.12, 0.0215, 0.55, 0.85, -1.3, -1.0, 0., 0., 0./)
        real(dl), dimension (9) :: maxparams=(/0.155, 0.0235, 0.85, 1.05,-0.7, 0.6,  0.1, 0.5, 1./)



	nmodels = 100

	do i = 1, nmodels
	   open(unit=11, file="NeutrinoGrowthRates/s-lhs.100.9_1", action="read")
           read(11,*)  omch2(i), ombh2(i), hubble(i), n_s(i), w0(i), wa(i), nu1(i), nu2(i), nu3(i)
	   omch2(i) = omch2(i)*(maxparams(1)-minparams(1))+minparams(1)
	   ombh2(i) = ombh2(i)*(maxparams(2)-minparams(2))+minparams(2)
           omch2(i) = omch2(i) - ombh2(i)
	   hubble(i) = hubble(i)*(maxparams(3)-minparams(3))+minparams(3)
	   n_s(i) = n_s(i)*(maxparams(4)-minparams(4))+minparams(4)
	   w0(i) = w0(i)*(maxparams(5)-minparams(5))+minparams(5)
	   wa(i) = wa(i)*(maxparams(6)-minparams(6))+minparams(6)
           nu1(i) = nu1(i)*(maxparams(7)-minparams(7))+minparams(7)
           nu2(i) = nu2(i)*(maxparams(8)-minparams(8))+minparams(8)
           nu3(i) = nu3(i)*(maxparams(9)-minparams(9))+minparams(9)
	   omnuh2(i) = (nu1(i)+nu2(i)+nu3(i))/93.1
!           write(*,*)  omch2(i), ombh2(i), sigma8(i), hubble(i), n_s(i), w0(i), wa(i), omnuh2(i), nu1(i), nu2(i), nu3(i)
        enddo

!	ombh2(1) = 0.0224_dl
!	omch2(1) = 0.1072_dl
!	omnuh2(1) = 0._dl
!	omk(1) = 0._dl
!	w0(1) = -1._dl
!	wa(1) = 0._dl
!	running(1) = 0._dl
!	n_s(1) = 0.97_dl
!	hubble(1) = 72._dl
!	sigma8(1) = 0.8_dl

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
          P%transfer%kmax          =  2
          P%transfer%k_per_logint  =  0 !!! this needs to be changed for high precision
          P%transfer%num_redshifts =  1

          transfer_interp_matterpower = .true.
!!!        transfer_power_var = Ini_read_int('transfer_power_var',transfer_power_var)       ???     what is this ???
          P%transfer%num_redshifts = itf
          P%transfer%redshifts(6)  = 0._dl
          P%transfer%redshifts(5)  = 0.48_dl
          P%transfer%redshifts(4)  = 0.49_dl
          P%transfer%redshifts(3)  = 0.5_dl
          P%transfer%redshifts(2)  = 0.51_dl
          P%transfer%redshifts(1)  = 0.52_dl

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
          P%InitPower%ScalarPowerAmp(1) = 2.1e-9





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
               &'ns = ', n_s(modelno)
          P%omegab = ombh2(modelno)/(P%H0/100)**2
          P%omegac = omch2(modelno)/(P%H0/100)**2
          P%omegan = omnuh2(modelno)/(P%H0/100)**2
!!$          P%omegan = 0._dl
          P%omegav = 1 - P%omegab-P%omegac - P%omegan

!!$          P%omegav = 1- P%omegab-P%omegac - P%omegan
          

          P%InitPower%an(1) = n_s(modelno)
          P%InitPower%n_run(1) = 0; 

          if (P%WantTensors) then
             P%InitPower%ant(1) = 0
             P%InitPower%rat(1) = 1
          end if

          nu_total_mass = nu1(modelno)+nu2(modelno)+nu3(modelno)

          if (P%omegan > 0) then
             P%Num_Nu_massless  = 0.046_dl!!! I don't know about this one...
             nmassive = 3._dl 
             P%Nu_mass_degeneracies(1)=1
             P%Nu_mass_degeneracies(2)=1
             P%Nu_mass_degeneracies(3)=1

             P%Nu_mass_fractions(1) =nu1(modelno)/nu_total_mass
             P%Nu_mass_fractions(2) =nu2(modelno)/nu_total_mass
             P%Nu_mass_fractions(3) =nu3(modelno)/nu_total_mass

             if (nu1(modelno)==0) then
                P%Num_Nu_massless = P%Num_Nu_massless + 1
                nmassive = nmassive -1 
                P%Nu_mass_degeneracies(3)=0
                P%Nu_mass_fractions(1) =nu2(modelno)/nu_total_mass
                P%Nu_mass_fractions(2) =nu3(modelno)/nu_total_mass
                P%Nu_mass_fractions(3) =0. 
             end if
             if (nu2(modelno)==0) then
                P%Num_Nu_massless = P%Num_Nu_massless + 1
                nmassive = nmassive -1 
                P%Nu_mass_degeneracies(3)=0
                P%Nu_mass_fractions(1) =nu1(modelno)/nu_total_mass
                P%Nu_mass_fractions(2) =nu3(modelno)/nu_total_mass
                P%Nu_mass_fractions(3) =0. 
                if (nu1(modelno)==0) then
                   P%Nu_mass_degeneracies(2)=0
                   P%Nu_mass_fractions(1) =nu3(modelno)/nu_total_mass
                   P%Nu_mass_fractions(2) =0.
                   P%Nu_mass_fractions(3) =0. 
                endif
             end if
             if (nu3(modelno)==0) then
                P%Num_Nu_massless = P%Num_Nu_massless + 1
                nmassive = nmassive -1 
                P%Nu_mass_degeneracies(3)=0
                P%Nu_mass_fractions(1) =nu1(modelno)/nu_total_mass
                P%Nu_mass_fractions(2) =nu2(modelno)/nu_total_mass
                P%Nu_mass_fractions(3) =0. 
                if (nu1(modelno)==0) then
                   P%Nu_mass_degeneracies(2)=0 
                   P%Nu_mass_fractions(1) =nu2(modelno)/nu_total_mass
                   P%Nu_mass_fractions(2) =0.
                endif
                if (nu2(modelno)==0) then
                   P%Nu_mass_degeneracies(2)=0
                   P%Nu_mass_fractions(1) =nu1(modelno)/nu_total_mass
                   P%Nu_mass_fractions(2) =0.
                endif
             end if

             !Store fractional numbers in the massless total
             P%Num_Nu_massive   = int(nmassive+1e-6)
             P%Num_Nu_massless  = P%Num_Nu_massless + nmassive-P%Num_Nu_massive
             P%nu_mass_splittings = .true.
             P%Nu_mass_eigenstates = P%Num_Nu_massive
             if (P%Nu_mass_eigenstates > max_nu) stop 'too many mass eigenstates'

             write(*,*) nu_total_mass, P%Num_Nu_massless, P%Num_Nu_massive, sum(P%Nu_mass_degeneracies(1:P%Nu_mass_eigenstates))
!!$             P%Nu_mass_fractions(:)=1
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


          do i=1,P%Transfer%num_redshifts
             call Transfer_GetMatterPower(MT,outpower(i,:), i, 1, minkh, dlnkh, npoints)
!!$             outpower(i,:) = outpower(i,:)*sigma8(modelno)**2.d0/MT%sigma_8(i,1)**2.d0
          enddo

          lnko=0
          dsig8o=0
          sig8=0
          sig8o=0
          sigr8=8._dl
          if (modelno < 10) then
             write (PowerFileName, "(A25,I1,A4)") "NeutrinoGrowthRates/f_M00", modelno, ".dat"
          else if (modelno < 100) then
             write (PowerFileName, "(A24,I2,A4)") "NeutrinoGrowthRates/f_M0", modelno, ".dat"
          else 
             write (PowerFileName, "(A23,I3,A4)") "NeutrinoGrowthRates/f_M", modelno, ".dat"
          endif
          open(unit=10, file=PowerFileName, action='write')


!!$          do i=1, npoints
!!$             k = exp(log(minkh) + dlnkh*(i-1))
!!$             x= k *sigr8
!!$             win =3._dl*(sin(x)-x*cos(x))/x**3._dl
!!$             lnk=k!log(k)
!!$             if (i==1) then
!!$                dlnk=0.5_dl 
!!$             else
!!$                dlnk=lnk-lnko
!!$             end if
!!$             dsig8=(win*lnk)**2*outpower(i)
!!$             sig8=sig8+(dsig8+dsig8o)*dlnk/2
!!$             dsig8o=dsig8
!!$             lnko=lnk
!!$             write(10,*) k, outpower(6,i)
!!$          end do
!!$          sig8 = sig8/acos(-1.)**2._dl/2._dl
!!$          write(*,*) "New sigma8 = ", sqrt(sig8), "sigma8 wanted = ", sigma8(modelno)
!!$          MT%sigma_8(1,1) = sqrt(sig8)
!!$          close(unit=10)

       

        scalef = 1./(P%Transfer%redshifts(3)+1.)
!!$        write(*,*) scalef
        growthrate(:) = -sqrt(outpower(1,:)/outpower(6,:))+8.*sqrt(outpower(2,:)/outpower(6,:))&
             &-8.*sqrt(outpower(4,:)/outpower(6,:))+sqrt(outpower(5,:)/outpower(6,:))
        growthrate(:) = growthrate(:)/12./(P%Transfer%redshifts(2)-P%Transfer%redshifts(3))*(-1./scalef**2) !!
        growthrate(:) = growthrate(:)*scalef/sqrt(outpower(3,:)/outpower(6,:)) !! converts dD/dz to dD/dlna

!!$        open(unit=20,file='NeutrinoGrowthRates/testpk.dat')
        do i=2,npoints-1
           write(10,*) minkh*exp((i-1)*dlnkh), outpower(6,i), growthrate(i)!(growthrate(i-1)+growthrate(i)+growthrate(i+1))/3.d0
        enddo
        close(unit=10)
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

