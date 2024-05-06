module global
    implicit none
    save
    integer , parameter :: dp = selected_real_kind(15 , 30)
    real(kind = dp)     :: pi = 4.0_dp*ATAN(1.0_dp)
end module global

program FT
    use global
    implicit none
    real(kind = dp)                  :: del_t = 0.01_dp , n_in = 0.0_dp , k_in = 0.0_dp , a_factor = 4.8_dp
    real(kind = dp)                  :: period_1
    complex(kind = dp)               :: W_const
    integer                          :: i = 0 , j = 0 , N_val = 10000 , m = 0 , istat = 0 , file_unit = 20
    complex(kind = dp) , allocatable :: W_matrix(: , :) , FT_vector(:) , h_vector(:) , f_vector(:)
    logical                          :: lexist
    character(len = 20)              :: filename    = 'testnorm.dat'
    character(len = 20)              :: window_case = 'none'
    character(len = 20)              :: pulse_case  = 'sine_signal'
    period_1 = 1.0_dp                               !Defining the period so that the period with relation to del_t so that it is within range

    allocate(FT_vector(N_val))                              !Allocating each of the required arrays here using N_val
    if (istat /= 0) STOP 'failure allocating FT_vector'     !Error codes in the event that allocation fails
    allocate(h_vector(N_val))
    if (istat /= 0) STOP 'failure allocating h_vector'
    allocate(f_vector(N_val))
    if (istat /= 0) STOP 'failure allocating f_vector'
    allocate(W_matrix(N_val , N_val))
    if (istat /= 0) STOP 'failure allocating W_matrix'
    W_matrix = 0.0_dp
    W_const = EXP(2.0_dp*pi*cmplx(0.0_dp , 1.0_dp , kind = dp)/real(N_val , kind = dp))
                                                            !Defining W_const as a global variable
    if(mod(N_val , 2) /= 0) STOP 'failure of N_val'         !Stopping program if N is not mathematical

!--------------------------------------------------------------------------------------!
!Nested do loop iterates i, j to N_val in 1's. i, j are then used to calculate n_in and!
!k_in respectively, they correspond to the n, k used in the mathematical explanation of!
!DFT.                                                                                  !
!n_in is effectively iterated betwen ±N/2 and k_in is effectively iterated between 0,  !
!N-1, this is done through converting the real(i),real(j) to the required value each   !
!iteration. f_vector is also used to calculate the range of frequencies over which the !
!DFT tests.                                                                            !
!--------------------------------------------------------------------------------------!

    do i = 1 , N_val , 1
        do j = 1 , N_val , 1
            n_in            = real(i , kind = dp) - 0.5_dp*N_val
            k_in            = real(j , kind = dp) - 1       !i,j vals are converted to reals n_in, k_in for W_matrix
            W_matrix(i , j) = W_const**(real(n_in , kind = dp)*real(k_in , kind = dp))
            f_vector(i)     = n_in/(del_t*real(N_val , kind = dp))
        end do
    end do

!--------------------------------------------------------------------------------------!
!k do loop iterates between 0, N_val-1 which is then used to calculate the required    !
!vals for t_k (t_sum), which are then used to calculate the discretised signal to be   !
!read by the DFT.                                                                      !
!--------------------------------------------------------------------------------------!

    call windowing_subroutine(del_t , N_val , h_vector , a_factor , window_case , pulse_case)

    FT_vector = matmul(W_matrix , h_vector)

    inquire(file=filename, exist=lexist)
        if (lexist) then                                    !using lexist logical to determine how to treat file
            print*,("old data file overwritten")            !checking to ensure that no file is ballooning in size
            open (unit = file_unit , file = filename , status = "replace" , action = "write" , position = "append" , iostat = istat)
            if (istat /= 0) stop "error opening test file"
        else
            open (unit = file_unit , file = filename , status = "new" , action = "write" , position = "append" , iostat = istat)
            if (istat /= 0) stop "error opening test file"
        end if
    do m = 1 , N_val , 1
        write (unit = file_unit , fmt = * , iostat = istat) f_vector(m)%re , abs(conjg(FT_vector(m)))
                                                            !writing/plotting only reals
    end do
    if (istat /= 0) stop "error writing test file"
    close (unit = file_unit , iostat = istat)
    if (istat /= 0) stop "error closing test file"

!--------------------------------------------------------------------------------------!
!windowing_subroutine takes in constants such as del_t, N, h_vector which are defined  !
!or otherwise calculated outside of the subroutine. I have passed them into the subrou-!
!tine for generality. I also pass in the window_selector and pulse_selector case state-!
!ments. This means that the signal and windowing of the signal can be defined with the !
!other constants of the system without necessitating changing many variables or functi-!
!ons within the program.                                                               !
!--------------------------------------------------------------------------------------!

contains
    subroutine windowing_subroutine(del_t_in , N_in , h_vector_in , length_factor , window_selector ,&
    &pulse_selector)
        implicit none
        real(kind = dp) , intent(IN)                     :: del_t_in , length_factor
        integer , intent(IN)                             :: N_in
        complex(kind = dp) , allocatable , intent(INOUT) :: h_vector_in(:)
        character(len = 20) , intent(IN)                 :: window_selector , pulse_selector
        integer                                          :: k = 0 , istat = 0 , file_unit2 = 20
        real(kind = dp)                                  :: t_sum = 0.0_dp , bell_width = 0.0_dp , amp = 1.0_dp , length
        character(len = 20)                              :: filename2 = 'gauwintest.dat'

        length = real(N_in , kind = dp)/length_factor       !The length of the signal is defined as a fraction of N here, using length_factor input

        inquire(file=filename2, exist=lexist)
        if (lexist) then                                    !using lexits logical to determine how to treat file
            print*,("old data file overwritten")            !checking to ensure that no file is ballooning in size
            open (unit = file_unit2 , file = filename2 , status = "replace" , action = "write" ,&
            &position = "append" , iostat = istat)
            if (istat /= 0) stop "error opening test file"
        else
            open (unit = file_unit2 , file = filename2 , status = "new" , action = "write" ,&
            &position = "append" , iostat = istat)
            if (istat /= 0) stop "error opening test file"
        end if

        select case(window_selector)                        !case statement for the triangle window
            case('triangle')
                bell_width = real(N_in , kind = dp)*del_t_in!iterating to calculate N signal calculations
                do k = 0 , (N_in - 1) , 1
                    t_sum = del_t_in*real(k , kind = dp)    !converting Nth signal calculation to a time value
                    if((k + 1) <= N_in/2) then              !calculating lhs of triangular window
                        amp = t_sum*(2.0_dp/bell_width)
                    else                                    !calculating rhs of triangular window
                        amp = t_sum*(-2.0_dp/bell_width) + 2.0_dp
                    end if
                    write (unit = file_unit2 , fmt = * , iostat = istat) t_sum , amp
                    if (istat /= 0) stop "error writing test file"
                    call signal_selector(t_sum , length , amp , h_vector_in , k , pulse_selector)
                end do                                      !calling subroutine to signal, as requested by pulse_selector string
            case('cosine_bell')                             !case statement for cosine bell window
                bell_width = real(N_in , kind = dp)*del_t_in
                do k = 0 , (N_in - 1) , 1
                    t_sum = del_t_in*real(k , kind = dp)
                    amp   = -0.5_dp*(COS(2.0_dp*pi*t_sum/bell_width) - 1.0_dp)
                    write (unit = file_unit2 , fmt = * , iostat = istat) t_sum , amp
                    if (istat /= 0) stop "error writing test file"
                    call signal_selector(t_sum , length , amp , h_vector_in , k , pulse_selector)
                end do
            case('gaussian_bell')                           !case for Gaussian bell window
                bell_width = 0.32                           !defining bell width variable (sigma)
                do k = 0 , (N_in - 1) , 1
                    t_sum = del_t_in*real(k , kind = dp)
                    amp   = EXP(-0.5_dp*((k - N_in*0.5_dp)/(bell_width*0.5_dp*N_in))**2.0_dp)
                    write (unit = file_unit2 , fmt = * , iostat = istat) t_sum , amp
                    if (istat /= 0) stop "error writing test file"
                    call signal_selector(t_sum , length , amp , h_vector_in , k , pulse_selector)
                end do
            case('none')                                    !case for no window required
                bell_width = 0.32
                amp = 1.0_dp
                do k = 0 , (N_in - 1) , 1
                    t_sum = del_t_in*real(k , kind = dp)
                    call signal_selector(t_sum , length , amp , h_vector_in , k , pulse_selector)
                end do
            case default                                    !run ends if no valid case is found
                STOP 'invalid window type'
        end select
        close (unit = file_unit2 , iostat = istat)
        if (istat /= 0) stop "error closing test file"
    end subroutine windowing_subroutine

!--------------------------------------------------------------------------------------!
!signal_selector is called to calculate signal at N different values with intervals of !
!delta t in windowing_subroutine.                                                      !
!signal_selector takes in t_sum_in, amp_in, k_in as constants which vary with k in its !
!do loop (see windowing_subroutine). signal_vector is what the calculated values are   !
!written to (hence INOUT). signal_type is string used in case statements and determine !
!which signal values are written to signal_vector.                                     !
!--------------------------------------------------------------------------------------!

    subroutine signal_selector(t_sum_in , L , amp_in , signal_vector , k_in , signal_type)
        implicit none
        real(kind = dp) , intent(IN)                     :: t_sum_in , amp_in , L
        complex(kind = dp) , allocatable , intent(INOUT) :: signal_vector(:)
        character(len = 20) , intent(IN)                 :: signal_type
        integer , intent(IN)                             :: k_in
        real(kind = dp)                                  :: sigma

        sigma = 0.23_dp

        select case(signal_type)
            case('sine_signal')                             !this is the perfect sine wave signal initially used
                signal_vector(k_in + 1) = cmplx(amp_in*SIN(2.0_dp*pi*t_sum_in/period_1) , 0.0_dp , kind = dp)
            case('rectangle_signal')                        !case for rectangular signal
                if((k_in + 1) <= L) then                    !defining none-zero part of rectangular signal
                    signal_vector(k_in + 1) = cmplx(amp_in , 0.0_dp , kind = dp)
                else                                        !defining null part of rectangular signal
                    signal_vector(k_in + 1) = cmplx(0.0_dp , 0.0_dp , kind = dp)
                end if
            case('triangle_signal')                         !case statement for triangular signal
                if((k_in + 1) <= L/2) then                  !defining lhs of triangular pulse
                    signal_vector(k_in + 1) = cmplx((2.0_dp*amp_in)/L , 0.0_dp , kind = dp)
                else if((k_in + 1) <= L) then               !defining rhs of triangular pulse
                    signal_vector(k_in + 1) = cmplx((-2.0_dp*amp_in)/L + 2.0_dp*amp_in , 0.0_dp , kind = dp)
                else                                        !defining null part of triangular pulse
                    signal_vector(k_in + 1) = cmplx(0.0_dp , 0.0_dp , kind = dp)
                end if
            case('gaussian_signal')                         !case statement for Gaussian pulse signal
                signal_vector(k_in + 1) = &                 !Gaussian distribution equation- please see Elog.
                &cmplx(amp_in*EXP(-0.5_dp*((k_in - L*0.5_dp)/(sigma*0.5_dp*L))**2.0_dp) , 0.0_dp , kind = dp)
            case default                                    !Kills program if there is an invalid signal case.
                STOP 'invalid signal case'
        end select
    end subroutine signal_selector
end program FT
