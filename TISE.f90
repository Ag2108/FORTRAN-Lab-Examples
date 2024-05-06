program TISE
    implicit none

!--------------------------------------------------------------------------------------!
!Initially, I define the global parameters and vairables. This means that I do not need!
!to define them for each function and subroutine that I subsequently use.              !
!These global variables include the initial starting value (t_init), independent       !
!variable stepsize (del_t), the length of the system (L), and the potential of the well!
!of the system (V_o).                                                                  !
!--------------------------------------------------------------------------------------!

    integer , parameter :: dp = selected_real_kind(15 , 300)
    real(kind = dp) :: y_start = 10.0_dp**(-12.0_dp) , z_start = 10_dp**(-12.0_dp) , &
    &del_t = 0.001_dp , t_init = 5.1_dp , L = 8.15_dp , Eigen_RK4 , Eigen_Numerov , Etest_1 , Etest_2 , n = 1.0_dp!defining the universal parameters
    real(kind = dp) :: V, V_o
    real(kind = dp) :: pi = 4.0_dp*ATAN(1.0_dp)
    !character (len = 20) :: filename_two = "test2.dat"
    integer :: istat, m , level

!--------------------------------------------------------------------------------------!
!Here, m is the multiplier for which its product with del_t makes an integer. This     !
!means that I can use integer addition when computing the position of the algorithms in!
!their respective subroutines.                                                         !
!--------------------------------------------------------------------------------------!

    level = int(n)
    m = int(1/del_t)
    V_o = 10000.0_dp !initially using a large potential to simulate an infinite sqare well.

!--------------------------------------------------------------------------------------!
!Here, the program uses the Etest_1 and Etest_2 variables within the energy_calculator !
!function (within contains) as starting points for the shooting method.                !
!The .true. and .false. dummy variable inputs dictate to the delta_calculator function !
!(within energy_calculator) whether to use the RK4 or Numerov evaluation when computing!
!the energy eigenvalue of the system.                                                  !
!--------------------------------------------------------------------------------------!

    Etest_1 = (n*pi/L)**2.0_dp - 0.04_dp!defining the initial guesses for the energy eigenvalue of the system
    Etest_2 = (n*pi/L)**2.0_dp + 0.06_dp
    Eigen_RK4 = energy_calculator(Etest_1 , Etest_2 , .true.)
    Eigen_Numerov = energy_calculator(Etest_1 , Etest_2 , .false.)!calling the function defined below
    print*, Eigen_RK4 , Eigen_Numerov!printing the eigenvalue calculated for the system

    contains

!--------------------------------------------------------------------------------------!
!I used my RK4 subroutine (developed for the Duffing oscillator lab) with some light   !
!modifications. I added in filename and file_unit dummy variables to the subroutine so !
!that these could be varied from one iteration of the subroutine to the other. This    !
!allowed both the Numerov and RK4 evaluations of psi to be plotted independently.      !
!I also used end_grad and DLS_out intent(OUT) variables to calculate the difference in !
!gradient and the DLS of the system at 0, as required for the shooting method.         !
!--------------------------------------------------------------------------------------!

    subroutine RK4(y_init , z_init , t_init , del_t , filename , file_unit , E_n , t_final , end_grad , DLS_out)
        implicit none
        real(kind = dp) , intent(IN) :: y_init , z_init , del_t , t_init , E_n , t_final
        real(kind = dp) , intent(OUT) :: end_grad , DLS_out!definging what we want the subroutine to output once it has finished the other calculations
        real(kind = dp) :: y_n , z_n , k_1y , k_2y , k_1z , k_2z , t_sum , &
        &k_3y , k_3z , k_4y , k_4z
        integer :: big_ind , file_unit , big_final
        logical :: lexist
        character(len = 7) , intent(IN) :: filename
        t_sum = t_init

!--------------------------------------------------------------------------------------!
!big_ind is the integer counter recycled for both RK4 and Numerov methods.             !
!--------------------------------------------------------------------------------------!

        big_ind = int(t_init*m)
        big_final = int(t_final*m)
        y_n = y_init !setting the initial parameters
        z_n = z_init

!--------------------------------------------------------------------------------------!
!Uses the inquire with the logical variable lesxist to check for old files with the    !
!same name as filename. Replaces the file if there is, makes a new file if there is not!
!This prevents the data from multiple runs from being plotted on the same file.        !
!--------------------------------------------------------------------------------------!

        inquire(file=filename, exist=lexist)
        if (lexist) then
            !print*,("old data file overwritten")!checking to ensure that no file is ballooning in size
            open (unit = file_unit , file = filename , status = "replace" , action = "write" , position = "append" , iostat = istat)
            if (istat /= 0) stop "error opening test file"!opening file for RK4 analysis
        else
            open (unit = file_unit , file = filename , status = "new" , action = "write" , position = "append" , iostat = istat)
            if (istat /= 0) stop "error opening test file"!opening file for RK4 analysis
        end if
        write (unit = file_unit , fmt = * , iostat = istat) t_sum , y_n , z_n , z_n/y_n!writing the initial 0th and 1st derivatives of the function to the data file

!--------------------------------------------------------------------------------------!
!Uses big_ind integer to check when the code has hit 0.                                !
!Until then, the code iterates over the k_yval and k_zval functions to evaluate the new!
!values of psi and the position for each iteration using the RK4 method.               !
!--------------------------------------------------------------------------------------!

        do while(big_ind /= big_final)!stops when it hits 0
            k_1y = k_yval(z_n , 0.0_dp)*del_t!k values are calculated using previous k values and the values at that iteration of the code
            k_1z = k_zval(y_n , t_sum , E_n)*del_t
            k_2y = k_yval(z_n , k_1z)*del_t
            k_2z = k_zval(y_n + (k_1y/2.0_dp) , t_sum + (del_t/2.0_dp) , E_n)*del_t
            k_3y = k_yval(z_n , k_2z)*del_t
            k_3z = k_zval(y_n + (k_2y/2.0_dp) , t_sum + (del_t/2.0_dp) , E_n)*del_t
            k_4y = k_yval(z_n , 2.0_dp*k_3z)*del_t
            k_4z = k_zval(y_n + k_3y , t_sum + del_t , E_n)*del_t
            z_n = z_n + ((k_1z + k_4z)/6.0_dp) + ((k_2z + k_3z)/3.0_dp)
            y_n = y_n + ((k_1y + k_4y)/6.0_dp) + ((k_2y + k_3y)/3.0_dp)

!--------------------------------------------------------------------------------------!
!Checks which way the code should be iterating, and moves it closer to 0. This ensures !
!that the do while() loop above us fulfilled.                                          !
!t_sum (real used in calculations) is then updated from big_ind using m.               !
!--------------------------------------------------------------------------------------!

            if(big_ind <= big_final) then !if statement ensures that both iterations converge towards 0
                big_ind = big_ind + 1
            else
                big_ind = big_ind - 1
            end if
            t_sum = real(big_ind , kind = dp)/real(m , kind = dp) !counting the sum of time using integers
            write (unit = file_unit , fmt = * , iostat = istat) t_sum , y_n , z_n , z_n/y_n! , COS(omega*t_sum)!values are written to RK4 analysis file
        end do
        close (unit = file_unit , iostat = istat)
        if (istat /= 0) stop "error closing test file"
        end_grad = z_n
        DLS_out = z_n/y_n
    end subroutine RK4

!--------------------------------------------------------------------------------------!
!k_zval function checks the position of the iteration using the independent variable   !
!input, t_in. The potential well is defined here.                                      !
!The values from k_val are used in both the Numerov and RK4 methods.                   !
!--------------------------------------------------------------------------------------!

    real(kind = dp) function k_zval(y_in , t_in , E_in)!function calculating the intermediary k values for the RK4 algorithm
        implicit none
        real(kind = dp) , intent(IN) :: y_in , t_in , E_in
        if(abs(t_in) >= (L/2.0_dp)) then!defining the discontinuous potential function
            V = V_o
        else
            V = 0.0_dp
        end if
        k_zval = -1.0_dp*(E_in - V)*y_in !constants are defined within the program
    end function k_zval

!--------------------------------------------------------------------------------------!
!Defines k_yval as used in RK4                                                         !
!--------------------------------------------------------------------------------------!

    real(kind = dp) function k_yval(z_in , k_in)
        implicit none
        real(kind = dp) , intent(IN) :: z_in , k_in
        k_yval = z_in + (k_in/2.0_dp)
    end function k_yval

!--------------------------------------------------------------------------------------!
!delta_calculator function uses the present boundaries to evaluate a single step of the!
!shooting method. The logical dummy variable switch determines if the RK4 algorithm    !
!(when .true.) or the Numerov method (otherwise) are used for the shooting method here.!
!Outputs are output using an array.                                                    !
!--------------------------------------------------------------------------------------!

    function delta_calculator(Etest_left , Etest_right , switch)
        implicit none
        real(kind = dp) , dimension(2) :: delta_calculator
        real(kind = dp) , intent(IN) :: Etest_left , Etest_right
        logical , intent(IN) :: switch
        real(kind = dp) :: grad_left1 , grad_right1 , grad_left2 , grad_right2 , DLS_left1 , DLS_left2 , DLS_right1 ,&
        &DLS_right2

!--------------------------------------------------------------------------------------!
!Calls RK4 for each side of the origin for both energies. Uses the intent(OUT)         !
!variables to evaluate the DSI and graient difference for those 2 energies.            !
!Checks for even n by using level variable. This changes how the eigenvalue is found   !
!--------------------------------------------------------------------------------------!

        if(mod(level , 2) == 0)then
            if(switch .eqv. .true.) then
                call RK4(y_start , z_start , -1.0_dp*t_init , del_t , "lhs.dat" , 30 , Etest_left , -L/(2.0_dp*n) , grad_left1 , &
                &DLS_left1)
                call RK4(-1.0_dp*y_start , z_start , t_init , -1.0_dp*del_t , "rhs.dat" , 40 , Etest_left , -L/(2.0_dp*n) ,&
                &grad_right1 , DLS_left2)
                call RK4(y_start , z_start , -1.0_dp*t_init , del_t , "lhs.dat" , 30 , Etest_right , -L/(2.0_dp*n) ,&
                & grad_left2 , DLS_right1)
                call RK4(-1.0_dp*y_start , z_start , t_init , -1.0_dp*del_t , "rhs.dat" , 40 , Etest_right , -L/(2.0_dp*n) ,&
                &grad_right2 , DLS_right2)
                print*, "RK4"
            else
                call Numerov(y_start , z_start , -1.0_dp*t_init , del_t , "lhn.dat" , 20 , Etest_left , -L/(2.0_dp*n) ,&
                &grad_left1 , DLS_left1)
                call Numerov(-1.0_dp*y_start , z_start , t_init , -1.0_dp*del_t , "rhn_n.dat" , 50 , Etest_left , -L/(2.0_dp*n) ,&
                &grad_right1 , DLS_left2)
                call Numerov(y_start , z_start , -1.0_dp*t_init , del_t , "lhn.dat" , 20 , Etest_right , -L/(2.0_dp*n) ,&
                &grad_left2 , DLS_right1)
                call Numerov(-1.0_dp*y_start , z_start , t_init , -1.0_dp*del_t , "rhn.dat" , 50 , Etest_right , -L/(2.0_dp*n) ,&
                &grad_right2 , DLS_right2)
                print*, "Numerov"
            end if
        else
            if(switch .eqv. .true.) then
                call RK4(y_start , z_start , -1.0_dp*t_init , del_t , "lhs.dat" , 30 , Etest_left , 0.0_dp , grad_left1 , &
                &DLS_left1)
                call RK4(y_start , -1.0_dp*z_start , t_init , -1.0_dp*del_t , "rhs.dat" , 40 , Etest_left , 0.0_dp , grad_right1&
                &, DLS_left2)
                call RK4(y_start , z_start , -1.0_dp*t_init , del_t , "lhs.dat" , 30 , Etest_right , 0.0_dp , grad_left2 , &
                &DLS_right1)
                call RK4(y_start , -1.0_dp*z_start , t_init , -1.0_dp*del_t , "rhs.dat" , 40 , Etest_right , 0.0_dp , grad_right2 ,&
                &DLS_right2)
                print*, "RK4"
            else
                call Numerov(y_start , z_start , -1.0_dp*t_init , del_t , "lhn.dat" , 20 , Etest_left , 0.0_dp , grad_left1 ,&
                &DLS_left1)
                call Numerov(y_start , -1.0_dp*z_start , t_init , -1.0_dp*del_t , "rhn_n.dat" , 50 , Etest_left , 0.0_dp ,&
                &grad_right1 , DLS_left2)
                call Numerov(y_start , z_start , -1.0_dp*t_init , del_t , "lhn.dat" , 20 , Etest_right , 0.0_dp , grad_left2 ,&
                &DLS_right1)
                call Numerov(y_start , -1.0_dp*z_start , t_init , -1.0_dp*del_t , "rhn.dat" , 50 , Etest_right , 0.0_dp ,&
                &grad_right2 , DLS_right2)
                print*, "Numerov"
            end if
        end if
        delta_calculator(1) = (grad_left1 - grad_right1)*(grad_left2 - grad_right2)
        delta_calculator(2) = DLS_left1 - DLS_left2
    end function delta_calculator

!--------------------------------------------------------------------------------------!
!Shooting method is iterated within energy_calculator function.                        !
!Given the use of the switch logical, both Numerov and RK4 methods use the same        !
!shooting algorithm.                                                                   !
!--------------------------------------------------------------------------------------!

    real(kind = dp) function energy_calculator(E_init1 , E_init2 , switch)
        implicit none
        real(kind = dp) , intent(IN) :: E_init1 , E_init2
        real(kind = dp) :: E_mid , E_1 , E_2
        real(kind = dp) , dimension(2) :: delta_prod
        logical , intent(IN) :: switch
        E_1 = E_init1
        E_2 = E_init2
        delta_prod = delta_calculator(E_1 , E_2 , switch)

!--------------------------------------------------------------------------------------!
!Shooting algorithm uses the delta_prod array output from the delta_calculator function!
!to determine which side of the midpoint it will find the solution and when the        !
!solution has been found.                                                              !
!--------------------------------------------------------------------------------------!

        if(delta_prod(1) > 0.0_dp) then
            stop "energy is not in range given, choose new energy boundaries."
        else
            do while(abs(delta_prod(2)) > 10.0_dp**(-4.0_dp))
                E_mid = 0.5_dp*(E_1 + E_2)
                delta_prod = delta_calculator(E_mid , E_2 , switch)
                if(delta_prod(1) < 0.0_dp) then
                    E_1 = E_mid
                else
                    E_2 = E_mid
                end if
                delta_prod = delta_calculator(E_1 , E_2 , switch)
            end do
        energy_calculator = E_mid!outputting the energy eigenvalue calculated by the function
        end if
    end function energy_calculator

!--------------------------------------------------------------------------------------!
!Numerov method is defined in the subroutine below.                                    !
!Many of the same routines are shared between the RK4 and Numerov subroutine.          !
!The same dummy variables are input into the Numerov subroutine as the RK4 subroutine. !
!RK4 functions are used to initialise G_n using the RK4 method.                        !
!--------------------------------------------------------------------------------------!

    subroutine Numerov(y_init , z_init , t_init , del_t , filename , file_unit , E_n , t_final , end_grad , DLS_out)
        implicit none
        real(kind = dp) , intent(IN) :: y_init , z_init , del_t , t_init , E_n , t_final
        real(kind = dp) , intent(OUT) :: end_grad , DLS_out!definging what we want the subroutine to output once it has finished the other calculations
        real(kind = dp) :: y_n , z_n , k_1y , k_2y , k_1z , k_2z , t_sum , &
        &k_3y , k_3z , k_4y , g_old , g_n , g_new , y_old
        integer :: big_ind , file_unit , big_final
        character(len = 7) , intent(IN) :: filename
        logical :: lexist
        t_sum = t_init
        big_ind = int(t_init*m)
        big_final = int(t_final*m)
        y_n = y_init !setting the initial parameters
        z_n = z_init
        g_old = G_function(y_n , del_t , t_sum , E_n)
        inquire(file=filename, exist=lexist)!checking if the file already exists and overriding it if it does
        if (lexist) then
!            print*,("old data file overwritten")!checking to ensure that no file is ballooning in size
            open (unit = file_unit , file = filename , status = "replace" , action = "write" , position = "append" , iostat = istat)
            if (istat /= 0) stop "error opening test file"!opening file for RK4 analysis
        else
            open (unit = file_unit , file = filename , status = "new" , action = "write" , position = "append" , iostat = istat)
            if (istat /= 0) stop "error opening test file"!opening file for RK4 analysis
        end if
        write (unit = file_unit , fmt = * , iostat = istat) t_sum , y_n , z_n , z_n/y_n!this part is simply copied from the above
        k_1y = k_yval(z_n , 0.0_dp)*del_t!Using RK4 to obtain the first iteration of the Numerov method
        k_1z = k_zval(y_n , t_sum , E_n)*del_t
        k_2y = k_yval(z_n , k_1z)*del_t
        k_2z = k_zval(y_n + (k_1y/2.0_dp) , t_sum + (del_t/2.0_dp) , E_n)*del_t
        k_3y = k_yval(z_n , k_2z)*del_t
        k_3z = k_zval(y_n + (k_2y/2.0_dp) , t_sum + (del_t/2.0_dp) , E_n)*del_t
        k_4y = k_yval(z_n , 2.0_dp*k_3z)*del_t
        y_n = y_n + ((k_1y + k_4y)/6.0_dp) + ((k_2y + k_3y)/3.0_dp)
        if(big_ind <= big_final) then !if statement ensures that both iterations converge towards 0
            big_ind = big_ind + 1
        else
            big_ind = big_ind - 1
        end if
        t_sum = real(big_ind , kind = dp)/real(m , kind = dp) !counting the sum of time using integers
        write (unit = file_unit , fmt = * , iostat = istat) t_sum , y_n , z_n , z_n/y_n

!--------------------------------------------------------------------------------------!
!G_function is used to initially evaluate g_n from the initial conditions and y_n from !
!RK4.                                                                                  !
!g_n is then iterated using the formula for the Numerov method and y_n is subsequently !
!obtained from the Numerov formula function to be plotted as psi at that point.        !
!The gradient is calculated by calculating the final change in y_n over the last step. !
!--------------------------------------------------------------------------------------!

        g_n = G_function(y_n , del_t , t_sum , E_n)
        do while(big_ind /= big_final)
            g_new = 2.0_dp*g_n - g_old + del_t**2.0_dp*k_zval(y_n , t_sum , E_n)
            g_old = g_n
            g_n = g_new
            y_old = y_n
            if(abs(t_sum) >= (L/2.0_dp)) then!defining the discontinuous potential function
                y_n = ((1.0_dp - ((del_t**2)/12.0_dp)*(E_n - V_o))**(-1.0_dp))*g_n
            else
                y_n = ((1.0_dp - ((del_t**2)/12.0_dp)*(E_n))**(-1.0_dp))*g_n
            end if
            if(big_ind <= big_final) then
                big_ind = big_ind + 1
            else
                big_ind = big_ind - 1
            end if
            t_sum = real(big_ind , kind = dp)/real(m , kind = dp) !counting the sum of time using integers
            write (unit = file_unit , fmt = * , iostat = istat) t_sum , y_n , z_n , z_n/y_n
        end do
        close (unit = file_unit , iostat = istat)
        if (istat /= 0) stop "error closing test file"
        end_grad = (y_n - y_old)/del_t
        DLS_out = end_grad/y_n
    end subroutine Numerov

!--------------------------------------------------------------------------------------!
!G_function simply calculates g_n from the initial condition and y_n from RK4.         !
!--------------------------------------------------------------------------------------!

    real(kind = dp) function G_function(y_in , del_t , t_in , E_in)
        implicit none
        real(kind = dp) , intent(IN) :: y_in , del_t , t_in , E_in
        G_function = (y_in + ((del_t**2.0_dp)/12.0_dp)*k_zval(y_in , t_in , E_in))
    end function G_function

end program TISE
