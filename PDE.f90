!--------------------------------------------------------------------------------------!
!Modules set up to define constants and algorithms used in previous experiments, in    !
!particular the random number generator.                                               !
!--------------------------------------------------------------------------------------!

module global
    implicit none
    save
    integer , parameter :: dp = selected_real_kind(15 , 30)
    real(kind = dp) :: pi = 4.0_dp*ATAN(1.0_dp)
    contains

!--------------------------------------------------------------------------------------!
!Defined the identity matrix in global module for simplicity.                          !
!Takes the 'number' to define the size of the identity matrix. This means that the sub-!
!routine can be used for the calculation of any identity matrix                        !
!--------------------------------------------------------------------------------------!

    subroutine identity(number , I_matrix)
        implicit none
        integer , intent(IN) :: number
        real(kind = dp) , allocatable , intent(OUT) :: I_matrix(: , :)
        integer :: k = 0
        allocate(I_matrix(number , number))
        I_matrix = 0.0_dp
        do k = 1 , number
            I_matrix(k , k) = 1.0_dp
        end do
    end subroutine
end module global

!--------------------------------------------------------------------------------------!
!Module invert uses code written and provided by Prof. Matt Probert and calls LAPACK   !
!subroutines from the openblas library.                                                !
!--------------------------------------------------------------------------------------!

module invert
    use global
    implicit none
    save
    contains
    subroutine invert_matrix(matrix)
        ! Invert the supplied matrix using LAPACK
        implicit none
        real(kind=dp), dimension(:,:), intent(inout) :: matrix
        integer :: N, LWORK, IERR
        integer, dimension(:), allocatable :: IPIV
        real(kind=dp), dimension(:), allocatable :: WORK
        if (size(matrix,1) /= size(matrix,2)) STOP "Matrix is not square"
        N = size(matrix,1)
        allocate(IPIV(N),stat=IERR)
        if (IERR/=0) STOP "Failed to allocate IPIV"
        LWORK = N**2
        allocate(WORK(LWORK),stat=IERR)
        if (IERR/=0) STOP "Failed to allocate WORK"
        call dgetrf(N,N,matrix,N,IPIV,IERR)
        if (IERR/=0) STOP "Error in dgetrf: Matrix is singular"
        call dgetri(N,matrix,N,IPIV,WORK,LWORK,IERR)
        if (IERR/=0) STOP "Error in dgetri: Matrix is singular"
    end subroutine
end module invert

module random
    use global
    implicit none
    save
    contains
    subroutine rng(z_in , z_out)
        implicit none
        real(kind = dp) , intent(IN) :: z_in
        real(kind = dp) , intent(OUT) :: z_out
        real(kind = dp) :: A = 137.0_dp , B = 15887.0_dp , M = 714025.0_dp
        z_out = mod(z_in*A + B , M)
    end subroutine rng
end module random

program PDE
    use global
    use random
    use invert
    implicit none
    real(kind = dp) :: del_t_in = 0.1 , alpha_in = 10.0_dp**(-6.0_dp) , length_in = 1.0_dp , diff_f , diff_b
    integer :: N_in = 100 , file_unit1 = 20 , file_unit2 = 30 , t_max_in = 100000 , big_t , file_unit3 = 40 , istat , &
    &global_case , file_unit4 = 50
    character(len = 20) :: filename1 , filename2 , filename3 , boundary_case1 , boundary_case2 , filename4
    logical :: llexist
    filename1 = 'FTCS.dat'
    filename2 = 'BTCS.dat'
    filename3 = 'heatdiff.dat'
    filename4 = 'lambdadiff.dat'
    boundary_case1 = 'periodic'
    boundary_case2 = 'dirichlet'
    global_case = 3

!--------------------------------------------------------------------------------------!
!inquire tests to ensure that the file can be written and overwrites the old file of   !
!the same name in the event that it can't be.                                          !
!Logical llexist is used for this, this same algorithm is also used within the         !
!energy_calc subroutine.                                                               !
!Each case corresponds to a different boundary condition, and whether timestep/ lambda !
!is being varied over iterations.                                                      !
!--------------------------------------------------------------------------------------!

    select case(global_case)
    case(1)
        inquire(file=filename3, exist=llexist)
        if (llexist) then
            !print*,("old data file overwritten")!checking to ensure that no file is ballooning in size
            open (unit = file_unit3 , file = filename3 , status = "replace" , action = "write" , position = "append" ,&
            &iostat = istat)
            if (istat /= 0) stop "error opening test file"
        else
            open (unit = file_unit3 , file = filename3 , status = "new" , action = "write" , position = "append" ,&
            &iostat = istat)
            if (istat /= 0) stop "error opening test file"
        end if

!--------------------------------------------------------------------------------------!
!do loop iterates big_t over the required range, big_t is then converted to the real   !
!del_t which is used within the heat_calc subroutine as the timestep.                  !
!diff_f, diff_b are intent(OUT) diff from the heat_calc subroutine, calculates the     !
!difference between the initial and final sums of the heat scalars in x_steps matrix   !
!--------------------------------------------------------------------------------------!

        do big_t = 10 , 100
            del_t_in = real(big_t , kind = dp)/100.0_dp
            call heat_calc(alpha_in , del_t_in , length_in , N_in , t_max_in , file_unit1 , filename1 , .false. , &
            &0.0_dp , 0.0_dp , boundary_case1 , diff_f)
            call heat_calc(alpha_in , del_t_in , length_in , N_in , t_max_in , file_unit2 , filename2 , .true. , &
            &0.0_dp , 0.0_dp , boundary_case1 , diff_b)
            write (unit = file_unit3 , fmt = * , iostat = istat) del_t_in , diff_f , diff_b
        end do
        close (unit = file_unit3 , iostat = istat)
    case(2)
        call heat_calc(alpha_in , del_t_in , length_in , N_in , t_max_in , file_unit1 , filename1 , .false. , &
        &0.0_dp , 0.0_dp , boundary_case1 , diff_f)
        call heat_calc(alpha_in , del_t_in , length_in , N_in , t_max_in , file_unit2 , filename2 , .true. , &
        &0.0_dp , 0.0_dp , boundary_case1 , diff_b)
    case(3)
        call heat_calc(alpha_in , del_t_in , length_in , N_in , t_max_in , file_unit1 , filename1 , .false. , &
        &0.0_dp , 0.0_dp , boundary_case2 , diff_f)
        call heat_calc(alpha_in , del_t_in , length_in , N_in , t_max_in , file_unit2 , filename2 , .true. , &
        &0.0_dp , 0.0_dp , boundary_case2 , diff_b)
        print*, del_t_in , alpha_in , N_in , del_t_in*N_in
    case(4)
        inquire(file=filename4, exist=llexist)
        if (llexist) then
            !print*,("old data file overwritten")!checking to ensure that no file is ballooning in size
            open (unit = file_unit4 , file = filename4 , status = "replace" , action = "write" , position = "append" ,&
            &iostat = istat)
            if (istat /= 0) stop "error opening test file"
        else
            open (unit = file_unit4 , file = filename4 , status = "new" , action = "write" , position = "append" ,&
            &iostat = istat)
            if (istat /= 0) stop "error opening test file"
        end if
        do big_t = 100 , 550
            print*, big_t
            del_t_in = real(big_t , kind = dp)/1000.0_dp
            call heat_calc(alpha_in , del_t_in , length_in , N_in , t_max_in , file_unit1 , filename1 , .false. , &
            &0.0_dp , 100.0_dp , boundary_case2 , diff_f)
            call heat_calc(alpha_in , del_t_in , length_in , N_in , t_max_in , file_unit2 , filename2 , .true. , &
            &0.0_dp , 100.0_dp , boundary_case2 , diff_b)
            write (unit = file_unit4 , fmt = * , iostat = istat) del_t_in , diff_f , diff_b
        end do
        close (unit = file_unit4 , iostat = istat)
    case default
        STOP 'invalid'
    end select

    contains
    subroutine heat_calc(alpha , del_t , length , N , t_max , file_unit , filename , switch , V_1 , V_2 , &
    &boundary_type , diff)
        implicit none
        real(kind = dp) , intent(IN) :: alpha , del_t , length , V_1 , V_2
        real(kind = dp) , intent(OUT) :: diff
        integer , intent(IN) :: N , file_unit , t_max
        character(len = 20) , intent(IN) :: filename , boundary_type
        logical , intent(IN) :: switch
        real(kind = dp) :: rand_old = 100.0_dp/714025.0_dp , rand_new , t_sum = 0.0_dp , lambda , init_sum
        real(kind = dp) , allocatable :: x_steps(: , :) , lambda_matrix(: , :) , x_steps_new(: , :) , identity_matrix(: , :) , &
        &v_array(: , :)
        integer :: i = 0 , j = 0 , file_unit5 = 60 , k = 0
        logical :: lexist
        character(len = 20) :: filename5 = 'threedplt.csv'

!--------------------------------------------------------------------------------------!
!Allocating the allocatable arrays as defined above by using the constant N, this keeps!
!the code general.                                                                     !
!--------------------------------------------------------------------------------------!

        allocate(x_steps(N , 1))
        x_steps = 100.0_dp
        allocate(x_steps_new(N , 1))
        allocate(lambda_matrix(N , N))
        lambda_matrix = 0.0_dp
        allocate(identity_matrix(N , N))
        allocate(v_array(N , 1))
        v_array = 0.0_dp
        lambda = alpha*del_t*(real(N , kind = dp)/length)**2.0_dp
        if(alpha/(lambda*(length/real(N , kind = dp))) > 6363.0_dp) STOP 'physically implausible'

!--------------------------------------------------------------------------------------!
!inquire checks if the file already exists and overwrites it if one does.              !
!--------------------------------------------------------------------------------------!

        inquire(file=filename, exist=lexist)
        if (lexist) then
            !print*,("old data file overwritten")!checking to ensure that no file is ballooning in size
            open (unit = file_unit , file = filename , status = "replace" , action = "write" , position = "append" , iostat = istat)
            if (istat /= 0) stop "error opening test file"
        else
            open (unit = file_unit , file = filename , status = "new" , action = "write" , position = "append" , iostat = istat)
            if (istat /= 0) stop "error opening test file"
        end if

!--------------------------------------------------------------------------------------!
!Setting the lambda array and instantiating the x_steps array with initial, normalised !
!random variables.                                                                     !
!--------------------------------------------------------------------------------------!

        do i = 1 , N
            call rng(rand_old , rand_new)
            rand_old = rand_new
!            x_steps(i , 1) = rand_old/714025.0_dp
!            print*, rand_old/714025.0_dp
            lambda_matrix(i , i) = 1.0_dp - 2.0_dp*lambda
            if(i /= 1) then
                lambda_matrix(i - 1 , i) = lambda
            end if
            if(i /= N) then
                lambda_matrix(i + 1 , i) = lambda
            end if
        end do
        init_sum = sum(x_steps)

        select case(boundary_type)
        case('periodic')
            lambda_matrix(1 , N) = lambda
            lambda_matrix(N , 1) = lambda

!--------------------------------------------------------------------------------------!
!the switch if statement checks if the system should be running the FTCS or the BTCS   !
!algorithm --> if switch is true, then the code changes the lambda matrix to perform   !
!the BTCS algorithm.                                                                   !
!--------------------------------------------------------------------------------------!

            if(switch .eqv. .true.)then
                call identity(N , identity_matrix)
                lambda_matrix = 2.0_dp*identity_matrix - lambda_matrix
                call invert_matrix(lambda_matrix)
            end if

!--------------------------------------------------------------------------------------!
!Multiplying the matricies to get the subsequent timestep.                             !
!j do loop also writes current time and current sum to file to be plotted and checked. !
!--------------------------------------------------------------------------------------!

            do j = 1 , t_max
                x_steps_new = matmul(lambda_matrix , x_steps)
                x_steps = x_steps_new
                t_sum = j*del_t
                write (unit = file_unit , fmt = * , iostat = istat) t_sum , sum(x_steps)*(length/real(N , kind = dp))
            end do
            close (unit = file_unit , iostat = istat)
            diff = abs(sum(x_steps) - init_sum)

!--------------------------------------------------------------------------------------!
!Case select used to determine the boundary conditions.                                !
!V_1, V_2 are defined within v_array, then code is iterated depending on the type of   !
!algorithm used. This is defined using the switch boolian.                             !
!--------------------------------------------------------------------------------------!

        case('dirichlet')
            v_array(1 , 1) = V_1
            v_array(N , 1) = V_2
            inquire(file=filename5, exist=llexist)
            if (llexist) then
            !print*,("old data file overwritten")!checking to ensure that no file is ballooning in size
            open (unit = file_unit5 , file = filename5 , status = "replace" , action = "write" , position = "append" ,&
            &iostat = istat)
            if (istat /= 0) stop "error opening test file"
            else
            open (unit = file_unit5 , file = filename5 , status = "new" , action = "write" , position = "append" ,&
            &iostat = istat)
            if (istat /= 0) stop "error opening test file"
            end if
            if(switch .eqv. .true.)then
                call identity(N , identity_matrix)
                lambda_matrix = 2.0_dp*identity_matrix - lambda_matrix
                call invert_matrix(lambda_matrix)
                do j = 1 , t_max
                    x_steps_new = matmul(lambda_matrix , x_steps) + matmul(lambda_matrix , (lambda*v_array))
                    x_steps = x_steps_new
                    t_sum = j*del_t
                    write (unit = file_unit , fmt = * , iostat = istat) t_sum , sum(x_steps)*(length/real(N , kind = dp))
                    do k = 1 , N
                        write (unit = file_unit5 , fmt = * , iostat = istat) t_sum , ' , ' ,&
                        &real(k , kind = dp)/real(N , kind = dp) , ' , ' , x_steps(k , 1)
                    end do
                end do
            else
                do j = 1 , t_max
                    x_steps_new = matmul(lambda_matrix , x_steps) + lambda*v_array
                    x_steps = x_steps_new
                    t_sum = j*del_t
                    write (unit = file_unit , fmt = * , iostat = istat) t_sum , sum(x_steps)*(length/real(N , kind = dp))
                    do k = 1 , N
                        write (unit = file_unit5 , fmt = * , iostat = istat) t_sum , real(k , kind = dp)/real(N , kind = dp) , &
                        &x_steps(k , 1)
                    end do
                end do
            end if
            close (unit = file_unit , iostat = istat)
            close (unit = file_unit5 , iostat = istat)
            diff = abs(sum(x_steps) - init_sum)
        case default
            STOP 'invalid boundary type'
        end select
    end subroutine heat_calc
end program PDE
