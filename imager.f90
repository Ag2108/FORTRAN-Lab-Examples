module global
    implicit none
    save
    integer , parameter :: dp = selected_real_kind(15 , 30)
    real(kind = dp)     :: pi = 4.0_dp*ATAN(1.0_dp)
end module global

module fft_module
    use global
    implicit none
    save

    integer , parameter :: i64 = selected_int_kind(18)

    integer , parameter :: fftw_forward = -1
    integer , parameter :: fftw_backward = +1
    integer , parameter :: fftw_destroy_input = 1
    integer , parameter :: fftw_unaligned = 2
    integer , parameter :: fftw_conserve_memory = 4
    integer , parameter :: fftw_exhaustive = 8
    integer , parameter :: fftw_preserve_input = 16
    integer , parameter :: fftw_patient = 32
    integer , parameter :: fftw_estimate = 64

end module fft_module

!--------------------------------------------------------------------------------------!
!Please note that I have heavily ammended my PBM read/write code such that it includes !
!error messages and 17-long lines for less verbose files.                              !
!The code for these ammendments was provided by Prof. Matt Probert and accessed on     !
!13/04/24.                                                                             !
!--------------------------------------------------------------------------------------!

module pbm_module
    use global
    use fft_module
    implicit none
    save
    contains

!--------------------------------------------------------------------------------------!
!pbm_readwrite subroutine takes in an array of data for the pbm data to write to or to !
!write to pbm, file_unit to ensure that different file unit is used for read and       !
!written pbm, max_greys for infor for pbm writing, pbm_case is an integer used to defi-!
!ne if the subroutine is reading or writing to a pbm file and filename_mod to define   !
!if the name of the file being manipulated.                                            !
!--------------------------------------------------------------------------------------!

    subroutine pbm_readwrite(pixels_in , file_unit , max_greys , pbm_case , filename_mod)

    integer , intent(IN)                  :: file_unit , pbm_case
    integer , intent(INOUT)               :: max_greys
    character(len = 20) , intent(IN)      :: filename_mod
    integer , allocatable , intent(INOUT) :: pixels_in(: , :)
    integer                               :: i , j , k , x , y , Nx , Ny
    integer                               :: istat , pixels_max , pixels_min
    integer , allocatable                 :: image_array(: , :)
    character(len = 20)                   :: magic_number

!--------------------------------------------------------------------------------------!
!case(1) makes the subroutine read a pbm file to an array, pixels_in.                  !
!Checks if the correct magic number, valid length and valid max greys is in the file   !
!being read.                                                                           !
!--------------------------------------------------------------------------------------!

    select case(pbm_case)
        case(1)
            open(file = filename_mod , unit = file_unit , status = 'old')
            read(file_unit , *) magic_number
            if (magic_number /= 'P2') STOP 'file in invalid format'                 !ensures that 17-long lines are valid

            read(file_unit , *) Nx , Ny
            if(Nx < 1 .or. Ny < 1) STOP 'invalid grid size'                         !protects against invalid runsizes

            allocate(pixels_in(1:Nx , 1:Ny) , stat = istat)
            if(istat /= 0) STOP 'unable to allocate pixels array'

            read(file_unit , *) max_greys
            if(max_greys <= 0 .or. max_greys > 256) STOP 'invalid max_greys'        !prevents invalid max_greys

!--------------------------------------------------------------------------------------!
!do loop here ensures that line length written to pbm is 17, this makes the file less  !
!verbose.                                                                              !
!--------------------------------------------------------------------------------------!

            do j = 1 , Ny
                do i = 1 , Nx - 17 , 17                                             !do loop changed from previous to allow for length 17 lines.
                    read(file_unit , *) (pixels_in(i + k - 1 , j) , k = 1 , 17)
                end do
                read(file_unit , *) (pixels_in(k , j) , k = i , Nx)
            end do

            close(unit = file_unit)
            if (istat /= 0) stop 'error closing file'

!--------------------------------------------------------------------------------------!
!case 2 writes pixels_in to a pbm file, this code was taken from my DLA experiment and !
!modified for improvement as explained above.                                          !
!--------------------------------------------------------------------------------------!

        case(2)
            Nx = size(pixels_in , dim = 1)
            Ny = size(pixels_in , dim = 2)
            allocate(image_array(1:Nx , 1:Ny) , stat = istat)
            if(istat /= 0) STOP 'error appending image_array'

            pixels_max = maxval(pixels_in)
            pixels_min = minval(pixels_in)

            y=lbound(pixels_in , dim = 2)
            do j=1 , Ny
                x=lbound(pixels_in , dim = 1)
                do i=1 , Nx
                    image_array(i,j)=int((pixels_in(x,y)-pixels_min)*max_greys/pixels_max)!want pixels: min->max mapped onto image_array: 0->max_greys
                    x=x+1
                end do
                y=y+1
            end do

            open(file = filename_mod , unit = file_unit , status = 'unknown' , iostat = istat)
            if (istat /= 0) stop 'error opening file'
            write(file_unit , 11) 'P2'
            write(file_unit , 12) Nx , Ny
            write(file_unit , 13) max_greys

            do j = 1 , Ny
                do i = 1 , Nx - 17 , 17                                             !do loop changed from previous to allow for length 17 lines.
                    write(file_unit , *) (image_array(i + k - 1 , j) , k = 1 , 17)
                end do
                write(file_unit , *) (image_array(k , j) , k = i , Nx)
            end do

            close(unit = file_unit)
            if (istat /= 0) stop 'error closing file'

            11 format(a2)
            12 format(i3 , 1x , i3)
            13 format(i5)

        case default
            STOP 'invalid case'                                                     !only valid cases are considered
        end select
    end subroutine pbm_readwrite
end module pbm_module

program imaging
    use global
    use pbm_module
    use fft_module
    implicit none
    character(len = 20) :: filename1 = 'clown.pgm'                                  !defining variables for subroutines
    character(len = 20) :: filename2 = 'blursign100.pgm'
    character(len = 20) :: filename3 = 'blurred.pgm'
    character(len = 20) :: edge_type = 'fft'
    character(len = 20) :: fft_case  = 'blur'
    real (kind = dp)    :: radius    = 30.0_dp , ampli = 100.0_dp , length = 21.0_dp

    select case(edge_type)
        case('sobel')
            call sobel_subroutine(filename1 , filename2)                             !calling subroutine from below
        case('fft')
            call fft_subroutine(filename3 , filename2 , radius , length , ampli , fft_case)
        case default
            STOP 'base case failure'
    end select

!--------------------------------------------------------------------------------------!
!sobel_subroutine is the main subroutine for doing Sobel edge detection on the image.  !
!filename1_in, filename2_in define the names of read and write files respectively.     !
!pixels_array (taken from the pgm) and gradient_array (passed to the pgm) are defined  !
!and allocated here.                                                                   !
!The result from the sobel operation is appended to the gradient_array.                !
!--------------------------------------------------------------------------------------!

    contains
    subroutine sobel_subroutine(filename1_in , filename2_in)
        implicit none
        character(len = 20) , intent(IN) :: filename1_in , filename2_in
        integer                          :: file_unit1 = 20 , file_unit2 = 30 , i = 0 , j = 0 , max_greys_out
        integer                          :: sobel_sumx = 0 , sobel_sumy = 0 , size_x , size_y , istat
        integer , allocatable            :: pixels_array(: , :) , gradient_array(: , :)

        call pbm_readwrite(pixels_array , file_unit1 , max_greys_out , 1 , filename1_in)

        size_x = size(pixels_array , dim = 1)                                       !defining pixels to be the correct size
        size_y = size(pixels_array , dim = 2)
        allocate(gradient_array(1:size_x , 1:size_y) , stat = istat)
        if(istat /= 0) STOP 'error appending gradient_array'
        gradient_array = 0

        do i = 2 , (size_x - 1) , 1
            do j = 2 , (size_y - 1) , 1
                call x_sobel(pixels_array , sobel_sumx , i , j)                     !running the sobel operator over pixels
                call y_sobel(pixels_array , sobel_sumy , i , j)
                gradient_array(i , j) = int((real(sobel_sumx , kind = dp)**2.0_dp + &
                &real(sobel_sumy , kind = dp)**2.0_dp)**0.5_dp)                     !calculating and appending gradient
            end do
        end do

        call pbm_readwrite(gradient_array , file_unit2 , max_greys_out , 2 , filename2_in)

    end subroutine sobel_subroutine

!--------------------------------------------------------------------------------------!
!x_sobel, y_sobel are defined in different subroutines but operated within the same do !
!loop wihtin sobel_subroutine.                                                         !
!x_sobel, y_sobel share the same input and output given that they operate on the same  !
!do loop and differ only in how they manipulate the data.                              !
!pixels_array_in defines the array of data being manipulated                           !
!sobel_sum is the output from the subroutine and is used to calculate the gradient.    !
!i_in, j_in are used to define which part of the pixels array is beinf manipulated.    !
!--------------------------------------------------------------------------------------!

    subroutine x_sobel(pixels_array_in , sobel_sum , i_in , j_in)
        implicit none
        integer , intent(IN)   :: i_in , j_in
        integer , intent(IN)   :: pixels_array_in(: , :)
        integer , dimension(6) :: sobel_array = 0
        integer , intent(OUT)  :: sobel_sum

        sobel_array(1) = -1*pixels_array_in(i_in - 1 , j_in + 1)                    !multiplies according to GX
        sobel_array(2) = 1*pixels_array_in(i_in + 1 , j_in + 1)
        sobel_array(3) = -2*pixels_array_in(i_in - 1 , j_in)
        sobel_array(4) = 2*pixels_array_in(i_in + 1 , j_in)
        sobel_array(5) = -1*pixels_array_in(i_in - 1 , j_in - 1)
        sobel_array(6) = 1*pixels_array_in(i_in + 1 , j_in - 1)

        sobel_sum = sum(sobel_array)
    end subroutine x_sobel

    subroutine y_sobel(pixels_array_in , sobel_sum , i_in , j_in)
        implicit none
        integer , intent(IN)   :: i_in , j_in
        integer , intent(IN)   :: pixels_array_in(: , :)
        integer , dimension(6) :: sobel_array = 0
        integer , intent(OUT)  :: sobel_sum

        sobel_array(1) = -1*pixels_array_in(i_in - 1 , j_in + 1)                    !multiplies according to GY
        sobel_array(2) = -2*pixels_array_in(i_in , j_in + 1)
        sobel_array(3) = -1*pixels_array_in(i_in + 1 , j_in + 1)
        sobel_array(4) = 2*pixels_array_in(i_in , j_in - 1)
        sobel_array(5) = 1*pixels_array_in(i_in - 1 , j_in - 1)
        sobel_array(6) = 1*pixels_array_in(i_in + 1 , j_in - 1)

        sobel_sum = sum(sobel_array)
    end subroutine y_sobel

!--------------------------------------------------------------------------------------!
!fft_subroutine calls from fftw3 library and takes filenames, mask radius defines the  !
!cut-off value for the mask when fft_subroutine is used for edge detection.            !
!blur_length defines the length of the rectangular convolution signal used when        !
!subroutine is used for deconvolution.                                                 !
!amp defines how large the rectangular convolution signal is.                          !
!case dictates if the fft_subroutine is used for deconvolution (blur) or edge detection!
!(edge).                                                                               !
!--------------------------------------------------------------------------------------!

    subroutine fft_subroutine(filename_in , filename_out , fft_radius , blur_length , amp , case)
        implicit none
        integer(kind = i64)                :: plan_forward, plan_backward
        integer                            :: file_unit_in = 30 , file_unit_out = 40 , max_greys , istat , i , j
        integer                            :: Nx , Ny
        character(len = 20) , intent(IN)   :: filename_in , filename_out , case
        real (kind = dp) , intent(IN)      :: fft_radius , blur_length , amp
        real (kind = dp) , allocatable     :: image(: , :)                      !image etc defined here
        complex (kind = dp) , allocatable  :: cimage(: , :)
        integer , allocatable              :: image_int(: , :)
        real (kind = dp) , allocatable     :: blur(:)
        complex (kind = dp) , allocatable  :: cblur(:)

        call pbm_readwrite(image_int , file_unit_in , max_greys , 1 , filename_in) !image_int is allocated here

        Nx = size(image_int , dim = 1)
        Ny = size(image_int , dim = 2)

        allocate(image(1:Nx , 1:Ny) , stat = istat)                             !image_int used to allocate image, cimage
        if(istat /= 0) STOP 'image failure'
        allocate(cimage(1:Nx/2 + 1 , 1:Ny) , stat = istat)
        if(istat /= 0) STOP 'cimage failure'

        do i = 1 , Nx , 1
            do j = 1 , Ny , 1
                image(i , j) = real(image_int(i , j) , kind = dp)               !pbm_readwrite produces integer array, but real array required for fftw3
            end do
        end do

        call dfftw_plan_dft_r2c_2d(plan_forward , Nx , Ny , image , cimage , FFTW_ESTIMATE)
        call dfftw_execute(plan_forward)                                        !converting image into Fourier space
        call dfftw_destroy_plan(plan_forward)

        select case(case)                                                       !dictating to do edges or deconvolution
            case('edge')                                                        !edge detection
                do i = 0 , Nx/2 - 1
                    do j = 0 , Ny - 1
                        if((real(i , kind = dp)**2.0_dp + real(j , kind = dp)**2.0_dp)**0.5_dp <= fft_radius .or. &
                        & (real(i , kind = dp)**2.0_dp + real(Ny - j + 1 , kind = dp)**2.0_dp)**0.5_dp <= fft_radius)then
                            cimage(i + 1 , j + 1) = 0.0_dp*cimage(i + 1 , j + 1)!implementation of the filter mask function
                        end if
                    end do
                end do
            case('blur')
                allocate(blur(1:Nx) , stat = istat)                             !allocating the blur function
                if(istat /= 0) STOP 'blur failure'
                allocate(cblur(1:Nx/2 + 1) , stat = istat)
                if(istat /= 0) STOP 'cblur failure'
                blur = 0.0_dp

                do i = 1 , int(blur_length) , 1                                 !defining the rectangular convolution signal
                    blur(i) = amp
                end do

                call dfftw_plan_dft_r2c_1d(plan_forward , Nx , blur , cblur , FFTW_ESTIMATE)
                call dfftw_execute(plan_forward)                                !passing rectangular signal to Fourier space
                call dfftw_destroy_plan(plan_forward)

                do j = 1 , Ny , 1                                               !deconvolution here
                    do i = 1 , Nx/2 + 1
                        if(abs(cblur(i)) >= 1E-7) then                          !ensures that no div. by 0
                            cimage(i , j) = cimage(i , j)/cblur(i)
                        else
                            cimage(i , j) = cimage(i , j)
                        end if
                    end do
                end do
            case default
                STOP 'blur case failure'
            end select

        image = 0.0_dp

        call dfftw_plan_dft_c2r_2d(plan_backward , Nx , Ny , cimage , image , FFTW_ESTIMATE)
        call dfftw_execute(plan_backward)
        call dfftw_destroy_plan(plan_backward)                                  !converts cimage back to real space (image)

        do i = 1 , Nx , 1
            do j = 1 , Ny , 1
                image_int(i , j) = int(abs(image(i , j)))                       !converts image to int for pgm write
            end do
        end do

        image_int = image_int/Nx

        call pbm_readwrite(image_int , file_unit_out , max_greys , 2 , filename_out)

    end subroutine fft_subroutine

end program imaging
