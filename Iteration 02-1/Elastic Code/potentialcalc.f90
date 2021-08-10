program potentialcalc
    implicit none
    integer i, j
    real pi
    real, dimension (182, 2) :: Input
    real, dimension (182, 2) :: Output
    pi = 4 * atan(1.0)
    open (unit = 1, file = "charge_enclosed.txt", status = 'old' )
        do i = 1, 182
            read (1,*,end = 10) (Input(i, j) ,  j = 1, 2)
        end do
    10 close (unit = 1)
    do i = 1, 182
        Output(i, 1) = Input(i, 1)
        Output(i, 2) = (54.0 - Input(i, 2)) / Input(i, 1)
    end do
    open(unit = 2 , file="calculated_potentials.txt")
    do i = 1, 182
        write (2, *) Output(i, 1), Output (i, 2)
    end do
    close (unit = 2)
end program potentialcalc