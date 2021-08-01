program charge
    implicit none
    
    integer i, j
    real pi, var
    real, dimension (182, 2) :: Input
    real, dimension (182, 2) :: Output
    pi = 4 * atan(1.0)
    open (unit = 1, file = "density.dat", status = 'old' )
        do i = 1, 182
            read (1,*,end = 10) (Input(i, j) ,  j = 1, 2)
        end do
    10 close (unit = 1)
    do i = 1, 182
        Output(i, 1) = Input(i, 1)
    end do
    Output(1, 2) = (Input(1, 1) ** 3) * Input(1, 2) * 4 * pi
    do i = 2, 182
        var = (Input(i, 1) ** 2) * 4 * pi * (Input(i, 1) - Input(i-1, 1))
        var = var * (Input(i, 2) + Input(i-1, 2)) / 2
        Output(i, 2) = Output(i-1, 2) + var
    end do
    open(unit = 2 , file="charge_enclosed.txt")
    do i = 1, 182
        write (2, *) Output(i, 1), Output (i, 2)
    end do
    close (unit = 2)
end program charge