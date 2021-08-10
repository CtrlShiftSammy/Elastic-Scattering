program divider
    implicit none
    integer i, j
    real, dimension (182, 2) :: Input
    open (unit = 1, file = "charge_enclosed.txt", status = 'old' )
    open(unit = 2 , file="charge_enclosed_modified.txt")
    do i = 1, 182
        read (1,*,end = 10) (Input(i, j) ,  j = 1, 2)
        Input(i, 2) = 54.0 - Input(i, 2) * 52.0 / 56.0
        write (2, *) Input(i, 1), Input(i, 2)
    end do
    10 close (unit = 1)
    close (unit = 2)
end program divider