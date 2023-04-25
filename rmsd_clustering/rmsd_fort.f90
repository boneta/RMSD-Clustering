
real(8) function rmsd(coord0, coord1, natoms)

    ! Calculate the RMSD between two sets of coordinates

    implicit none
    integer(4), intent(in) :: natoms
    real(8), intent(in) :: coord0(natoms, 3), coord1(natoms, 3)
    integer(4) :: i, j

    rmsd = 0.D0
    do i = 1, natoms
        do j = 1, 3
            rmsd = rmsd + (coord0(i, j) - coord1(i, j))**2
        end do
    end do

    rmsd = sqrt(rmsd / natoms)

end function rmsd

real(8) function det_3x3(matrix)

    ! Calculate the determinant of a 3x3 matrix

    implicit none
    real(8), intent(in) :: matrix(3, 3)

    det_3x3 = matrix(1,1) * (matrix(2,2) * matrix(3,3) - matrix(2,3) * matrix(3,2)) &
            - matrix(1,2) * (matrix(2,1) * matrix(3,3) - matrix(2,3) * matrix(3,1)) &
            + matrix(1,3) * (matrix(2,1) * matrix(3,2) - matrix(2,2) * matrix(3,1))

end function det_3x3

subroutine svd_3x3(matrix, u, s, vt)

    ! SVD decomposition of a 3x3 matrix
    ! matrix = u * s * vt

    implicit none
    real(8), intent(in) :: matrix(3,3)
    real(8), intent(out) :: u(3,3), s(3), vt(3,3)
    integer(4) :: info
    real(8) :: work(1000)

    ! query the optimal workspace
    call DGESVD( 'A', 'A', 3, 3, matrix, 3, s, u, 3, vt, 3, work, -1, info)
    ! compute SVD
    call DGESVD( 'A', 'A', 3, 3, matrix, 3, s, u, 3, vt, 3, work, int(work(1)), info)
    ! check for convergence
    if ( info .ne. 0 ) then
        print*, 'SVD FAILED'
        stop 1
    end if

end subroutine svd_3x3

subroutine kabsch(coord0, coord1, rot_matrix, natoms)

    ! Calculate the rotation matrix that minimizes the RMSD between two sets of coordinates

    implicit none
    integer(4), intent(in) :: natoms
    real(8), intent(in) :: coord0(natoms, 3), coord1(natoms, 3)
    real(8), intent(out) :: rot_matrix(3,3)
    real(8), dimension(3,3) :: covar(3,3), u(3,3), s(3), vt(3,3)
    real(8), external :: det_3x3

    covar = matmul(transpose(coord0), coord1)
    call svd_3x3(covar, u, s, vt)
    if ( det_3x3(u) * det_3x3(vt) .lt. 0 ) then
        s(3) = -s(3)
        u(:,3) = -u(:,3)
    end if
    rot_matrix = matmul(u, vt)

end subroutine kabsch

real(8) function rmsd_kabsch(coord0, coord1, natoms)

    ! Calculate the RMSD between two sets of coordinates after Kabsch alignment

    implicit none
    integer(4), intent(in) :: natoms
    real(8), intent(in) :: coord0(natoms, 3), coord1(natoms, 3)
    real(8) :: coord0_rot(natoms, 3), rot_matrix(3,3)
    real(8), external :: rmsd

    ! calculate rotation matrix
    call kabsch(coord0, coord1, rot_matrix, natoms)

    ! rotate coord1
    coord0_rot = matmul(coord0, rot_matrix)

    ! calculate RMSD
    rmsd_kabsch = rmsd(coord0_rot, coord1, natoms)

end function rmsd_kabsch

subroutine rmsd_matrix(molecs, align, matrix, natoms, nmolecs)

    ! Calculate the RMSD matrix between a set of molecules

    implicit none
    integer(4), intent(in) :: natoms, nmolecs
    real(8), intent(in) :: molecs(nmolecs, natoms, 3)
    logical, intent(in) :: align
    real(8), intent(out) :: matrix(nmolecs, nmolecs)
    real(8), dimension(natoms, 3) :: coord0, coord1
    real(8), external :: rmsd, rmsd_kabsch
    integer(4) :: i, j

    matrix = 0.D0
    !$OMP PARALLEL DO SHARED(matrix) PRIVATE(i, j, coord0, coord1)
    do i = 1, nmolecs
        coord0 = molecs(i, :, :)
        do j = i+1, nmolecs
            coord1 = molecs(j, :, :)
            if ( align ) then
                matrix(i, j) = rmsd_kabsch(coord0, coord1, natoms)
            else
                matrix(i, j) = rmsd(coord0, coord1, natoms)
            end if
            matrix(j, i) = matrix(i, j)
        end do
    end do
    !$OMP END PARALLEL DO

end subroutine rmsd_matrix
