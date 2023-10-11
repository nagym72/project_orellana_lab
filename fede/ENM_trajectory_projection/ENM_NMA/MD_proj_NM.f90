program MD_proj_NM
	implicit none

	integer :: n_CA, i, j
	character (len = 50) :: frame_name, ref_NM_name
	character (len = 10) :: n_CA_str
	double precision, allocatable :: frame_coord_matrix(:)
	double precision, allocatable :: NM_vectors(:,:)
	character (len = 4) :: label
	integer :: num_atom, num_res
	character (len = 3) :: atom_type
	character (len = 3) :: res_type
	character (len = 1) :: chain
	double precision :: x, y, z
	double precision, allocatable :: red_frame_conf1_coord_matrix(:)
	double precision, allocatable :: ref_coord_matrix(:)
	double precision :: proj(10)


	call getarg(1, n_CA_str)
	call getarg(2, frame_name)
	call getarg(3, ref_NM_name)

	read (n_CA_str, *) n_CA

	
	allocate(NM_vectors(3*n_CA,10))
400 format(10f10.5)
	open(11,file="NMs.txt",status='old')
	do i = 1, 3*n_CA
		read(11,400) NM_vectors(i,1), NM_vectors(i,2), NM_vectors(i,3), NM_vectors(i,4), NM_vectors(i,5), NM_vectors(i,6), &
		NM_vectors(i,7), NM_vectors(i,8), NM_vectors(i,9), NM_vectors(i,10)
	end do
	close(11)

	!!! PROJECT MD FRAME CONFORMATION ON NM SPACE !!!

	allocate(frame_coord_matrix(3*n_CA))
555 format(a4,i7,2x,a3,1x,a3,1x,a1,i4,4x,3f8.3)
	open(12,file=trim(frame_name),status = 'old')
	do i = 1, n_CA
		read(12,555) label, num_atom, atom_type, res_type, chain, num_res, x, y, z
		frame_coord_matrix(3*i-2) = x
		frame_coord_matrix(3*i-1) = y
		frame_coord_matrix(3*i) = z
	end do
	close(12)

	allocate(ref_coord_matrix(3*n_CA))
	open(13,file=trim(ref_NM_name),status = 'old')
	do i = 1, n_CA
		read(13,555) label, num_atom, atom_type, res_type, chain, num_res, x, y, z
		ref_coord_matrix(3*i-2) = x
		ref_coord_matrix(3*i-1) = y
		ref_coord_matrix(3*i) = z
	end do
	close(13)

	allocate(red_frame_conf1_coord_matrix(3*n_CA))
	do i = 1, 3*n_CA
		red_frame_conf1_coord_matrix(i) = frame_coord_matrix(i) - ref_coord_matrix(i)
	end do

	proj = 0.0d0
	do i = 1, 3*n_CA
		do j = 1, 10
			proj(j) = proj(j) + NM_vectors(i,j)*red_frame_conf1_coord_matrix(i)
		end do
	end do  
	

300 format(10f9.2)
	open(14,file='proj_'//trim(frame_name)//'.txt',status='new')
	write(14,300) proj(1), proj(2), proj(3), proj(4), proj(5), proj(6), proj(7), proj(8), proj(9), proj(10)
	close(14)

end program MD_proj_NM