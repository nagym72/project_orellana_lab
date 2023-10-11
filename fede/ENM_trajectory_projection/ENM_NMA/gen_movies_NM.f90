program gen_movies_NM
	implicit none

	interface
		subroutine get_CA(pdb_id, chain, x_coord, y_coord, z_coord, res_mass, res_cha, res_type)
			character (len = 20), intent(in) :: chain
			character (len = 20), intent(in) :: pdb_id
			double precision, intent(out), allocatable :: x_coord(:), y_coord(:), z_coord(:), res_mass(:)
			character (len = 1), intent(out), allocatable :: res_cha(:)
			character (len = 3), intent(out), allocatable :: res_type(:)
		end subroutine get_CA
	end interface		
	
	character (len = 20) :: protein_id, chains_protein
	integer :: n_CA, info, n, i, num_modes, k
	character (len = 10) :: scale_ampl_str
	double precision, allocatable :: X_CA(:), Y_CA(:), Z_CA(:), CA_mass(:)
	character, allocatable :: residue_chain(:)
	character (len = 3), allocatable :: res_name(:)
	logical :: exist
	double precision, allocatable :: nms_to_plot(:,:)
	double precision :: scale_ampl
	integer :: num_frames
	character (len = 2) :: auxstr

	double precision, parameter :: PI = 3.1415

	!!! READ INPUTS !!!
 	call getarg(1, protein_id)
 	call getarg(2, chains_protein)
 	call getarg(3, scale_ampl_str)

 	read(scale_ampl_str, *) scale_ampl


 	num_modes = 10
 	num_frames = 10

	!!!! READ THE PDB AND C-ALPHA ATOMS (COORDS, MASSES) !!!!

	call get_CA(protein_id, chains_protein, X_CA, Y_CA, Z_CA, CA_mass, residue_chain, res_name)
	n_CA = size(CA_mass)


	!!! GENERATE MOVIES OF THE MODES !!!

	allocate(nms_to_plot(3*n_CA,num_modes))

200 format(10f10.5)
	open(11,file = "NMs.txt",status = 'old')	
	do i = 1, 3*n_CA
		read(11,200) nms_to_plot(i,1), nms_to_plot(i,2), nms_to_plot(i,3), nms_to_plot(i,4), nms_to_plot(i,5),&
		nms_to_plot(i,6), nms_to_plot(i,7), nms_to_plot(i,8), nms_to_plot(i,9), nms_to_plot(i,10)
	end do
	close(11)

100 format(a4,i7,2x,a3,a1,a3,1x,a1,i4,4x,3f8.3)
	do n = 1, num_modes
		write(auxstr,'(i2.2)') n
		inquire(file = trim(protein_id)//"_"//trim(chains_protein)//"_NM_"//trim(auxstr)//".pdb", exist=exist)
		if (exist) then
			open(12,file = trim(protein_id)//"_"//trim(chains_protein)//"_NM_"//trim(auxstr)//".pdb",status = 'old')
			close(12,status='delete')
		end if
		open(13,file = trim(protein_id)//"_"//trim(chains_protein)//"_NM_"//trim(auxstr)//".pdb",status = 'new')
		do k = 1, 2*num_frames + 3
			write(13,'(a5)') "MODEL"
			do i = 1, n_CA
				write(13,100) "ATOM", i, "CA ", " ", res_name(i), residue_chain(i), i, &
				X_CA(i)+scale_ampl*nms_to_plot(3*i-2,n)*sin(PI*(k-1)/(num_frames+1)), &
				Y_CA(i)+scale_ampl*nms_to_plot(3*i-1,n)*sin(PI*(k-1)/(num_frames+1)), &
				Z_CA(i)+scale_ampl*nms_to_plot(3*i,n)*sin(PI*(k-1)/(num_frames+1))
			end do
			write(13,'(a3)') "TER"
			write(13,'(a6)') "ENDMDL"
		end do
		close(13)
	end do
	
end program gen_movies_NM


subroutine get_CA(pdb_id, chain, x_coord, y_coord, z_coord, res_mass, res_cha, res_type)
	implicit none
	
	character (len = 1) :: alt_loc, chain_id
	character (len = 20), intent(in) :: pdb_id, chain
	character (len = 3) :: atom_type, res_name
	character (len = 4) :: label
	character (len = 80) :: string
	integer :: pdb_iostat
	integer :: num_atoms, num_res, n, atom_num, res_num, num_chains, i
	double precision :: x, y, z, occup, b_factor
	double precision, intent(out), allocatable :: x_coord(:), y_coord(:), z_coord(:), res_mass(:)
	character (len = 1), intent(out), allocatable :: res_cha(:)
	character (len = 3), intent(out), allocatable :: res_type(:)
	logical :: exist


	open(11,file = trim(pdb_id)//".pdb", status = 'old', iostat = pdb_iostat)
	if (pdb_iostat .ne. 0) then
		print *, "I couldn't open the file "//trim(pdb_id)//".pdb "//" or the file does not exist!"
		stop
	end if

	inquire(file = trim(pdb_id)//"_ATOM.pdb", exist=exist)
	if (exist) then
		open(12,file = trim(pdb_id)//"_ATOM.pdb",status = 'old')
		close(12,status='delete')
	end if
	open(13,file = trim(pdb_id)//"_ATOM.pdb",status = 'new')

	num_atoms = 0
	do
	    read(11,'(a80)',end = 10) string
	    if (trim(string(1:4)) == 'ATOM') then
	    	num_atoms = num_atoms + 1
	    	write (13,'(a80)') string
	    end if
	end do
10  close(11)
	close(13)

	! format :ATOM    821  CA  LEU A 109      28.830  71.526  48.547  1.00 31.28           C  
20  format(a4,i7,2x,a3,a1,a3,1x,a1,i4,4x,3f8.3,2x,f4.2,f6.2)

	num_res = 0
	num_chains = len(trim(chain))
	do i = 1, num_chains
		open(14,file = trim(pdb_id)//"_ATOM.pdb",status = 'old')
		do n = 1, num_atoms
			read(14,20) label, atom_num, atom_type, alt_loc, res_name, chain_id, res_num, x, y, z, occup, b_factor
			if (chain_id == chain(i:i)) then 
				if (atom_type(1:2) == 'CA') then
					if (alt_loc == 'A' .or. alt_loc == ' ') then
						num_res = num_res + 1
					end if
				end if
			end if
		end do
		close(14)
	end do
	
	allocate (x_coord(num_res))
	allocate (y_coord(num_res))
	allocate (z_coord(num_res))
	allocate (res_mass(num_res))
	allocate (res_cha(num_res))
	allocate (res_type(num_res))

	num_res = 0
	do i = 1, num_chains
		open(15,file = trim(pdb_id)//"_ATOM.pdb",status = 'old')
		do n = 1, num_atoms
			read(15,20) label, atom_num, atom_type, alt_loc, res_name, chain_id, res_num, x, y, z, occup, b_factor
			if (chain_id == chain(i:i)) then 
				if (atom_type(1:2) == 'CA') then
					if (alt_loc == 'A' .or. alt_loc == ' ') then
						num_res = num_res + 1
						x_coord(num_res) = x
						y_coord(num_res) = y
						z_coord(num_res) = z
						res_cha(num_res) = chain_id
						res_type(num_res) = res_name
						if (res_name == 'ALA') then
							res_mass(num_res) = 71
						else if (res_name == 'ARG') then
							res_mass(num_res) = 156
						else if (res_name == 'ASN') then
							res_mass(num_res) = 114
						else if (res_name == 'ASP') then
							res_mass(num_res) = 115
						else if (res_name == 'CYS') then
							res_mass(num_res) = 103
						else if (res_name == 'GLU') then
							res_mass(num_res) = 129
						else if (res_name == 'GLN') then
							res_mass(num_res) = 128
						else if (res_name == 'GLY') then
							res_mass(num_res) = 57
						else if (res_name == 'HIS') then
							res_mass(num_res) = 137
						else if (res_name == 'ILE') then
							res_mass(num_res) = 113
						else if (res_name == 'LEU') then
							res_mass(num_res) = 113
						else if (res_name == 'LYS') then
							res_mass(num_res) = 128
						else if (res_name == 'MET') then
							res_mass(num_res) = 131
						else if (res_name == 'PHE') then
							res_mass(num_res) = 147
						else if (res_name == 'PRO') then
							res_mass(num_res) = 97
						else if (res_name == 'SER') then
							res_mass(num_res) = 87
						else if (res_name == 'THR') then
							res_mass(num_res) = 101
						else if (res_name == 'TRP') then
							res_mass(num_res) = 186
						else if (res_name == 'TYR') then
							res_mass(num_res) = 163
						else if (res_name == 'VAL') then
							res_mass(num_res) = 99
						else
							res_mass(num_res) = 100
						end if
					end if
				end if
			end if
		end do
		close(15)
	end do

	open(16,file = trim(pdb_id)//"_ATOM.pdb",status = 'old')
	close(16,status='delete')

end subroutine get_CA