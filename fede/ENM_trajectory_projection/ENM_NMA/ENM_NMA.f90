program ENM_NMA
	implicit none

	interface
		subroutine get_CA(pdb_id, chain, x_coord, y_coord, z_coord, res_mass, res_cha, res_number)
			character (len = 100), intent(in) :: chain
			character (len = 100), intent(in) :: pdb_id
			double precision, intent(out), allocatable :: x_coord(:), y_coord(:), z_coord(:), res_mass(:)
			character (len = 1), intent(out), allocatable :: res_cha(:)
			integer, intent(out), allocatable :: res_number(:)
		end subroutine get_CA
	end interface	
	
	character (len = 100) :: protein_id, chains_protein
	integer :: n_CA, info, n, i, j, k, cutoff_chains
	character (len = 10) :: cutoff_chains_str
	double precision, allocatable :: X_CA(:), Y_CA(:), Z_CA(:), CA_mass(:)
	character, allocatable :: residue_chain(:)
	integer, allocatable :: residue_number(:)
	double precision, allocatable :: spring_matrix(:,:), hessian(:,:), mass_w_hessian(:,:), mass_w_hessian_sacr(:,:)
	double precision, allocatable :: evals_real(:), evals_imag(:), dummy(:,:), evecs(:,:), wrk(:)
	double precision, allocatable :: sorted_evals(:), sorted_evecs(:,:), frequencies_GHz(:)
	double precision, allocatable :: coord_matrix(:,:)
	logical :: exist
	character (len = 4) :: pdb_name, atom_label
	character (len = 100) :: chain_exp
	integer :: iostat, count, num_atom_l, num_res_l
	double precision :: x_l, y_l, z_l
	character (len = 1) :: alt_loc_l, chain_l
	character (len = 3) :: res_type_l
	character (len = 3) :: atom_type_l

	double precision, parameter :: mass_conv_u_to_kg = 0.16605, hessian_conv_kcalmolA2_to_Nm = 0.6950
	double precision, parameter :: PI = 3.1415, kB = 0.001380649

	!Scaling of masses: 1 uma = 0.16605E-26 kg. The E-26 term is applied after the eigendecomp. to avoid numerical errors
	!Scaling of stiffnesses: 1 kcal/(molÅ^2) = 0.6950 N/m
	!kB is expressed in kg*Å^2/(K*s^2)

	external dgeev


	!!! READ INPUTS !!!
 	call getarg(1, protein_id)
 	call getarg(2, chains_protein)
 	call getarg(3, cutoff_chains_str)

 	read(cutoff_chains_str, *) cutoff_chains

	!!!! READ THE PDB AND C-ALPHA ATOMS (COORDS, EXP B_FACTORS, MASSES) !!!!
	call get_CA(protein_id, chains_protein, X_CA, Y_CA, Z_CA, CA_mass, residue_chain, residue_number)
	n_CA = size(CA_mass)
	CA_mass = CA_mass*mass_conv_u_to_kg


	!!!! BUILD THE edENM !!!!
	allocate (spring_matrix(n_CA,n_CA))

	call get_spring_matrix(n_CA,X_CA,Y_CA,Z_CA,residue_chain,cutoff_chains,spring_matrix,residue_number)
	spring_matrix = spring_matrix*hessian_conv_kcalmolA2_to_Nm

	allocate (hessian(3*n_CA,3*n_CA))
	call get_hessian(n_CA,X_CA,Y_CA,Z_CA,spring_matrix,hessian)


	!!!! GENERATE MASS-WEIGHTED HESSIAN AND EXTRACT FREQUENCIES AND NORMAL MODES !!!!
	allocate (mass_w_hessian(3*n_CA,3*n_CA))
	allocate (mass_w_hessian_sacr(3*n_CA,3*n_CA))
	call get_mass_w_hessian(n_CA,CA_mass,hessian,mass_w_hessian)

	allocate (evals_real(3*n_CA))
	allocate (evals_imag(3*n_CA))
	allocate (evecs(3*n_CA,3*n_CA))
	allocate (wrk(8*3*n_CA))

	mass_w_hessian_sacr = mass_w_hessian
	call dgeev('N','V',3*n_CA,mass_w_hessian_sacr,3*n_CA,evals_real,evals_imag,dummy,1,evecs,3*n_CA,wrk,8*3*n_CA,info)


	!!!! SORT EIGENVALUES IN ASCENDING ORDER !!!!
	allocate (sorted_evals(3*n_CA))
	allocate (sorted_evecs(3*n_CA,3*n_CA))

	call sort_eigs_increasing(n_CA,evals_real,evecs,sorted_evals,sorted_evecs)

	allocate (frequencies_GHz(3*n_CA-6))
	do n = 1, 3*n_CA - 6
		frequencies_GHz(n) = (sqrt(sorted_evals(n+6))/(2*PI))*(10**4)
	end do

	inquire(file = "NM_frequencies_GHz.txt", exist = exist)
	if (exist) then
		open(21,file = "NM_frequencies_GHz.txt",status = 'old')
		close(21,status='delete')
	end if
	open(22,file = "NM_frequencies_GHz.txt",status = 'new')	
	do n = 1, 3*n_CA - 6
		write(22,'(f7.1)') frequencies_GHz(n)
	end do
	close(22)

	inquire(file = "NMs.txt", exist = exist)
	if (exist) then
		open(23,file = "NMs.txt",status = 'old')
		close(23,status='delete')
	end if
	open(24,file = "NMs.txt",status = 'new')	
400 format(10f10.5)
	do i = 1, 3*n_CA
		write(24,400) sorted_evecs(i,7), sorted_evecs(i,8), sorted_evecs(i,9), sorted_evecs(i,10),&
		sorted_evecs(i,11), sorted_evecs(i,12), sorted_evecs(i,13), sorted_evecs(i,14),&
		sorted_evecs(i,15), sorted_evecs(i,16)
	end do
	close(24)

end program ENM_NMA


subroutine get_CA(pdb_id, chain, x_coord, y_coord, z_coord, res_mass, res_cha, res_number)
	implicit none
	
	character (len = 1) :: alt_loc, chain_id
	character (len = 100), intent(in) :: pdb_id, chain
	character (len = 3) :: atom_type, res_type
	character (len = 4) :: label
	character (len = 80) :: string
	integer :: pdb_iostat
	integer :: num_atoms, num_res, n, atom_num, res_num, num_chains, i
	double precision :: x, y, z, occup, b_factor
	double precision, intent(out), allocatable :: x_coord(:), y_coord(:), z_coord(:), res_mass(:)
	character (len = 1), intent(out), allocatable :: res_cha(:)
	integer, intent(out), allocatable :: res_number(:)
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
			read(14,20) label, atom_num, atom_type, alt_loc, res_type, chain_id, res_num, x, y, z, occup, b_factor
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
	allocate (res_number(num_res))

	num_res = 0
	do i = 1, num_chains
		open(15,file = trim(pdb_id)//"_ATOM.pdb",status = 'old')
		do n = 1, num_atoms
			read(15,20) label, atom_num, atom_type, alt_loc, res_type, chain_id, res_num, x, y, z, occup, b_factor
			if (chain_id == chain(i:i)) then 
				if (atom_type(1:2) == 'CA') then
					if (alt_loc == 'A' .or. alt_loc == ' ') then
						num_res = num_res + 1
						x_coord(num_res) = x
						y_coord(num_res) = y
						z_coord(num_res) = z
						res_cha(num_res) = chain_id
						res_number(num_res) = res_num
						if (res_type == 'ALA') then
							res_mass(num_res) = 71
						else if (res_type == 'ARG') then
							res_mass(num_res) = 156
						else if (res_type == 'ASN') then
							res_mass(num_res) = 114
						else if (res_type == 'ASP') then
							res_mass(num_res) = 115
						else if (res_type == 'CYS') then
							res_mass(num_res) = 103
						else if (res_type == 'GLU') then
							res_mass(num_res) = 129
						else if (res_type == 'GLN') then
							res_mass(num_res) = 128
						else if (res_type == 'GLY') then
							res_mass(num_res) = 57
						else if (res_type == 'HIS') then
							res_mass(num_res) = 137
						else if (res_type == 'ILE') then
							res_mass(num_res) = 113
						else if (res_type == 'LEU') then
							res_mass(num_res) = 113
						else if (res_type == 'LYS') then
							res_mass(num_res) = 128
						else if (res_type == 'MET') then
							res_mass(num_res) = 131
						else if (res_type == 'PHE') then
							res_mass(num_res) = 147
						else if (res_type == 'PRO') then
							res_mass(num_res) = 97
						else if (res_type == 'SER') then
							res_mass(num_res) = 87
						else if (res_type == 'THR') then
							res_mass(num_res) = 101
						else if (res_type == 'TRP') then
							res_mass(num_res) = 186
						else if (res_type == 'TYR') then
							res_mass(num_res) = 163
						else if (res_type == 'VAL') then
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


subroutine get_spring_matrix(n,x_coord,y_coord,z_coord,res_chain_label,cutoff_chains,spr_matrix,res_numb)
	implicit none

	integer, intent(in) :: n, cutoff_chains
	double precision, intent(in) :: x_coord(n), y_coord(n), z_coord(n)
	double precision, intent(out) :: spr_matrix(n,n)
	character (len = 1), intent(in) :: res_chain_label(n)
	integer, intent(in) :: res_numb(n)

	integer :: i, j, edENM_cutoff
	integer, parameter :: edENM_M = 3, edENM_seq_exp = 2, edENM_cart_exp = 6
	double precision, parameter :: edENM_seq_force_constant = 60, edENM_cart_force_constant = 6
	double precision :: dist
	integer :: seq_dist

			
	edENM_cutoff = nint(2.9*log(real(n)) - 2.9)
	if (edENM_cutoff >= 20) then
		edENM_cutoff = 20
	else if (edENM_cutoff <= 8) then
		edENM_cutoff = 8
	end if

	spr_matrix = 0.0d0
	do i = 1, n - 1
		do j = i + 1, n
			if (res_chain_label(i) == res_chain_label(j)) then
				seq_dist = abs(res_numb(j) - res_numb(i))
				if (seq_dist <= edENM_M) then
					spr_matrix(i,j) = edENM_seq_force_constant/dble((j-i)**edENM_seq_exp)
				else
					dist = sqrt((x_coord(i)-x_coord(j))**2 + (y_coord(i)-y_coord(j))**2 &
					+ (z_coord(i)-z_coord(j))**2)
					if (dist <= dble(edENM_cutoff)) then
						spr_matrix(i,j) = (edENM_cart_force_constant/dist)**edENM_cart_exp
					end if					
				end if
			else
				dist = sqrt((x_coord(i)-x_coord(j))**2 + (y_coord(i)-y_coord(j))**2 &
				+ (z_coord(i)-z_coord(j))**2)
				if (dist <= dble(cutoff_chains)) then
					spr_matrix(i,j) = (edENM_cart_force_constant/dist)**edENM_cart_exp
				end if
			end if
		end do
	end do
end subroutine get_spring_matrix

subroutine get_hessian(n,x,y,z,spring_matrix,hessian_matrix)
	implicit none

	integer, intent(in) :: n
	double precision, intent(in) :: x(n), y(n), z(n), spring_matrix(n,n)
	double precision, intent(out) :: hessian_matrix(3*n,3*n)

	integer :: i, j
	double precision :: dist, spring_constant

	hessian_matrix = 0.0d0
	do i = 1, n
		do j = 1, n
			if (j.ne.i) then
				dist = sqrt((x(i) - x(j))**2 + (y(i) - y(j))**2 + (z(i) - z(j))**2)
				
				if (j > i) then
					spring_constant = spring_matrix(i,j)
				else
					spring_constant = spring_matrix(j,i)
				end if

				hessian_matrix(3*i-2,3*j-2) = -spring_constant*(x(i)-x(j))*(x(i)-x(j))/(dist**2)
				hessian_matrix(3*i-2,3*j-1) = -spring_constant*(x(i)-x(j))*(y(i)-y(j))/(dist**2)
				hessian_matrix(3*i-2,3*j) = -spring_constant*(x(i)-x(j))*(z(i)-z(j))/(dist**2)
				hessian_matrix(3*i-1,3*j-2) = -spring_constant*(y(i)-y(j))*(x(i)-x(j))/(dist**2)
				hessian_matrix(3*i-1,3*j-1) = -spring_constant*(y(i)-y(j))*(y(i)-y(j))/(dist**2)
				hessian_matrix(3*i-1,3*j) = -spring_constant*(y(i)-y(j))*(z(i)-z(j))/(dist**2)
				hessian_matrix(3*i,3*j-2) = -spring_constant*(z(i)-z(j))*(x(i)-x(j))/(dist**2)
				hessian_matrix(3*i,3*j-1) = -spring_constant*(z(i)-z(j))*(y(i)-y(j))/(dist**2)
				hessian_matrix(3*i,3*j) = -spring_constant*(z(i)-z(j))*(z(i)-z(j))/(dist**2)

				hessian_matrix(3*i-2,3*i-2) = hessian_matrix(3*i-2,3*i-2) - hessian_matrix(3*i-2,3*j-2)
				hessian_matrix(3*i-2,3*i-1) = hessian_matrix(3*i-2,3*i-1) - hessian_matrix(3*i-2,3*j-1)
				hessian_matrix(3*i-2,3*i) = hessian_matrix(3*i-2,3*i) - hessian_matrix(3*i-2,3*j)
				hessian_matrix(3*i-1,3*i-2) = hessian_matrix(3*i-1,3*i-2) - hessian_matrix(3*i-1,3*j-2)
				hessian_matrix(3*i-1,3*i-1) = hessian_matrix(3*i-1,3*i-1) - hessian_matrix(3*i-1,3*j-1)
				hessian_matrix(3*i-1,3*i) = hessian_matrix(3*i-1,3*i) - hessian_matrix(3*i-1,3*j)
				hessian_matrix(3*i,3*i-2) = hessian_matrix(3*i,3*i-2) - hessian_matrix(3*i,3*j-2)
				hessian_matrix(3*i,3*i-1) = hessian_matrix(3*i,3*i-1) - hessian_matrix(3*i,3*j-1)
				hessian_matrix(3*i,3*i) = hessian_matrix(3*i,3*i) - hessian_matrix(3*i,3*j)
			end if
		end do
	end do
end subroutine get_hessian

subroutine get_mass_w_hessian(n,mass_res,hessian,mass_w_hessian)
	implicit none

	integer, intent(in) :: n
	double precision, intent(in) :: mass_res(n), hessian(3*n,3*n)
	double precision, intent(out) :: mass_w_hessian(3*n,3*n)

	integer :: i, j

	do i = 1, n
		do j = 1, n
			mass_w_hessian(3*i-2,3*j-2) = hessian(3*i-2,3*j-2)/sqrt(mass_res(i)*mass_res(j))
			mass_w_hessian(3*i-2,3*j-1) = hessian(3*i-2,3*j-1)/sqrt(mass_res(i)*mass_res(j))
			mass_w_hessian(3*i-2,3*j) = hessian(3*i-2,3*j)/sqrt(mass_res(i)*mass_res(j))
			mass_w_hessian(3*i-1,3*j-2) = hessian(3*i-1,3*j-2)/sqrt(mass_res(i)*mass_res(j))
			mass_w_hessian(3*i-1,3*j-1) = hessian(3*i-1,3*j-1)/sqrt(mass_res(i)*mass_res(j))
			mass_w_hessian(3*i-1,3*j) = hessian(3*i-1,3*j)/sqrt(mass_res(i)*mass_res(j))
			mass_w_hessian(3*i,3*j-2) = hessian(3*i,3*j-2)/sqrt(mass_res(i)*mass_res(j))
			mass_w_hessian(3*i,3*j-1) = hessian(3*i,3*j-1)/sqrt(mass_res(i)*mass_res(j))
			mass_w_hessian(3*i,3*j) = hessian(3*i,3*j)/sqrt(mass_res(i)*mass_res(j))
		end do
	end do
end subroutine get_mass_w_hessian

subroutine sort_eigs_increasing(n,evals,evecs,sorted_evals,sorted_evecs) 
	implicit none
	
	integer :: i, j
	integer, intent(in) :: n
	integer :: sort_index(3*n)
	double precision, intent(in) :: evals(3*n), evecs(3*n,3*n)
	double precision, intent(out) :: sorted_evals(3*n), sorted_evecs(3*n,3*n)
	double precision :: min_value, sacr_evals(3*n)


	!Define sort_indexes vector - evals in ascending order!

	sacr_evals = evals
	do i = 1,3*n
		min_value = maxval(sacr_evals)
		do j = 1,3*n
			if (sacr_evals(j) <= min_value) then
				min_value = sacr_evals(j)
				sort_index(i) = j
			end if
		end do
		sacr_evals(sort_index(i)) = maxval(sacr_evals) + 1D+0
	end do

	!Sort evals and evecs!

	do i = 1,3*n
		sorted_evals(i) = evals(sort_index(i))
		do j = 1,3*n
			sorted_evecs(j,i) = evecs(j,sort_index(i))
		end do
	end do
end subroutine sort_eigs_increasing