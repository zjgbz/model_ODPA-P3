import numpy as np

# Defining the function used for the generating header of the data file.
# The input of this function includes the number of chains, the number of monomers
# in each chain, the number of bond types, the number of atom types,
# the number of angle types, the number of dihedral types, the size of
# the box, and the mass of each kind of atom.

def print_header(polymer_name, monomer_num, chain_num, atom_types_num, bond_types_num, angle_types_num, dihedral_types_num, xlo,xhi,ylo,yhi,zlo,zhi, mass):
	
	# Cacluating the number of bonds, angles and dihedrals based on the number of monomers and chains.

	if atom_types_num == 0:
		atoms_num = 0
	else:
		atoms_num = monomer_num * chain_num
	if bond_types_num == 0:
		bonds_num = 0
	else:
		bonds_num = (monomer_num - 1) * chain_num
	if angle_types_num == 0:
		angles_num = 0
	else:
		angles_num = (monomer_num - 2) * chain_num
	if dihedral_types_num == 0:
		dihedrals_num = 0
	else:
		dihedrals_num = (monomer_num - 3) * chain_num

	print("#Models for %s\n"%polymer_name)

	print("\t%d\tatoms"%atoms_num)
	print("\t%d\tbonds"%bonds_num)
	print("\t%d\tangles"%angles_num)
	print("\t%d\tdihedrals\n"%dihedrals_num)

	print("\t%d\tatom types"%atom_types_num)
	print("\t%d\tbond types"%bond_types_num)
	print("\t%d\tangle types"%angle_types_num)
	print("\t%d\tdihedral types\n"%dihedral_types_num)

	print("\t%.4f\t%.4f\t xlo xhi"%(xlo, xhi))
	print("\t%.4f\t%.4f\t ylo yhi"%(ylo, yhi))
	print("\t%.4f\t%.4f\t zlo zhi\n"%(zlo, zhi))

	print("Masses\n")
	print("\t%d\t%.2f\n"%(atom_types_num,mass))
	return 0

# The file of configuration of the polymers contains four parts: atoms, bonds, angles
# and dihedral. The content of any of them could be printed by the function below.

def print_content(monomer_list, title):
	if title == "Atoms":
		info_col = 3
	else:
		info_col = 2
	print("%s\n"%title)
	line_num = np.size(monomer_list,0)
	element_num = np.size(monomer_list,1)
	for line_i in range(line_num):
		for element_i in range(info_col):
			print("%d\t"%monomer_list[line_i,element_i],end='')
		for element_i in range(info_col,element_num):
			if title == "Atoms":
				print("%.5f\t"%monomer_list[line_i,element_i],end='')
			else:
				print("%d\t"%monomer_list[line_i,element_i],end='')
		print("\n",end='')
	print("\n",end='')
	return 0

# This function is used to generate the content of the atom part of the configuration file.
# By this file, a bunch of uniformly distributed polymer chains will be created based on the size
# of the box and the number of polymer chains and the number monomers in each chain.

def generate_particles(xlo,xhi,ylo,yhi,zlo,zhi,x_redundant,y_redundant,z_redundant,x_interval,y_interval,z_interval,monomer_num,x_chain,y_chain):
	total_monomer_num = monomer_num * chain_num
	monomer_list = np.zeros((total_monomer_num, 9))
	monomer_list[:,0] = range(1, total_monomer_num + 1)
	monomer_list[:,2] = np.ones(total_monomer_num)
	start_i = 0
	chain_i = 1
	for x_chain_i in range(x_chain):
		start_x = xlo + x_chain_i * x_interval + x_redundant
		for y_chain_i in range(y_chain):
			start_y = ylo + y_chain_i * y_interval + y_redundant
			for monomer_i in range(monomer_num):
				monomer_list[monomer_i + start_i,3] = start_x
				monomer_list[monomer_i + start_i,4] = start_y
				monomer_list[monomer_i + start_i,5] = zlo + monomer_i * z_interval + z_redundant
				monomer_list[monomer_i + start_i,1] = chain_i
			start_i = start_i + monomer_num
			chain_i = chain_i + 1
	return monomer_list

# The function here is used to generate all other connection include bonds, angles and dihedrals.

def generate_connections(monomer_num,x_chain,y_chain,connection_type):
	if connection_type == "Bonds":
		involve_particles = 2
	elif connection_type == "Angles":
		involve_particles = 3
	elif connection_type == "Dihedrals":
		involve_particles = 4
	chain_num = x_chain * y_chain
	connection_num = (monomer_num - involve_particles + 1) * chain_num
	connection_per_chain = monomer_num - involve_particles + 1
	connection_list = np.zeros((connection_num,2 + involve_particles))
	connection_list[:,0] = range(1,connection_num + 1)
	connection_list[:,1] = np.ones(connection_num)
	start_i = 0 - (involve_particles - 1) + 1
	for chain_i in range(chain_num):
		start_i = start_i + (involve_particles - 1)
		for connection_i in range(connection_per_chain):
			line_i = chain_i * connection_per_chain + connection_i
			for particle_i in range(2,2 + involve_particles):
				connection_list[line_i,particle_i] = start_i
				start_i = start_i + 1
			start_i = start_i - (involve_particles - 1)
	return (connection_list,connection_type)

# This function is used to print the correspoinding title of certain block of content.

def print_particles_connections(atom_types_num,bond_types_num,angle_types_num,dihedral_types_num,monomer_list,bond_boundle,angle_boundle,dihedral_boundle):
	if atom_types_num == 1:
		print_content(monomer_list,"Atoms")
	if bond_types_num == 1:
		print_content(bond_boundle[0],bond_boundle[1])
	if angle_types_num == 1:
		print_content(angle_boundle[0],angle_boundle[1])
	if dihedral_types_num == 1:
		print_content(dihedral_boundle[0],dihedral_boundle[1])
	return 0

# This function can be used to compute the proper box size based on the LJ coefficient
# and the size of each monomer.
def box_size(x_num,y_num,z_num,diameter,chain_coeff):
	x_size = (max(x_num,y_num,z_num) + 1) * (diameter * chain_coeff)
	y_size = x_size
	z_size = x_size
	xlo = - x_size / 2
	xhi = x_size / 2
	ylo = xlo
	yhi = xhi
	zlo = xlo
	zhi = xhi
	return x_size,y_size,z_size,xlo,xhi,ylo,yhi,zlo,zhi

if __name__ == '__main__':
	# Initialization the parameters.
	polymer_name = "ODPA-P3"
	monomer_num = 85
	x_chain = 50
	y_chain = 50
	chain_num = x_chain * y_chain
	diameter = 2 ** (-1/6) # diameter of the monomer.
	chain_coeff = 4.15310957874
	equil_len = 1.53 # stable position of the particle.
	monomoer_coeff = equil_len / diameter

	atoms_num = monomer_num * chain_num
	bonds_num = (monomer_num - 1) * chain_num
	angles_num = (monomer_num - 2) * chain_num
	dihedrals_num = (monomer_num - 3) * chain_num
	atom_types_num = 1
	bond_types_num = 1
	angle_types_num = 1
	dihedral_types_num = 1
	mass = 14.01
	x_size,y_size,z_size,xlo,xhi,ylo,yhi,zlo,zhi = box_size(x_chain,y_chain,monomer_num,diameter,chain_coeff)
	
	x_redundant = (x_size - (x_chain - 1) * diameter * chain_coeff) / 2
	y_redundant = x_redundant
	z_redundant = (z_size - (monomer_num - 1) * diameter * monomoer_coeff) / 2
	x_interval = diameter * chain_coeff
	y_interval = diameter * chain_coeff
	z_interval = diameter * monomoer_coeff
	
	# Print the header of the configuration file based on the initialization.
	print_header(polymer_name, monomer_num, chain_num, atom_types_num, bond_types_num, angle_types_num, dihedral_types_num, xlo,xhi,ylo,yhi,zlo,zhi, mass)

	# Compute and print the content of each part of content.
	monomer_list = generate_particles(xlo,xhi,ylo,yhi,zlo,zhi,x_redundant,y_redundant,z_redundant,x_interval,y_interval,z_interval,monomer_num,x_chain,y_chain)
	bond_boundle = generate_connections(monomer_num,x_chain,y_chain,"Bonds")
	angle_boundle = generate_connections(monomer_num,x_chain,y_chain,"Angles")
	dihedral_boundle = generate_connections(monomer_num,x_chain,y_chain,"Dihedrals")
	print_particles_connections(atom_types_num,bond_types_num,angle_types_num,dihedral_types_num,monomer_list,bond_boundle,angle_boundle,dihedral_boundle)
