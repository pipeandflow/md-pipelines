def bench_bosons_single(nbosons, boson_positions, socket_id, workdir):
    assert len(boson_positions) == nbosons

    boson_positions_config = str(flatten(boson_positions))

    boson_labels = str(["E"] * nbosons)
    boson_masses = str(["1.0"] * nbosons)
    bosons_list = str(list(range(nbosons)))

    config = ipi_config(boson_positions_config, boson_masses, boson_labels, bosons_list, socket_id, workdir)
    return config

def lammps_config(num_bosons, num_beads, num_steps, seed):
    num_beads = float(num_beads)
    return f"""
#Units and dims
units electron
dimension 2
boundary p p p 
atom_style atomic
atom_modify sort 0 0.0 map yes

#User-defined parameters

#Time step is in femtoseconds
timestep 1.0
#Define temperature in Kelvin 
variable Temp equal 5.8
#Define force constant
variable k equal 1.21647924e-8  # Hartree / Bohr ^2
#Number of beads
variable Nbeads equal {num_beads}
variable ibead uloop ${{Nbeads}} pad
variable seed equal {seed}

#Create box and atoms. All distances are in Bohr
region box block -1500 1500 -1500 1500 -1500 1500
create_box 1 box

variable a loop {num_bosons}
label loop
create_atoms 1 single 0.0 0.0 0.0  # X Y Z coordinates  are in units of Bohr
next a
jump SELF loop

#Mass is in amu. This is the electron mass in amu multiplied by 0.84
mass 1 0.00054858

variable gauss_height equal 0.1895659039436153
variable gauss_sigma equal 33.67207052340757

#Non-interacting particles
pair_style      none

#Initialize velocities
velocity all create ${{Temp}} ${{seed}}${{ibead}} mom yes rot yes dist gaussian

#Define variables for force constants
variable fx atom -v_k*x
variable fy atom -v_k*y
variable harm2d atom 0.5*v_k*(x^2+y^2)
variable trapvir atom -0.5*(x*v_fx+y*v_fy)/v_Nbeads
compute trapvir all reduce sum v_trapvir
compute trap all reduce sum v_harm2d

#Add harmonic external force
#fix harm all addforce v_fx v_fy v_fz energy v_harm3d
fix harm all addforce v_fx v_fy 0.0 energy v_harm2d
#Add harmonic potential energy to total energy and potential energy
fix_modify harm energy yes

dump xyz all xyz 1000 system_${{ibead}}.xyz

variable nhc equal 4
fix 12 all pimdb method nmpimd fmass 1.0 sp 1.0 temp ${{Temp}} nhc ${{nhc}}
variable newvir equal f_12[3]

#print thermodynamics info to log file every n steps
thermo_style custom step time temp ke pe c_trap etotal fmax fnorm v_newvir c_trapvir
thermo 100

fix 5 all enforce2d

run             {num_steps}
    """

if __name__ == '__main__':
    num_bosons = int(snakemake.params.num_bosons)
    num_beads = int(snakemake.params.num_beads)
    seed = int(snakemake.params.seed)
    num_steps = int(snakemake.params.num_steps)
    config = lammps_config(num_bosons=num_bosons, num_beads=num_beads, num_steps=num_steps, seed=seed)

    output_path = snakemake.output[0]
    with open(output_path, 'w') as f:
        f.write(config)
