hn_min=7
hn_max=40

depth=4
c=5.695

fn_template="TiS2_l40_PAD_Nu_p_%deV_Sqrt_l%d.dat"

for hn in range(hn_min, hn_max+1):
    fn=fn_template % (hn, depth)
    with open(fn, "w") as f:
        f.write("&Control\n")
        f.write("Calculation PAD\n")
        f.write("Log_file TiS2_l40_PAD_Nu_p_%deV_Sqrt_l%d.log\n" % (hn, depth))
        f.write("Console_log true\n")
        f.write("Output_file TiS2_l40_PAD_Nu_p_%deV_Sqrt_l%d.hdf5\n" % (hn, depth))
        f.write("/\n")
        f.write("&PAD\n")
        f.write("Input_file /path/to/hdf5/file\n")
        f.write("E_min -2\n")
        f.write("E_max 0.5\n")
        f.write("E_pixel 0.01\n")
        f.write("dE 0.1\n")
        f.write("Initial_state PAO\n")
        f.write("Final_state FP_PAO\n")
        f.write("FPFS_energy_step 0.05\n")
        f.write("Excitation_energy %d\n" % hn)
        f.write("VPS_file /path/to/Pseudopotentials.hdf5\n")
        f.write("Polarization Linear\n")
        f.write("Theta 40\n")
        f.write("Phi 0\n")
        f.write("Ignore_nonlocal true\n")
        f.write("Atomic_orbitals_file /path/to/PAO_and_AO_after_opt.hdf5\n")
        f.write("Weighting True\n")
        f.write("Weighting_shape Sqrt\n")
        f.write("Weighting_axis 0 0 1\n")
        f.write("Weighting_origin 244.885\n")
        f.write("Weighting_width %.3f\n" % (-depth*c,))
        f.write("FPFS_Numerov true\n")
        f.write("/\n")

with open("SPADExp.sh", "w") as f:
    f.write("#!/bin/sh\n")
    f.write("\n")
    f.write("#SBATCH -p F1cpu\n")
    f.write("#SBATCH -N 1\n")
    f.write("#SBATCH -n 1\n")
    f.write("#SBATCH -c 128\n")
    f.write("#SBATCH --mail-type=ALL\n")           
    f.write("#SBATCH --mail-user=your@mail.address\n")
    f.write("\n")
    f.write("ulimit -s unlimited\n")
    for hn in range(hn_min, hn_max+1):
        fn=fn_template % (hn, depth)
        f.write("srun ../../SPADExp.o < %s\n" % fn)
