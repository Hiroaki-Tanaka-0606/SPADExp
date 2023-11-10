hn_list=[25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40]
kR_list=[2.688,2.636,2.587,2.540,2.496,2.454,2.414,2.376,2.340,2.305,2.272,2.240,2.210,2.181,2.152,2.125]
hn_max=40

depth=4
c=5.695

fn_template="TiS2_l40_PAD_Bulk_p_%deV_Sqrt_l%d.dat"

for i, hn in enumerate(hn_list):
    kR=kR_list[i]
    fn=fn_template % (hn, depth)
    with open(fn, "w") as f:
        f.write("&Control\n")
        f.write("Calculation PAD\n")
        f.write("Log_file TiS2_l40_PAD_Bulk_p_%deV_Sqrt_l%d.log\n" % (hn, depth))
        f.write("Console_log true\n")
        f.write("Output_file TiS2_l40_PAD_Bulk_p_%deV_Sqrt_l%d.hdf5\n" % (hn, depth))
        f.write("/\n")
        f.write("&PAD\n")
        f.write("Input_file /path/to/hdf5/file\n")
        f.write("E_min -1\n")
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
        f.write("FPFS_bulk 22.78 239.19 38\n")
        f.write("FPFS_kRange %.3f\n" % kR)
        f.write("/\n")

shell_template="SPADExp_%d.sh"
with open("SPADExp_le.sh", "w") as f:
    f.write("#!/bin/sh\n")
    f.write("\n")
    f.write("#SBATCH -p L1cpu\n")
    f.write("#SBATCH -N 1\n")
    f.write("#SBATCH -n 1\n")
    f.write("#SBATCH -c 128\n")
    f.write("#SBATCH --mail-type=ALL\n")           
    f.write("#SBATCH --mail-user=your@mail.address\n")
    f.write("\n")
    f.write("ulimit -s unlimited\n")
    for hn in reversed(hn_list):
        fn=fn_template % (hn, depth)
        f.write("srun ../../SPADExp.o < %s\n" % fn)
