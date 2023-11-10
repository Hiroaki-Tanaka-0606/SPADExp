hn=7

depth=4
c=5.695

fn_template="TiS2_l40_PAD_Nu_%ddeg_%deV_Sqrt_l%d.dat"

angles=[0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180]
thetas=[90.0000,96.4086,102.7000,108.7472,114.4045,119.4987,
        123.8258,127.1586,129.2735,130.0000,129.2735,127.1586,
        123.8258,119.4987,114.4045,108.7472,102.7000,96.4086,90.0000]
phis=[90.0000,97.6926,105.5794,113.8587,122.7324,132.3941,
      142.9955,154.5862,167.0375,180.0000,192.9625,205.4138,
      217.0045,227.6059,237.2676,246.1413,254.4206,262.3074,270.0000]

for i, angle in enumerate(angles):
    theta=thetas[i]
    phi=phis[i]
    fn=fn_template % (angle, hn, depth)
    with open(fn, "w") as f:
        f.write("&Control\n")
        f.write("Calculation PAD\n")
        f.write("Log_file TiS2_l40_PAD_Nu_%ddeg_%deV_Sqrt_l%d.log\n" % (angle, hn, depth))
        f.write("Console_log true\n")
        f.write("Output_file TiS2_l40_PAD_Nu_%ddeg_%deV_Sqrt_l%d.hdf5\n" % (angle, hn, depth))
        f.write("/\n")
        f.write("&PAD\n")
        f.write("Input_file /path/to/hdf5/file\n")
        f.write("E_min -1\n")
        f.write("E_max 0.5\n")
        f.write("E_pixel 0.002\n")
        f.write("dE 0.06\n")
        f.write("Initial_state PAO\n")
        f.write("Final_state FP_PAO\n")
        f.write("FPFS_energy_step 0.05\n")
        f.write("Excitation_energy %d\n" % hn)
        f.write("VPS_file /path/to/Pseudopotentials.hdf5\n")
        f.write("Polarization Linear\n")
        f.write("Theta %.4f\n" % theta)
        f.write("Phi %.4f\n" % phi)
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
    for angle in angles:
        fn=fn_template % (angle, hn, depth)
        f.write("srun ../../SPADExp.o < %s\n" % fn)
