# Config for SPADExp_GUI and SPADExp_Viewer

# For Gui
fontFamilies=["Segoe UI", "Yu Gothic UI"]
# font sizes are in unit of pixel
fontSize_normal=16
fontSize_large=24
ContentsMargins=[5,5,5,5]

tickLength=-30

sigma_max=5

Eh=27.2114 # (eV)
au_ang=0.529177 # (Ang)

# Pen for vLine and hLine
pen1=(0, 255, 0)
pen2=(255, 0, 0)
pen3=(255, 255, 0)
gridAlpha=50 # 100 =max

plot3D_minWidth=400
plot3D_minHeight=400

PAO_and_AO="/path/to/PAO_and_AO_after_opt.hdf5"
pen_AO=(0, 0, 255)
pen_PAO=(0, 128, 255)
pen_finalp1=(255, 128, 0)
pen_finalm1=(128, 255, 0)

# for realSpace drawing
elements_file="/path/to/VESTA/VESTA-win64/elements.ini"
radius_index=0
not_found_element=96
radius_coeff=0.3
reciprocal_coeff=1
pen_kx=(1, 0, 0)
pen_ky=(0, 0.5, 1)
reciprocal_axis_width=3
polarization_length=5
pen_pol=(0, 1, 0)
polarization_width=10
Weighting_offset=0.2
