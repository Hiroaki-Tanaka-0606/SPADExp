# Config for OpenMX_orbitals


from pyqtgraph.Qt import QtGui, QtCore, QtWidgets

# Working directory
workingDirectory="/home/hiroaki/OpenMX/openmx3.9/DFT_DATA19"

# Element symbols
el_symbol=[ \
      "E" , "H" , "He", "Li", "Be", "B" , "C" , "N" , "O" , "F" , "Ne",\
            "Na", "Mg", "Al", "Si", "P" , "S" , "Cl", "Ar", "K" , "Ca",\
            "Sc", "Ti", "V" , "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",\
            "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y" , "Zr",\
            "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",\
            "Sb", "Te", "I" , "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",\
            "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",\
            "Lu", "Hf", "Ta", "W" , "Re", "Os", "Ir", "Pt", "Au", "Hg",\
            "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",\
            "Pa", "U" , "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",\
            "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",\
            "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"]

azimuthal=["s", "p", "d", "f", "g"]

analysisTypes=["---- select one from the following ----",\
               "PAO: Before and after optimization",\
               "PAO: Before optimization and from VPS",\
               "PAO: After optimization and from VPS",
               "PAO and AO from VPS",
               "PAO and AO after optimization"]

OrthNorm_threshold=0.01
Contraction_threshold=0.01

search_UpperE_coeff=2.0
num_of_partition_coeff=2.0

invert_criterion=-0.2

# (modified) adpack
adpack="/home/hiroaki/adpack2.2/work/adpack"

# For Gui
fontFamilies=["Segoe UI", "Yu Gothic UI"]
# For table
unselected_item=QtGui.QBrush(QtGui.QColor(255,255,255))
selected_item=QtGui.QBrush(QtGui.QColor(0, 170, 255))
error_item=QtGui.QBrush(QtGui.QColor(255, 0, 153))
# For graph
pen_after=(0, 170, 255)
pen_before=(0, 0, 255)
pen_fromVPS=(128, 212, 255)
pen_reproduced=(41, 204, 41)
# font sizes are in unit of pixel
fontSize_normal=16
fontSize_large=24
ContentsMargins=[5,5,5,5]
