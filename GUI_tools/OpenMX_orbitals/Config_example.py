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

# hdf5 file
hdf5File="PAO_and_AO_after_opt.hdf5"

# Z = 59, 61, 63, 64, 65, 68, 69, 70, 85, (86,) 87-118 do not exist in the database

database_list=[\
      [0,  "Kr10.0", "E"],\
      [1,  "H6.0",   "H_CA19"],   [1,  "H6.0",   "H_PBE19"],\
      [2,  "He8.0",  "He_CA19"],  [2,  "He8.0",  "He_PBE19"],\

      [3,  "Li8.0",  "Li_CA19"],  [3,  "Li8.0",  "Li_PBE19"],\
      [4,  "Be7.0",  "Be_CA19"],  [4,  "Be7.0",  "Be_PBE19"],\
      [5,  "B7.0",   "B_CA19"],   [5,  "B7.0",   "B_PBE19"],\
      [6,  "C6.0",   "C_CA19"],   [6,  "C6.0",   "C_PBE19"],\
      [7,  "N6.0",   "N_CA19"],   [7,  "N6.0",   "N_PBE19"],\
      [8,  "O6.0",   "O_CA19"],   [8,  "O6.0",   "O_PBE19"],\
      [9,  "F6.0",   "F_CA19"],   [9,  "F6.0",   "F_PBE19"],\
      [10, "Ne9.0",  "Ne_CA19"],  [10, "Ne9.0",  "Ne_PBE19"],\

      [11, "Na9.0",  "Na_CA19"],  [11, "Na9.0",  "Na_PBE19"],\
      [12, "Mg9.0",  "Mg_CA19"],  [12, "Mg9.0",  "Mg_PBE19"],\
      [13, "Al7.0",  "Al_CA19"],  [13, "Al7.0",  "Al_PBE19"],\
      [14, "Si7.0",  "Si_CA19"],  [14, "Si7.0",  "Si_PBE19"],\
      [15, "P7.0",   "P_CA19"],   [15, "P7.0",   "P_PBE19"],\
      [16, "S7.0",   "S_CA19"],   [16, "S7.0",   "S_PBE19"],\
      [17, "Cl7.0",  "Cl_CA19"],  [17, "Cl7.0",  "Cl_PBE19"],\
      [18, "Ar9.0",  "Ar_CA19"],  [18, "Ar9.0",  "Ar_PBE19"],\

      [10, "K10.0",  "K_CA19"],   [10, "K10.0",  "K_PBE19"],\
      [20, "Ca9.0",  "Ca_CA19"],  [20, "Ca9.0",  "Ca_PBE19"],\
      [21, "Sc9.0",  "Sc_CA19"],  [21, "Sc9.0",  "Sc_PBE19"],\
      [22, "Ti7.0",  "Ti_CA19"],  [22, "Ti7.0",  "Ti_PBE19"],\
      [23, "V6.0",   "V_CA19"],   [23, "V6.0",   "V_PBE19"],\
      [24, "Cr6.0",  "Cr_CA19"],  [24, "Cr6.0",  "Cr_PBE19"],\
      [25, "Mn6.0",  "Mb_CA19"],  [25, "Mn6.0",  "Mb_PBE19"],\
      [26, "Fe5.5H", "Fe_CA19H"], [26, "Fe5.5H", "Fe_PBE19H"],\
      [26, "Fe6.0S", "Fe_CA19S"], [26, "Fe6.0S", "Fe_PBE19S"],\
      [27, "Co6.0H", "Co_CA19H"], [27, "Co6.0H", "Co_PBE19H"],\
      [27, "Co6.0S", "Co_CA19S"], [27, "Co6.0S", "Co_PBE19S"],\
      [28, "Ni6.0H", "Ni_CA19H"], [28, "Ni6.0H", "Ni_PBE19H"],\
      [28, "Ni6.0S", "Ni_CA19S"], [28, "Ni6.0S", "Ni_PBE19S"],\
      [29, "Cu6.0H", "Cu_CA19H"], [29, "Cu6.0H", "Cu_PBE19H"],\
      [29, "Cu6.0S", "Cu_CA19S"], [29, "Cu6.0S", "Cu_PBE19S"],\
      [30, "Zn6.0H", "Zn_CA19H"], [30, "Zn6.0H", "Zn_PBE19H"],\
      [30, "Zn6.0S", "Zn_CA19S"], [30, "Zn6.0S", "Zn_PBE19S"],\
      [31, "Ga7.0",  "Ga_CA19"],  [31, "Ga7.0",  "Ga_PBE19"],\
      [32, "Ge7.0",  "Ge_CA19"],  [32, "Ge7.0",  "Ge_PBE19"],\
      [33, "As7.0",  "As_CA19"],  [33, "As7.0",  "As_PBE19"],\
      [34, "Se7.0",  "Se_CA19"],  [34, "Se7.0",  "Se_PBE19"],\
      [35, "Br7.0",  "Br_CA19"],  [35, "Br7.0",  "Br_PBE19"],\
      [36, "Kr10.0", "Kr_CA19"],  [36, "Kr10.0", "Kr_PBE19"],\

      [37, "Rb11.0", "Rb_CA19"],  [37, "Rb11.0", "Rb_PBE19"],\
      [38, "Sr10.0", "Sr_CA19"],  [38, "Sr10.0", "Sr_PBE19"],\
      [39, "Y10.0",  "Y1_CA19"],  [39, "Y10.0",  "Y1_PBE19"],\
      [40, "Zr7.0",  "Zr_CA19"],  [40, "Zr7.0",  "Zr_PBE19"],\
      [41, "Nb7.0",  "Nb_CA19"],  [41, "Nb7.0",  "Nb_PBE19"],\
      [42, "Mo7.0",  "Mo_CA19"],  [42, "Mo7.0",  "Mo_PBE19"],\
      [43, "Tc7.0",  "Tc_CA19"],  [43, "Tc7.0",  "Tc_PBE19"],\
      [44, "Ru7.0",  "Ru_CA19"],  [44, "Ru7.0",  "Ru_PBE19"],\
      [45, "Rh7.0",  "Rh_CA19"],  [45, "Rh7.0",  "Rh_PBE19"],\
      [46, "Pd7.0",  "Pd_CA19"],  [46, "Pd7.0",  "Pd_PBE19"],\
      [47, "Ag7.0",  "Ag_CA19"],  [47, "Ag7.0",  "Ag_PBE19"],\
      [48, "Cd7.0",  "Cd_CA19"],  [48, "Cd7.0",  "Cd_PBE19"],\
      [49, "In7.0",  "In_CA19"],  [49, "In7.0",  "In_PBE19"],\
      [50, "Sn7.0",  "Sn_CA19"],  [50, "Sn7.0",  "Sn_PBE19"],\
      [51, "Sb7.0",  "Sb_CA19"],  [51, "Sb7.0",  "Sb_PBE19"],\
      [52, "Te7.0",  "Te_CA19"],  [52, "Te7.0",  "Te_PBE19"],\
      [53, "I7.0",   "I_CA19"],   [53, "I7.0",   "I_PBE19"],\
      [54, "Xe11.0", "Xe_CA19"],  [54, "Xe11.0", "Xe_PBE19"],\

      [55, "Cs12.0", "Cs_CA19"],  [55, "Cs12.0", "Cs_PBE19"],\
      [56, "Ba10.0", "Ba_CA19"],  [56, "Ba10.0", "Ba_PBE19"],\
      [57, "La8.0",  "La_CA19"],  [57, "La8.0",  "La_PBE19"],\
      [58, "Ce8.0",  "Ce_CA19"],  [58, "Ce8.0",  "Ce_PBE19"],\
      
      [60, "Nd8.0",  "Nd_CA19"],  [60, "Nd8.0",  "Nd_PBE19"],\
      
      [62, "Sm8.0",  "Sm_CA19"],  [62, "Sm8.0",  "Sm_PBE19"],\
      
      
      
      [66, "Dy8.0",  "Dy_CA19"],  [66, "Dy8.0",  "Dy_PBE19"],\
      [67, "Ho8.0",  "Ho_CA19"],  [67, "Ho8.0",  "Ho_PBE19"],\
  
  
  
      [71, "Lu8.0",  "Lu_CA19"],  [71, "Lu8.0",  "Lu_PBE19"],\
      [72, "Hf9.0",  "Hf_CA19"],  [72, "Hf9.0",  "Hf_PBE19"],\
      [73, "Ta7.0",  "Ta_CA19"],  [73, "Ta7.0",  "Ta_PBE19"],\
      [74, "W7.0",   "W_CA19"],   [74, "W7.0",   "W_PBE19"],\ 
      [75, "Re7.0",  "Re_CA19"],  [75, "Re7.0",  "Re_PBE19"],\
      [76, "Os7.0",  "Os_CA19"],  [76, "Os7.0",  "Os_PBE19"],\
      [77, "Ir7.0",  "Ir_CA19"],  [77, "Ir7.0",  "Ir_PBE19"],\
      [78, "Pt7.0",  "Pt_CA19"],  [78, "Pt7.0",  "Pt_PBE19"],\
      [79, "Au7.0",  "Au_CA19"],  [79, "Au7.0",  "Au_PBE19"],\
      [80, "Hg8.0",  "Hg_CA19"],  [80, "Hg8.0",  "Hg_PBE19"],\
      [81, "Tl8.0",  "Tl_CA19"],  [81, "Tl8.0",  "Tl_PBE19"],\
      [82, "Pb8.0",  "Pb_CA19"],  [82, "Pb8.0",  "Pb_PBE19"],\
      [83, "Bi8.0",  "Bi_CA19"],  [83, "Bi8.0",  "Bi_PBE19"],\
      [84, "Po8.0",  "Po_CA19"],  [84, "Po8.0",  "Po_PBE19"],\  
]


