# DFT_DATA10_mod.md

This note explains the modification of DFT_DATA19 to adjust our SPADExp package.

| Z | Element | PAO | VPS | radial.cutoff.pao | ocupied | others |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | H | H5.0 | H_CA19&rarr;H5.0_CA19, H_PBE19&rarr;H5.0_PBE19 | 6.0&rarr;5.0 | corrected | |
| 1 | H | H6.0 | H_CA19&rarr;H6.0_CA19, H_PBE19&rarr;H6.0_PBE19 | 6.0 | corrected | |
| 1 | H | H7.0 | H_CA19&rarr;H7.0_CA19, H_PBE19&rarr;H7.0_PBE19 | 6.0&rarr;7.0 | corrected | |
| 2 | He | He8.0 | He_CA19&rarr;He8.0_CA19, He_PBE19&rarr;He8.0_PBE19 | 8.0 | no error | |
| 2 | He | He10.0 | He_CA19&rarr;He10.0_CA19, He_PBE19&rarr;He10.0_PBE19 | 8.0&rarr;10.0 | no error | |
| 3 | Li | Li8.0 | Li_CA19, Li_PBE19 | 8.0 | no error | |
| 4 | Be | Be8.0 | Be_CA19, Be_PBE19 | 7.0 | no error | |
| 5 | B | B7.0 | B_CA19, B_PBE19 | 6.0&rarr;7.0 | no error | |
| 6 | C | C6.0 | C_CA19, C_PBE19 | 6.0 | corrected | |
| 7 | N | N6.0 | N_CA19, N_PBE19 | 5.0&rarr;6.0 | no error | |
| 8 | O | O6.0 | O_CA19, O_PBE19 | 6.0 | corrected | |
| 9 | F | F6.0 | F_CA19, F_PBE19 | 4.5&rarr;6.0 | no error | |
| 10 | Ne | Ne9.0 | Ne_CA19, Ne_PBE19 | 9.0 | no error | |
| 11 | Na | Na9.0 | Na_CA19, Na_PBE19 | 9.0 | no error | |
| 12 | Mg | Mg9.0 | Mg_CA19, Mg_PBE19 | 7.0&rarr;9.0 | no error | |
| 13 | Al | Al7.0 | Al_CA19, Al_PBE19 | 8.0&rarr;7.0 | no error | System.Name Al_XX19S&rarr;Al_XX19 |
| 14 | Si | Si7.0 | Si_CA19, Si_PBE19 | 7.5&rarr;7.0 | corrected | |
| 15 | P | P7.0 | P_CA19, P_PBE19 | 7.0 | no error | |
| 16 | S | S7.0 | S_CA19, S_PBE19 | 7.0 | no error | |
| 17 | Cl | Cl7.0 | Cl_CA19, Cl_PBE19 | 4.5&rarr;7.0 | no error | |
| 18 | Ar | Ar9.0 | Ar_CA19, Ar_PBE19 | 9.0 | no error | |
| 19 | K | K10.0 | K_CA19, K_PBE19 | 9.0&rarr;10.0 | no error | |
| 20 | Ca | Ca9.0 | Ca_CA19, Ca_PBE19 | 9.0 | no error | |
| 21 | Sc | Sc9.0 | Sc_CA19, Sc_PBE19 | 7.0&rarr;9.0 | no error | |
| 22 | Ti | Ti7.0 | Ti_CA19, Ti_PBE19 | 7.0 | no error | |
| 23 | V | V6.0 | V_CA19, V_PBE19 | 6.0 | no error | |
| 24 | Cr | Cr6.0 | Cr_CA19, Cr_PBE19 | 6.0 | no error | |
| 25 | Mn | Mn6.0 | Mn_CA19, Mn_PBE19 | 7.0&rarr;6.0 | no error | num.pao 5&rarr;15, search.UpperE 20&rarr;60|
| 26 | Fe(H) | Fe5.5H | Fe_CA19H, Fe_PBE19H | 6.0&rarr;5.5 | no error | One of **Input file** rows are removed in Fe5.5H.pao |
| 26 | Fe(S) | Fe6.0S | Fe_CA19S, Fe_PBE19S | 15.0&rarr;6.0 | no error | maxL.pao 3&rarr;4, num.pao 5&rarr;15, System.Name Fe_XX19&rarr;Fe_XX19S |
| 27 | Co(H) | Co6.0H | Co_CA19H, Co_PBE19H | 6.0 | no error | |
| 27 | Co(S) | Co6.0S | Co_CA19S, Co_PBE19S | 6.0 | no error | |
| 28 | Ni(H) | Ni6.0H | Ni_CA19H, Ni_PBE19H | 6.0 | no error | |
| 28 | Ni(S) | Ni6.0S | Ni_CA19S, Ni_PBE19S | 6.0 | no error | |
| 29 | Cu(H) | Cu6.0H | Cu_CA19H, Cu_PBE19H | 6.0 | no error | System.Name Cu_XX19&rarr;Cu_XX19H |
| 29 | Cu(S) | Cu6.0S | Cu_CA19S, Cu_PBE19S | 6.0 | no error | |
| 30 | Zn(H) | Zn6.0H | Zn_CA19H, Zn_PBE19H | 6.0 | no error | |
| 30 | Zn(S) | Zn6.0S | Zn_CA19S, Zn_PBE19S | 6.0 | no error | System.Name Zn_XX19&rarr;Zn_XX19S |
| 31 | Ga | Ga7.0 | Ga_CA19, Ga_PBE19 | 7.0 | no error | |
| 32 | Ge | Ge7.0 | Ge_CA19, Ge_PBE19 | 7.0 | no error | |
| 33 | As | As7.0 | As_CA19, As_PBE19 | 7.0 | no error | One of **Input file** rows are removed in As7.0.pao |
| 34 | Se | Se7.0 | Se_CA19, Se_PBE19 | 7.0 | no error | |
| 35 | Br | Br7.0 | Br_CA19, Br_PBE19 | 7.0 | no error | |
| 36 | Kr | Kr10.0 | Kr_CA19, Kr_PBE19 | 10.0 | no error | |
| 37 | Rb | Rb11.0 | Rb_CA19, Rb_PBE19 | 11.0 | no error | |
| 38 | Sr | Sr10.0 | Sr_CA19, Sr_PBE19 | 10.0 | no error | |
| 39 | Y | Y10.0 | Y_CA19, Y_PBE19 | 11.0&rarr;10.0 | no error | num.pao 5&rarr;15, search.UpperE 30&rarr;60 |
| 40 | Zr | Zr7.0 | Zr_CA19, Zr_PBE19 | 9.0&rarr;7.0 | no error | search.UpperE 30&rarr;50 |
| 41 | Nb | Nb7.0 | Nb_CA19, Nb_PBE19 | 11.0&rarr;7.0 | no error | search.UpperE 30&rarr;50 |
| 42 | Mo | Mo7.0 | Mo_CA19, Mo_PBE19 | 7.0 | no error | |
| 43 | Tc | Tc7.0 | Tc_CA19, Tc_PBE19 | 7.0 | no error | |
| 44 | Ru | Ru7.0 | Ru_CA19, Ru_PBE19 | 7.0 | no error | |
| 45 | Rh | Rh7.0 | Rh_CA19, Rh_PBE19 | 7.0 | no error | |
| 46 | Pd | Pd7.0 | Pd_CA19, Pd_PBE19 | 7.0 | no error | |
| 47 | Ag | Ag7.0 | Ag_CA19, Ag_PBE19 | 7.0 | no error | maxL.pao 3&rarr;4 |
| 48 | Cd | Cd7.0 | Cd_CA19, Cd_PBE19 | 7.0 | no error | |
| 49 | In | In7.0 | In_CA19, In_PBE19 | 7.0 | no error | |
| 50 | Sn | Sn7.0 | Sn_CA19, Sn_PBE19 | 7.0 | no error | |
| 51 | Sb | Sb7.0 | Sb_CA19, Sb_PBE19 | 7.0 | no error | |
| 52 | Te | Te7.0 | Te_CA19, Te_PBE19 | 7.0 | no error | |
| 53 | I | I7.0 | I_CA19, I_PBE19 | 7.0 | no error | |
| 54 | Xe | Xe11.0 | Xe_CA19, Xe_PBE19 | 7.0&rarr;11.0 | no error | |
| 55 | Cs | Cs12.0 | Cs_CA19, Cs_PBE19 | 12.0 | no error | |
| 56 | Ba | Ba10.0 | Ba_CA19, Ba_PBE19 | 10.0 | no error | |
| 57 | La | La10.0 | La_CA19, La_PBE19 | 8.0 | no error | |
| 58 | Ce | Ce8.0 | Ce_CA19, Ce_PBE19 | 8.0 | no error | |
| 60 | Nd | Nd8.0 | Nd_CA19, Nd_PBE19 | 8.0 | no error | |
| 62 | Sm | Sm8.0 | Sm_CA19, Sm_PBE19 | 8.0 | no error | One of **Input file** rows are removed in Sm8.0.pao |
| 66 | Dy | Dy8.0 | Dy_CA19, Dy_PBE19 | 8.0 | no error | |
| 67 | Ho | Ho8.0 | Ho_CA19, Ho_PBE19 | 8.0 | no error | |
| 71 | Lu | Lu8.0 | Lu_CA19, Lu_PBE19 | 8.0 | no error | |
| 72 | Hf | Hf9.0 | Hf_CA19, Hf_PBE19 | 7.0&rarr;9.0 | no error | search.UpperE 7&rarr;50 |
| 73 | Ta | Ta7.0 | Ta_CA19, Ta_PBE19 | 7.0 | no error | |
| 74 | W | W7.0 | W_CA19, W_PBE19 | 7.0 | no error | |
| 75 | Re | Re7.0 | Re_CA19, Re_PBE19 | 7.0 | no error | |
| 76 | Os | Os7.0 | Os_CA19, Os_PBE19 | 7.0 | no error | Os_PBE19.pao cannot be generated due to some errors (not yet resolved) |
| 77 | Ir | Ir7.0 | Ir_CA19, Ir_PBE19 | 7.0 | no error | Ir_PBE19.pao cannot be generated due to some errors (not yet resolved) |
| 78 | Pt | Pt7.0 | Pt_CA19, Pt_PBE19 | 7.0 | no error | search.UpperE 7&rarr;60 |
| 79 | Au | Au7.0 | Au_CA19, Au_PBE19 | 7.0 | no error | Au_CA19.pao and Au_PBE19.pao cannot be generated due to some errors (not yet resolved) |
| 80 | Hg | Hg8.0 | Hg_CA19, Hg_PBE19 | 8.0 | no error | |
| 81 | Tl | Tl8.0 | Tl_CA19, Tl_PBE19 | 8.0 | no error | |
| 82 | Pb | Pb8.0 | Pb_CA19, Pb_PBE19 | 8.0 | no error | |
| 83 | Bi | Bi8.0 | Bi_CA19, Bi_PBE19 | 8.0 | no error | |

