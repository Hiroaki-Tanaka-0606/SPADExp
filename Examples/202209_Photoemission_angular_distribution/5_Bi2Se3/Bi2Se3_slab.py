# Output fractional coordinates for Bi2Se3 slab

BiX=[0.00000, 0.66667, 0.33333, 0.00000, 0.66667, 0.33333]
BiY=[0.00000, 0.33333, 0.66667, 0.00000, 0.33333, 0.66667]
BiZ=[0.40046, 0.73379, 1.06713, 0.59954, 0.93287, 0.26621]

SeX=[0.00000, 0.66667, 0.33333, 0.00000, 0.66667, 0.33333, 0.00000, 0.66667, 0.33333]
SeY=[0.00000, 0.33333, 0.66667, 0.00000, 0.33333, 0.66667, 0.00000, 0.33333, 0.66667]
SeZ=[0.20970, 0.54303, 0.87637, 0.79030, 1.12363, 0.45697, 1.00000, 0.33333, 0.66667]

index=1
num_layer=3
num_space=1
for i in range(num_space, num_layer+num_space):
    for j in range(0, 6):
        print(("{0:-3d} Bi {1:10.8f} {2:10.8f} {3:10.8f} 7.5 7.5 45 0 45 0 0 off").format(index, (BiZ[j]+i-0.2)/(num_layer+num_space*2), BiX[j]+0.05, BiY[j]+0.05))
        index+=1
    for j in range(0, 9):
        print(("{0:-3d} Se {1:10.8f} {2:10.8f} {3:10.8f} 3.0 3.0 45 0 45 0 0 off").format(index, (SeZ[j]+i-0.2)/(num_layer+num_space*2), SeX[j]+0.05, SeY[j]+0.05))
        index+=1
