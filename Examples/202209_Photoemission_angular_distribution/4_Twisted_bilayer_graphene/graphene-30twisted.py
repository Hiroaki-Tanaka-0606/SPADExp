# output fractional coordinates of 30deg-twisted bilayer graphene

x=[0.6000, 0.9333, 0.1000, 0.4333]
y=[0.1000, 0.1000, 0.6000, 0.6000]

x2=[0.1000, 0.4333, 0.6000, 0.9333]

# case 1: 4x7 unit cell
index=1
for i in range(0, 4):
    for j in range(0, 7):
        for k in range(0, 4):
            print(("{0:-3d} C {1:6.4f} {2:6.4f} {3:6.4f} 2.0 2.0").format(index, (x[k]+i)/4, (y[k]+j)/7, 0.5))
            index+=1
            print(("{0:-3d} C {1:6.4f} {2:6.4f} {3:6.4f} 2.0 2.0").format(index, (y[k]+j)/7, (x[k]+i)/4, 0.6677))
            index+=1
            
print("")
print("")
            
# case 2: 4x7 unit cell, top layer is shifted by (a/2, 0)
index=1
for i in range(0, 4):
    for j in range(0, 7):
        for k in range(0, 4):
            print(("{0:-3d} C {1:6.4f} {2:6.4f} {3:6.4f} 2.0 2.0").format(index, (x[k]+i)/4, (y[k]+j)/7, 0.5))
            index+=1
            print(("{0:-3d} C {1:6.4f} {2:6.4f} {3:6.4f} 2.0 2.0").format(index, (y[k]+j)/7, (x2[k]+i)/4, 0.6677))
            index+=1

print("")
print("")
            
# case 3: 7x12 unit cell
index=1
for i in range(0,7):
    for j in range(0, 12):
        for k in range(0, 4):
            print(("{0:-3d} C {1:6.4f} {2:6.4f} {3:6.4f} 2.0 2.0").format(index, (x[k]+i)/7, (y[k]+j)/12, 0.5))
            index+=1
            print(("{0:-3d} C {1:6.4f} {2:6.4f} {3:6.4f} 2.0 2.0").format(index, (y[k]+j)/12, (x[k]+i)/7, 0.6677))
            index+=1
