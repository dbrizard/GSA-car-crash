# LS-DYNA simplified car model for global sensitivity analysis
## Model description and characteritics
The car is 1397mm wide, 1397mm high and 4318mm long. 
The total mass is around 500kg (depending on the thicknesses of the parts).

Initial velocity along x axis is -15640.0 mm/s (56km/h).

### List of parts
1. body
2. bumper
3. hood
4. grill
5. longeron (front)
6. longeron (rear)
7. roof
8. post
9. rigid cylinder
10. rigid plane

The hood is attached to the grill by 1 node and to the body by 2 nodes 
(`*CONSTRAINED_GENERALIZED_WELD_SPOT` without failure).

### Materials


## Origin of model
This simplified car model was taken from LS-OPT exemples 
("small car pole crash", `car6_crash.k`), see:

Stander, N., Roux, W., & Basudhar, A. (2018). LS-OPT training class. Optimization and robust design. Tutorial problems. Livermore Software Technology Corporation.



## Modifications
The mesh of the front end of the car was refined, which allows -in particular- 
the local buckling of the front longerons:

- bumper, grill, hood and front longeron were refined;
- front part of body was refined;


## Position of the post
Y location of `*RIGIDWALL_GEOMETRIC_CYLINDER_ID` :
- center: 660.0;
- left:   209.752 (centered on the axis of the left longeron);
- right: 1187.31  (centered on the axis of the right longeron).
