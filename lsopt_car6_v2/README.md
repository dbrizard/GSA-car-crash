# LS-DYNA simplified car model for global sensitivity analysis
## Model description and characteritics
The car is 1397mm wide, 1397mm high and 4318mm long. 
The total mass is around 500kg (depending on the thicknesses of the parts).

Initial velocity along x axis is -15640.0 mm/s (56km/h). 50ms simulation time
Initial velocity along x axis is -8333.0 mm/s (30km/h). 80ms simulation time


### List of parts
1. body
2. bumper
3. hood
4. grill
5. longeron (front)
6. longeron (rear)
7. roof
8. post (rigid cylinder)
9. ground (rigid plane)

The hood is attached to the grill by 1 node and to the body by 2 nodes 
(`*CONSTRAINED_GENERALIZED_WELD_SPOT` without failure).

### Materials

### Model outputs

The final velocity is taken as the last *X-rigid body velocity* of the body (part #1) in MATSUM. 
Indeed, GLSTAT *global_x_velocity* is far from the final velocity.
Order of magnitude: 6000

The maximum force transmitted to the post (RWFORC, *x_force*). 
The post is RigidWall #2 (RigidWall #1 is the ground).
Order of magnitude: 5e5

GLSTAT, *internal_energy*. 
Order of magnitude: 5e7

The final *x_displacement* of the car, taken at node 605 (rear body).
Order of magnitude: 500


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

## Accelerate computations

- reduce DT for D3PLOT file writing does not change 
- NCPU=4 seems better than NCPU=8....
