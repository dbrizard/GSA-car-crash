

# [ ] Refine mesh
Mesh > EleEdit > Split/Merge

# [ ] Reduce output file size
* remove writting of d3dump: __d=nodump__ in expression
* no runrsf file: remove __*DATABASE_BINARY_RUNRSF__ keyword

# Check
* [ ] single vs double precision
* [ ] termination criterion (zero velocity?)

# LOG
## car6_crash.k
* original model from LS-OPT example

## car6_crash_v2.k
* model rewritten by ls-prepost (no more freeformat)
* 4CPU, 1s
* ~19.7MB

## car6_crash_v21.k
* main_v21.k: main file with parameters, calls the following files:
* car6_v21.k: rest of keywords

## car6_crash_v211.k
* coarse mesh
* try to get local square tube buckling
  - en passant de 5 à 2mm l'épaisseur du longeron avant, on voir apparaître du flanbement local, mais les mailles sont trop longues...
  - le longeron raffiné passe à travers le poteau...

## car6_crash_v3.k
* refined mesh
* 4CPU, 9s
* 130MB !!

## car6_crash_v22.k
* half refined mesh. Rails not refined!
* 4CPU, 3s
* 47MB
* [ ] check duplicate nodes

## car6_crash_v221.k
* half refined mesh (see v22)
* refined rails (to get local buckling)
* [x]  duplicate nodes removed
* /!\: bcp hourglass sur longeron

## car6_crash_v222.k
* shorter elements in the length of the right rail
* refined connection front/rear rail. 
* [x] il reste 2 éléments à diviser...
* [ ] better hourglass control is required ?
  - HE/TE=2.32e7/4.83e8 at the end
  - HE/IE=2.32e7/1.07e8 at the end
  - highest HE: body > bumper > hood
  - change value of QM on *HOURGLASS ?
  - IHQ=0 means 1=2=3 viscous type.   4=5 is stiffness
  - IHQ=0 et QM=0.05 sur les 3 premières parts => HE=2.01e7
  - IHQ=4 et QM=0.05 sur les 3 premières parts => HE=9.67e6
  - IHQ=4 et QM=0.05 sur toutes les parts => max HE=4e6


# POST
Y location of *RIGIDWALL_GEOMETRIC_CYLINDER_ID :
- center: 660.0
- left:   209.752
- right: 1187.31

# MODEL OUTPUT
## dimensions
1397 x 4318 x 1397
 t o t a l  m a s s          = 0.48067436E+01

capot attaché en trois points

initial velocity: -15640.0 vx

## Nodout
* 167: bumper, left
* 184: grill, center
* 432: roof, center, front
* 444: roof, center, rear 

## Rwforc
* forces:
  - 1: planar
  - 2: geometric cylinder
* transducer: cycle, x_force, y_force, z_force



# MODEL MODIFICATION
## longeron
Déformation locale avec initiateur




# notes CG 01/12/2023
- hide rigidwall: assembly and part select / SelPart menus

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

