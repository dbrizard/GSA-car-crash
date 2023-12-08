
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




# MODEL OUTPUT
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
- [ ] Déformation locale avec initiateur

