Model v223 from v222
Adding variation of yield stress

# 2024-03-27
Contact problem between front rail and rigid cylinder:
the front rail penetrates the cylinder
- not possible to modify contact properties with cylinder rigidwall
- [ ] check `CONTACT_AUTOMATIC_SINGLE_SURFACE_ID`
  - slave segment part set ID 1 contains only 1, 2, 3, 4; not 5 (front rail!)
  - [x] part 5 added to part set 1 => much better tube crushing
- [ ] check model outputs (forces and displacements)
- [ ] reconsider endtime: 0.080 => 0.200s


# MODEL OUTPUT
## Nodout
* 167: bumper, left
* 184: grill, center
* 432: roof, center, front
* 444: roof, center, rear
* 605: body, center, rear
