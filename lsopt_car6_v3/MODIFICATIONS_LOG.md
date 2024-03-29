Model v223 from v222
Adding variation of yield stress

# 2024-03-27
Contact problem between front rail and rigid cylinder:
the front rail penetrates the cylinder
- not possible to modify contact properties with cylinder rigidwall
- [ ] check `CONTACT_AUTOMATIC_SINGLE_SURFACE_ID`
  - slave segment part set ID 1 contains only 1, 2, 3, 4; not 5 (front rail!)
  - [x] part 5 added to part set 1 => much better tube crushing
- [x] reconsider endtime: 0.080 => 0.200s
- [ ] check model outputs for base design
  - 'fmax': 534474.5 N = 5e5 N
  - 'dmax': -1158 mm
  - 'vfin': -4242 mm/s  (initial 8333 mm/s = 30km/h)
  - 'IE': 94043328 N.mm = 94e6 N.mm


# MODEL OUTPUT
## Nodout
* 167: bumper, left
* 184: grill, center
* 432: roof, center, front
* 444: roof, center, rear
* 605: body, center, rear

## Units
MASS  LENGTH   TIME  FORCE  STRESS  ENERGY  DENSITY  YOUNG's   35 mph     GRAVITY
                                                               56.33 kph
tonne mm       s     N      MPa     N-mm    7.83e-9  2.07e+05  1.56e+04   9.806e+03
