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

## base design output 2024-05-03
{'fmax': 772836.8125,
 'dmax': -641.872802734375,
 'vfin': -710.947021484375,
 'IE': array([1.20394832e+08])}

{'fmax': 775579.75,
 'dmax': -643.6537475585938,
 'vfin': -716.1944580078125,
 'IE': array([1.2032208e+08])}

{'fmax': 772913.375,
 'dmax': -642.1265869140625,
 'vfin': -693.0195922851562,
 'IE': array([1.20379344e+08])}

{'fmax': 770100.4375,
 'dmax': -643.90771484375,
 'vfin': -730.5367431640625,
 'IE': array([1.2035144e+08])}

{'fmax': 773196.375,
 'dmax': -642.2013549804688,
 'vfin': -706.6123657226562,
 'IE': array([1.20424064e+08])}

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

# Conputation time
## NCPU=1
 Total CPU time     =        26 seconds (   0 hours  0 minutes 26 seconds)
 Elapsed time      26 seconds for   31280 cycles using  1 SMP thread

## NCPU=4
 Total CPU time     =        66 seconds (   0 hours  1 minutes  6 seconds)
 Elapsed time      16 seconds for   34121 cycles using  4 SMP threads

## NCPU=8
 Total CPU time     =       169 seconds (   0 hours  2 minutes 49 seconds)
 Elapsed time      21 seconds for   33707 cycles using  8 SMP threads
