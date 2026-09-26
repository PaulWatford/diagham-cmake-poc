# FQHECylinderConfiningPotentialCoefficients — manual

Source: DiagHam wiki page `FQHECylinderConfiningPotentialCoefficients`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHECylinderConfiningPotentialCoefficients](../FQHEOnCylinder/FQHECylinderConfiningPotentialCoefficients.md)
[FQHECylinderConfiningPotentialCoefficients](FQHECylinderConfiningPotentialCoefficients.md) compute the one-body matrix for a polynomial confining potential on the cylinder geometry. Right now, this code is limited to the lowest Landau level. 


## y-Translation invariant potential

The confining potential is defined by $V_L(x)+V_R(x)$ where $x$ is the coordinate along the cylinder axis, the origin $x=0$ being centered. Here we have

*$V_R(x)= V_{0,R} \left(x-x_{0,R}\right)^{\alpha_R}$ if $x > x_{0,R}$ , $V_R(x)=0$ for $x < x_{0,R}$.
*$V_L(x)= V_{0,L} \left(x_{0,L}-x\right)^{\alpha_L}$ if $x < x_{0,L}$ , $V_L(x)=0$ for $x > x_{0,L}$.

The potential has no $y$ dependence, i.e. we assume it is translation invariant along the cylinder perimeter. The various parameters of $V_R(x)$ can be set through the following options

-  --confining-rightpower is used to define $\alpha_R$.
-  --confining-rightoffset is used to define $x_{0,R}$.
-  --confining-rightstrength is used to define $V_{0,R}$.

Similar options are available for $V_L(x)$.

A typical usage of this code is 
*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnCylinder/FQHECylinderConfiningPotentialCoefficients -s 24 --cylinder-perimeter 10.0 --confining-rightpower 2 --spinful --error 1e-7  --confining-rightoffset 4.5*

Here --error 1e-7 discards matrix elements that are lower than 1e-7. --spinful allows to generate a file compatible with spinful system. Since -confining-leftstrength is set to zero by default, we only have a confining potential to the right of the cylinder. The geometry of the cylinder is  set by -s (number of lux quanta) and --cylinder-perimeter (the cylinder perimeter). The program creates a files confining_cylinder_perimeter_10.000000_2s_24_alphar_2_x0r_4.500000_v0r_1.000000_alphal_2_x0l_0.000000_v0l_0.000000.dat that should contain

```text
# confining potential defined by :
# right alpha = 2, right V0 = 1, right shift = 4.5
# left alpha = 2, left V0 = 0, left shift = 0
# on a cylinder with perimeter L=10 and N_phi=24
OneBodyPotentialUpUp     = 0 0 0 0 0 0 0 0 0 0 0 0 0 1.3505752684037e-06 9.8927615977435e-05 0.0031684867734634 0.043495858744239 0.24578538577285 0.51848034728227 0.31121507153451 0.71242220021481 1.8256539068685 3.6791924321347 6.3153304987157 9.740519660135
OneBodyPotentialDownDown = 0 0 0 0 0 0 0 0 0 0 0 0 0 1.3505752684037e-06 9.8927615977435e-05 0.0031684867734634 0.043495858744239 0.24578538577285 0.51848034728227 0.31121507153451 0.71242220021481 1.8256539068685 3.6791924321347 6.3153304987157 9.740519660135
```

[FQHECylinderConfiningPotentialCoefficients](FQHECylinderConfiningPotentialCoefficients.md) can also be used to generate a confining potential directly defined in momentum/orbital space i.e. $V_m=  V_0 \left(m-m_0\right)^{\alpha}$ for the non-zero contribution. This mode is activated when using the --confining-momentum option. Note that $m_0$ is defined as $m_0=\frac{L_y}{2 \pi l_B^2}x_{0}$. For example, we can create a linear potential in momentum space to the right of a cylinder with the following command

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnCylinder/FQHECylinderConfiningPotentialCoefficients --confining-rightpower 1 --cylinder-perimeter 8.0 --nbr-flux 26 --spinful --confining-momentum*

The output file confining_momentum_cylinder_perimeter_8.000000_2s_26_alphar_1_x0r_0.000000_v0r_1.000000_alphal_2_x0l_0.000000_v0l_0.000000.dat should look like

```text
# confining potential in momentum space defined by :
# right alpha = 1, right V0 = 1, right shift = 0
# left alpha = 2, left V0 = 0, left shift = 0
# on a cylinder with perimeter L=8 and N_phi=26
OneBodyPotentialUpUp     = 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 2 3 4 5 6 7 8 9 10 11 12 13
OneBodyPotentialDownDown = 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 2 3 4 5 6 7 8 9 10 11 12 13
```


If we want to get a symmetric linear confining potential in momentum space, we just have to run 

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnCylinder/FQHECylinderConfiningPotentialCoefficients --confining-rightpower 1 --cylinder-perimeter 8.0 --nbr-flux 21 --spinful --confining-leftpower 1 --confining-leftstrength 1 --confining-momentum*

The output file should be

```text
# confining potential in momentum space defined by :
# right alpha = 1, right V0 = 1, right shift = 0
# left alpha = 1, left V0 = 1, left shift = 0
# on a cylinder with perimeter L=8 and N_phi=21
OneBodyPotentialUpUp     = 10.5 9.5 8.5 7.5 6.5 5.5 4.5 3.5 2.5 1.5 0.5 0.5 1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5 9.5 10.5
OneBodyPotentialDownDown = 10.5 9.5 8.5 7.5 6.5 5.5 4.5 3.5 2.5 1.5 0.5 0.5 1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5 9.5 10.5
```


The orbital index can also be shifted to mimic a flux insertion along the cylinder axis using the --flux-insertion option. The convention is that the flux insertion shifts the orbital index from $m$ to $m+\Phi/\Phi_0$ where $\Phi_0$ is the flux quantum. The value of $\Phi$ (in $\Phi_0$ units) is set by the option --flux-insertion . For example, 

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnCylinder/FQHECylinderConfiningPotentialCoefficients --confining-rightpower 1 --cylinder-perimeter 8.0 --nbr-flux 21 --spinful --confining-leftpower 1 --confining-leftstrength 1 --confining-momentum --flux-insertion 0.25* 

will create a file that should contain confining_momentum_cylinder_perimeter_8.000000_2s_21_alphar_1_x0r_0.000000_v0r_1.000000_alphal_1_x0l_0.000000_v0l_1.000000_flux_0.250000.dat that should contain 

```text
# confining potential in momentum space defined by :
# right alpha = 1, right V0 = 1, right shift = 0
# left alpha = 1, left V0 = 1, left shift = 0
# on a cylinder with perimeter L=8 and N_phi=21
OneBodyPotentialUpUp     = 10.25 9.25 8.25 7.25 6.25 5.25 4.25 3.25 2.25 1.25 0.25 0.75 1.75 2.75 3.75 4.75 5.75 6.75 7.75 8.75 9.75 10.75
OneBodyPotentialDownDown = 10.25 9.25 8.25 7.25 6.25 5.25 4.25 3.25 2.25 1.25 0.25 0.75 1.75 2.75 3.75 4.75 5.75 6.75 7.75 8.75 9.75 10.75
```


## Finite y extension potential

The confining potential can also be defined with the a finite extension along the cylinder perimeter. The new right and left confining potentials are defined as

*$V_R(x,y)= V_{0,R} \left(x-x_{0,R}\right)^{\alpha_R}$ if $x > x_{0,R}$ and $y < |l/2|$ , $V_R(x)=0$ {\rm} for $x < x_{0,R}$ or $y > |l/2|$.
*$V_L(x,y)= V_{0,L} \left(x_{0,L}-x\right)^{\alpha_L}$ if $x < x_{0,L}$ and $y < |l/2|$ , $V_L(x)=0$ {\rm} for $x > x_{0,L}$ or $y > |l/2|$.

The finite extension $l$ is set bye the option --y-extension. Since the finite extension leads to off diagonal hopping terms in the momentum basis. A truncation on the hopping term range is applied via the option --max-momentumtransfer. If this truncation parameter is set to a negative value, no truncation is performed. Due to these hopping terms, the output file has a different format. If we run the command 

*\$PATHTODIAGHAM/FQHECylinderConfiningPotentialCoefficients --cylinder-perimeter 8.0 -s 18 --spinful --confining-momentum --confining-rightpower 1 --confining-leftstrength 1.0 --confining-leftpower 1.0 --y-extension 7.0 --max-momentumtransfer 1*

we generate a potential file confining_momentum_cylinder_perimeter_8.000000_2s_18_l_7.000000_alphar_1_x0r_0.000000_v0r_1.000000_alphal_1_x0l_0.000000_v0l_1.000000.dat that contains

```text
# confining potential in momentum space defined by :
# right alpha = 1, right V0 = 1, right shift = 0
# left alpha = 1, left V0 = 1, left shift = 0
# with a finite extension 7 along the cylinder perimeter
# on a cylinder with perimeter L=8 and N_phi=18
0 0 9 9
0 1 0.88743192034094 0.88743192034094
1 0 0.88743192034094 0.88743192034094
1 1 8 8
1 2 0.78302816500672 0.78302816500672
2 1 0.78302816500672 0.78302816500672
2 2 7 7
2 3 0.67862440967249 0.67862440967249
3 2 0.67862440967249 0.67862440967249
...
```

The first column is the creation index, the second column is the annihilation index, the third column is the corresponding matrix element. If the option --spinful is used, a fourth column is added, which is just a copy of the third one.

[FQHECylinderConfiningPotentialCoefficients](FQHECylinderConfiningPotentialCoefficients.md) can also be used to generate pairing matrix elements with a superconducting region that can have a finite extension along the cylinder axis. For that purpose we can use the following command

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnCylinder/FQHECylinderConfiningPotentialCoefficients --cylinder-perimeter 8.0 -s 18  --confining-momentum --confining-rightpower  0 --confining-rightstrength 2.0 --confining-leftstrength 2.0 --confining-leftpower 0 --y-extension 7 --max-momentumtransfer 1 --confining-rightoffset 2.6 --confining-leftoffset -2.6 --confining-phase --confining-leftphase 0.25 --confining-rightphase 0.0*

We have a few new options in this example. --confining-phase allows to add Josephson phases to the left and right superconducting region. These phases are set by -confining-leftphase and --confining-rightphase where the phases should be provided in $pi$ units. The above command line will generate a file confining_momentum_cylinder_perimeter_8.000000_2s_18_l_7.000000_alphar_0_x0r_2.600000_v0r_2.000000_phir_0.000000_alphal_0_x0l_-2.600000_v0l_2.000000_phil_0.250000.dat (note that the file name contains the information about the Josephson phases via _phir_ and _phil_). It should look like

```text
# confining potential in momentum space defined by :
# right alpha = 0, right V0 = 2, right shift = 2.6
# left alpha = 0, left V0 = 2, left shift = -2.6
# with a finite extension 7 along the cylinder perimeter
# on a cylinder with perimeter L=8 and N_phi=18
0 0 2 0.25
0 1 0.20880751066846 0.25
1 0 0.20880751066846 0.25
1 1 2 0.25
1 2 0.20880751066846 0.25
2 1 0.20880751066846 0.25
2 2 2 0.25
2 3 0.20880751066846 0.25
3 2 0.20880751066846 0.25
3 3 2 0.25
3 4 0.20880751066846 0.25
4 3 0.20880751066846 0.25
4 4 2 0.25
4 5 0.20880751066846 0.25
5 4 0.20880751066846 0.25
5 5 2 0.25
5 6 0.20880751066846 0.25
6 5 0.20880751066846 0.25
6 6 0 0
6 7 0 0
7 6 0 0
7 7 0 0
7 8 0 0
8 7 0 0
8 8 0 0
8 9 0 0
9 8 0 0
9 9 0 0
9 10 0 0
10 9 0 0
10 10 0 0
10 11 0 0
11 10 0 0
11 11 0 0
11 12 0 0
12 11 0 0
12 12 0 0
12 13 0.20880751066846 0
13 12 0.20880751066846 0
13 13 2 0
13 14 0.20880751066846 0
14 13 0.20880751066846 0
14 14 2 0
14 15 0.20880751066846 0
15 14 0.20880751066846 0
15 15 2 0
15 16 0.20880751066846 0
16 15 0.20880751066846 0
16 16 2 0
16 17 0.20880751066846 0
17 16 0.20880751066846 0
17 17 2 0
17 18 0.20880751066846 0
18 17 0.20880751066846 0
18 18 2 0
```

Note that the right region where the confining potential is non zero can be shifted along the cylinder perimeter using the --y-rightshift option. In that case, $V_R(x,y)$ is non zero if $-l/2 + \delta < y < l/2 + \delta$ where $\delta$ is set by --y-rightshift and an additional phase is induced.
