# FQHESphereJackGenerator — manual

Source: DiagHam wiki page `FQHESphereJackGenerator`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHESphereJackGenerator](../FQHEOnSphere/FQHESphereJackGenerator.md)

> **Options in this manual that the current program does not list**: `--add-index`, `--haldane`, `--max-range`, `--min-range`, `--use-gmp`. They may have been renamed, removed, or depend on a build option not enabled here (GMP, MPI). Trust `--help`.
FQHESphereJackGenerator allows to generate decomposition of Jack polynomials on symmetric monomial basis or decomposition of Jack polynomial time Van Der Monde on the unnormalized Slater basis.

## Basic usage

FQHESphereJackGenerator needs a text file that describes the polynomials including the number of particles (NbrParticles), the number of flux quanta (LzMax) and the root partition (ReferenceState). Such a text file looks like

```text
NbrParticles=8
LzMax=14
ReferenceState=1 0 1 0 1 0 1 0 1 0 1 0 1 0 1
```

write this file with the name laughlin_n_8_2s_14.dat and then run from the directory where you save this file

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereJackGenerator -a -2 -t bosons_unnormalized_haldane_laughlin_n_8_2s_14_lz_0.0.vec.txt --reference-file laughlin_n_8_2s_14.dat*

\$PATHTODIAGHAM should be replace by the path to the DiagHam directory. The -a option allows to set the alpha parameter of the Jack polynomial. Here it is set to -2 for the Laughlin state. It should be -3 for the Moore-Read state, -(k+1) for the Read-Rezayi k state, or -(k+1)/(r-1) for any (k,r) clustered state (using the bosonic value of r). The -t option sets the name of the file where the Jack polynomial decomposition will be stored. Using the -t option will write down the decomposition in a ASCII file, i.e. a file that can be open with any text editor. If you want to process the result with another DiagHam program, you should better store the result in binary mode using the -o option

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereJackGenerator -a -2 -o bosons_unnormalized_haldane_laughlin_n_8_2s_14_lz_0.0.vec --fermion --reference-file laughlin_n_8_2s_14.dat*

The binary vectors *(wiki page binary vectors)* can easily be manipulated using the tools provided with DiagHam. By default, the Jack generator expresses the result in the unnormalized basis and you convert it to the sphere geometry using [FQHESphereUnnormalizeState](FQHESphereUnnormalizeState.md). Actually, you can get directly the sphere geometry normalization by adding the -n option.

Fermionic counterparts of Jack polynomials (i.e. Jack polynomial times a Van der Monde determinant) can be computed in a similar way. The reference file should contains the fermionic root configuration, for example the Laughlin state at 1/3 should look like 

```text
NbrParticles=8
LzMax=21
ReferenceState=1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1
```

then when running the Jack generator, you should add the --fermion option

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereJackGenerator -a -2 -o fermions_unnormalized_haldane_laughlin_n_8_2s_21_lz_0.0.vec --reference-file laughlin_n_8_2s_21.dat*

assuming the reference file has been named laughlin_n_8_2s_21.dat. Notice that the alpha value is the same than the bosonic case.

## Advanced usage

**check singularity** :

The recursion formula used in the Jack generator may sometimes lead to 0/0 (this is discussed in [this paper](http://arxiv.org/abs/0902.4320)). States such as Laughlin, Moore-Read or the Read Rezayi series,... are safe. For some other clustered (k,r) states, this might happen. In that case the coefficients in a Jack expansion might just be crap. That why for these kin of states, it's safe to check for singular coefficients using the --check-singularity option. For example, the k=3 r=4 clustered for 12 particles using the following reference file named nonunitary_3_4_n_12_2s_12.dat

```text
NbrParticles=12
LzMax=12
ReferenceState= 3 0 0 0 3 0 0 0 3 0 0 0 3
```


will leads to 6 singular coefficients

'\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereJackGenerator --reference-file nonunitary_3_4_n_12_2s_12.dat -a -1.3333333333333333 -o test.vec --check-singularity

```text
memory requested for Hilbert space = 13Mo
partitions that may lead to singular coefficients :
[12,12,9,9,9,9,2,2,2,2,2,2] = 0 0 6 0 0 0 0 0 0 4 0 0 2
[11,11,11,11,8,8,2,2,2,2,2,2] = 0 0 6 0 0 0 0 0 2 0 0 4 0
[11,11,5,5,5,5,5,5,5,5,5,5] = 0 0 0 0 0 10 0 0 0 0 0 2 0
[10,10,10,10,10,10,4,4,1,1,1,1] = 0 4 0 0 2 0 0 0 0 0 6 0 0
[10,10,10,10,10,10,3,3,3,3] = 2 0 0 4 0 0 0 0 0 0 6 0 0
[7,7,7,7,7,7,7,7,7,7,1,1] = 0 2 0 0 0 0 0 10 0 0 0 0 0
```


**rational numbers** 

The Jack generator provides basic support of rational numbers. This allows to compute Jacks (only bosonic one, no huge mode) relying only on rational numbers. Typical usage is

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereJackGenerator --numerator-alpha -3 --denominator-alpha 2 --rational -t bosons_unnormalized_haldane_gaffnian_n_12_2s_15_lz_0.0.vec.txt --reference-file gaffnian_n_12_2s_15.dat*

--rational option should be turned on. The alpha coefficient is given through its numerator (--numerator-alpha option) and its denominator (--denominator-alpha option). It should be notice that no test on integer overflow is done so this option should be used with caution. In case of a singular coefficient is obtained at some point of the calculation, the code will stop and only output the Jack coefficients that have been computed, filling the other with zeroes.

The --check-singularity can also be used in rational mode. In that case, the singular coefficients are detected without any numerical error.

By default, the precision of rational depends on the machine/OS integer type (32bits or 64bits). Beware that no test is done to check integer overflow. The option --use-longlong allows to use larger integer numbers. If the GMP support is available, the --use-longlong option is replaced by the --use-gmp option that allows to use rational with arbitrary precision.

Singular coefficients may be fixed using the --use-symbolic option. This mode is highly CPU and memory demanding and should only be used if singular coefficients are present. It also requires arbitrary precision rational numbers even for moderate system sizes.

## Huge mode

In this mode, the Jack generator can handle Hilbert spaces that cannot fit into the computer memory or for Hilbert spaces dimension larger than 2^31.

;create a text file describing the state (bosonic N=20 Moore-Read state)

```text
NbrParticles=20
LzMax=18
ReferenceState=2 0 2 0 2 0 2 0 2 0 2 0 2 0 2 0 2 0 2
```


;compute and save the squeezed Hilbert space :
: time build/FQHE/src/Programs/FQHEOnSphere/FQHESphereJackGenerator -o bosons_unnormalized_haldane_pfaffian_n_20_2s_18_lz_0.0.vec -a -3 --reference-file pfaffian_n_20_2s_18.dat --symmetrized-basis --huge-basis --memory 1000 --file-size 1 --save-hilbert pfaffian_n_20_2s_18.hil

;compute the state itself
: time build/FQHE/src/Programs/FQHEOnSphere/FQHESphereJackGenerator -o bosons_unnormalized_haldane_pfaffian_n_20_2s_18_lz_0.0.vec -a -3 --huge-vector 2000 --reference-file pfaffian_n_20_2s_18.dat --symmetrized-basis --huge-basis --memory 1 --file-size 1 --load-hilbert pfaffian_n_20_2s_18.hil -S --processors 2 --huge-blocks 1000



Here is the typical output on a dual AMD Opteron 2376 with g++ 4.3.4

;computing and saving the squeezed Hilbert space :

total dimension = 386905330  
total dimension = 386905330  (7798)  
total requested disk space = 77 Mbytes splitted in 289 files  
largest file = 865 KBytes (110780 elements)  
disk storage is useless, keep everything in memory  
Haldane space dimension = 193181910  
final Hilbert space dimension 193181910  
Nbr Saved states = 193181910  
  
real    7m51.997s  
user    7m45.777s  
sys     0m3.788s  

;computing the state itself

TotalCount = 330236  
Factorization ratio = 0.225442%  
look up size 1048576  
Hilbert space factorization done in 16.4077s  
dimension = 193181910  
using disk storage for the Hilbert space and output state, some options might be disabled  
vector does not require temporary disk storage  
number of precalculation blocks = 12102  
183054853 / 193181910 (94%)  
  
real    157m52.731s  
user    230m30.528s  
sys     10m10.810s  ''

## Examples of large states

Typical examples of large state are the Laughlin state 1/3 (or 1/2) for more than 14 particles. 

Starting with the Laughlin 1/3 for N=15 particles, we use the root configuration file 

```text
NbrParticles = 15
LzMax = 42
ReferenceState = 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1
```

As indicated in the previous section, the Hilbert space generation and stats generation should be done in two consecutive steps, using e.g. the two following command lines:

'\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereJackGenerator -a -2  --fermion --huge-basis --memory 1000 --file-size 1 --save-hilbert laughlin3_n_15_2s_42.hil --reference-file root_laughlin3_n_15.txt -o fermions_haldane_unnormalized_laughlin3_n_15_2s_42_lz_0.0.vec'


'\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereJackGenerator -a -2  --fermion --huge-basis --memory 1000 --file-size 1 --load-hilbert laughlin3_n_15_2s_42.hil  --huge-blocks 4000 -S --processors 4 --reference-file root_laughlin3_n_15.txt -o fermions_haldane_unnormalized_laughlin3_n_15_2s_42_lz_0.0.vec'


On a dual Intel Xeon 6242, here is the total running time:

```text
total dimension = 1521967986
total dimension = 1521967986  (72486)
total requested disk space = 75 Mbytes splitted in 402 files 
largest file = 906 KBytes (116013 elements)
disk storage is useless, keep everything in memory
Haldane space dimension = 671906876
final Hilbert space dimension 671906876
Nbr Saved states = 671906876
TotalCount = 426068
Factorization ratio = 0.0981771%
look up size 1048576
using disk storage for the Hilbert space and output state, some options might be disabled
number of precalculation blocks = 38707
670134292 / 671906876 (99%)

real    1088m49.036s
user    793m17.042s
sys     686m52.188s
```


We can check that the state has been correctly generated by direct inspection of its component, either the first ones

'\$PATHTODIAGHAM/build/src/Programs/VectorBinary2Ascii -i fermions_haldane_unnormalized_laughlin3_n_15_2s_42_lz_0.0.vec --max-range 4 -s --add-index'
```text
0 1
1 -3
2 -3
3 6
```

or the last ones (which are more prone to fine accuracy issues)


'\$PATHTODIAGHAM/build/src/Programs/VectorBinary2Ascii -i fermions_haldane_unnormalized_laughlin3_n_15_2s_42_lz_0.0.vec --min-range -4 -s --add-index'
```text
671906872 -2.1582980275106e+15
671906873 7.1152682225625e+14
671906874 2.9884126534762e+15
671906875 -6.1902833536294e+15
```

Finally, the conversion to, e.g., the sphere normalization can be performed as follow

'\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereUnnormalizeState -i fermions_haldane_unnormalized_laughlin3_n_15_2s_42_lz_0.0.vec --haldane --reference-file root_laughlin3_n_15.txt --load-hilbert laughlin3_n_15_2s_42.hil -o fermions_haldane_laughlin3_n_15_2s_42_lz_0.0.vec --normalize --huge-basis'

```text
TotalCount = 426068
Factorization ratio = 0.0981771%
look up size 1048576
671906876 67190687676 (99%)           

real    2m20.657s
user    2m12.428s
sys     0m8.086s
```


For sake of comparison, a similar calculation for the bosonic Laughlin 1/2 state gives (on a dual AMD EPYC 7543)

'\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereJackGenerator -a -2   --huge-basis --memory 1000 --file-size 1 --save-hilbert laughlin2_n_15_2s_28.hil --reference-file root_laughlin2_n_15.txt -o bosons_haldane_unnormalized_laughlin2_n_15_2s_28_lz_0.0.vec'

'\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereJackGenerator -a -2   --huge-basis --memory 1000 --file-size 1 --load-hilbert laughlin2_n_15_2s_28.hil  --huge-blocks 4000 -S --processors 4 --reference-file root_laughlin2_n_15.txt -o bosons_haldane_unnormalized_laughlin2_n_15_2s_28_lz_0.0.vec'

```text
total dimension = 1521967986
total dimension = 1521967986  (72486)
total requested disk space = 75 Mbytes splitted in 402 files 
largest file = 906 KBytes (116013 elements)
disk storage is useless, keep everything in memory
Haldane space dimension = 671906876
final Hilbert space dimension 671906876
Nbr Saved states = 671906876
TotalCount = 426068
Factorization ratio = 0.0981771%
look up size 1048576
Hilbert space factorization done in 18.4612s
dimension = 671906876
using disk storage for the Hilbert space and output state, some options might be disabled
number of precalculation blocks = 57393
590631364 / 671906876 (87%)

real    2511m9.512s
user    1229m40.935s
sys     1743m8.276s
```


'\$PATHTODIAGHAM/build/src/Programs/VectorBinary2Ascii -i bosons_haldane_unnormalized_laughlin2_n_15_2s_28_lz_0.0.vec --min-range -4 -s --add-index'

```text
671906872 -12454041600
671906873 -12454041600
671906874 87178291200
671906875 -1307674368000
```

'\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereUnnormalizeState -i bosons_haldane_unnormalized_laughlin2_n_15_2s_28_lz_0.0.vec --haldane --reference-file root_laughlin2_n_15.txt  --load-hilbert laughlin2_n_15_2s_28.hil -o bosons_haldane_laughlin2_n_15_2s_28_lz_0.0.vec --normalize --huge-basis'

```text
TotalCount = 426068
Factorization ratio = 0.0981771%
look up size 1048576
Hilbert space factorization done in 16.2783s
671906876 67190687676 (99%)

real    10m45.511s
user    10m32.892s
sys     0m11.006s
```


Ramping up, we can also generate the N=16 Laughlin 1/2 state (on a dual Intel Xeon 6242)

'\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereJackGenerator -a -2   --huge-basis --memory 100000 --file-size 100 --save-hilbert laughlin2_n_16_2s_30.hil --reference-file root_laughlin2_n_16.txt -o bosons_haldane_unnormalized_laughlin2_n_16_2s_30_lz_0.0.vec --symmetrized-basis'

```text
total dimension = 9020077206
total dimension = 9020077206  (3797)
total requested disk space = 5492 Mbytes splitted in 237 files
largest file = 90654 KBytes (11603734 elements)
disk storage is useless, keep everything in memory
Haldane space dimension = 3846342253
final Hilbert space dimension 3846342253
Nbr Saved states = 3846342253

real    153m51.430s
user    152m55.547s
sys     0m22.618s
```

Generating the state itself requires up to 60Gb of RAM

'\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereJackGenerator -a -2   --huge-basis --memory 100000 --file-size 100 --load-hilbert laughlin2_n_16_2s_30.hil  --huge-blocks 40000 -S --processors 4 --reference-file root_laughlin2_n_16.txt -o bosons_haldane_unnormalized_laughlin2_n_16_2s_30_lz_0.0.vec --symmetrized-basis'

```text
memory requested for Hilbert space = 44029Mo
Hilbert space factorization done in 70.8907s
dimension = 3846342253
3846342144 / 3846342253 (99%)

real    3996m14.633s
user    3951m5.114s
sys     7m35.480s
```

A quick inspection of the last coefficients reveals

'\$PATHTODIAGHAM/build/src/Programs/VectorBinary2Ascii -i bosons_haldane_unnormalized_laughlin2_n_16_2s_30_lz_0.0.vec --min-range -4 -s  --add-index'

```text
3846342249 174356582400
3846342250 174356582400
3846342251 -1307674368000
3846342252 20922789888000
```

(note that up a sign, the last coefficient is always N!).
