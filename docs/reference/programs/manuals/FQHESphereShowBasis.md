# FQHESphereShowBasis — manual

Source: DiagHam wiki page `FQHESphereShowBasis`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHESphereShowBasis](../FQHEOnSphere/FQHESphereShowBasis.md)
FQHESphereShowBasis displays the n-body basis used in DiagHam. For example,

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereShowBasis -p 7 -l 8 -z 0 --boson*

shows the n-body basis for 7 bosons with 8 flux quanta and a total Lz value of zero. Display style uses occupation numbers, leftmost occupation being the corresponding to the one body state with the smallest Lz value). Typical output will look like :

```text
memory requested for Hilbert space = 2Mo                                                                                                          
3 0 0 0 1 0 0 0 3                                                                                                                                 
2 1 0 1 0 0 0 0 3                                                                                                                                 
2 0 2 0 0 0 0 0 3                                                                                                                                 
1 2 1 0 0 0 0 0 3                                                                                                                                 
0 4 0 0 0 0 0 0 3
```

For fermions, you can just remove the --boson i.e. 

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereShowBasis -p 6 -l 15 -z 0*

```text
memory requested for Hilbert space = 4Mo                                                                                           
1 1 1 0 0 0 0 0 0 0 0 0 0 1 1 1                                                                                                    
1 1 0 1 0 0 0 0 0 0 0 0 1 0 1 1                                                                                                    
1 1 0 0 1 0 0 0 0 0 0 1 0 0 1 1                                                                                                    
1 0 1 1 0 0 0 0 0 0 0 1 0 0 1 1                                                                                                    
1 1 0 0 0 1 0 0 0 0 1 0 0 0 1 1                                                                                                    
1 0 1 0 1 0 0 0 0 0 1 0 0 0 1 1                                                                                                    
0 1 1 1 0 0 0 0 0 0 1 0 0 0 1 1
```

This code also handles spinful cases. For example, to consider the case of SU(2) spinfule fermions with a pairing term (i.e. without particle conservation), we just have to type

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereShowBasis --fermion --su2-spin --use-pairing  -l 2 -z 0 -s 0*

```text
Hilbert space dimension = 8
memory requested for Hilbert space = 96
memory requested for lookup table = 12
X X X 
X 0 X 
d X u 
d 0 u 
u X d 
u 0 d 
0 X 0 
0 0 0
```

## Displaying components of a vector in the Fock basis

The code above allows to display all the components of a vector in the given Fock basis. This can be done simply by using the option "--state <vector_filename>" (if the vector is complex, one should also use the option "-c" or "--complex-vector").
