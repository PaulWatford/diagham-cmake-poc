# Write a spin Hilbert space with a discrete symmetry

Purpose: how to derive a new spin-chain Hilbert-space class that adds a discrete symmetry (here: k → −k inversion in the k = 0 sector) to translation symmetry.
Source: upstream wiki "Create_a_new_spin_Hilbert_space_with_a_discrete_symmetry" (as of 2026-09-24). Changed: every class and method it cites was checked to exist at r4493 (`Spin1_2ChainWithTranslations`, `FermionOnSphereSymmetricBasis::GetCanonicalState`, `FermionOnSphereInvertTable`, `FindCanonicalForm`, `FindNumberTranslation`, `CompatibilityWithMomentum`, `RescalingFactors`, `NbrStateInOrbit`, `SmiSpj`); a pointer to classes that now implement exactly this pattern is added at the end. The wiki's bit-diagram image (`invertbits.png`) was not archived.

This example writes a new Hilbert space for spin-1/2 chains including
translations in the k = 0 sector and the discrete symmetry that maps k to
−k. It is based on the existing class `Spin1_2ChainWithTranslations`
(`Spin/src/HilbertSpace/`); call the new class
`Spin1_2ChainWithTranslationsKInversionSymmetry`.

## The inversion method

The first method to implement takes a configuration and returns its image
under the inversion. Such code already exists for the FQHE Lz → −Lz
symmetry: `FermionOnSphereSymmetricBasis::GetCanonicalState` (inline, in
`FQHE/src/HilbertSpace/FermionOnSphereSymmetricBasis.h`):

```cpp
// get canonical expression of a given state
//
// initialState = state that has to be converted to its canonical expression
// return value = corresponding canonical state

inline unsigned long FermionOnSphereSymmetricBasis::GetCanonicalState (unsigned long initialState)
{
  initialState <<= this->InvertShift;
#ifdef __64_BITS__
  unsigned long TmpState = FermionOnSphereInvertTable[initialState & 0xff] << 56;
  TmpState |= FermionOnSphereInvertTable[(initialState >> 8) & 0xff] << 48;
  TmpState |= FermionOnSphereInvertTable[(initialState >> 16) & 0xff] << 40;
  TmpState |= FermionOnSphereInvertTable[(initialState >> 24) & 0xff] << 32;
  TmpState |= FermionOnSphereInvertTable[(initialState >> 32) & 0xff] << 24;
  TmpState |= FermionOnSphereInvertTable[(initialState >> 40) & 0xff] << 16;
  TmpState |= FermionOnSphereInvertTable[(initialState >> 48) & 0xff] << 8;
  TmpState |= FermionOnSphereInvertTable[initialState >> 56];
#else
  unsigned long TmpState = FermionOnSphereInvertTable[initialState & 0xff] << 24;
  TmpState |= FermionOnSphereInvertTable[(initialState >> 8) & 0xff] << 16;
  TmpState |= FermionOnSphereInvertTable[(initialState >> 16) & 0xff] << 8;
  TmpState |= FermionOnSphereInvertTable[initialState >> 24];
#endif
  initialState >>= this->InvertShift;
  TmpState >>= this->InvertUnshift;
  if (TmpState < initialState)
    return TmpState;
  else
    return initialState;
}
```

To apply the inversion to a state description you first shift the
configuration so that it is centred on bit 16 (32-bit mode) or 32 (64-bit
mode); the inversion is then applied byte by byte, every possible 8-bit
configuration having been precalculated in `FermionOnSphereInvertTable`,
and the result assembled at the right place; finally the state is
unshifted back to the start of the word. Implement the same thing in the
new class as `GetKInversionCanonicalState`.

## Generating the Hilbert space

Inherit from `Spin1_2ChainWithTranslations`. The only methods to rewrite
are the constructors, assignment, and the basic operators such as
`SmiSpj`.

In `Spin1_2ChainWithTranslations` the Hilbert space is generated in the
constructor `Spin1_2ChainWithTranslations(int chainLength, int momentum,
int translationStep, int sz, int memorySize, int memorySlice)`, in this
part:

```cpp
this->StateDescription = new unsigned long [this->EvaluateHilbertSpaceDimension(this->ChainLength, this->Sz)];
long TmpHilbertSpaceDimension = this->GenerateStates(0l, this->ChainLength - 1, (this->ChainLength + this->Sz) >> 1);
this->LargeHilbertSpaceDimension = 0l;
unsigned long TmpState;
unsigned long TmpState2;
int NbrTranslation;
int CurrentNbrStateInOrbit;
unsigned long DicardFlag = ~0x0ul;
for (long i = 0l; i < TmpHilbertSpaceDimension; ++i)
  {
    TmpState = this->StateDescription[i];
    TmpState2 = this->FindCanonicalForm(TmpState, NbrTranslation);
    if (TmpState2 == TmpState)
      {
        CurrentNbrStateInOrbit = this->FindNumberTranslation(TmpState2);
        if (this->CompatibilityWithMomentum[CurrentNbrStateInOrbit] == true)
          {
            ++this->LargeHilbertSpaceDimension;
          }
        else
          {
            this->StateDescription[i] = DicardFlag;
          }
      }
    else
      {
        this->StateDescription[i] = DicardFlag;
      }
  }
```

The first two lines generate the Hilbert space with only the total-Sz
constraint. For each state its canonical state under translation is
computed — the translated state with the largest integer representation —
and only canonical states are kept; a non-canonical state is set to
`DicardFlag`. To add a discrete symmetry, add a similar test before
incrementing the dimension, using `GetKInversionCanonicalState`. Beware
that the order in which the Hilbert space is built defines how a state's
canonicity is checked: first the canonical state under translation, then
under inversion. In particular the canonical state is **not** the one with
the largest integer representation under both transformations at once.

The constructor takes one extra parameter, the parity sector, stored in a
new member — as a `double` with value +1.0 or −1.0, for a reason that
becomes clear below.

## Implementing operators

The operators are trickier. Take `SmiSpj`; in `Spin1_2ChainWithTranslations`:

```cpp
unsigned long tmpState = this->StateDescription[state];
unsigned long State = tmpState;
unsigned long tmpState2 = tmpState;
tmpState >>= i;
tmpState &= 0x1ul;
if (i != j)
  {
    tmpState2 >>= j;
    tmpState2 &= 0x1ul;
    tmpState2 <<= 1;
    tmpState2 |= tmpState;
    if (tmpState2 == 0x1ul)
      {
        State = this->FindCanonicalForm((State | (0x1ul << j)) & ~(0x1ul << i), nbrTranslation, i);
        if (this->CompatibilityWithMomentum[i] == false)
          return this->HilbertSpaceDimension;
        j = this->FindStateIndex(State);
        coefficient = this->RescalingFactors[this->NbrStateInOrbit[state]][i];
        return j;
      }
    else
      {
        coefficient = 0.0;
        return this->HilbertSpaceDimension;
      }
  }
if (tmpState == 0)
  {
    coefficient = -0.25;
    return state;
  }
return this->HilbertSpaceDimension;
```

If the operator acting on the input `state` gives a non-zero result, check
whether the result is canonical under the inversion; if not, multiply
`coefficient` by the sector parity. The difficult part is that the
operator may connect a state that is invariant under the inversion to one
that is not, or the reverse: check the invariance of the input and of the
output state, and apply the scaling factor √2 (invariant → non-invariant)
or 1/√2 (non-invariant → invariant). Do these checks only once the result
is known to be non-zero.

## Real implementations of this pattern

Since this page was written, upstream added Hilbert spaces that implement
translations plus inversion (and Sz) symmetry exactly this way; read them
alongside the recipe:
`Spin/src/HilbertSpace/PairHoppingP1AsSpin1ChainWithTranslationsAndInversionSzSymmetry.{h,cc}`
and `PairHoppingGenericPAsSpin1ChainWithTranslationsAndInversionSzSymmetry.{h,cc}`
(with `…Long` variants for chains beyond 64 sites).
