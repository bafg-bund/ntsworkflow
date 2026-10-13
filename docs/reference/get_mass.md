# Calculate the mass/charge (m/z) from molecular formula and adduct form

If adduct is not given or \[M\] assumes a \[M+H\]+ or \[M-H\]- adduct if
charge is not 0. Adduct must be one of "\[M+H\]+", "\[M-H\]-", "\[M\]+",
"\[M\]-", "\[M\]".

## Usage

``` r
get_mass(formula, charge, adduct = "[M]")
```

## Arguments

- formula:

  Molecular formula as \`character\`, insert a space before atomic
  numbers

- charge:

  Ion charge as a positive or negative \`numeric\`

- adduct:

  Adduct type as a \`character\`, default is no adduct (\`\[M\]\`). See
  details.

## Value

m/z as numeric

## Details

Uses a modification of the \`OrgMassSpecR\` Function by Nathan Dodder.
The modification allows for the addition of isotopes 2H, 13C 15N and
37Cl (a whitespace charachter must precede the isotope in order to
separate it from the other elements). The adduct must be one of
\`\[M+H\]\`, \`\[M-H\]-\`, \`\[M\]+\`, \`\[M\]-\`, \`\[M\]\`.

## Examples

``` r
# [M+H]+ Mass of Carbamazepine-13C15N
get_mass(formula = "C14 13CH12 15NNO", charge = 1, adduct = "[M+H]+")
#> [1] 239.1026
```
