# Run GenForm with various Commands

This is a wrapper function for GenForm, requiring MS, MS/MS, ion
settings and directories. Original function Copyright Emma Schymanski
(schymane/GenFormR). Modified to work on CentOS.

## Usage

``` r
RunGenForm(ms_file, msms_file, ion, GenFormDir, ResultsDir,ppm=5, acc=15,
elements="CHNOPS",ff="", dbe=TRUE, oei=TRUE, mz=0, exist=TRUE, an_loss=TRUE,an_loss_limit=500,
sort=FALSE, cleanMSMS=FALSE)
```

## Arguments

- ms_file:

  File name of a plain text, two column file (e.g. `MS.txt`) containing
  space-separated mz and intensity values of the MS1 spectrum (full
  scan). The most intense peak is taken as "m" unless this is defined
  with `"mz"`.

- msms_file:

  A plain text, two column file (e.g. `MSMS.txt`) containing
  space-separated mz and intensity values of the tandem (MS2) spectrum,
  i.e. the fragments. While GenForm can be run without MSMS, this
  function requires its use at this stage.

- ion:

  Ion setting to use. Options: `c("-e","+e","-H","+H","+Na")`.

- GenFormDir:

  Directory containing `GenForm.exe`.

- ResultsDir:

  Directory to save the GenForm results. Existing files will be
  overwritten.

- ppm:

  Accuracy setting for MS1 (default `5` ppm).

- acc:

  Accuracy setting for MS2 (default `15` ppm).

- elements:

  Selection of elements to consider. Default `"CHNOPS"`, implemented are
  `"CHBrClFINOPSSi"`.

- ff:

  A "fuzzy formula" input to offer more control beyond `elements`. This
  should be in format `C0-6H0-20O1-3`.

- dbe:

  Default `TRUE` writes double bond equivalents to output; `FALSE`
  suppresses this output.

- oei:

  Default `TRUE` allows odd electron ions for explaining MS/MS peaks,
  `FALSE` means only even electron ions are considered.

- mz:

  If mz is \>0, this is the mass GenForm will use to calculate the
  formulas, considering the `ion` setting. Default is the most intense
  peak in `ms_file`.

- exist:

  Default `TRUE` activates the existance filter (allowing only formulae
  for which at least one structural formula exists). `FALSE` is not
  recommended for most use cases.

- an_loss:

  Activate the annotation of MSMS fragments with subformulas and losses.

- an_loss_limit:

  Define the cut-off limit for annotating subformulas and losses. For
  many peaks and formulas, this becomes time consuming and produces
  massive files of limited use.

- sort:

  Default `FALSE` will produce default output. If `TRUE`, generates an
  additional sorted file sorted by combined match value.

- cleanMSMS:

  Default `FALSE`. If `TRUE`, produces an additional cleaned MSMS file
  containing only peaks with subformula annotation.

## Value

Several files containing various GenForm outputs.

## References

GenForm on Source Forge: https://sourceforge.net/projects/GenForm/

## Author

Emma Schymanski (R wrapper, `emma.schymanski@uni.lu`) and Markus
Meringer (GenForm)
