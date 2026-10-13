# Add additional peak and sample ID information to peaklists

takes a list of peaklists and adds sample ID and peak ID information

## Usage

``` r
new_peaklist(peaklist, datenList)
```

## Arguments

- peaklist:

  \`list\` of peak tables

- datenList:

  \`list\` of \`xcms::xcmsRAW\` objects

## Value

New peak table \`list\` in the same form as before but with additional
columns sample_id and peak_id_all. If these were already present, they
are over-written

## Details

This needs to be run before alignment
