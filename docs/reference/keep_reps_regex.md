# Replicate feature based on regular expression

use a regular expression on the filename to gather the groups of samples
in the batch into sets. Then apply the keepReps function on each of
these sets in turn. This will be quite slow because we are looping
through the sets and each set loops through each row. But it is the
easiest way for now.

## Usage

``` r
keep_reps_regex(alignment, samples, sampleList, regexp, least)
```

## Arguments

- alignment:

  alignment table (matrix, grouped)

- samples:

  Sample IDs which should be considered (e.g. only samples and no
  blanks)

- sampleList:

  Sample list dataframe with columns ID, File etc.

- regexp:

  string length 1 regular expression by which to group the samples into
  replicates this is done removing the matching string with
  stringr::str_replace(filename, regexp, "\1"). The part of the regexp
  in brackets is kept (\1) and so should be the same in all replicates.

- least:

  How many times does the feature need to appear. Must be less than the
  number of replicates in each group. If it is higher, will be changed
  to the number of replicates in the group with a message.

## Value

alignment table with rows removed where the number of aligned features
in each replicate set is less than 'least'.
