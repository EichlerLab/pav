# PAV changes

## 1.2

* PAV merges using sequence context by default. Sequences are only merged if
  they are 80% identical by alignment parameters. Merging is done with SV-Pop
  parameters "match=0.8,2,-1,-4,-0.25,500000,9" (80% identical, match=2, mismatch=-1,
  gap-open=-4, gap-extend=-0.25, align up to 500,000 bp (Jaccard similarity if larger),
  Jaccard k-mer size 9). Parameter is tunable with option "merge_align", which may be
  "true" (default), "false" (no sequence matching), or a string that replaces everything
  after "match=".
