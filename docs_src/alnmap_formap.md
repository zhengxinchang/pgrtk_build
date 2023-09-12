
The command line tool `pgr-alnmap` can generate alingment map between a set of assembled contigs
and a reference file. It creates a number of files that are useful for downstram analysis.
This document describe the format of the generated files.

## `*.alnmap` files

A `*.alnmap` files contain the blocks that mapped from the assembly contigs (also called "query") to
the reference fiel (also called "target") from the PGR-TK-WGS's whole genome alignment code. (We will
describe the algorithm behind it in a different documents.) The file contains chains of alignments Each
chains is a set of blocks. 

Each chain of the align block has an integer id which is the first field of the line (also called as a 
record.) The second field is the record types:

- "B": the record represents the begin of the chain of the align blocks.
- "E": the record represents the end of the chain of the align blocks.
- "M": the record represents a full alignment (no variant) block.
- "M_D": the record represents a full alignment (no variant) block, however, two or more query are mapped to the same target block
- "V": the record represetns there are variants between the query and the target, the variant information are appended.
- "V_D": the record represents there are variants between the query and the target, the variant information are appended. There are other query  blocks mapping to the same target block.
- "S": the record represents potential structral variants between the query and the target block.
- "S_D": the record represents potential structral variants between the query and the target block. There are other query  blocks mapping to the same target block.

All records share the following common 9 fields seperated by `tab`:

`aligned_chain_id, block_type, target_name, target_start, target_end, query_name, query_start, query_end, query_strand`

The "B" and "V/V_D" blocks contains additional fields.


The follow command generate the unqiuely mapped blocks:

```
cat hg19_to_grch38.alnmap | awk '$2 == "V" || $2 =="M" || $2 == "S" ' | cut  -f1-9 | sort  -k3,3 -k4,4n -u >  hg19_to_grch38_uniq_blocks.alnmap

```

and the duplicate mapped blocks:

```
cat hg19_to_grch38.alnmap | awk '$2 == "V_D" || $2 =="M_D" || $2 == "S_D" ' | cut  -f1-9 | sort  -k3,3 -k4,4n -u >  hg19_to_grch38_dup_blocks.alnmap
```

For the "V" and "V_D", six addition fields are appended:
(`variant_position_in_the_target_block`, `variant_position_in_the_query_block`, `variant_position_in_the_target_sequence`, `variant_type`, `ref_seq`, `variant_seq`)

For the "B" and "E" blocks, two additional fields are appended: `query_sequence_length`, `the_alingment_orientation_of_the_contig`.

For the "S" and "S_D" blocks, two addional fileds are appended: `the_alingment_orientation_of_the_contig`, `sv_candidate_ type`.




