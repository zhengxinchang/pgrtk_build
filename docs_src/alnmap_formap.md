
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
- "M_O": the record represents a full alignment (no variant) block, however, two or more query are mapped to the same target block from overlapped alignment chains
- "V": the record represetns there are variants between the query and the target, the variant information are appended.
- "V_D": the record represents there are variants between the query and the target, the variant information are appended. There are other query  blocks mapping to the same target block.
- "V_O": the record represents there are variants between the query and the target, the variant information are appended. There are other query  blocks mapping to the same target block from overlapped alignment chains.
- "S": the record represents potential structral variants between the query and the target block.
- "S_D": the record represents potential structral variants between the query and the target block. There are other query  blocks mapping to the same target block.
- "S_O": the record represents potential structral variants between the query and the target block. There are other query  blocks mapping to the same target blocke

All records share the following common 9 fields seperated by `tab`:

`aligned_chain_id, block_type, target_name, target_start, target_end, query_name, query_start, query_end, query_strand`


The follow command generate the unqiuely mapped blocks:

```
cat grch38_to_chm13.alnmap | awk '$2 == "V" || $2 =="M" || $2 == "S" ' | cut  -f1-9 | sort  -k3,3 -k4,4n -u >  grch38_to_chm13_unique_blocks.alnmap

```

and the duplicate mapped blocks:

```
cat grch38_to_chm13.alnmap | awk '$2 == "V_D" || $2 =="M_D" || $2 == "S_D" ' | cut  -f1-9 | sort  -k3,3 -k4,4n -u >  grch38_to_chm13_dup_blocks.alnmap
```

For the "V", "V_D" and "V_O" records, six addition fields are appended:
(`variant_position_in_the_target_block`, `variant_position_in_the_query_block`, `variant_position_in_the_target_sequence`, `variant_type`, `ref_seq`, `variant_seq`)

For the "B" and "E" blocks, two additional fields are appended: `query_sequence_length`, `the_alingment_orientation_of_the_contig`.

For the "S", "S_D" and "S_O" blocks, two addional fileds are appended: `the_alingment_orientation_of_the_contig`, `sv_candidate_ type`.




