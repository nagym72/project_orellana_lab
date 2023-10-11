package require saltbr
package require hbonds

animate delete  beg 0 end 1 skip 0 0

set all_atoms [atomselect top all]
saltbr -sel $all_atoms -outdir ./salt_bridges
hbonds -sel1 $all_atoms -writefile yes -type all

exit
