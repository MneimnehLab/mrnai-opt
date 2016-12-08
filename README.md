# mrnai-opt
Dynamic programming algorithm to predict the optimal structure of multiple interacting RNA strands.

To install, goto ```mrnai-opt```, and write ```make all```

To run any program, you have to provide sequence information (sequences, their names, and their parity) in a file. All sequences should be provided in 5'->3'. For example, suppose the file ```yeast_trunc``` contains
```
NNNNGUAUGUNNNN&NNNNACAGAGAUGAUCAGCNNNN&NNNNGCUUAGAUCAAGUGUAGUANNNN&NNNNUACUAACACCNNNN
I1&U6&U2&I2
1,3|0,2
```
in the first three lines. The first line contains all sequences, separated by ```&```. The second line has their names. The third line contains all even sequences' numbers followed by a `|`, followed by all odd sequences' numbers.
Your file can contain other metadata as long as the first three lines are in this format.

To find the optimal structure, run
```
./run.sh inputs/yeast_trunc
```
assuming `yeast_trunc` is in the folder `inputs`. Running this script generates several output files. Most of them are intermediary files that contain output from RNAup, and have names such as `default-X_Y_itemized.out`, where `X` is an even RNA and `Y` is an odd RNA. A file `partial_matrix.out` will be created by the DP algorithm's program. This file contains the dynamic programming matrix, W(i_1, i_2, ..., i_m) for all indices. This file is needed by the suboptimal-enumeration algorithm's implementation. Finally, a file `compiled_weights.out` will be created, which contains data from all `default-X_Y_itemized.out` files, taking care of orientation information where necessary. This can be used with sampling programs.

