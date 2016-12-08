# mrnai-opt
Dynamic programming algorithm to predict the optimal structure of multiple interacting RNA strands.

To install, goto ```mrnai-opt```, and write ```make all```

To run any program, you have to provide sequence information (sequences, their names, and their parity) to stdin. Since this will be needed for multiple programs, we recommend saving it in a file. All sequences should be provided in 5'->3'. For example, suppose the file ```yeast_trunc``` contains
```
NNNNGUAUGUNNNN&NNNNACAGAGAUGAUCAGCNNNN&NNNNGCUUAGAUCAAGUGUAGUANNNN&NNNNUACUAACACCNNNN
I1&U6&U2&I2
1,3|0,2
```
in the first three lines. The first line contains all sequences, separated by ```&```. The second line has their names. The third line contains all even sequences' numbers followed by a `|`, followed by all odd sequences' numbers.
Your file can contain other metadata as long as the first three lines are in this format. 

Note that RNAs are layered in the order specified in line 1. To find the optimal structure from scratch, run
```
./run.sh inputs/yeast_trunc
```
assuming `yeast_trunc` is in the folder `inputs`. Results are written to stdout. The core program that is run here is `findone`. To see the full list of arguments that can be passed, execute `build/findone -h`.

Running this script generates several output files in the folder `output`. Most of them are intermediary files that contain output from RNAup, and have names such as `default-X_Y_itemized.out`, where `X` is an even RNA and `Y` is an odd RNA. A file `partial_matrix.out` will be created by the DP algorithm's program. This file contains the dynamic programming matrix, W(i_1, i_2, ..., i_m) for all indices. This file is needed by the suboptimal-enumeration algorithm's implementation. Finally, a file `compiled_weights.out` will be created, which contains data from all `default-X_Y_itemized.out` files, taking care of orientation information where necessary. This can be used with sampling programs.

### Heuristically Converging to a Permutation
The script `./run_convg` can be used with the same arguments to heuristically converge to an optimal ordering of RNAs. The initial ordering doesn't matter, since the program randomizes the seed configuration every time it is run.

### All Permutations
To take a brute force approach, generate structures for all permuations of RNAs. Use the script `./run_allperms.sh`.

### subopt-enumeration
The program to generate all suboptimal solutions lying within the range [OPT, OPT(1+epsilon)] is also available in this package. As mentioned above, this program needs the partial dynamic programming matrix in order to run. To run everything from scratch (e.g. with epsilon = 0.1), use
```
./run_subopt.sh inputs/yeast_trunc 0.1
```
This will write all suboptimal structures and their energies to a file `subopt_structs.out`. If you already have the partial matrix, or prefer writing to another file, run the file `subopt-enumeration/subopt_enum.py` using the corresponding flags. Note that `subopt_enum.py` needs `&`-separated sequences from stdin.
