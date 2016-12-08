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
assuming `yeast_trunc` is in the folder `inputs`.

