# mrnai-opt
Dynamic programming algorithm to predict the optimal structure of multiple interacting RNA strands.

To install, goto $PWD, and write

```
make
```

To run any program, you have to provide sequence information (sequences, their names, and their parity) in a file. All sequences should be provided in 5'->3'. For example, suppose the file ```yeast``` contains
```
NNNNGUAUGUNNNN&NNNNACAGAGAUGAUCAGCNNNN&NNNNGCUUAGAUCAAGUGUAGUANNNN&NNNNUACUAACACCNNNN
I1&U6&U2&I2
1,3|0,2
```
