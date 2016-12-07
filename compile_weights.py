#!/usr/bin/python

import sys

def main():
    inp = sys.stdin.readlines()
    sequences = inp[0].strip().split("&")
    names = inp[1].strip().split("&")

    # convert "a,b,c|d,e,f" to even = [a,b,c], odd = [d,e,f]
    even, odd = map(lambda x: map(int, x.split(",")), inp[2].strip().split("|"))


    N = len(sequences) - 1  # len-1 if no wrap around
    if '-w' in sys.argv:
        N = len(sequences)
        
    for level in range(N):
        # file formats are "default-EVEN_ODD_fpGfunc.out"

        if level in even:

            # not much work required
            # look at "default-<level>_<level+1>_fpGfunc.out"
            fname = "output/default-%d_%d_fpGfunc.out" % (level, (level+1)%len(sequences))
            with open(fname) as f:
                for line in f:
                    print '%d\t%s' % (level, line.strip())

        else:
            # look at "default-<level+1>_<level>_fpGfunc.out"
            fname = "output/default-%d_%d_fpGfunc.out" % ((level+1)%len(sequences), level)
            with open(fname) as f:
                for line in f:
                    data = line.strip().split()

                    # replace rna indices
                    a,b,c,d = data[:4]
                    data[0:2] = c,d
                    data[2:4] = a,b

                    # replace unpaired energies
                    data[9:11] = (data[9:11])[::-1]
                    
                    print '%d\t%s' % (level, '\t'.join(data))




if __name__ == '__main__':
    main()
