# source: https://github.com/ekimb/rust-mdbg/blob/f155058b6db7a375fb11940715f3d25bb15fcd73/utils/gfa_break_loops.py

# MIT License
# 
# Copyright (c) 2021 Barış Ekim & Rayan Chikhi & Bonnie Berger
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# this script artificially breaks loop nodes by cutting one of the extremity edges
# e.g. 
# L x + y - 10M
# L x + y + 10M
# will remove one of the two lines
# also removes 
# L x + y + 10M
# L x + y + 10M
# (it doesn't care about orientation)

# new: and also let's make it remove x<->x loops

seen_edges = set()
import sys, traceback
for line in open(sys.argv[1]):
    if not line.startswith('L'): 
        print(line.strip())
        continue
    try:
        e_in, e_out = line.split()[1], line.split()[3]
        e_tuple = tuple(sorted([e_in, e_out]))
        to_remove = e_tuple in seen_edges
        if e_in == e_out: # self loop
            to_remove = True
        seen_edges.add(e_tuple)
        if not to_remove:
            print(line.strip())
    except:
        print("at line:",line.strip())
        traceback.print_exc(file=sys.stdout)
        exit(1)
