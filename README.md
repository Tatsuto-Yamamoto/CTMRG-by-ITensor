# CTMRG-by-ITensor

This the source code for CTMRG(Corner Transfer Matrix Renormalization Group) method(https://arxiv.org/abs/cond-mat/9507087v5).

CTMRG is one of the DMRG methods, and it calculates the thermodynamical quantity of the two-dimensional classical spin system efficiently.

For tensor contraction, ITensor(http://itensor.org/index.html) is called.

To compile this, you have to rewrite "Makefile", as
LIBRARY_DIR=(Pass to ITensor)/itensor
... and to type "make" to compile.

In this source code, we calculate for the Ising model, but you can use it for other models by rewriting it just a little.

This source code was mainly written by qope(https://github.com/qope).
