OPENQASM 3.0;
include "stdgates.inc";
qubit[2] qqqq;
h qqqq[0];
h qqqq[1];
U(pi/2, pi/4, pi/8) qqqq[0];
cx qqqq[0], qqqq[1];
