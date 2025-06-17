# SLPConstructor
A package for Constructing Subsystem Lifted Product Codes from a Base Matrix

In our GNARSIL paper,  https://ieeexplore.ieee.org/document/10821364, we introduced the concept of Subsystem Lifted Product Codes or SLP codes. A More in-depth view SLP codes can be found here:

https://wp.optics.arizona.edu/opti646/wp-content/uploads/sites/55/2022/11/topological_quantum_computing.pdf


These are exciting codes due to their QLDPC nature*. They can be seen as natural extensions of Subsystem Hypergraph Product (SHP) codes.

Our package is very simple to use:

  1. Build a Base Matrix using np arrays and Galois Polynomial Objects. The Field should always be set to GF2 inside of each object instance.

  2. call SLP with the base matrix and Lift size L as shown.

  3. Have Fun!!

Please Cite the GNARSIL paper if you use this package for publication!!!

* The gauges are QLDPC, Stabilizers are not.
