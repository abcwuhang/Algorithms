# Find Convex Hull Based on Stern-Brocot Tree
Given a monotonic convex/concave function f(x), this algorithm shows how to obtain the (upper) convex hull of the integer lattice points under f(x) for <a href="https://www.codecogs.com/eqnedit.php?latex=x1&space;\leq&space;x&space;\lt&space;x2" target="_blank"><img src="https://latex.codecogs.com/gif.latex?x1&space;\leq&space;x&space;\lt&space;x2" title="x1 \leq x \lt x2" /></a>.

Note that you should try to make the interval <a href="https://www.codecogs.com/eqnedit.php?latex=x1&space;\leq&space;x&space;\lt&space;x2" target="_blank"><img src="https://latex.codecogs.com/gif.latex?x1&space;\leq&space;x&space;\lt&space;x2" title="x1 \leq x \lt x2" /></a> as short as possible to reduce the time complexity. Let the number of vertices of the convex hull be <a href="https://www.codecogs.com/eqnedit.php?latex=m" target="_blank"><img src="https://latex.codecogs.com/gif.latex?m" title="m" /></a>. Typically, the algorithm runs in <a href="https://www.codecogs.com/eqnedit.php?latex=O(m)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?O(m)" title="O(m)" /></a>.

Some examples are added to further validate the algorithm. You may enable the "print_detail" flag to see the convex hull.

A frequent usage of this algorithm is to count the number of integer lattice points under f(x) together with the Pick's theorem.

Moreover, other weighted-lattice-points counting can also be solved based on this algorithm.

### Updates:
[DIVCNT1](https://www.spoj.com/problems/DIVCNT1/) is solved by the above algorithm. The source code is also uploaded.

[AFS3](https://www.spoj.com/problems/AFS3/) is uploaded. However, it is too slow to be accepted. 
