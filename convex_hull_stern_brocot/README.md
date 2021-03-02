# Find Convex Hull Based on Stern-Brocot Tree
Given a monotonic convex/concave function f(x), this algorithm shows how to obtain the (upper) convex hull of the integer lattice points under f(x) for $x1 \le x \lt x2$.

Note that you should try to make the interval $x1 \le x \lt x2$ as short as possible to reduce the time complexity. Let the number of vertices of the convex hull be $m$. Typically, the algorithm runs in $O(m)$.

Some examples are added to further validate the algorithm. You may enable the 'print_debug' flag to see the convex hull.

A frequent usage of this algorithm is to count the number of integer lattice points under f(x) together with the Pick's theorem.

Moreover, other weighted-lattice-points counting can also be solved based on this algorithm.

Updates:
DIVCNT1 is solved by the above algorithm. The source code is also uploaded.
