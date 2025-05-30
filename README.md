# structure_theory_python

Some structures and techniques from Algebraic Structure Theory of Sequential Machines (1966), in python.

The book deals with the properties of the SP lattice, which had been described by earlier researchers in systems (Ashby 1956). Hartmanis and Stearns independently developed the concept for engineering applications while working for General Electric, and showed how it could be used to synthesize circuits realizing functions under different constraints with respect to components, error, etc. This finite state work is much less famous than their development of computational complexity theory.

They present one algorithm for enumerating the SP lattice of a machine, which naturally suffers a worst-case exponential running time since the size of the lattice is worst-case exponential in the size of the machine. The algorithm still improves over naive enumeration as the demonstration in `plot_sp.py` shows.
