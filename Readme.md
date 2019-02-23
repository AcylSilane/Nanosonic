# ¯\\\_(ツ)\_/¯

The inspiration for this project was from: https://www.jasondavies.com/graph-music/

So basically, this works by first building an adjacency matrix from the atoms in a nanoparticle (Cu is just assumed by default. I'll probably change this later). This provides a translation and rotation-invariant representation of the NP.

From there, the eigenvalues of the adjacency matrix are calculated to generate what's known as a graph spectrum (http://mathworld.wolfram.com/GraphSpectrum.html).

The graph spectrum is then normalized between the notes C0 and B8, duplicate notes are removed, and they're played through a sine tone generator. Try it out and see what you get!
