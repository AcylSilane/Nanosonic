
import ase.cluster
import ase.neighborlist
import numpy as np
import pyaudio


# Set up globals for defaults
DEFAULT_ELEMENTS = ("Cu", "Cu")
DEFAULT_RADIUS = 2.8


def buildAdjacencyMatrix(atoms_object: "ase.Atoms",
                         radius_dictionary: "dict" = {DEFAULT_ELEMENTS: DEFAULT_RADIUS}) -> "np.ndarray":
    """
    Sparse matrix representation from an ase atoms object.
    Args:
    atoms_object (ase.Atoms): An ASE atoms object representing the system of interest
    radius_dictionary (dict): A dictionary with the atom-atom radii at-which a bond is considered a
                              bond. If no dict is supplied, Cu-Cu bonds of max-len 2.8 are assumed.
    Returns:
    np.ndarray : A numpy array representing the sparse matrix of the ase object
    """
    # Construct the list of bonds
    sources, destinations = ase.neighborlist.neighbor_list("ij", atoms_object, radius_dictionary)
    # Generate the matrix
    adjacency_matrix = np.zeros([len(atoms_object), len(atoms_object)])
    for bond in zip(sources, destinations):
        adjacency_matrix[bond[0], bond[1]] += 1
    return np.matrix(adjacency_matrix)

def normalize_and_scale(eigenvals, highest_tone = 7902.13, lowest_tone = 16.35):
    scale = (eigenvals - np.average(eigenvals)) / np.std(eigenvals)
    scale = np.real(np.round(scale - np.min(scale), 4))
    scale /= (np.max(scale) / (highest_tone-lowest_tone))
    scale += lowest_tone
    return scale

def sine(frequency, length, rate):
    # Sine generation code is from https://stackoverflow.com/questions/42192239/remove-control-clicking-sound-using-pyaudio-as-an-oscillator
    length = int(length * rate)
    factor = (float(frequency) * (np.math.pi * 2) / rate)
    return np.sin(np.arange(length) * factor)

def tonify(tone, stream):
    # Tonification code is from https://stackoverflow.com/questions/42192239/remove-control-clicking-sound-using-pyaudio-as-an-oscillator
    volume = 0.5
    duration = 0.05
    samples = [sine(tone, duration, 44100)]
    sample = np.concatenate(samples) * 0.25
    fade = 200
    fade_in = np.arange(0., 1., 1/fade)
    fade_out = np.arange(1., 0., -1/fade)

    sample[:fade] = np.multiply(sample[:fade], fade_in)
    sample[-fade:] = np.multiply(sample[-fade:], fade_out)

    stream.write(sample.astype(np.float32).tostring())

def soundify(tones, n_repeats = 1):
    pa = pyaudio.PyAudio()
    stream = pa.open(format=pyaudio.paFloat32,
                    channels=1,
                    rate=44100,
                    output=True)
    for i in range(0,n_repeats):
        for j in tones:
            tonify(j, stream)
    stream.close()
    pa.terminate()
def midify(tones):
    pass

if __name__ == "__main__":
    for i in range(2,5):
        x1 = ase.cluster.FaceCenteredCubic("Cu", np.identity(3), [i]*3)
        x2 = ase.cluster.Icosahedron("Cu", i)
        x3 = ase.cluster.Octahedron("Cu", i)
        for x in [x1, x2, x3]:
            print(x)
            adjmat = buildAdjacencyMatrix(x)
            y = np.linalg.eig(adjmat)[0]
            z=sorted(set(normalize_and_scale(y)))
            soundify(z)

