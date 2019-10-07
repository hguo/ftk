import numpy as np

def get_levelsets(filename):
    raw_level_sets = np.fromfile(filename, dtype=np.float32)

    print(len(raw_level_sets))

    index = 0
    level_sets = []
    ntimestep = 0

    while index < len(raw_level_sets):
        ncoord = int(raw_level_sets[index] / 4)
        index += 1

        level_set = []
        for i in range(ncoord):
            coord = np.zeros((4))

            coord[0] = raw_level_sets[index + 0]
            coord[1] = raw_level_sets[index + 1]

            coord[2] = raw_level_sets[index + 2]  # time
            coord[3] = raw_level_sets[index + 3]  # value

            index += 4

            level_set.append(coord)
            if coord[2]+1 > ntimestep:
                ntimestep = coord[2] + 1

        level_sets.append(level_set)

        # break

    ntimestep = int(ntimestep)

    return level_sets, ntimestep


def level_sets_to_3d(filename, shape, level_sets):
    volume = np.zeros(shape).astype('float32')

    nlevelset = len(level_sets)
    nele = 0
    max_nele_set = 0
    for i in range(nlevelset):
        nele_set = len(level_sets[i])

        nele += nele_set
        max_nele_set = max(max_nele_set, nele_set)

        level_set = level_sets[i]
        for coord in level_set:
            volume[int(coord[2])][int(coord[1])][int(coord[0])] = i + 1
            # volume[int(coord[2])][int(coord[1])][int(coord[0])] = coord[3]


    print("# of Connected Components: " + str(nlevelset))
    print("# of element: " + str(nele))
    print("Ratio of largest component: " + str(float(max_nele_set) / nele))

    volume.astype('float32').tofile(filename)


if __name__ == '__main__':

    # dataset = "synthetic/syn_512"
    dataset = "bout/bout_512"

    input_filename = "../data/" + dataset + ".level_set"
    level_sets, ntimestep = get_levelsets(input_filename)

    # Shape is in the order of time (t), height(h), width (w)
    if 'synthetic' in dataset:
        shape = (128, 128, 128)
    elif 'bout' in dataset:
        shape = (701, 880, 425)  # time (t) height(h) width (w)
    output_filename = "../data/" + dataset + ".level_set.raw"
    level_sets_to_3d(output_filename, shape, level_sets)
