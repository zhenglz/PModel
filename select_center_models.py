import mdtraj as mt
import sys
import sklearn
import numpy as np
from sklearn import cluster
import pandas as pd
import argparse

def load_pdb(fn, stride=10):
    traj = None
    for i, t in enumerate(mt.iterload(fn, stride=stride, chunck=10)):
        top = t.topology
        protein_indices = top.select("CA")
        t = t.atom_slice(protein_indices)

        if i == 0:
            traj = t
        else:
            traj = mt.join([traj, t])

    return traj[int(traj.n_frames * 0.8): ]

def distance_matrix(dat):
    dist = []
    for x in dat:
        for y in dat:
            d = rmse(x, y)
            dist.append(d)

    return np.array(dist).reshape((dat.shape[0], -1))

def rmse(x, y):
    return np.sqrt(np.mean(np.square(x - y)))


if __name__ == "__main__":
    d = """
    """

    parser = argparse.ArgumentParser(description=d)
    parser.add_argument("--fn", type=str, nargs="+",
                        help="Input, str, multiple. Input protein trajectory files.")
    parser.add_argument("--stride", type=int, default=10,
                        help="Input, int, optional, default is 10. The stride for trajectory reading. \n"
                             "N stride means only read the frames every n frames.")
    parser.add_argument("--n_clusters", type=int, default=5,
                        hep="Input, int, optional, default is 5. Output 5 clusters' centers. ")

    args, unkown = parser.parse_known_args()

    traj_fns = args.fn

    traj = None

    for i, fn in enumerate(traj_fns):
        if i == 0:
            traj = load_pdb(fn, args.stride)
        else:
            traj = traj.join(load_pdb(fn, args.stride))
        print("Load PDB progress: ", i, fn)

    traj = traj.superpose(traj[0])
    print("Number of frames used: ", traj.n_frames)

    xyz = traj.xyz.reshape((traj.n_frames, -1))
    pca = sklearn.decomposition.PCA(n_components=5)
    xyz_t = pca.fit_transform(xyz)
    print("PCA done")

    dist_matrix = distance_matrix(xyz_t)
    dist_matrix = pd.DataFrame(dist_matrix)
    print("Distance matrix done")
    print("Distance matrix shape ", dist_matrix.shape)
 
    clust = cluster.AgglomerativeClustering(n_clusters=args.n_clusters)
    clust.fit(xyz_t)
    labels = clust.labels_
    print("Group labels of the frames: ", labels)

    for i in range(args.n_clusters):
        # cluster i
        i_label = (labels == i)
        dist = dist_matrix[i_label]
        print(dist.shape, i_label.sum(), "Group ", i)
        i_pdb = dist.sum(axis=1).sort_values().index.values[0]
        traj[i_pdb].save_pdb("center_%d.pdb" % i, )

        print("Save cluster center structure ", i, " frame index in traj: ", i_pdb)
