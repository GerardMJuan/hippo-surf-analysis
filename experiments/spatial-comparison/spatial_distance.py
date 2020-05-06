"""
Compare similarity between two meshes. NEW VERSION.

This version computes the similarity between two meshes using cosine similarity, outputting a mean and a mesh vertex per vertex,
and do a permutation test. It saves the mesh on a given output dir, and returns the mean and p-value of the permutation test.

Example:
"""

import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity
import sys
import vtk
from scipy.stats import percentileofscore
from joblib import Parallel, delayed

def test_cosine_distance(mesh_1, mesh_2, out_file, avg_mesh, perm=False):
    """
    computes cosine distance.
    mesh_1 is the first mesh to compare.
    mesh_2 is the second mesh to compare.
    out_file is where the cosine mesh will be saved.
    avg_mesh is the average mesh used to draw the distances on it.
    perm: parameter to decide if the permutation test is done or not, False by default
    """

    # Load the two meshes
    df_values_mesh_1 = pd.read_csv(mesh_1, header=None)
    df_values_mesh_2 = pd.read_csv(mesh_2, header=None)

    # Remember that those meshes are organized as:
    # 0 Hotellings's T
    # 1 angle
    # 2 x component
    # 3 y component
    # 4 z component
    components_1 = df_values_mesh_1.iloc[:, [2, 3, 4]].to_numpy()
    components_2 = df_values_mesh_2.iloc[:, [2, 3, 4]].to_numpy()

    # Compute cosine similarity
    result_similarity = [cosine_similarity(x.reshape(1, -1), y.reshape(1, -1)) for (x, y) in zip(components_1, components_2)]

    # Get median of the cosine similarity
    mean_res =(np.mean(result_similarity) + 1.0)/2
    print('Mean similarity: ' + str(mean_res))

    # Save as resulting mesh
    t_values = result_similarity
    vtk_mesh = vtk.vtkDoubleArray()
    vtk_mesh.SetName("Similarity")
    for x in t_values:
        vtk_mesh.InsertNextValue(x)

    # Save it to disk
    rd = vtk.vtkPolyDataReader()
    rd.SetFileName(avg_mesh)
    rd.Update()
    mesh = rd.GetOutput()
    mesh.GetPointData().SetScalars(vtk_mesh)
    rd2 = vtk.vtkPolyDataWriter()
    rd2.NewInstance()
    rd2.SetInputData(mesh)
    rd2.SetFileName(out_file + '.vtk')
    rd2.Write()

    #Start permutation test
    N = 100000

    #Structures that save the permutation results
    #mean_dist = np.zeros(shape=N)
    #mesh_dist = np.zeros(shape=(N, components_1.shape[0]))

    ## Parallelizing
    mesh_dist = []
    for i in range(N):
        rd_point_1 = np.random.normal(size=3)
        rd_point_2 = np.random.normal(size=3)

        # we add 1 so that the distribution is over numbers > 0 so that later we can do percentileofscore.
        # We are just translating the distribution.
        res_sim = cosine_similarity(rd_point_1.reshape(1, -1), rd_point_2.reshape(1, -1)) + 1.0
        mesh_dist.append(res_sim)
    # mesh_dist = Parallel(n_jobs=-1)(delayed(permutation_par)(i) for i in range(N))
    # mesh_dist = np.array(0)

    # Score over the surface of the hippocampus
    # Compute percentileofscore over each vertex
    p_test_mesh = np.zeros(shape=components_1.shape[0])
    i = 0

    for ver in list(result_similarity):
        #TODO: need to do a sanity check here
        # We add one so that we make sure that instead of the range being -1 .. 1, we transfer it to 0 .. 2, so
        # that hte function percentileofscore works.
        p_test_mesh[i] = percentileofscore(mesh_dist, ver+1.0, kind='mean')
        p_test_mesh[i] = 1 - p_test_mesh[i] / 100
        i += 1

    # Save it as csv
    df_pval = pd.DataFrame({"values": p_test_mesh})
    df_pval.to_csv(out_file + '_pval.csv', sep='\t')

    #For paraview visualization purposes
    #Save the p_test_mesh onto disk:
    t_values = p_test_mesh
    vtk_mesh = vtk.vtkDoubleArray()
    vtk_mesh.SetName("Similarity p values")
    for x in list(t_values):
        vtk_mesh.InsertNextValue(x)

    # Save it to disk
    rd = vtk.vtkPolyDataReader()
    rd.SetFileName(avg_mesh)
    rd.Update()
    mesh = rd.GetOutput()
    mesh.GetPointData().SetScalars(vtk_mesh)
    rd2 = vtk.vtkPolyDataWriter()
    rd2.NewInstance()
    rd2.SetInputData(mesh)
    rd2.SetFileName(out_file + '_pval.vtk')
    rd2.Write()

    return mean_res, out_file + '_pval.csv'


def permutation_par(i):
    """
    Auxiliary function to paralelize code.

    i doesnt actually do anything lmaoooo
    """
    ##TODO1: AIXO S'HA DE PARALELITZAR D'ALGUNA MANERA. joblib?
    ##TODO2: i fer al cluster segurament
    rd_point_1 = np.random.normal(size=3)
    rd_point_2 = np.random.normal(size=3)

    # we add 1 so that the distribution is over numbers > 0 so that later we can do percentileofscore.
    # We are just translating the distribution.
    result_similarity = cosine_similarity(rd_point_1.reshape(1, -1), rd_point_2.reshape(1, -1)) + 1.0
    #print(np.array(result_similarity).shape)
    return result_similarity
    # save it
    # mean_dist[i] = mean_res
    # mesh_dist[i, :] = result_similarity




if __name__ == "__main__":

    np.random.seed(1714)

    # Path of mesh 1
    mesh_1 = sys.argv[1]

    # Path of mesh 2
    mesh_2 = sys.argv[2]

    # output file of out_file
    out_file = sys.argv[3]

    # average mesh
    avg_mesh = sys.argv[4]

    mean, pmesh = test_cosine_distance(mesh_1, mesh_2, out_file, avg_mesh, True)
    print(mean)