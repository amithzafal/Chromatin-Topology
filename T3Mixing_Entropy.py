import os
import sys
import numpy as np
from vtkReader import vtkReader
from scipy.spatial import cKDTree

class LocalMixingAnalysis:

    def __init__(self, outputDir, chrom1, chrom2, chrom3, initFrame, radius):
        self.reader1 = vtkReader(outputDir, chrom1, initFrame, readLiq=False, readPoly=True)
        self.reader2 = vtkReader(outputDir, chrom2, initFrame, readLiq=False, readPoly=True)
        self.reader3 = vtkReader(outputDir, chrom3, initFrame, readLiq=False, readPoly=True)

        self.radius = radius
        self.mixingFile = os.path.join(self.reader1.outputDir, f"r{initFrame}_localEntropy_r{radius}.res")
        if os.path.exists(self.mixingFile):
            print(f"File '{self.mixingFile}' already exists - aborting")
            sys.exit()

    def Compute(self):
        self.ProcessFrame()

    def ProcessFrame(self):
        data1 = next(self.reader1)
        data2 = next(self.reader2)
        data3 = next(self.reader3)
        
        pos1 = data1.polyPos
        pos2 = data2.polyPos
        pos3 = data3.polyPos
        pos_total = np.concatenate((pos1, pos2, pos3), axis=0)

        num_monomers = len(pos1)
        #total_monomers = len(pos_total)
        n_polymers = 3
        polymer_types = np.repeat(np.arange(n_polymers), num_monomers)

        tree = cKDTree(pos_total)
        neighbors = tree.query_ball_point(pos_total, self.radius)
    
        # Compute local compositions and mixing entropy
        entropies = []
        min_neighbors = 2
        for i, local_neighbors in enumerate(neighbors):
            if len(local_neighbors) > min_neighbors:
                # Get polymer types of local neighbors
                local_types = np.concatenate(([polymer_types[i]], polymer_types[local_neighbors]))
                #type_counts = Counter(local_types)
                
                n_total = len(local_types)
                composition = np.array([np.sum(local_types == j) for j in range(n_polymers)]) / n_total
            
                # Compute local mixing entropy
                non_zero = composition > 0
                local_entropy = -np.sum(composition[non_zero] * np.log(composition[non_zero]))
            
                entropies.append(local_entropy)
    
        # Compute average entropy
        norm_entropy = np.array(entropies) / np.log(n_polymers)

        self.avg_mixing_index = norm_entropy

    def Print(self):
        np.savetxt(self.mixingFile, self.avg_mixing_index)
        print(f"\033[1;32mPrinted local entropy to '{self.mixingFile}'\033[0m")

if __name__ == "__main__":
    if len(sys.argv) != 7:
        print("\033[1;31mUsage is %s outputDir chrom1 chrom2 chrom3 initFrame radius\033[0m" % sys.argv[0])
        sys.exit()

    outputDir = sys.argv[1]
    chrom1 = sys.argv[2]
    chrom2 = sys.argv[3]
    chrom3 = sys.argv[4]
    initFrame = int(sys.argv[5])
    radius = float(sys.argv[6])

    mixing_analyzer = LocalMixingAnalysis(outputDir, chrom1, chrom2, chrom3, initFrame=initFrame, radius=radius)
    mixing_analyzer.Compute()
    mixing_analyzer.Print()