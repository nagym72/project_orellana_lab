from sklearn.cluster import AffinityPropagation, MeanShift, AgglomerativeClustering, KMeans, DBSCAN, OPTICS # you can experiment with different algorithms. Usually Kmeans does the job well. Sometimes OPTICS performs better but for your case you shall go with Kmeans
# you have also the option downstream to do OPTICS in the code.

from sklearn.metrics import silhouette_score # this is the metric we need to judge how many clusters are optimal (see downstairs the function _obtain_optimal_clusters)
import numpy as np # this is default good to have you will need it later for all kind of metrics we compute.
import math # same as np useful always to be there.
import os
from matplotlib import pyplot as plt # plotting default





# this you might need to parse your pdb/mmcif and superpose them if you have multiple structures and then you can do all kind of distance measurements.
#from Bio.PDB import MMCIFParser, MMCIFIO, PDBParser, PDBIO, Select, Superimposer, Model, Structure, Chain, Residue, Atom, Selection, NeighborSearch 




class AutoClusterClass:
    """Works but needs documentation."""
    def __init__(self, work_dir, pca_result_dict, pca_expl_var, logging=True, clusters="auto"):

        self.work_dir = work_dir
        self.pca_result_dict = pca_result_dict
        self.pca_explained_var = pca_expl_var
        self.logging = logging
        self.num_clusters = clusters
        self.cluster_results = {}
        self.representative_strucs_per_cluster = {}

    def run_clustering(self, automated=True, OPTICS=False, optics_min_samples=None):
        if self.pca_result_dict:
            #then lets go and convert this to a np.array

            # we cluster for each range and result.
            for ranges, inner_dict in self.pca_result_dict.items():
                labels = np.array(list(inner_dict.keys())) # which we need later for identification of pdbs
                labels = [os.path.basename(x) for x in labels] #shorten it to pdb_code

                
                print(f"{ranges=}, {labels=}")
                X = np.array(list(inner_dict.values())) #pc1 and pc2 as tuples 

                max_samples = X.shape[0]
                if automated:
                    optimal_clust_num = self._obtain_optimal_clusters(X, max_samples)
                    if optimal_clust_num == 0:
                        print(f"No clusters found for {ranges=}")
                        print("inside optimal clust -")
                        print(labels)
                        self.cluster_results[ranges] = (labels, 0, X)
                        continue
                        
                    kmeans = KMeans(n_clusters=optimal_clust_num, random_state=7)    
                    clust_labels = kmeans.fit_predict(X)
                    
                    self.cluster_results[ranges] = (labels, clust_labels, X)

                if OPTICS:
                    if optics_min_samples:
                        optics = OPTICS(min_samples=optics_min_samples) # think how to do that properly
                        optics.fit(X)
                        clust_labels = optics.labels_
                        self.cluster_results[ranges] = (labels, clust_labels, X) #pdb codes , the predicted class , the PC1 PC2 array for plotting.
            
            return self.cluster_results


    def get_representative_structures(self):
        
        if self.cluster_results:
            for ranges, lbl_clust_X in self.cluster_results.items():
                labels, clust_labels, X = lbl_clust_X # unpack tuple.

                # now we get the labels which is the pdb code, clust_labels which is assigned class and X which is a 2 member arr.. PC1 and PC2
                tmp_dict = dict()
                
                if np.all(clust_labels) == None:
                    self.representative_strucs_per_cluster[ranges] = dict(enumerate(labels))
                    continue
                    
                for clust_label in np.unique(clust_labels):

                    cluster_idx = np.where(clust_labels == clust_label)[0]  #idx of those hits that have the same class label

                    cluster_X = X[cluster_idx]

                    # compute mean
                    mean_PC = np.mean(cluster_X, axis=0)
                    
                    distances = np.linalg.norm(cluster_X - mean_PC, axis=1)  # col wise

                    closest_idx = np.argmin(distances) # get idx of smallest dist.

                    representative_structure = labels[cluster_idx[closest_idx]]
                    
                    tmp_dict[clust_label+1] = representative_structure


                self.representative_strucs_per_cluster[ranges] = tmp_dict

    
    def plot_clusters(self, num_labels_to_show=None):
        if self.cluster_results:
            for ranges, lbl_clust_X in self.cluster_results.items():
    
                expl_variance = self.pca_explained_var[ranges]
                
                labels, clust_labels, X = lbl_clust_X  #unpack the tuple
    
                plt.figure(figsize=(10,8))
    
                # Plot each cluster with its own label and color
                unique_labels = np.unique(clust_labels)
                for label in unique_labels:
                    # Select data points belonging to the current cluster label
                    cluster_points = X[clust_labels == label]
                    # Plot these points with a specific label for the legend
                    plt.scatter(cluster_points[:, 0], cluster_points[:, 1], label=f'Cluster {label+1}')
    
                plt.title('Clustering')
                plt.xlabel('PC1')
                plt.ylabel('PC2')
    
                if num_labels_to_show:
                    for i, label in enumerate(labels):
                        if i % num_labels_to_show == 0:
                            plt.annotate(label, (X[i, 0], X[i, 1]), textcoords="offset points", xytext=(0,10), ha='center')
                
                plt.legend()
                prot_name = os.path.basename(self.work_dir)
                save_path = os.path.join(self.work_dir, f"{prot_name}_{ranges}_cluster_results.pdf")
    
                plt.title(f'{prot_name}')
                plt.xlabel(f'Principal Component 1 {np.round(float(expl_variance[0]),4)*100:.2f}%')
                plt.ylabel(f'Principal Component 2 {np.round(float(expl_variance[1]),4)*100:.2f}%')
                plt.savefig(save_path)
                plt.show()


    def _obtain_optimal_clusters(self, X, max_samples, visualize=False):
        
        #more than 10 conformers are unrealistic
        
        max_classes = min(12, max_samples)
                
        cluster_range = range(2, max_classes)
        best_num_clusters = 0
        best_silhouette_score = -1
        elbow_scores = []
        silhouette_scores = []
        for n_clusters in cluster_range:
            kmeans = KMeans(n_clusters=n_clusters, random_state=7)
            cluster_labels = kmeans.fit_predict(X)
            elbow_scores.append(kmeans.inertia_) # the lower the bett
            #silhouette score
            silhouette_avg = silhouette_score(X, cluster_labels)
            silhouette_scores.append(silhouette_avg)
            print(f"{silhouette_avg=}, {n_clusters=}")
            #if silhouette is better:
            if silhouette_avg > best_silhouette_score:
                best_silhouette_score = silhouette_avg
                best_num_clusters = n_clusters 

        #print(silhouette_scores)
        if visualize:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 6))
            # Elbow plot on the first subplot
            ax1.plot(cluster_range, elbow_scores, 'bo-')
            ax1.set_title('Elbow Method For Optimal Number of Clusters')
            ax1.set_xlabel('Number of Clusters')
            ax1.set_ylabel('Inertia')
            
            # Silhouette scores plot on the second subplot
            ax2.plot(cluster_range, silhouette_scores, 'bo-')
            ax2.set_title('Silhouette Score For Optimal Number of Clusters')
            ax2.set_xlabel('Number of Clusters')
            ax2.set_ylabel('Silhouette Score')
            
            # Display the subplots
            plt.tight_layout()  # Adjust subplots to fit into figure area.
            plt.show()

        return best_num_clusters
