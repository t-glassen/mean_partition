package meanPartition;

/**
 * Provides various methods for partition distributions
 * and consensus clustering.
 * 
 * @author Thomas J. Glassen
 */

import java.util.SplittableRandom;

public class PartitionDistribution {
	
	private int[][] partitions;
	private int[] initialMP;
	private SplittableRandom r;
    
	// constructors
	/////////////////////////////////////////////////////////////////////
	
	//create an empty instance
	public PartitionDistribution() {
		partitions = null;
		initialMP = null;
		initRNG();
	}

	public PartitionDistribution(int[][] partitions){
		setNewPartitions(partitions);
		initRNG();
	}
	
	//generates K random partitions of size N with C clusters
	public PartitionDistribution(int K, int N, int C) {
		initRNG();
		setNewPartitions(K,N,C);
	}
	
	// initialization routines
	/////////////////////////////////////////////////////////////////////
	
	private void initRNG() {
		long seed = System.currentTimeMillis();
		System.out.println("[PartitionDistribution] seed is: " + seed);		
		r = new SplittableRandom(seed);
	}
	
	public void setNewPartitions(int[][] partitions) {
		this.partitions = new int[partitions.length][partitions[0].length];
		for(int i = 0; i < partitions.length;i++)
			this.partitions[i] = partitions[i].clone();
		
		initialMP = partitions[0].clone();
	}
	
	//generates K random partitions of size N with C clusters
	public void setNewPartitions(int K, int N, int C) {
		partitions = new int[K][N];		
		for(int i = 0; i < K; i++)
			partitions[i] = generateRandomPartition(N,C);
		
		initialMP = partitions[0].clone();
	}

	//sets the starting point for the local search procedures
	//(the initial mean partition)
	public void setInitialMP(int[] partition) {
		initialMP = partition.clone();
	}
	
	// random partition generation
	/////////////////////////////////////////////////////////////////////
	
	//generates a random partition of length N with C clusters
	public int[] generateRandomPartition(int N, int C) {
		int[] partition = new int[N];
		int numOfClusters = 0;
		while(numOfClusters != C) {
			numOfClusters = 0;
			int[] clusterSizes = new int[C];
			for(int i = 0; i < N; i++) {
				partition[i] = r.nextInt(C);
				if(++clusterSizes[partition[i]] == 1)
					numOfClusters++;
			}
		}
		return partition;
	}	
	
	// routines for gathering local search procedure runtime information
	/////////////////////////////////////////////////////////////////////
	
	private int[] numOfIterations;
	private int currentRun;
	
	public void prepareArrayToKeepIterations(int numOfRuns) {
		numOfIterations = new int[numOfRuns];
		currentRun = 0;
	}
	
	public int[] getNumOfIterations() {
		return numOfIterations;
	}
		
	// routines for calculating sums or averages of similarity to a given partition
	/////////////////////////////////////////////////////////////////////
	
	/*
	 * returns the sum of the specified similarity of every available partition to p
	 * runtime: O(K * (N + C^3))
	 */
	public double getSumOfSimilarities(int[] p, int measure){
		double sumOfSimilarities = 0;
		for(int i=0;i<partitions.length;i++)
			sumOfSimilarities += getSimilarity(p,partitions[i], measure);
		
		return sumOfSimilarities;
	}
	
	/*
	 * returns the average similarity to p
	 * runtime: O(K * (N + C^3))
	 */
	public double getAverageSimilarity(int[] p, int measure) {
		return getSumOfSimilarities(p, measure) / (double) partitions.length;
	}
	
	/*
	 * returns the similarity between two partitions calculated via a specified measure
	 * runtime: O(N + C^3)
	 */
	public double getSimilarity(int[] p1, int[] p2, int measure) {
		if(measure == 0)
			return getPartitionDistance(p1,p2);
		else
			return (p1.length > 1 ? 1 - getPartitionDistance(p1,p2) / (double) (p1.length - 1) : 1);
	}
	
	/*	Calculation of the partition distance via a reduction of the problem 
	 *  to a linear sum assignment problem according to Konovalov, Litow & Bajema
	 *  (Partition-distance via the assignment problem; 2005, S. 2464)
	 *  with the faster reduction of Glassen
	 * (Psychologisch orientierte Kategorisierung in der kognitiven Robotik mit dem Hierarchischen Dirichlet Prozess; 2018)
	 *  runtime: O(N + C^3)
	 */ 
	public int getPartitionDistance(int[] partition1, int[] partition2){
				
		int numOfClustersP1 = 0, numOfClustersP2 = 0;
		int[] clusterIndicesP1 = new int[partition1.length + 1], clusterIndicesP2 = new int[partition2.length + 1];
		
		//count number of clusters in both partitions and assign new cluster numbers starting from 1
		for(int x=0;x<partition1.length;x++){
			if (clusterIndicesP1[partition1[x]] == 0){
				numOfClustersP1++;
				clusterIndicesP1[partition1[x]] = numOfClustersP1;
			}
			if (clusterIndicesP2[partition2[x]] == 0){
				numOfClustersP2++;
				clusterIndicesP2[partition2[x]] = numOfClustersP2;
			}
		}
		
		//switch partitions if partiton 1 has more clusters 
		if(numOfClustersP2 < numOfClustersP1){
			int[] tmp = clusterIndicesP1;
			clusterIndicesP1 = clusterIndicesP2;
			clusterIndicesP2 = tmp;
			int tmp2 = numOfClustersP1;
			numOfClustersP1 = numOfClustersP2;
			numOfClustersP2 = tmp2;
		    tmp = partition1;
		    partition1 = partition2;
		    partition2 = tmp;
		}
		
		double[][] costMatrix = new double[numOfClustersP1][numOfClustersP2];
		int[] clusterSizes = new int[numOfClustersP1];
		
		//count the size of the clusters in partition 1 and decrement every cell in the cost matrix
		//where an item is shared between the corresponding clusters
		for(int x=0;x<partition1.length;x++){						
			int i = clusterIndicesP1[partition1[x]] - 1;			
			clusterSizes[i]++;
			costMatrix[i][clusterIndicesP2[partition2[x]] - 1]--;
		}
		
		//add the cluster size of cluster i to every cell in row i
		for(int i=0;i<numOfClustersP1;i++)
			for(int j=0;j<numOfClustersP2;j++)
				costMatrix[i][j] += clusterSizes[i];		
		
		int jobIndexForWorker[] = (new HungarianAlgorithm(costMatrix)).execute();
		
		int overallCosts = 0;
		for(int worker=0;worker < jobIndexForWorker.length;worker++)
			overallCosts += costMatrix[worker][jobIndexForWorker[worker]];
		
		return overallCosts;
	}
		
	/* 
	 * The heuristic method of Huelsenbeck & Andolfatto 
	 * (Inference of Population Structure Under a Dirichlet Process Model; 2007, S.1792)
	 * using the faster reduction of Glassen
	 * (Psychologisch orientierte Kategorisierung in der kognitiven Robotik mit dem Hierarchischen Dirichlet Prozess; 2018)
	 * runtime: O(K*N^2*C + K*N*C*(N + C^3)) per cycle
	 */		
	public int[] getMeanPartition1(){		
		double bestDistanceSum;
		
		DynamicPartition mp = new DynamicPartition(initialMP);
		bestDistanceSum = getSumOfSimilarities(mp.getPartition(), 0);
		
		boolean success = true;
		int numOfIter = 0;
	    while(success){
			numOfIter++;
	        success = false;
	    	for(int obj=0;obj<initialMP.length;obj++){
	    		int[] clusterIndicesMP = mp.getAscendingClusterNrs();
	    		int bestObjClusterIndex = mp.getClusterIndexOfObj(obj), oldObjClusterIndex = bestObjClusterIndex;
	    		int sizeOfClusterOfObj = mp.getSizeOfClusterOfObj(obj);
	    		for(int clusterIndex=0;clusterIndex<clusterIndicesMP.length;clusterIndex++){
	    			if (clusterIndicesMP[clusterIndex] == oldObjClusterIndex){
						if (sizeOfClusterOfObj == 1)
							continue;
						mp.movObjToNew(obj);
	    			}
	    			else
	    				mp.movObjTo(obj,clusterIndicesMP[clusterIndex]); 
 	                
	                double tmpDistanceSum = getSumOfSimilarities(mp.getPartition(), 0); //one can trade space for time via keeping the cost matrizes and updating them constantly  

	                if(tmpDistanceSum < bestDistanceSum){
	                	bestDistanceSum = tmpDistanceSum;
	                	bestObjClusterIndex = mp.getClusterIndexOfObj(obj);
	                	success = true;
	                }
	    		}
	    		mp.movObjTo(obj, bestObjClusterIndex);	    		
	    	}
	    }
		System.out.println("iter (0): " + numOfIter);
		numOfIterations[currentRun++] = numOfIter;
	    return mp.getPartition();
	}
	
	/*
	 * version based on the first improvement of Glassen, von Oertzen & Konovalov (Finding the mean in a partition distribution; unpublished)
	 * runtime: O(K*N*C*(C^2)) per cycle
	 */
	public int[] getMeanPartition2(){
		DynamicPartitionDistance[] partDist = new DynamicPartitionDistance[partitions.length];
		int bestDistanceSum = 0;
		
		for(int p=0;p<partitions.length;p++){
			partDist[p] = new DynamicPartitionDistance(partitions[p], initialMP);
			bestDistanceSum += partDist[p].get(); 
		}
		boolean flag = true;
		int numOfIter = 0;
		while(flag){
			numOfIter++;
			flag = false;
			for(int obj=0;obj<partitions[0].length;obj++){  			
	    		for(int p = 0;p < partitions.length;p++)
	    			partDist[p].saveState(obj);
	    		int[] clusterIndices = partDist[0].getAscendingClusterNrs();	    		
	    		int oldObjClusterNr = partDist[0].getClusterNrOfObj(obj);
	    		int sizeOfClusterOfObj = partDist[0].getSizeOfClusterOfObj(obj);	    		
	    		for(int clusterIndex=0;clusterIndex<clusterIndices.length;clusterIndex++){
	    			int distanceSum = 0;			
	    			if (clusterIndices[clusterIndex] == oldObjClusterNr){ //create new cluster
	    				if (sizeOfClusterOfObj == 1) //object is already in its own cluster
	    					continue;
	    				for(int p = 0;p < partitions.length;p++)
	    					distanceSum += partDist[p].movObjToNew(obj);
	    			}
	    			else { //choose existing cluster
	    				for(int p = 0;p < partitions.length;p++)
	    					distanceSum += partDist[p].movObjTo(obj, clusterIndices[clusterIndex]);
	    			}
	    			
	    			if(distanceSum < bestDistanceSum) {
	    				bestDistanceSum = distanceSum;
	    				flag = true;
	    				for(int p = 0;p < partitions.length;p++)
		    				partDist[p].saveState(obj);
	    			}
	    		}
	    		for(int p = 0;p < partitions.length;p++)
	    			partDist[p].restoreState();
			}
		}
		System.out.println("iter (0): " + numOfIter);
		numOfIterations[currentRun++] = numOfIter;
		return partDist[0].getPartition();
	}
}
