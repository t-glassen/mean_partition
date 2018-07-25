package meanPartition;

import java.util.Arrays;

/**
 * Class for dynamically calculating the distance between two partitions.
 * The algorithms use the implementation of the Hungarian method of Kevin L. Stern (2012).
 * 
 * @author Thomas J. Glassen
 */

/* Copyright (c) 2012 Kevin L. Stern
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/* Stern's original description of his implementation:
* "An implementation of the Hungarian algorithm for solving the assignment
* problem. An instance of the assignment problem consists of a number of
* workers along with a number of jobs and a cost matrix which gives the cost of
* assigning the i'th worker to the j'th job at position (i, j). The goal is to
* find an assignment of workers to jobs so that no job is assigned more than
* one worker and so that no worker is assigned to more than one job in such a
* manner so as to minimize the total cost of completing the jobs.
* <p>
* 
* An assignment for a cost matrix that has more workers than jobs will
* necessarily include unassigned workers, indicated by an assignment value of
* -1; in no other circumstance will there be unassigned workers. Similarly, an
* assignment for a cost matrix that has more jobs than workers will necessarily
* include unassigned jobs; in no other circumstance will there be unassigned
* jobs. For completeness, an assignment for a square cost matrix will give
* exactly one unique worker to each job.
* <p>
* 
* This version of the Hungarian algorithm runs in time O(n^3), where n is the
* maximum among the number of workers and the number of jobs."
*/

public class DynamicPartitionDistance {
  private int maxExpectedNumOfClusters;
  private double[][] costMatrix;
  private int rows, cols, dim;
  private double[] labelByWorker, labelByJob;
  private int[] minSlackWorkerByJob;
  private double[] minSlackValueByJob;
  private int[] matchJobByWorker, matchWorkerByJob;
  private int[] parentWorkerByCommittedJob;
  private boolean[] committedWorkers;
  
  private int[][] partitions;
  private int[] numOfClusters;
  private int[][] clusterSizes;
  
  private int[] finalRowDecrease;
  private int finalDecrease;
  
  private double initialNewWorkerLabel;
  
  private int firstEmptyCluster;
  private int backup_firstEmptyCluster;
  
  private int currentCosts;
  private int backup_currentCosts;
  
  private double[] backup_labelByWorker, backup_labelByJob, backup_rowCostMatrix;
  private int backup_rowIndex, backup_objIndex;
  private int backup_clusterNr;
  private int backup_finalDecrease,backup_finalRowDecrease;
  private int backup_numOfClusters;
  private int[] backup_clusterSizes;
  private int[] backup_matchJobByWorker, backup_matchWorkerByJob;
  private int backup_cols;
  private double backup_initialNewWorkerLabel;
  
  
  public void saveState(int index){
	  backup_rowIndex = partitions[0][index] - 1;
	  backup_objIndex = index;
	  backup_cols = cols; 
	  for(int i = 0;i< cols + 1;i++){ //there is at most one additional cluster till the call of restoreBackup
		  backup_labelByWorker[i] = labelByWorker[i];
		  backup_labelByJob[i] = labelByJob[i];
		  backup_rowCostMatrix[i] = costMatrix[backup_rowIndex][i];		  
		  backup_matchJobByWorker[i] = matchJobByWorker[i];
		  backup_matchWorkerByJob[i] = matchWorkerByJob[i];
		  backup_clusterSizes[i] = clusterSizes[1][i];
	  }
	  backup_clusterNr = partitions[1][index];
	  backup_finalDecrease = finalDecrease;
	  backup_finalRowDecrease = finalRowDecrease[backup_rowIndex];
	  backup_numOfClusters = numOfClusters[1];
	  backup_currentCosts = currentCosts;
	  backup_firstEmptyCluster = firstEmptyCluster;
	  backup_initialNewWorkerLabel = initialNewWorkerLabel;
  }
  
  public void restoreState(){
	  cols = backup_cols;
	  dim = Math.max(rows, cols);
	  for(int i = 0;i< cols + 1;i++){
		  labelByWorker[i] = backup_labelByWorker[i];
		  labelByJob[i] = backup_labelByJob[i];
		  costMatrix[backup_rowIndex][i] = backup_rowCostMatrix[i];
		  matchJobByWorker[i] = backup_matchJobByWorker[i];
		  matchWorkerByJob[i] = backup_matchWorkerByJob[i];
		  clusterSizes[1][i] = backup_clusterSizes[i];
	  }
	  partitions[1][backup_objIndex] = backup_clusterNr;
	  finalDecrease = backup_finalDecrease;
	  finalRowDecrease[backup_rowIndex] = backup_finalRowDecrease;
	  numOfClusters[1] = backup_numOfClusters;
	  currentCosts = backup_currentCosts;
	  firstEmptyCluster = backup_firstEmptyCluster;
	  initialNewWorkerLabel = backup_initialNewWorkerLabel;
  }
  
  private void increaseSizeOfArrays() {
	  
		//set the maximum expected number of clusters. Keep this low to avoid heap size errors
		maxExpectedNumOfClusters = Math.min(Math.max(numOfClusters[0], numOfClusters[1] + 1) * 2, partitions[1].length + 1);
	  
	  	double[] temp = new double[maxExpectedNumOfClusters];
	  	int[] temp2 = new int[maxExpectedNumOfClusters];
	  	boolean[] temp3 = new boolean[maxExpectedNumOfClusters];
	  	double[][] newCostMatrix = new double[maxExpectedNumOfClusters][maxExpectedNumOfClusters];
	  	int[][] newClusterSizes = new int[2][maxExpectedNumOfClusters];
	  	
	  	for(int i = 0; i < costMatrix.length;i++) {
	  		newClusterSizes[0][i] = clusterSizes[0][i];
	  		newClusterSizes[1][i] = clusterSizes[1][i];
	  		temp2[i] = backup_clusterSizes[i];
	  		for(int j = 0; j < costMatrix[0].length;j++)
	  			newCostMatrix[i][j] = costMatrix[i][j];
	  	}
	  	costMatrix = newCostMatrix; clusterSizes = newClusterSizes; backup_clusterSizes = temp2.clone();		
	  	
	  	System.arraycopy(backup_rowCostMatrix, 0, temp, 0, backup_rowCostMatrix.length);
	  	backup_rowCostMatrix = temp.clone();
	  	System.arraycopy(labelByWorker, 0, temp, 0, labelByWorker.length);
	  	labelByWorker = temp.clone();	  	
	  	System.arraycopy(backup_labelByWorker, 0, temp, 0, backup_labelByWorker.length);
	  	backup_labelByWorker = temp.clone();	   
	  	System.arraycopy(labelByJob, 0, temp, 0, labelByJob.length);
	  	labelByJob = temp.clone();
	  	System.arraycopy(backup_labelByJob, 0, temp, 0, backup_labelByJob.length);
	  	backup_labelByJob = temp.clone();

	  	minSlackWorkerByJob = temp2.clone();
	  	minSlackValueByJob = temp.clone();
	  	committedWorkers = temp3.clone();
	  	parentWorkerByCommittedJob = temp2.clone();
	  	
	  	Arrays.fill(temp2, -1);
	  	System.arraycopy(matchJobByWorker, 0, temp2, 0, matchJobByWorker.length);
	  	matchJobByWorker = temp2.clone();
	  	System.arraycopy(backup_matchJobByWorker, 0, temp2, 0, backup_matchJobByWorker.length);
	  	backup_matchJobByWorker = temp2.clone();
	  	System.arraycopy(matchWorkerByJob, 0, temp2, 0, matchWorkerByJob.length);
	  	matchWorkerByJob = temp2.clone();
	  	System.arraycopy(backup_matchWorkerByJob, 0, temp2, 0, backup_matchWorkerByJob.length);
	  	backup_matchWorkerByJob = temp2.clone();	  
  }
  
  public DynamicPartitionDistance(int[] partition1,int[] partition2){
		
		//build cost matrix
  		buildCostMatrix(partition1, partition2);
  		
  		firstEmptyCluster = numOfClusters[1] + 1;
  		
  		finalDecrease = 0;
  		initialNewWorkerLabel = 0;
  	      	    
  	    finalRowDecrease = new int[numOfClusters[0]];
  	    backup_rowCostMatrix = new double[maxExpectedNumOfClusters];
  	    labelByWorker = new double[maxExpectedNumOfClusters];
  	    backup_labelByWorker = labelByWorker.clone();
  	    labelByJob = new double[maxExpectedNumOfClusters];
  	    backup_labelByJob = labelByJob.clone();
  	    minSlackWorkerByJob = new int[maxExpectedNumOfClusters];
  	    minSlackValueByJob = new double[maxExpectedNumOfClusters];
  	    committedWorkers = new boolean[maxExpectedNumOfClusters];
  	    parentWorkerByCommittedJob = new int[maxExpectedNumOfClusters];
  	    matchJobByWorker = new int[maxExpectedNumOfClusters];
  	    backup_matchJobByWorker = matchJobByWorker.clone();
  	    Arrays.fill(matchJobByWorker, -1);
  	    matchWorkerByJob = new int[maxExpectedNumOfClusters];
  	    backup_matchWorkerByJob = matchWorkerByJob.clone();
  	    Arrays.fill(matchWorkerByJob, -1);
   
	}
	
	//this function awaits the mean partition as second partition
	private void buildCostMatrix(int[] partition1, int[] partition2){
		
	    int[][] part = {partition1.clone(),partition2.clone()};
	    partitions = part;	    
	    int[][] clusterIndices = new int[2][partition1.length + 1];
		numOfClusters = new int[2];
						
		//count number of clusters in partition i and assign new cluster numbers starting from 1
		for(int i=0;i<2;i++)
			for(int x=0;x<partitions[0].length;x++){
				if (clusterIndices[i][partitions[i][x]] == 0){
					numOfClusters[i]++;
					clusterIndices[i][partitions[i][x]] = numOfClusters[i];				
				}		
				partitions[i][x] = clusterIndices[i][partitions[i][x]];			
			}
		
		//set the maximum expected number of clusters. Keep this low to avoid heap size errors
		maxExpectedNumOfClusters = Math.min(Math.max(numOfClusters[0], numOfClusters[1] + 1) * 2, partitions[1].length + 1);
				
		costMatrix = new double[maxExpectedNumOfClusters][maxExpectedNumOfClusters];
		rows = numOfClusters[0];
	    dim = Math.max(rows, numOfClusters[1]);
	    cols = dim;
	    
		clusterSizes = new int[2][maxExpectedNumOfClusters];
		backup_clusterSizes = clusterSizes[1].clone();		
		
		//count the size of the clusters in both partitions and decrement every cell in the cost matrix
		//where an item is shared between the corresponding clusters
		for(int x=0;x<partitions[0].length;x++){						
			int i = partitions[0][x] - 1;
			int j = partitions[1][x] - 1;
			clusterSizes[0][i]++;
			clusterSizes[1][j]++;
			costMatrix[i][j]--;
		}
		
		//add the cluster size of cluster i to every cell in row i
		for(int i=0;i<rows;i++)
			for(int j=0;j<cols;j++)
				costMatrix[i][j] += clusterSizes[0][i];
	}
	  
	  /**
	   * Move object in partition 2 with index index to cluster k
	   * 
	   * @return the partition distance 
	   */
	  public int movObjTo(int index, int k){
		  	  
		  int i = partitions[0][index] - 1;
		  int j = partitions[1][index] - 1;
		  
		  if(j == k - 1)
			  return currentCosts;
		  
		  //update cluster sizes
		  clusterSizes[1][j]--;
		  clusterSizes[1][k-1]++;
		  
		  //update number of clusters
		  if(clusterSizes[1][k-1] == 1){ //a new cluster is needed
			  numOfClusters[1]++;
			  
			  //if necessary update the relevant dimensions of the matrix and its costs
			  //Runtime: O(C)
			  if(k > cols){
				  dim = cols = k;
				  
				  //set labels of new worker and job and match them
				  labelByWorker[cols-1] = initialNewWorkerLabel;
				  labelByJob[cols-1] = initialNewWorkerLabel * -1;
				  match(cols-1,cols-1); //match last row and last column directly to keep needed reconnections as minimal as possible
				  			  
				  //update the costs for the new cluster column
				  //(a new higher k is always an increment by one to cols)			  
				  for(int m=0;m<rows;m++)
					  costMatrix[m][cols - 1] = clusterSizes[0][m] + finalRowDecrease[m];			  
			  }
		  }
		  if(clusterSizes[1][j] == 0){ //a cluster became empty
			  numOfClusters[1]--;
			  
			  //keep cluster nr if j + 1 is lower than firstEmptyCluster
			  if(firstEmptyCluster > j + 1)
			  	firstEmptyCluster = j + 1;
			  
			  //if necessary update cols by searching for the next non empty cluster below cols - 1
			  //Runtime: O(C)
			  if(j == cols - 1){			  
				  for(int m = cols - 1;m >= rows;m--){
					 if(m == cols - 1){
						  unmatchWorker(cols-1); //the squared cost matrix loses its last row -> free worker and matched job
					  	  unmatchJob(cols-1);
					  } else if(labelByWorker[m] + labelByJob[m] == costMatrix[m][m]) {
						  match(matchWorkerByJob[m],matchJobByWorker[m]); //cross connection to keep needed reconnections as minimal as possible
						  match(m,m);
					  }					
					  cols = m;
					  if(clusterSizes[1][m-1] > 0)
						  break;
				  }
				  dim = cols;
			  }
		  }
		  
		  partitions[1][index] = k;
		  
		  //update costs after moving an object. Runtime: O(C)
		  unmatchWorker(i);	  
		  costMatrix[i][j]++;     				
		  for(int m=0;m<this.dim;m++)
			  if(m != k - 1)
				  costMatrix[i][m]++;
		  finalRowDecrease[i]++;	// keep track of the additional costs for new clusters in that row
		  finalDecrease++;

		  //reconnect the disconnected worker. Runtime: O(C^2)
		  currentCosts = calcDist(proceed()) - finalDecrease;	  	  
		  return currentCosts;
	  }
	  
	  /**
	   * Move object in partition 2 with index index to new cluster 
	   * 
	   * @return the partition distance 
	   */
	  public int movObjToNew(int index){			  

		  //get the next empty cluster. Runtime: O(C)
		  int emptyClusterNr = firstEmptyCluster;
		  firstEmptyCluster++;
		  
		  while(firstEmptyCluster <= clusterSizes[1].length && clusterSizes[1][firstEmptyCluster - 1] > 0)
			  firstEmptyCluster++;
		  
		  //increase size of arrays if maxExpectedNumOfClusters is nearly reached
		  if(firstEmptyCluster == maxExpectedNumOfClusters)
			  increaseSizeOfArrays();
		  
		  return movObjTo(index, emptyClusterNr);
	  }
	  
	  /**
	   * Provides the cluster index of an object in partition 2
	   * 
	   * @return the cluster index of obj 
	   */
	  public int getClusterNrOfObj(int index){
		  return partitions[1][index];
	  }
	  
	  public int getSizeOfClusterOfObj(int index){
		  return clusterSizes[1][partitions[1][index] - 1];
	  }
	  
	  public int[] getPartition(){
		  return partitions[1];
	  }
	  
	  /**
	   * calculates the unique ascending cluster indices of partition 2
	   * 
	   * @return the ascending cluster indices
	   */
	  public int[] getAscendingClusterNrs(){
		  int counter = numOfClusters[1];
		  int[] clusterNrs = new int[counter];
		  for(int i=0;counter > 0;i++)
			  if(clusterSizes[1][i] > 0){
				  clusterNrs[clusterNrs.length - counter] = i+1;
				  counter--;
			  }
		  return clusterNrs;
	  }
	  	  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The subsequent code is mainly from Kevin L. Stern (2012) with some modifications 

  /**
   * Compute an initial feasible solution by assigning zero labels to the
   * workers and by assigning to each job a label equal to the minimum cost
   * among its incident edges.
   */
  protected void computeInitialFeasibleSolution() {
    for (int j = 0; j < dim; j++) {
      labelByJob[j] = Double.POSITIVE_INFINITY;
    }
    for (int w = 0; w < dim; w++) {
      for (int j = 0; j < dim; j++) {
        if (costMatrix[w][j] < labelByJob[j]) {
          labelByJob[j] = costMatrix[w][j];
        }
      }
    }
  }
  
  protected int calcDist(int[] jobIndexForWorker){
	  int DistanceSum = 0;
	  for(int worker=0;worker < jobIndexForWorker.length;worker++)
			DistanceSum += costMatrix[worker][jobIndexForWorker[worker]];

	  return DistanceSum;
  }

  /**
   * Execute the algorithm.
   * 
   * @return the minimum cost matching of workers to jobs based upon the
   *         provided cost matrix. A matching value of -1 indicates that the
   *         corresponding worker is unassigned.
   */ 
  public int get(){	  
	  greedyMatch();
	  currentCosts = calcDist(proceed()) - finalDecrease;
	  return currentCosts;
  }
  
  /**
   * Continue execution of the algorithm.
   * 
   * @return the minimum cost 
   */
  private int[] proceed() {
	  int w = fetchUnmatchedWorker();
	  while (w < dim) {
		  initializePhase(w);
		  executePhase();
	      w = fetchUnmatchedWorker();
	  }
	  int[] result = Arrays.copyOf(matchJobByWorker, rows);
	  for (w = 0; w < result.length; w++) {
		  if (result[w] >= cols) {
			  result[w] = -1;
		  }
	  }
	  return result;
  }
  
  /**
   * Execute a single phase of the algorithm. A phase of the Hungarian algorithm
   * consists of building a set of committed workers and a set of committed jobs
   * from a root unmatched worker by following alternating unmatched/matched
   * zero-slack edges. If an unmatched job is encountered, then an augmenting
   * path has been found and the matching is grown. If the connected zero-slack
   * edges have been exhausted, the labels of committed workers are increased by
   * the minimum slack among committed workers and non-committed jobs to create
   * more zero-slack edges (the labels of committed jobs are simultaneously
   * decreased by the same amount in order to maintain a feasible labeling).
   * <p>
   * 
   * The runtime of a single phase of the algorithm is O(n^2), where n is the
   * dimension of the internal square cost matrix, since each edge is visited at
   * most once and since increasing the labeling is accomplished in time O(n) by
   * maintaining the minimum slack values among non-committed jobs. When a phase
   * completes, the matching will have increased in size.
   */
  protected void executePhase() {
    while (true) {
      int minSlackWorker = -1, minSlackJob = -1;
      double minSlackValue = Double.POSITIVE_INFINITY;
      for (int j = 0; j < dim; j++) {
        if (parentWorkerByCommittedJob[j] == -1) {
          if (minSlackValueByJob[j] < minSlackValue) {
            minSlackValue = minSlackValueByJob[j];
            minSlackWorker = minSlackWorkerByJob[j];
            minSlackJob = j;
          }
        }
      }
      if (minSlackValue > 0) {
        updateLabeling(minSlackValue);
      }
      parentWorkerByCommittedJob[minSlackJob] = minSlackWorker;
      if (matchWorkerByJob[minSlackJob] == -1) {
        /*
         * An augmenting path has been found.
         */
        int committedJob = minSlackJob;
        int parentWorker = parentWorkerByCommittedJob[committedJob];
        while (true) {
          int temp = matchJobByWorker[parentWorker];
          match(parentWorker, committedJob);
          committedJob = temp;
          if (committedJob == -1) {
            break;
          }
          parentWorker = parentWorkerByCommittedJob[committedJob];
        }
        return;
      } else {
        /*
         * Update slack values since we increased the size of the committed
         * workers set.
         */
        int worker = matchWorkerByJob[minSlackJob];
        committedWorkers[worker] = true;
        for (int j = 0; j < dim; j++) {
          if (parentWorkerByCommittedJob[j] == -1) {
            double slack = costMatrix[worker][j] - labelByWorker[worker]
                - labelByJob[j];
            if (minSlackValueByJob[j] > slack) {
              minSlackValueByJob[j] = slack;
              minSlackWorkerByJob[j] = worker;
            }
          }
        }
      }
    }
  }

  /**
   * 
   * @return the first unmatched worker or {@link #dim} if none.
   */
  protected int fetchUnmatchedWorker() {
    int w;
    for (w = 0; w < dim; w++) {
      if (matchJobByWorker[w] == -1) {
        break;
      }
    }
    return w;
  }

  /**
   * Find a valid matching by greedily selecting among zero-cost matchings. This
   * is a heuristic to jump-start the augmentation algorithm.
   */
  protected void greedyMatch() {
    for (int w = 0; w < dim; w++) {
      for (int j = 0; j < dim; j++) {
        if (matchJobByWorker[w] == -1 && matchWorkerByJob[j] == -1
            && costMatrix[w][j] - labelByWorker[w] - labelByJob[j] == 0) {
          match(w, j);
        }
      }
    }
  }

  /**
   * Initialize the next phase of the algorithm by clearing the committed
   * workers and jobs sets and by initializing the slack arrays to the values
   * corresponding to the specified root worker.
   * 
   * @param w
   *          the worker at which to root the next phase.
   */
  protected void initializePhase(int w) {
    Arrays.fill(committedWorkers, false);
    Arrays.fill(parentWorkerByCommittedJob, -1);
    committedWorkers[w] = true;
    for (int j = 0; j < dim; j++) {    	
      minSlackValueByJob[j] = costMatrix[w][j] - labelByWorker[w]
          - labelByJob[j];
      minSlackWorkerByJob[j] = w;
    }
  }

  /**
   * Helper method to record a matching between worker w and job j.
   */
  protected void match(int w, int j) {
    matchJobByWorker[w] = j;
    matchWorkerByJob[j] = w;
  }

  /**
   * Helper method to umatch worker w and job j.
   */
  protected void unmatchWorker(int w) {
	  if(matchJobByWorker[w] != - 1)
		  matchWorkerByJob[matchJobByWorker[w]] = -1;
	  matchJobByWorker[w] = -1;
  }
  
  /**
   * Helper method to umatch job j and worker w.
   */
  protected void unmatchJob(int j) {
	  if(matchWorkerByJob[j] != - 1)
		  matchJobByWorker[matchWorkerByJob[j]] = -1;
	matchWorkerByJob[j] = - 1;
  }
  
  /**
   * Reduce the cost matrix by subtracting the smallest element of each row from
   * all elements of the row as well as the smallest element of each column from
   * all elements of the column. Note that an optimal assignment for a reduced
   * cost matrix is optimal for the original cost matrix.
   */
  protected void reduce() {
    for (int w = 0; w < dim; w++) {
      double min = Double.POSITIVE_INFINITY;
      for (int j = 0; j < dim; j++) {
        if (costMatrix[w][j] < min) {
          min = costMatrix[w][j];
        }
      }
      for (int j = 0; j < dim; j++) {
        costMatrix[w][j] -= min;
      }
    }
    double[] min = new double[dim];
    for (int j = 0; j < dim; j++) {
      min[j] = Double.POSITIVE_INFINITY;
    }
    for (int w = 0; w < dim; w++) {
      for (int j = 0; j < dim; j++) {
        if (costMatrix[w][j] < min[j]) {
          min[j] = costMatrix[w][j];
        }
      }
    }
    for (int w = 0; w < dim; w++) {
      for (int j = 0; j < dim; j++) {
        costMatrix[w][j] -= min[j];
      }
    }
  }

  /**
   * Update labels with the specified slack by adding the slack value for
   * committed workers and by subtracting the slack value for committed jobs. In
   * addition, update the minimum slack values appropriately.
   */
  protected void updateLabeling(double slack) {
    for (int w = 0; w < dim; w++) {
      if (committedWorkers[w]) {
        labelByWorker[w] += slack;
        
        //keep highest difference between worker label and the sum of its cluster size and its row decrease
        double critDiff = labelByWorker[w] - clusterSizes[0][w];
        if(w < finalRowDecrease.length)
        		critDiff -= finalRowDecrease[w];
        if(initialNewWorkerLabel < critDiff)
        	initialNewWorkerLabel = critDiff;
      }
    }
    for (int j = 0; j < dim; j++) {
      if (parentWorkerByCommittedJob[j] != -1) {
        labelByJob[j] -= slack;
      } else {
        minSlackValueByJob[j] -= slack;
      }
    }
  }
}
