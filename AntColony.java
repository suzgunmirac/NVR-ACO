package aco_sba;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.StringTokenizer;

public class AntColony {

    int peakCount, aaCount, antCount, unavailableValue, alpha, beta, gamma,
        binaryDistance[][], cnstrntPkCnt, nextPeakSelect, constForInfeasible,
        pheromoneUpdateCount, localSearchVersion, localSearchStop, pheromoneReset,
        totalScoreCalcCoef, candidateListValue, iterationCount, elitistAntCount, 
        bestUpdate,                     // iteration number for best assignment update
        scoreConstraints[];             // number of 10e9 values per peak
    double evaporationConst, pheromoneInitialVal, infeasibleScore, localPhrmnUpdtValue,
           combinedScore[][], heuristic[][], pheromone[][], assignmentDelta,
           pheromoneMin, pheromoneMax, heuristicBetaPower[][];
    final double EPSILON = .00001;
    final int RESET_COUNT = 5;
    String paramsFile = "Parameters.txt", binaryFile, combinedScoreFile, noeListFile, resultFile, matlabFile;
    ArrayList<Integer> noeList[], candidateList[], candidateListAA[];
    boolean localPheromoneUpdate, paralelRun, infeasibleSelect, unavailableSelect, 
            outputMatlab, keyBestSet;

    class AssignmentInfo {
        Assignment assignment[];
        double totalScore;
        double accuracy;

        public AssignmentInfo(Assignment[] a, double score, double acc) {
            assignment = a;
            totalScore = score;
            accuracy = acc;
        }
    }

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        try {
            AntColony app = new AntColony();
            int runCount = 0, alpha1 = 0, alpha2 = 0, beta1 = 0, beta2 = 0;
            
            if (args.length > 5) {
                app.paramsFile = args[0];
                runCount = Integer.parseInt(args[1]);
                alpha1 = Integer.parseInt(args[2]);
                alpha2 = Integer.parseInt(args[3]);
                beta1 = Integer.parseInt(args[4]);
                beta2 = Integer.parseInt(args[5]);
            }
            
            app.mainProgram();
            String resFile = app.resultFile,
                   matFile = app.matlabFile;
            if (runCount > 1) {
                for (int i = 1; i <= runCount; i++)
                for (app.alpha = alpha1; app.alpha < alpha2; app.alpha++)
                for (app.beta = beta1; app.beta < beta2; app.beta++) {
                    app.initializePheromoneTrail();
                    app.resultFile = resFile + i + "_A" + app.alpha + "_B" + app.beta + ".txt";
                    app.matlabFile = matFile + i + ".txt";
                    app.runAnts();
                }
            } else {
                for (app.alpha = alpha1; app.alpha < alpha2; app.alpha++)
                for (app.beta = beta1; app.beta < beta2; app.beta++) {
                    app.initializePheromoneTrail();
                    app.resultFile = resFile + "_A" + app.alpha + "_B" + app.beta + ".txt";
                    app.matlabFile = matFile + ".txt";
                    app.runAnts();
                }
            }
        } catch (Exception e) { e.printStackTrace(System.out); }
    }

    private void mainProgram() {
        try {
            readParameters();
            
            // parameters for pheromone reinitializing after some iteration. if -1 no pheromone reinitializing
            if (pheromoneReset != -1)
                pheromoneReset = iterationCount * pheromoneReset / 100;

            // reading combined score file
            if (readCombinedScore() < 0) {
                System.out.println("Error in reading file " + combinedScoreFile);
            } else {
                infeasibleScore = getInfeasibleScore();
                if (infeasibleScore == -1.0) 
                    throw new Exception("infeasible score calculation failed");
                else
                    System.out.println("Infeasible score is " + infeasibleScore);
                setHeuristicValues();
                scoreConstraints = new int[peakCount];
                for (int i = 0; i < peakCount; i++)
                for (int j = 0; j < aaCount; j++)
                    if (combinedScore[i][j] == unavailableValue) scoreConstraints[i]++;
            }

            if (readBinaryDistance() < 0) {
                System.out.println("Error in reading file " + binaryFile);
            }
            if (readNoeList() < 0) {
                System.out.println("Error in reading file " + noeListFile);
            }

            if (pheromoneInitialVal <= 0) {     // set initial pheromone value according to "parameters.txt"
                if (pheromoneInitialVal == -1) {
                    pheromoneInitialVal = 1 / ((1 - evaporationConst) * peakCount * getCombinedScoreAverage());
                }
                if (pheromoneInitialVal == -2) {
                    pheromoneInitialVal = 1 / getCombinedScoreAverage();
                }
            }

            pheromone = new double[peakCount][aaCount];
            initializePheromoneTrail();         // setting pheromone values to initial value

            int i, j;
            int[] temp = new int[aaCount];
            double[] tempScore = new double[aaCount];
            candidateList = new ArrayList[peakCount];
            for (i = 0; i < peakCount; i++) {
                for (j = 0; j < aaCount; j++) {
                    temp[j] = j;
                    tempScore[j] = combinedScore[i][j];
                }
                quicksort(tempScore, temp, 0, aaCount - 1);
                candidateList[i] = new ArrayList<Integer>(aaCount);
                for (j = 0; j < aaCount; j++) {
                    candidateList[i].add(temp[j]);
                }
            }

            ArrayList<Integer> cand;
            candidateListAA = new ArrayList[aaCount];
            for (i = 0; i < aaCount; i++) {
                candidateListAA[i] = new ArrayList<Integer>(peakCount);
            }
            for (i = 0; i < peakCount; i++) {
                cand = candidateList[i];
                for (j = 0; j < cand.size(); j++) {
                    candidateListAA[cand.get(j)].add(i);
                }
            }

            combinedScoreAnalysis();
            binaryAnalysis();
            noeListAnalysis();
            System.out.println("Diagonal score is " + calculateDiagonalScore());
        } catch (Exception e) { e.printStackTrace(System.out); }
    }
    
    private void readParameters() {
        try {
            FileReader fr;
            BufferedReader br;
            fr = new FileReader(paramsFile);
            br = new BufferedReader(fr);
            String s, str;
            StringTokenizer st;
            while ((s = br.readLine()) != null) {
               st = new StringTokenizer(s);
               str = st.nextToken();
               if (str.compareTo("peak") == 0) {
                   peakCount = Integer.parseInt(st.nextToken());
                   System.out.println(str + " = " + peakCount);
               }
               if (str.compareTo("aminoacid") == 0) {
                   aaCount = Integer.parseInt(st.nextToken());
                   System.out.println(str + " = " + aaCount);
               }
               if (str.compareTo("ant") == 0) {
                   antCount = Integer.parseInt(st.nextToken());
                   System.out.println(str + " = " + antCount);
               }
               if (str.compareTo("evaporation") == 0) {
                   evaporationConst = Double.parseDouble(st.nextToken());
                   System.out.println(str + " = " + evaporationConst);
               }
               if (str.compareTo("iteration") == 0) {
                   iterationCount = Integer.parseInt(st.nextToken());
                   System.out.println(str + " = " + iterationCount);
               }
               if (str.compareTo("paralelRun") == 0) {
                   paralelRun = Boolean.parseBoolean(st.nextToken());
                   System.out.println(str + " = " + paralelRun);
               }
               if (str.compareTo("unavailableSelect") == 0) {
                   unavailableSelect = Boolean.parseBoolean(st.nextToken());
                   System.out.println(str + " = " + unavailableSelect);
               }
               if (str.compareTo("unavailableValue") == 0) {
                   unavailableValue = Integer.parseInt(st.nextToken());
                   System.out.println(str + " = " + unavailableValue);
               }
               if (str.compareTo("infeasibleSelect") == 0) {
                   infeasibleSelect = Boolean.parseBoolean(st.nextToken());
                   System.out.println(str + " = " + infeasibleSelect);
               }
               if (str.compareTo("binaryFile") == 0) {
                   binaryFile = st.nextToken();
                   System.out.println(str + " = " + binaryFile);
               }
               if (str.compareTo("combinedScoreFile") == 0) {
                   combinedScoreFile = st.nextToken();
                   System.out.println(str + " = " + combinedScoreFile);
               }
               if (str.compareTo("noeListFile") == 0) {
                   noeListFile = st.nextToken();
                   System.out.println(str + " = " + noeListFile);
               }
               if (str.compareTo("resultFile") == 0) {
                   resultFile = st.nextToken();
                   System.out.println(str + " = " + resultFile);
               }
               if (str.compareTo("localPhrmnUpdate") == 0) {
                   localPheromoneUpdate = Boolean.parseBoolean(st.nextToken());
                   System.out.println(str + " = " + localPheromoneUpdate);
               }
               if (str.compareTo("localUpdateVal") == 0) {
                   localPhrmnUpdtValue = Double.parseDouble(st.nextToken());
                   System.out.println(str + " = " + localPhrmnUpdtValue);
               }
               if (str.compareTo("cnstrntPeakCount") == 0) {
                   cnstrntPkCnt = Integer.parseInt(st.nextToken());
                   System.out.println(str + " = " + cnstrntPkCnt);
               }
               if (str.compareTo("nextPeakSelect") == 0) {
                   nextPeakSelect = Integer.parseInt(st.nextToken());
                   System.out.println(str + " = " + nextPeakSelect);
               }
               if (str.compareTo("constInfeasible") == 0) {
                   constForInfeasible = Integer.parseInt(st.nextToken());
                   System.out.println(str + " = " + constForInfeasible);
               }
               if (str.compareTo("phrmnUpdateCount") == 0) {
                   pheromoneUpdateCount = (int) Math.round(peakCount * Double.parseDouble(st.nextToken()) / 100);
                   System.out.println(str + " = " + pheromoneUpdateCount);
               }
               if (str.compareTo("elitistAntCount") == 0) {
                   elitistAntCount = peakCount * Integer.parseInt(st.nextToken()) / 100;
                   System.out.println(str + " = " + elitistAntCount);
               }               
               if (str.compareTo("pheromoneInitial") == 0) {
                   pheromoneInitialVal = Double.parseDouble(st.nextToken());
                   System.out.println(str + " = " + pheromoneInitialVal);
               }
               if (str.compareTo("localSearch") == 0) {
                   localSearchVersion = Integer.parseInt(st.nextToken());
                   System.out.println(str + " = " + localSearchVersion);
               }
               if (str.compareTo("localSearchStop") == 0) {
                   localSearchStop = Integer.parseInt(st.nextToken());
                   System.out.println(str + " = " + localSearchStop);
               }
               if (str.compareTo("keyBestSet") == 0) {
                   keyBestSet = Boolean.parseBoolean(st.nextToken());
                   assignmentDelta = Double.parseDouble(st.nextToken()) / 100;
                   System.out.println(str + " = " + assignmentDelta);
               }
               if (str.compareTo("pheromoneReset") == 0) {
                   pheromoneReset = Integer.parseInt(st.nextToken());
                   System.out.println(str + " = " + pheromoneReset);
               }
               if (str.compareTo("totalScoreCalc") == 0) {
                   totalScoreCalcCoef = Integer.parseInt(st.nextToken());
                   System.out.println(str + " = " + totalScoreCalcCoef);
               }
               if (str.compareTo("candidateList") == 0) {
                   candidateListValue = Integer.parseInt(st.nextToken());
                   System.out.println(str + " = " + candidateListValue);
               }
               if (str.compareTo("pheromoneMinMax") == 0) {
                   if (Boolean.parseBoolean(st.nextToken()) == true) {
                       pheromoneMin = Double.parseDouble(st.nextToken());
                       pheromoneMax = Double.parseDouble(st.nextToken());
                   } else {
                       pheromoneMax = Double.MAX_VALUE;
                       pheromoneMin = Double.NEGATIVE_INFINITY;
                       st.nextToken();
                       st.nextToken();
                   }
                   System.out.println(str + " = " + pheromoneMin + " - " + pheromoneMax);
               }
               if (str.compareTo("outputMatlab") == 0) {
                   outputMatlab = Boolean.parseBoolean(st.nextToken());
                   matlabFile = st.nextToken();
                   System.out.println(str + " = " + outputMatlab + "; " + matlabFile);
               }
               if (str.equals("gamma")) {
                   gamma = Integer.parseInt(st.nextToken());
               }
           }
           br.close();
           fr.close();
        } catch (Exception exc) { exc.printStackTrace(System.out); }
    }

    private int readCombinedScore() {
        try {
            FileReader fr;
            BufferedReader br;
            fr = new FileReader(combinedScoreFile);
            br = new BufferedReader(fr);
            String s;
            StringTokenizer st;
            combinedScore = new double[peakCount][aaCount];
            int i = 0, j;
            while ((s = br.readLine()) != null) {
                st = new StringTokenizer(s);
                for (j=0; j<aaCount; j++)
                    combinedScore[i][j] = Double.parseDouble(st.nextToken());
                i++;
           }

            System.out.println(combinedScoreFile + " is succesfully migrated.");
            br.close();
            fr.close();
            return 1;

        } catch (Exception exc) { exc.printStackTrace(System.out); return -1; }
    }
    
    private double getInfeasibleScore() {
        try {
            double res = 0;
            int i, j;
            for (i = 0; i < peakCount; i++)
            for (j = 0; j < aaCount; j++)
                if ((combinedScore[i][j] != 1000000000) && (combinedScore[i][j] > res))
                    res = combinedScore[i][j];
            return res * constForInfeasible;
        } catch (Exception e) { 
            e.printStackTrace(System.out);
            return -1.0;
        }
    }
    
    private void setHeuristicValues() {
        try {
            int i, j;
            heuristic = new double[peakCount][aaCount];
            heuristicBetaPower = new double[peakCount][aaCount];
            for (i = 0; i < peakCount; i++)
            for (j = 0; j < aaCount; j++) {
                heuristic[i][j] = 1 / combinedScore[i][j];
                heuristicBetaPower[i][j] = Math.pow(1 / combinedScore[i][j], beta);
            }
        } catch (Exception e) { e.printStackTrace(System.out); }
    }

    public void quicksort(double[] array, int[] index, int left, int right) {
        int pivot, leftIdx = left, rightIdx = right;
        double k;
        if (right - left > 0) {
            pivot = (left + right) / 2;
            while (leftIdx <= pivot && rightIdx >= pivot) {
                while (array[leftIdx] < array[pivot] && leftIdx <= pivot) leftIdx++;
                while (array[rightIdx] > array[pivot] && rightIdx >= pivot) rightIdx--;
                k = array[leftIdx]; array[leftIdx] = array[rightIdx];   array[rightIdx] = k;
                k = index[leftIdx]; index[leftIdx] = index[rightIdx];   index[rightIdx] = (int)k;
                leftIdx++;  rightIdx--;
                if (leftIdx - 1 == pivot) pivot = ++rightIdx;
                else if (rightIdx + 1 == pivot) pivot = --leftIdx;
            }
            quicksort(array, index, left, pivot - 1);
            quicksort(array, index, pivot + 1, right);
        }
    }

    private void runAnts() {
        FileWriter fr;
        StringBuffer s;
        double minTotalScore = Double.MAX_VALUE, worstGoodTotalScore = -1, bestCheck = -1, 
               stDevDelAvg, stDevDelAvgLast = 1, lastTotalScore = -1;
        int i, j, k, l, count, bestSoFar = -1, iterationBest[] = elitistAntCount > 0 ? new int[elitistAntCount] : null, 
                bestCheckCount = 0, bestAssignmentIteration = -1, improvementCount = 0, improvementAvg = 0, improvementDelta = 0,
                stopOnStDevAvg = 0,                      // parameter to stop when StDev/AVG <= threshold for predefined number of iteration
                stopOnNonImprovement = 20, // iterationCount/100*5;  // parameter to stop when there is not improvement in best total score for predefined number of iteration
                resetCount = 0,                 // parameter to stop whole program after resetting. counting all nonimprovement resets
                wholeResetCount = 0,            // all resets in program 
                bestAssignmentResetNo = -1,
                stopIteration;
        long    bestAssignmentTime = -1;        // time and number of the reset of the best assignment
                
        boolean bestSoFarUpdate = false;        // indicates whether best-so-far update is used
        ArrayList<AssignmentInfo> goodAssignmentList = new ArrayList<AssignmentInfo>();     // k-best solution set
        ArrayList<Double> goodSetValues = new ArrayList<Double>();
        AssignmentInfo bestAssignment = null;
        Ant ant[] = new Ant[antCount],
            bestAnt = new Ant(peakCount, combinedScore, binaryDistance, noeList, infeasibleScore, alpha, gamma, //beta,
                    pheromone, heuristicBetaPower, unavailableSelect, unavailableValue, infeasibleSelect,
                    localPhrmnUpdtValue, cnstrntPkCnt, nextPeakSelect, localSearchStop, totalScoreCalcCoef,
                    candidateListValue * aaCount / 100, candidateList, candidateListAA, scoreConstraints);
        bestAnt.stopWorking = false;
        bestAnt.curStep = peakCount;
            
        for (i = 0; i < antCount; i++) {
            ant[i] = new Ant(i, combinedScore, binaryDistance, noeList, infeasibleScore, alpha, gamma, //beta,
                    pheromone, heuristicBetaPower, unavailableSelect, unavailableValue, infeasibleSelect,
                    localPhrmnUpdtValue, cnstrntPkCnt, nextPeakSelect, localSearchStop, totalScoreCalcCoef,
                    candidateListValue * aaCount / 100, candidateList, candidateListAA, scoreConstraints);
        }
        
//        System.out.println(ant[0].whatScore());
//        if (true) return;
        
        try {
            fr = new FileWriter(resultFile);
            fr.write("MIN\tMAX\tAVG\tSTDev\tDiagMIN\tDiagMAX\tDiagAVG\tDiagSTDev\tIterBest\tBestSoFar\tStDev/AVG\tBest Ants\n");
            fr.close();
        } catch (Exception e) { e.printStackTrace(System.out); }

        long tTotal = System.nanoTime(), tIter = System.nanoTime(), tDecision = 0, 
             tNextPeak = 0, tGoodSet = 0, tUpdateBestAnalyse = 0, time, tLS = 0;
        System.out.println("Test: " + resultFile);
        
        try {
            for (i = 0; i < iterationCount; i++) {                      // main iteration loop
            System.out.println("\nITERATION " + i);
                improvementDelta++;

                for (j = 0; j < antCount; j++) ant[j].resetAnt((i * antCount + j) % peakCount);

                if (paralelRun)
                    for (k = 0; k < peakCount; k++) {
                        time = System.nanoTime();
                        for (l = 0; l < antCount; l++)
                            if (!ant[l].stopWorking) ant[l].makeDecision();
                        tDecision += System.nanoTime() - time;

                        time = System.nanoTime();
                        for (l = 0; l < antCount; l++)
                            if (!ant[l].stopWorking) {
                                if (localPheromoneUpdate) ant[l].localUpdatePheromone();
                                ant[l].selectNextPeak();
                            }
                        tNextPeak += System.nanoTime() - time;
                    }
 //               System.out.print("\tMakeDecision");
                else
                    for (l = 0; l < antCount; l++)
                        for (k = 0; k < peakCount; k++)
                            if (!ant[l].stopWorking) {
                                ant[l].makeDecision();
                                if (!ant[l].stopWorking) {
                                    if (localPheromoneUpdate) ant[l].localUpdatePheromone();
                                    ant[l].selectNextPeak();             // sequential'da local update nerede ?
                                }
                            }
                pheremoneEvaporate();                                     // evaoprate pheromone
            System.out.print("\tEvaporate");

                count = 0;
                int pastePlace = -1;                                                // for sorting iteration best assignments
                for (j = 0; j < iterationBest.length; j++) iterationBest[j] = -1;   // reset iteration best ant list
                for (l = 0; l < antCount; l++) {                                    // analysis part
//            ant[l].printAssignResults();
                    if (ant[l].getTotalScore() == -1) continue;                     // if didn't find total assignment then continue;

                    // look if found minimum total score in the iteration
                    time = System.nanoTime();                                       // for tUpdateBestAnalyse
                    if (elitistAntCount > 0) {                                 // part for iteration best k ants pheromone updates
                        for (j = 0; j < elitistAntCount; j++)
                            if ((iterationBest[j] == -1) || (ant[l].getTotalScore() < ant[iterationBest[j]].getTotalScore())) {
                                pastePlace = j;
                                break;
                            }
                        if (j != elitistAntCount) {
                            if (iterationBest[j] != -1)
                                for (j = elitistAntCount - 1; j > pastePlace; j--) 
                                    iterationBest[j] = iterationBest[j - 1];
                            iterationBest[pastePlace] = l;
                        }
                        count++;
                    }

                // add assignment to good assignment list if total score is less than minimum + assignment delta
                if (keyBestSet)
                if (ant[l].getTotalScore() <= minTotalScore * (1 + assignmentDelta)) {
                    for (j=0; j<goodSetValues.size(); j++)
                        if (ant[l].getTotalScore() == goodSetValues.get(j)) break;
                    if (j == goodSetValues.size()) {
                        goodAssignmentList.add(new AssignmentInfo(ant[l].getCompleteAssignment(), ant[l].getTotalScore(), ant[l].accuracy()));
                        goodSetValues.add(ant[l].getTotalScore());
                    }
                }

                    // pheromone update part 1
                    if (pheromoneUpdateCount == 0) ant[l].updatePheromone(1);               // every ant which fonud total assignment
                    if ((pheromoneUpdateCount == -2) && (ant[l].getInfeasibleCount() == 0)) // every ant which found feasible assignment
                        ant[l].updatePheromone(1);
//                ant[l].printAssignResults();
                    tUpdateBestAnalyse += System.nanoTime() - time;
                }
            System.out.print("\tAnalyze ");

                time = System.nanoTime();
                if (localSearchVersion != -1)              // ELITIST LOCAL SEARCH
                    for (j = 0; j < elitistAntCount; j++) {
                        if (iterationBest[j] != -1) {
                    System.out.print(iterationBest[j] + ": " + ant[iterationBest[j]].getTotalScore());
                            switch (localSearchVersion) {                                       // make elitist local search
                                case 2:
                                    ant[iterationBest[j]].localSearchTwoOpt();
                                    break;
                                case 3:
                                    ant[iterationBest[j]].localSearchThreeOpt();
                                    break;
                            }
                    System.out.println(" --> " + ant[iterationBest[j]].getTotalScore());               // print elitist local search results
                        } else break;
                        
                        // add assignment to good assignment list if total score is less than minimum + assignment delta
                        if (keyBestSet)
                        if (ant[iterationBest[j]].getTotalScore() <= minTotalScore * (1 + assignmentDelta)) {
                            for (k = 0; k < goodSetValues.size(); k++)
                                if (ant[iterationBest[j]].getTotalScore() == goodSetValues.get(k)) break;
                            if (k == goodSetValues.size()) {
                                goodAssignmentList.add(new AssignmentInfo(ant[iterationBest[j]].getCompleteAssignment(), ant[iterationBest[j]].getTotalScore(), ant[iterationBest[j]].accuracy()));
                                goodSetValues.add(ant[iterationBest[j]].getTotalScore());
                            }
                        }
                    }
                tLS += System.nanoTime() - time;
            System.out.print("\tLS");

                for (j = 0; j < elitistAntCount - 1; j++)
                    for (k = j + 1; k < elitistAntCount; k++)
                        if (ant[iterationBest[j]].getTotalScore() > ant[iterationBest[k]].getTotalScore()) {
                            l = iterationBest[j];
                            iterationBest[j] = iterationBest[k];
                            iterationBest[k] = l;
                        }

//                time = System.nanoTime();                         // time for tGoodSet
                if (minTotalScore-EPSILON > ant[iterationBest[0]].getTotalScore()) {                 // if minimum total score fonud.
                    minTotalScore = ant[iterationBest[0]].getTotalScore();
                    bestSoFar = iterationBest[0];
                    bestAssignment = new AssignmentInfo(ant[iterationBest[0]].getCompleteAssignment(), 
                                            ant[iterationBest[0]].getTotalScore(), ant[iterationBest[0]].accuracy());
                    bestAnt.peakAssignment = bestAssignment.assignment;
                    bestAnt.totalScore = bestAssignment.totalScore;
                    bestAnt.setUpdateScore();
                    bestAssignmentIteration = i;
                    bestAssignmentResetNo = wholeResetCount;
                    bestAssignmentTime = System.nanoTime() - tTotal;
                            
                    improvementCount++;
                    improvementAvg += improvementDelta;
                    improvementDelta = 0;
                    
                    // removing k-best assignments which are higher thn new bound value
                    if (keyBestSet)
                        cleanGoodAssignmentSet(goodAssignmentList, goodSetValues, minTotalScore * (1 + assignmentDelta));
                }
                time = System.nanoTime();

                double worstScore = 0;                         // finding worst score of elitist ants. (for relatively weighted update)
                for (j = pheromoneUpdateCount-1; j >= 0; j--)
                    if (iterationBest[j] != -1) {
                        worstScore = ant[iterationBest[j]].getTotalScore();
                        break;
                    }

                // pheromone update part
                if (pheromoneUpdateCount > 0)                            // iteration best k ants
                    for (j = 0; j < pheromoneUpdateCount; j++)
                        if (iterationBest[j] != -1)
//                            ant[iterationBest[j]].updatePheromone(pheromoneUpdateCount - j);              // using rank for pheromone update
                            ant[iterationBest[j]].updatePheromone( worstScore / ant[iterationBest[j]].getTotalScore() );        // using relative weight for pheromone update
                        else break;
//            System.out.println("\tPherUpdate");
//            if ((pheromoneUpdateCount == -1) && (bestSoFar != -1))              // burasi cok subheli. bunu deyismen lazim.
//                ant[bestSoFar].updatePheromone(1);       // best so far ant
//                if (bestSoFarUpdate)
//                    bestAnt.updatePheromone(pheromoneUpdateCount+1);                 // using rank for pheromone update
//                    bestAnt.updatePheromone( worstScore / bestAnt.getTotalScore() );        // using relative weight for pheromone update
 
                if (bestAssignment != null)               // this part is for resetting pheromone trails after some iteration
                    if (bestAssignment.totalScore != bestCheck) {
                        bestCheckCount = 0;
                        bestCheck = bestAssignment.totalScore;
                    } else
                        bestCheckCount++;
                else
                    bestCheckCount++;

                if (bestCheckCount == pheromoneReset) {
                    bestCheckCount = 0;
                    initializePheromoneTrail();
                }

//            System.out.println("Minimum total score after iteration " + i + " = " + minTotalScore);

//        for (j=0; j<peakCount; j++) System.out.println("Pheromone analysis: Average: " + getPheromoneAverage(pheromone[j]) +
//                "; StDev: " + getPheromoneStDev(pheromone[j]) + "; Median: " + getPheromoneMedian(pheromone[j]));
                if (iterationBest != null && iterationBest[0] != -1)
                    stDevDelAvg = analyzePheromone(ant[iterationBest[0]].getTotalScore(), bestAssignment.totalScore, iterationBest, ant);
                else if (bestAssignment != null)
                    stDevDelAvg = analyzePheromone(-1.0, bestAssignment.totalScore, null, null);
                else
                    stDevDelAvg = analyzePheromone(-1.0, -1.0, null, null);
//            System.out.print("StDev/AVG = " + stDevDelAvg + ".\t");
//            System.out.print(Math.abs(stDevDelAvg/stDevDelAvgLast-1) + " -> ");
                if (Math.abs(stDevDelAvg / stDevDelAvgLast - 1) <= 0.001) { // 5 iterasyon boyunca stdev /avg orani < 0.001 olursa programi durdur
                    stopOnStDevAvg++;                                       // bunlar daha parametrik hale getirilmedi.
                    if ((stopOnStDevAvg > 5) && (bestCheckCount > stopOnNonImprovement)) {
                        pheromoneInitialVal = 1 / ((1-evaporationConst) * bestAnt.getTotalScore());   // it depends on the method we use
                        initializePheromoneTrail();                      
                        wholeResetCount++;
                        for (j=0; j<peakCount; j++) if (!(bestAssignment.assignment[j].feasibility)) break;
                        bestSoFarUpdate = true;             // make best-so-far update
                        bestAnt.printAssignResults();       // print current best result.
                        if ((j==peakCount) && (lastTotalScore-EPSILON <= bestAssignment.totalScore)) {
                            resetCount++;
                        } else if (j==peakCount)
                            resetCount = 1;
                        else 
                            resetCount = 0;
                        if (resetCount == RESET_COUNT)
                            break;                  // if best feasible result is repeated stop the program
                        lastTotalScore = bestAssignment.totalScore;
                    }
                } else
                    stopOnStDevAvg = 0;
                stDevDelAvgLast = stDevDelAvg;

                tUpdateBestAnalyse += System.nanoTime() - time;

                // check pheromone mini max
                if (pheromoneMin != Double.NEGATIVE_INFINITY)
                    for (j = 0; j < peakCount; j++)
                    for (l = 0; l < aaCount; l++)
                        if      (pheromone[j][l] > pheromoneMax) pheromone[j][l] = pheromoneMax;
                        else if (pheromone[j][l] < pheromoneMin) pheromone[j][l] = pheromoneMin;

//            System.out.println("Iteration: " + (System.nanoTime() - tIter) / 1000000);
                tIter = System.nanoTime();
//            System.out.println();
            }
        } catch (Exception e) { e.printStackTrace(System.out); }
        
        stopIteration = i;
        improvementAvg /= improvementCount;
        tTotal = System.nanoTime() - tTotal;

//        for (i = 0; i < antCount; i++) ant[i].sayTime(resultFile);

        // printing results

        // printing best assignment
        try {
            if (outputMatlab) {
                fr = new FileWriter(matlabFile);
                for (i=0; i<peakCount; i++) {
                    count = bestAssignment.assignment[i].aminoAcid;                    
                    for (j=0; j<count; j++) fr.write("0 ");
                    fr.write("1 ");
                    for (j=count+1; j<aaCount; j++) fr.write("0 ");
                    fr.write(i + " " + count + " ");
                    if (i==count)
                        fr.write("1\n");
                    else 
                        fr.write("0\n");
                }
                fr.close();
            }
                        
            fr = new FileWriter(resultFile, true);
            
            // printing good assignment list
            fr.write("\n");
            if (keyBestSet)
            for (i = 0; i < goodAssignmentList.size(); i++) {
                count = 0;
                for (j = 0; j < peakCount; j++) {
//                    fr.write("<" + j + "-" + goodAssignmentList.get(i).assignment[j].aminoAcid + ", " + goodAssignmentList.get(i).assignment[j].feasibility + "> ");
                    fr.write("<" + goodAssignmentList.get(i).assignment[j].aminoAcid + ", " + goodAssignmentList.get(i).assignment[j].feasibility + "> ");
                    if (!goodAssignmentList.get(i).assignment[j].feasibility) count++;
                }
                fr.write("\tTotal score: " + goodAssignmentList.get(i).totalScore + "; Number of infeasible assignments: " + count + "; Accuracy: " + +goodAssignmentList.get(i).accuracy + "\n");
            }
                        
            count = 0;
            if (bestAssignment != null) {
                fr.write("\n\nPrinting best assignment:");
                for (j = 0; j < peakCount; j++) {
                    fr.write("<" + j + "-" + bestAssignment.assignment[j].aminoAcid + ", " + bestAssignment.assignment[j].feasibility + "> ");
                    if (!bestAssignment.assignment[j].feasibility) {
                        count++;
                    }
                }
                fr.write("\nNumber of infeasible assignments: " + count);
                fr.write("\nMinimum total score = " + bestAssignment.totalScore);
                fr.write("\nBest assignments iteration: " + bestAssignmentIteration);
                fr.write("\nAccuracy = " + bestAssignment.accuracy);
                fr.write("\nBest assignment time = " + bestAssignmentTime/1000000000);
                fr.write("\nBest assignment reset no = " + bestAssignmentResetNo);
                fr.write("\nPenalty score = " + infeasibleScore);
                fr.write("\nGood assignment set size: " + goodAssignmentList.size());
                fr.write("\nImprovement count: " + improvementCount);
                fr.write("\nAverage improvement rate: " + improvementAvg);
            } else {
                fr.write("No total assignments found.");
            }
            fr.write("\nGood assignment set size: " + goodAssignmentList.size());
            fr.write("\nMake decision: " + tDecision / 1000000);
            fr.write("\nNext peak select: " + tNextPeak / 1000000);
            fr.write("\nLocal search: " + tLS / 1000000);
            fr.write("\nGood assignment set operations: " + tGoodSet / 1000000);
            fr.write("\nPheromone update, analyse & best assignment operations : " + tUpdateBestAnalyse / 1000000);
            fr.write("\nProblem is done in " + tTotal / 1000000 + " milliseconds.");
            fr.write("\n\n" + bestAssignment.totalScore + "\t" + bestAssignment.accuracy + "\t" +
                    (bestAssignmentIteration+1) + "\t" + stopIteration + "\t" + tTotal / 1000000000 + "\t" + wholeResetCount);
            fr.close();
        } catch (Exception e) { e.printStackTrace(System.out); }
    }

    private int readBinaryDistance() {     // reading binary distance values from binary file to binaryDistance matrix
        FileReader fr;
        BufferedReader br;
        try {
            fr = new FileReader(binaryFile);
            br = new BufferedReader(fr);
            String s;
            StringTokenizer st;
            binaryDistance = new int[aaCount][aaCount];
            int i = 0, j;

            while ((s = br.readLine()) != null) {
                st = new StringTokenizer(s);
                for (j=0; j<aaCount; j++)
                    binaryDistance[i][j] = Integer.parseInt(st.nextToken());
                i++;
           }

            System.out.println(binaryFile + " is succesfully migrated.");
            br.close();
            fr.close();
            return 1;

        } catch (Exception exc) { exc.printStackTrace(System.out); return -1; }
    }

    private int readNoeList() {     // reading noe values from file to noe arrays
        FileReader fr;                              // "-1": peak numbering is from 1 in file
        BufferedReader br;                          // so, peak=peak-1, noePeak=noePeak-1
        try {
            fr = new FileReader(noeListFile);
            br = new BufferedReader(fr);
            String s;
            StringTokenizer st;
            int peak, row = 0, i, count, noePeak;
            noeList = new ArrayList[peakCount];
            for (i=0; i<peakCount; i++) noeList[i] = new ArrayList<Integer>();

            while ((s = br.readLine()) != null) {
                st = new StringTokenizer(s);
                count = Integer.parseInt(st.nextToken());                   // number of noe peaks
                peak = Integer.parseInt(st.nextToken()) - 1;                // number of peak
                for (i = 0; i < count; i++) {
                    noePeak = Integer.parseInt(st.nextToken())-1;
                    if (noeList[peak].indexOf(noePeak) == -1) {
                        noeList[peak].add(noePeak);
                        noeList[noePeak].add(peak);
                    }
                }
            }

            System.out.println(noeListFile + " is succesfully migrated.");
            br.close();
            fr.close();
            return 1;

        } catch (Exception exc) { exc.printStackTrace(System.out); return -1; }
    }

    private void pheremoneEvaporate() {
        int i,j;
        for (i=0; i<peakCount; i++)
        for (j=0; j<aaCount; j++)
            pheromone[i][j] *= evaporationConst;
    }

    private void initializePheromoneTrail() {
        for (int i = 0; i < peakCount; i++)
            for (int j = 0; j < aaCount; j++)
                pheromone[i][j] = pheromoneInitialVal;
        System.out.println("Pheromone trails initialized to value " + pheromoneInitialVal);
    }

    private int getAvailableCount() {
        int count = 0, i, j;
        for (i = 0; i < peakCount; i++) {
            for (j = 0; j < aaCount; j++) {
                if (combinedScore[i][j] != unavailableValue) {
                    count++;
                }
            }
        }
        return count;
    }

    private int getBinaryCloseCount() {
        int i, j, res = 0;
        for (i = 0; i < aaCount; i++) {
            for (j = 0; j < aaCount; j++) {
                if (binaryDistance[i][j] == 1) {
                    res++;
                }
            }
        }
        return res;
    }

    private int getNoeCount() {
        int count = 0, i;
        for (i = 0; i < peakCount; i++) {
            count += noeList[i].size();
        }
        return count;
    }

    private int getMinAvailableCount() {
        int min = peakCount + 1, count, i, j;
        for (i = 0; i < peakCount; i++) {
            count = 0;
            for (j = 0; j < aaCount; j++) {
                if (combinedScore[i][j] != unavailableValue) {
                    count++;
                }
            }
            if (min > count) {
                min = count;
            }
        }
        return min;
    }

    private int getMinNoeCount() {
        int min = peakCount + 1, i;
        for (i = 0; i < peakCount; i++) {
            if (min > noeList[i].size()) {
                min = noeList[i].size();
            }
        }
        return min;
    }

    private double getBinaryMinCloseCount() {
        int min = aaCount + 1, count, i, j;
        for (i = 0; i < aaCount; i++) {
            count = 0;
            for (j = 0; j < aaCount; j++) {
                if (binaryDistance[i][j] == 1) {
                    count++;
                }
            }
            if (min > count) {
                min = count;
            }
        }
        return min;
    }

    private int getMaxAvailableCount() {
        int max = -1, count, i, j;
        for (i = 0; i < peakCount; i++) {
            count = 0;
            for (j = 0; j < aaCount; j++) {
                if (combinedScore[i][j] != unavailableValue) {
                    count++;
                }
            }
            if (max < count) {
                max = count;
            }
        }
        return max;
    }

    private int getBinaryMaxCloseCount() {
        int max = -1, count, i, j;
        for (i = 0; i < aaCount; i++) {
            count = 0;
            for (j = 0; j < aaCount; j++) {
                if (binaryDistance[i][j] == 1) {
                    count++;
                }
            }
            if (max < count) {
                max = count;
            }
        }
        return max;
    }

    private int getMaxNoeCount() {
        int max = -1, i;
        for (i = 0; i < peakCount; i++) {
            if (max < noeList[i].size()) {
                max = noeList[i].size();
            }
        }
        return max;
    }

    private double getAverageAvailableCount() {
        return (double) getAvailableCount() / peakCount;
    }

    private double getBinaryAverageCloseCount() {
        return (double) getBinaryCloseCount() / aaCount;
    }

    private double getAverageNoeCount() {
        return (double) getNoeCount() / peakCount;
    }

    private int getMedianAvailableCount() {
        int count, i, j, m, n, list[] = new int[peakCount];
        for (i = 0; i < peakCount; i++) {
            list[i] = -1;
        }
        for (i = 0; i < peakCount; i++) {
            count = 0;
            for (j = 0; j < aaCount; j++) {
                if (combinedScore[i][j] != unavailableValue) {
                    count++;
                }
            }
            for (m = 0; m < peakCount; m++) {
                if (list[m] == -1) {
                    list[m] = count;
                    break;
                }
                if (list[m] > count) {
                    for (n = peakCount - 1; n > m; n--) {
                        list[n] = list[n - 1];
                    }
                    list[m] = count;
                    break;
                }
            }
        }
//        System.out.print("Median list: ");
        for (i = 0; i < peakCount; i++) {
//            System.out.print(list[i] + " ");
        }
//        System.out.println();
        return list[peakCount / 2 - 1];
    }

    private int getBinaryMedianCloseCount() {
        int count, i, j, m, n, list[] = new int[aaCount];
        for (i = 0; i < aaCount; i++) {
            list[i] = -1;
        }
        for (i = 0; i < aaCount; i++) {
            count = 0;
            for (j = 0; j < aaCount; j++) {
                if (binaryDistance[i][j] == 1) {
                    count++;
                }
            }
            for (m = 0; m < aaCount; m++) {
                if (list[m] == -1) {
                    list[m] = count;
                    break;
                }
                if (list[m] > count) {
                    for (n = aaCount - 1; n > m; n--) {
                        list[n] = list[n - 1];
                    }
                    list[m] = count;
                    break;
                }
            }
        }
//        System.out.print("Median list: ");
        for (i = 0; i < aaCount; i++) {
//            System.out.print(list[i] + " ");
        }
//        System.out.println();
        return list[aaCount / 2 - 1];
    }

    private int getMedianNoeCount() {
        int count, i, m, n, list[] = new int[peakCount];
        for (i = 0; i < peakCount; i++) {
            list[i] = -1;
        }
        for (i = 0; i < peakCount; i++) {
            count = noeList[i].size();
            for (m = 0; m < peakCount; m++) {
                if (list[m] == -1) {
                    list[m] = count;
                    break;
                }
                if (list[m] > count) {
                    for (n = peakCount - 1; n > m; n--) {
                        list[n] = list[n - 1];
                    }
                    list[m] = count;
                    break;
                }
            }
        }
//        System.out.print("Median list: ");
        for (i = 0; i < peakCount; i++) {
//            System.out.print(list[i] + " ");
        }
//        System.out.println();
        return list[peakCount / 2 - 1];
    }

    private double getPheromoneAverage(double pher[]) {
        double res = 0;
        for (int i = 0; i < aaCount; i++) {
            res += pher[i];
        }
        return res / aaCount;
    }

    private double getPheromoneMedian(double pher[]) {
        int place = -1, i, j;
        double list[] = new double[aaCount];
        for (i = 0; i < aaCount; i++) list[i] = -1;

        for (i = 0; i < aaCount; i++) {
            for (j = 0; j < aaCount; j++)
                if (list[j] == -1 || pher[i] <= list[j]) {
                    place = j;
                    break;
                }
            if (list[j] != -1)
                for (j = aaCount - 1; j > place; j--)
                    list[j] = list[j - 1];
            list[place] = pher[i];
        }

        return list[aaCount / 2 - 1];
    }

    private double getPheromoneStDev(double pher[]) {
        double avg = getPheromoneAverage(pher), res = 0;

        for (int i = 0; i < aaCount; i++)
            res += Math.pow(pher[i] - avg, 2.0);

        return Math.sqrt(res / aaCount);
    }

    private void combinedScoreAnalysis() {
        System.out.println();
        System.out.println("Available assignment count: " + getAvailableCount());
        System.out.println("Min avilable assignment count of peaks: " + getMinAvailableCount());
        System.out.println("Max available assignment conut of peaks: " + getMaxAvailableCount());
        System.out.println("Average available assignment conut of peaks: " + getAverageAvailableCount());
        System.out.println("Median available assignment count of peaks: " + getMedianAvailableCount());
        System.out.println();
        System.out.println("Mean of combined score matrix: " + getCombinedScoreAverage());
        System.out.println("Standart deviation of combined score matrix: " + getCombinedScoreStdev());
        System.out.println("'StDev/Avg' of combined score matrix: " + getCombinedScoreStdev() / getCombinedScoreAverage());
    }

    private void binaryAnalysis() {
        System.out.println();
        System.out.println("Close amino acid count: " + getBinaryCloseCount());
        System.out.println("Min close amino acid count of amino acids: " + getBinaryMinCloseCount());
        System.out.println("Max close amino acid conut of amino acids: " + getBinaryMaxCloseCount());
        System.out.println("Average close amino acid count of amino acids: " + getBinaryAverageCloseCount());
        System.out.println("Median close amino acid count of amiono acids: " + getBinaryMedianCloseCount());
    }

    private void noeListAnalysis() {
        System.out.println();
        System.out.println("NOE constraint count: " + getNoeCount());
        System.out.println("Min NOE constraint count of peaks: " + getMinNoeCount());
        System.out.println("Max NOE constraint conut of peaks: " + getMaxNoeCount());
        System.out.println("Average NOE constraint count of peaks: " + getAverageNoeCount());
        System.out.println("Median NOE constraint count of peaks: " + getMedianNoeCount());
    }

    private double analyzePheromone(double iterBest, double best, int[] iterBestAnts, Ant[] ants) {
        int i, j;
        double min = pheromone[0][0], max = min, avg = 0, stDev = 0,
               dMin = min, dMax = max, dAvg = 0, dStDev = 0, stDevDelAvg;
        for (i=0; i<peakCount; i++) {
            if (pheromone[i][i] > dMax) dMax = pheromone[i][i];
            if (pheromone[i][i] < dMin) dMin = pheromone[i][i];
            dAvg += pheromone[i][i];
            for (j=0; j<aaCount; j++) {
                if (pheromone[i][j] > max) max = pheromone[i][j];
                if (pheromone[i][j] < min) min = pheromone[i][j];
                avg += pheromone[i][j];
            }
        }
        dAvg /= peakCount;
        avg /= peakCount * aaCount;
        for (i=0; i<peakCount; i++) {
            dStDev += Math.pow(dAvg - pheromone[i][i], 2.0);
            for (j=0; j<aaCount; j++)
                stDev += Math.pow(avg - pheromone[i][j], 2.0);
        }
        dStDev = Math.sqrt(dStDev);
        stDev = Math.sqrt(stDev);
        stDevDelAvg = stDev / avg;
        try {
            DecimalFormat df = new DecimalFormat("0.000000"),
                          af = new DecimalFormat("0.00");
            FileWriter fr = new FileWriter(resultFile, true);
            String s = df.format(min)+"\t"+df.format(max)+"\t"+df.format(avg)+"\t"+df.format(stDev)+
                       "\t"+df.format(dMin)+"\t"+df.format(dMax)+"\t"+df.format(dAvg)+"\t"+df.format(dStDev)+
                       "\t"+df.format(iterBest)+"\t"+df.format(best)+"\t" + df.format(stDevDelAvg) + "\t";
            StringBuilder st = new StringBuilder(s);
            if (iterBestAnts != null)
                for (i=0; i<iterBestAnts.length; i++) {
                    s = df.format(ants[iterBestAnts[i]].totalScoreBeforeLS) + "->" + df.format(ants[iterBestAnts[i]].getTotalScore()) + ", " +
                        af.format(ants[iterBestAnts[i]].accuracyBeforeLS) + "->" + af.format(ants[iterBestAnts[i]].accuracy()) + "; ";
                    st.append(s);
                }
            fr.write(st.toString() + "\n");
            fr.close();
            return stDevDelAvg;
        } catch (Exception e) { e.printStackTrace(System.out); return -1; }
    }

    private double getCombinedScoreAverage() {  // average of the combined score matrix without 10e9's
        int i, j, count = 0;
        double avg = 0;
        for (i=0; i<peakCount; i++)
        for (j=0; j<aaCount; j++)
            if (combinedScore[i][j] != unavailableValue) {
                avg += combinedScore[i][j];
                count++;
            }
        return avg/count;
    }
    
    private double getCombinedScoreStdev() {
        double avg = getCombinedScoreAverage(), res = 0;
        int i, j, count = 0;
        
        for (i = 0; i < peakCount; i++)
        for (j = 0; j < aaCount; j++)
            if (combinedScore[i][j] != unavailableValue) {
                res += Math.pow(combinedScore[i][j] - avg, 2.0);
                count++;
            }
        return Math.sqrt(res / count);
    }

    private double calculateDiagonalScore() {
        double res = 0;
        for (int i=0; i<peakCount; i++)
            res += combinedScore[i][i];
        return res;
    }

    private void cleanGoodAssignmentSet( ArrayList<AssignmentInfo> set, ArrayList<Double> valSet, double limit ) {
        int i = 0;
        while (i < set.size()) {
            if (set.get(i).totalScore > limit)
                set.remove(i);
            else
                i++;
        }
        i = 0;
        while (i < valSet.size()) {
            if (valSet.get(i) > limit)
                valSet.remove(i);
            else
                i++;
        }
    }
}