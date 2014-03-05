package aco_sba;

import java.io.*;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Random;

/*
 * Bu class bir karincayi temsil eder. Bir karincanin yaptigi tum islemleri yapar. * 
 */

public class Ant {

    private int startPeak, curPeak, assignedAACount,
                number,              // ant's number. used without reason. 
                candSize,            // will display the candidate list search point (if 25% and peak=80 will be =20 then +1, +1, ...)
                availableAA[],       // availability list of AA's ; "1=available"; "0=infeasible"; "-1=unavailable";
                peakConstraints[],   // peak constraint count list (for most constraint peak selection)
                peakMemoryReverse[]; // reverse array of peakMemory
    public int curStep, 
               peakMemory[];// memory of the ant. peaks which are assigned sequentially

    public Assignment peakAssignment[];           // assignment of peak to AA with Feasibility; "i = peak"
    public boolean stopWorking = false;           // karincanin search'e devam edip etmedigini gosterir
    private final boolean infeasibleSelect,       // whether NOE infeasible AA's will be available
                          unavailableSelect;      // select unavailable 10e+9 values or not

    private final int peakCount,                // number of peaks
                      aaCount,                  // number of AA's
                      binaryDistance[][],       // binary distance matrix
                      alpha, gamma, //beta,     // standard parameters of ACO
                      constraintPeakCount,  // percent value for selection from the most constrained peaks
                      nextPeakSelect,       // number of the method for selection next peak
                      scoreConstraints[],   // number of 10e9 values per peak
                      totalScoreCalcCoef,   // a cooficient for total score calculation
                      localSearchStop,      // specifies when to stop local search. see "parameters.txt",
                      candSizeStart,        // initial number of candSize (percent * aaCount / 100)
                      unavailableValue;      // which value will be unavailable
                                             // (may be deleted if it always 10e+9)
    public double totalScore,              // totalScore - current sum of combined scores,
                  totalScoreBeforeLS,       // local search'e girmeden olan skor
                  accuracyBeforeLS,         // local search'e girmeden olan accuracy
                  updateScore,             // updateScore - pheromone update score (setUpdateScore)
                  pheromone[][],            // pheromne matrix
                  heuristicBetaPower[][];   // heuristic fonksiyon ^ beta. bu deger sabit oldugundan her seferinde 
                                            // hesaplanmamasi icin onceden hesaplanarak memory'de tutulur
    private final double infeasibleScore,   // infeasibleScore - score for 10e+9 assignments. penalty skoru
                         combinedScore[][],     // combined skor matris
                        localPheromoneUpdateValue;  // 
    private final ArrayList<Integer> noeList[],  // noe listesi. her peak icin bir ArrayList<Integer> turunde
                                    candidateList[],    // list of aa's to which a peak can be assigned
                                     candidateListAA[];             // list of peaks to which an aa can be assigned

    private long tProb = 0, tConstr = 0, tFind = 0, tNext = 0, tSort = 0;   // performans (zaman) hesabi icin degiskenler

    // ciktiya shu anki performans (sure, zaman) durumunu veren fonskiyon
    public void sayTime(String fileName) {
        try {
            FileWriter fr = new FileWriter(fileName, true);
            fr.write("Ant " + number + ": Constraints: " + tConstr + " Probabilities: "
                     + tProb + " Find: " + tFind + " NextPeakSelect: " + tNext + " Sort: " + tSort + "\n");
            fr.close();
        } catch (Exception e) {};
    }

    // Ant class'in konstruktoru.
    public Ant(int num, double[][] cmbnd, int[][] bnr, ArrayList<Integer>[] noe, double infsblScr, int a, int g, //int b, was for beta
               double[][] phrmn, double[][] hrstc, boolean unavailSlct, int unavailVal, boolean infsblSel,
               double lclPhrmnUpdtVal, int cnstrntPkCnt, int nxtPkSlct, int searchStop, int ttlScrCoef,
               int candListSize, ArrayList<Integer>[] candidate, ArrayList<Integer>[] candidateAA,
               int[] scrConstr) {
        number = num;
        binaryDistance = bnr;
        combinedScore = cmbnd;
        pheromone = phrmn;
        heuristicBetaPower = hrstc;
        noeList = noe;
        peakCount = noeList.length;
        aaCount = combinedScore[0].length;
        infeasibleScore = infsblScr;
        alpha = a;
        gamma = g;
        candSizeStart = candListSize;
        infeasibleSelect = infsblSel;
        unavailableSelect = unavailSlct;
        unavailableValue = unavailVal;
        localPheromoneUpdateValue = lclPhrmnUpdtVal;
        constraintPeakCount = peakCount * cnstrntPkCnt / 100;
        nextPeakSelect = nxtPkSlct;
        localSearchStop = searchStop;
        totalScoreCalcCoef = ttlScrCoef;
        candidateList = candidate;
        candidateListAA = candidateAA;
        scoreConstraints = scrConstr;
    }

    /* karincayi baslangic pozisyona dondurur. tum parametreleri resetler. 
    strt parametresi karincanin basta uzerinde olacagi peak'dir
    */
    public void resetAnt(int strt) {
        curPeak = startPeak = strt;
        stopWorking = false;
        totalScore = totalScoreBeforeLS = accuracyBeforeLS = 0;
        updateScore = 0;
        curStep = 0;
        candSize = candSizeStart;
        assignedAACount = 0;
        peakConstraints = new int[peakCount];
        peakAssignment = new Assignment[peakCount];
        peakMemory = new int[peakCount];
        peakMemoryReverse = new int[peakCount];
        availableAA = new int[aaCount];
        for (int i=0; i<aaCount; i++) availableAA[i] = 0;            // ESKIDEN 1 OLUYORDU. SIMDI 0
        for (int i=0; i<peakCount; i++) {
//            peakConstraints[i] = -1;
            peakConstraints[i] = unavailableSelect ? getScoreConstraints(i, new int[candidateList[i].size()], candidateList[i]) // new int array is passed for no reason
                                                   : getScoreConstraints(i, new int[candidateList[i].size()], candidateList[i]) * 2 ;
            peakAssignment[i] = new Assignment(-1, false);
            peakMemory[i] = -1;
            peakMemoryReverse[i] = -1;
        }
        peakMemory[0] = curPeak;
        peakMemoryReverse[curPeak] = 0;
    }

    public void updatePheromone( double weight ) {
        if (stopWorking) return;
        for (int i=0; i<peakCount; i++)
             pheromone[i][peakAssignment[i].aminoAcid] += (updateScore * weight);
    }

    public void localUpdatePheromone() {
        pheromone[curPeak][peakAssignment[curPeak].aminoAcid] *= localPheromoneUpdateValue;
    }

    public int getCurrentPeak() { return curPeak; }

    public int getCurrentAssignment() { return peakAssignment[curPeak].aminoAcid; }

    // Will find an AA and assign it to current peak
    // When an ant stopped working it will wait until other ants finish current
    //      iteration and will start at the next iteration
    public void makeDecision() {
        long time1 = System.nanoTime();
        Assignment assignment = findAminoAcid();
        tFind += System.nanoTime() - time1;
        if (assignment != null) {                                        // if something is found to assign
            peakAssignment[curPeak].aminoAcid = assignment.aminoAcid;
            peakAssignment[curPeak].feasibility = assignment.feasibility;
            // next contraint peak operations
            if (nextPeakSelect == 2) {
                int i, j;
                ArrayList<Integer> cand = candidateListAA[assignment.aminoAcid];
                peakConstraints[curPeak] = -1;
                for (i = 0; i < cand.size(); i++)
                    if (peakConstraints[cand.get(i)] != -1) {
                        peakConstraints[cand.get(i)] += 2;
                        if (combinedScore[cand.get(i)][assignment.aminoAcid] == unavailableValue)
                            peakConstraints[cand.get(i)] -= unavailableSelect ? 1 : 2;
                    }
                ArrayList<Integer> noe = noeList[curPeak];
                int noeCount = noe.size(), noePeak;
                for (i = 0; i < noeCount; i++) {
                    noePeak = noe.get(i);
                    if (peakConstraints[noePeak] == -1) continue;
                    cand = candidateList[noePeak];
                    for (j = 0; j < cand.size(); j++) {
                        // bu kisim tartisilir kisimdir. 1, 2 olayindan dolayi. yani bir peak secildiyi zaman 1 mi 2 mi yoxsa her ikisi mi?
                        if (binaryDistance[cand.get(j)][assignment.aminoAcid] == 0)
                            peakConstraints[noePeak] += infeasibleSelect ? 1 : 2;
//                        if (combinedScore[noePeak][cand.get(j)] == unavailableValue )
//                            peakConstraints[noePeak] -= unavailableSelect ? 1 : 2;
                    }
                }
            }

            if (candSize < aaCount) candSize++;             // adds one to cand to extend the candidate list
            assignedAACount++;
            availableAA[assignment.aminoAcid] = -1;
            peakMemory[curStep] = curPeak;
            peakMemoryReverse[curPeak] = curStep;
        } else                               // if nothing is found
            stopWorking = true;
    }

    public void selectNextPeak() {
        curStep++;
        if ( !completed() )                 // if it hasn't finished current circle
            switch (nextPeakSelect) {
                case 0: { curPeak = nextSequentialPeak(); break; }
                case 1: { curPeak = nextRandomPeak();     break; }
                case 2: { curPeak = nextConstraintPeak(); break; }
            }
        else {                               // if it finished it's current circle
            totalScore = calculateDeltaScore(peakAssignment);
            setUpdateScore();
        }
    }

    private Assignment findAminoAcid() {

        ArrayList<Integer> candidate = candidateList[curPeak];
        int i, highIndex = -1, length = candSize,                 // bunu (length'i) 10e+9 sonra hepsini kapsayacak sekilde yap.
            avail[] = new int[length];
        double values[] = new double[length],
               sum = 0, sum1 = 0, randValue;

        long time1 = System.nanoTime();
        for (i=0; i<length; i++) avail[i] = availableAA[candidate.get(i)];
        getScoreConstraints(curPeak, avail, candidate);
        if (assignedAACount > 0) getNoeConstraints(curPeak, avail, candidate);
        tConstr += System.nanoTime() - time1;

        // values may be determined in the ant class as a attribute
        time1 = System.nanoTime();
        for (i=0; i<length; i++) {                                // initialize values
            if ( avail[i] == -1 )
                values[i] = 0;
            else if ( avail[i] == 0 )
                values[i] = Math.pow(pheromone[curPeak][candidate.get(i)], alpha) * heuristicBetaPower[curPeak][candidate.get(i)];
            else
                values[i] = Math.pow(pheromone[curPeak][candidate.get(i)], alpha) * heuristicBetaPower[curPeak][candidate.get(i)] / Math.pow( avail[i] * infeasibleScore, gamma );

            if (Double.isInfinite(values[i])) values[i] = Double.MAX_VALUE;
            sum += values[i];
        }
        if (sum == 0) {                             // nothing to assign. all avail[i]'s are = -1 that is sum = 0
            tProb += System.nanoTime() - time1;     // can't be if feasibleSelect and unavailableSelect are TRUE
            return null;
        }

        if (number==1)
            number = 1;

        randValue = new Random().nextDouble();                      // selecting among values
        for (i=0; i<length; i++) {                                  // normalize values
            sum1 += values[i] /= sum;
            if ( sum1 > randValue ) { highIndex = i; break; }
        }
        tProb += System.nanoTime() - time1;

        if ( avail[highIndex] == 0 )          // seperate infeasible (false) and feasible (true) assignments with boolean value
            return new Assignment(candidate.get(highIndex), true);
        else
            return new Assignment(candidate.get(highIndex), false);
    }

    // for sequential peak assignment e.g. IF startPeak=5 THEN nextPeak=6 THEN nextPeak=7
    private int nextSequentialPeak() { return (curPeak + 1) % peakCount; }

    // random selection of next peak.
    // this code can be optimised !!!!!!!
    private int nextRandomPeak() {
        ArrayList<Integer> availablePeaks = new ArrayList<Integer>();
        for (int i=0; i<peakCount; i++)
            if (peakAssignment[i].aminoAcid == -1) availablePeaks.add(i);
        return availablePeaks.get( new Random().nextInt(availablePeaks.size()) );
    }

    // selection of the next peak is biased toward contrainted peaks (which has less AA's to assign )
    private int nextConstraintPeak() {
        int mostConstraintPeaks[] = new int[peakCount],
            constraints[] = new int[peakCount], i;

        for (i = 0; i < peakCount; i++) {
            mostConstraintPeaks[i] = i;
            constraints[i] = peakConstraints[i];
        }
        long time = System.nanoTime();
        quicksort(constraints, mostConstraintPeaks, 0, mostConstraintPeaks.length-1);
        tSort += System.nanoTime() - time;

        int randBound = constraintPeakCount;
        for (i=0; i<constraintPeakCount; i++)
            if (constraints[i]==-1) {
                randBound = i;
                break;
            }
        return mostConstraintPeaks[new Random().nextInt(randBound)];
    }

    public void setUpdateScore() { updateScore = 1 / Math.pow(totalScore, 1.0/totalScoreCalcCoef); }

    // in this procedure contraint is 10e9 value, so this procedure try to prevent 10e9 from being selected
    private int getScoreConstraints(int peak, int[] avail, ArrayList<Integer> cand) {
        int res = 0;
        for (int i = 0; i < avail.length; i++)
            if ( (combinedScore[peak][cand.get(i)] == unavailableValue) && (avail[i]!=-1) ) {
                if (unavailableSelect)
                    avail[i]++;
                else
                    avail[i] = -1;
                res++;
            }
        return res;
    }

    private void getNoeConstraints(int peak, int[] avail, ArrayList<Integer> cand) {
        int i, j, noeCount, noePeak, noeAA;
        ArrayList<Integer> noe = noeList[peak];
        noeCount = noe.size();
        for (i = 0; i < noeCount; i++) {
            noePeak = noe.get(i);
            noeAA = peakAssignment[noePeak].aminoAcid;
            if (noeAA != -1)
                for (j = 0; j < avail.length; j++)
                    if ( (binaryDistance[noeAA][cand.get(j)]==0) && (avail[j]!=-1) )
                        if (infeasibleSelect)
                            avail[j]++;
                        else
                            avail[j] = -1;
        }
    }

    public void printAssignResults() {
//        if (stopWorking) {
//            System.out.println("Ant "+number+" didn't find any complete assignments.");
//            return;
//        }
//        System.out.println("Printing assign results for ant " + number);
        for (int i=0; i<peakAssignment.length; i++)
            System.out.print("<"+i+"-"+peakAssignment[i].aminoAcid+", "+peakAssignment[i].feasibility+"> ");
//        int i;
//        for (i=0; i<peakMemory.length; i++)
//            if (peakMemory[i] != -1)
//                System.out.print("<"+peakMemory[i]+"-"+peakAssignment[peakMemory[i]].aminoAcid+", "+peakAssignment[peakMemory[i]].feasibility+"> ");
//            else break;
//        System.out.println("\n" + curStep);
        System.out.println("\nTotal score of the assignment = " + new DecimalFormat("0.0000").format(totalScore) );
        System.out.println("Accuracy = " + accuracy());
        System.out.println("No. of infeasible assignments: " + getInfeasibleCount());
        System.out.println();
    }

    public double getTotalScore() { 
        if ( !stopWorking ) return (double) Math.round(totalScore * 1000000) / 1000000;
        else return -1;
    }

    private boolean completed() { return curStep == peakCount; }

    private double setTotalScore( Assignment[] completeAssignment ) {
        double res = 0;
        for (int i = 0; i < peakCount; i++) {
            if ( ! completeAssignment[i].feasibility )                  // bu duruma gore degisiyor. BUNA 31/03/2011 TARIHINDE ANLASTIK
                res += infeasibleScore;
            else
                res += combinedScore[i][completeAssignment[i].aminoAcid];
//            res += combinedScore[i][completeAssignment[i].aminoAcid];  // eskiden yaptigimdir
//            if ( ! completeAssignment[i].feasibility )
//                res += infeasibleScore;
        }
        return res;
    }

    public void localSearchTwoOpt() {
        // consider LS stop conditions. becuase you changed this code. and now it works only with "stop when no improvement"

        totalScoreBeforeLS = totalScore;
        accuracyBeforeLS = accuracy();
        if ( !completed() || stopWorking ) return;

        ArrayList<Integer> swapList = new ArrayList<Integer>(), noe;
        int i, j, k, aa1, aa2, maxI, maxJ, pCount;
//        Assignment assignment[] = new Assignment[peakCount];
        double maxDelta, delta[][] = new double[peakCount][peakCount], tempScore;
        for (i=0; i<peakCount; i++) {                            // create new empty assignment
//            assignment[i] = new Assignment(-1, false);
            swapList.add(i);
        }

        System.out.println("\n" + number);
        while (true) {
//            for (k = 0; k < peakCount; k++) {                   // take current assignment
//                assignment[k].aminoAcid = peakAssignment[k].aminoAcid;
//                assignment[k].feasibility = peakAssignment[k].feasibility;
//            }

            pCount = swapList.size();
            for (i=0; i<pCount; i++) {                 // swapLIst.size() = count; 
                aa1 = peakAssignment[swapList.get(i)].aminoAcid;
                for (j = 0; j < peakCount; j++) {
                    if (j == swapList.get(i)) {         // don't swap the same peaks :)
                        delta[j][j] = Double.POSITIVE_INFINITY;
                        continue;
                    }
                    aa2 = peakAssignment[j].aminoAcid;
                    if (!unavailableSelect)
                        if ((combinedScore[swapList.get(i)][aa2] == unavailableValue) ||
                            (combinedScore[j][aa1] == unavailableValue)) {
                            delta[swapList.get(i)][j] = Double.POSITIVE_INFINITY;
                            continue;
                        }
                    peakAssignment[swapList.get(i)].aminoAcid = aa2;            // change amino acids
                    peakAssignment[j].aminoAcid = aa1;
                    
                    tempScore = 0;
                    tempScore -= peakAssignment[swapList.get(i)].feasibility ? combinedScore[swapList.get(i)][aa1] : infeasibleScore;
                    tempScore -= peakAssignment[j].feasibility ? combinedScore[j][aa2] : infeasibleScore;
                    delta[swapList.get(i)][j] = deltaScore(peakAssignment, swapList.get(i), j, tempScore);
                    delta[j][swapList.get(i)] = delta[swapList.get(i)][j];

                    peakAssignment[swapList.get(i)].aminoAcid = aa1;
                    peakAssignment[j].aminoAcid = aa2;
                }
            }

            maxDelta = delta[0][0]; maxI = 0; maxJ = 0;
            for (i=0; i<peakCount; i++)
            for (j=0; j<peakCount; j++)
                if (delta[i][j] < maxDelta) {
                    maxDelta = delta[i][j];
                    maxI = i;
                    maxJ = j;
                }

//            for (int qw=0; qw<peakCount; qw++) {
//                for (int qe=0; qe<peakCount; qe++)
//                    System.out.print(delta[qw][qe]+" ");
//                System.out.println();
//            }

            // if we gained imporvement
            if (maxDelta < -0.0001) {
                totalScore += maxDelta;
                aa1 = peakAssignment[maxI].aminoAcid;
                peakAssignment[maxI].aminoAcid = peakAssignment[maxJ].aminoAcid;
                peakAssignment[maxJ].aminoAcid = aa1;
                System.out.println(maxI + " <-> " + maxJ + ": " + maxDelta + "; ");
                calculateDeltaScore(peakAssignment);
                
                swapList = new ArrayList<Integer>();         // update swap list
                swapList.add(maxI);
                swapList.add(maxJ);
                noe = noeList[maxI];
                pCount = noe.size();
                for (i=0; i<pCount; i++)
                    if (swapList.indexOf(noe.get(i)) == -1) swapList.add(noe.get(i));
                noe = noeList[maxJ];
                pCount = noe.size();
                for (i=0; i<pCount; i++)
                    if (swapList.indexOf(noe.get(i)) == -1) swapList.add(noe.get(i));
                k = swapList.size();
                for (j=2; j<k; j++) {
                    noe = noeList[swapList.get(j)];
                    pCount = noe.size();
                    for (i=0; i<pCount; i++)
                        if (swapList.indexOf(noe.get(i)) == -1) swapList.add(noe.get(i));
                }
//                for (i=0; i<swapList.size(); i++)
//                    System.out.print(swapList.get(i)+", ");
//                System.out.println(" :" + totalScore);
            } else break;
        }
        setUpdateScore();
    }

    private double calculateDeltaScore( Assignment[] assignment ) {
        int aa, i, j, l, peakCur;
        ArrayList<Integer> noe;

        for (i=0; i<peakCount; i++) {
            aa = assignment[i].aminoAcid;
            if ( combinedScore[i][aa] == unavailableValue ) {
                assignment[i].feasibility = false;              // if 10E+9 is selected then make it feasiblity false
                continue;                                       // NOTE: if we proceeded to this point then we have unavilableSelect = true
            }
            assignment[i].feasibility = true;
            noe = noeList[i];                                   // look for noe violation
            l = noe.size();
            for (j=0; j<l; j++) {
                peakCur = noe.get(j);
                if (binaryDistance[aa][assignment[peakCur].aminoAcid] == 0) {
                    assignment[i].feasibility = false;
                    if ( !infeasibleSelect ) return Double.NEGATIVE_INFINITY;
                    break;
                }
            }
        }
        return setTotalScore(assignment);
    }

    public Assignment[] getCompleteAssignment() {
        if ( !completed() || stopWorking ) return null;

        Assignment res[] = new Assignment[peakCount];
        for (int i=0; i<peakCount; i++) {
            res[i] = new Assignment(peakAssignment[i].aminoAcid, peakAssignment[i].feasibility);
        }
        return res;
    }

    public int getInfeasibleCount() {
        if ( !completed() || stopWorking ) return -1;

        int res = 0;
        for (int i=0; i<peakCount; i++)
            if ( peakAssignment[i].feasibility == false ) res++;
        return res;
    }

    // not yet optimised for newer methods. also may run incorrectly
    public void localSearchThreeOpt() {                 // consider 10e9 values !!!!
        if ( !completed() || stopWorking ) return;

        ArrayList<Integer> swapList = new ArrayList<Integer>(), noe;
        int i, j, k, w, aa1, aa2, aa3, maxI, maxJ, maxK, maxC, iteration = 0, pCount;
        final int combination = 5;
//        Assignment assignment[] = new Assignment[peakCount];
        double maxDelta, delta[][][][] = new double[peakCount][peakCount][peakCount][combination];
        for (i=0; i<peakCount; i++) {
//            assignment[i] = new Assignment(-1, false);
            swapList.add(i);
        }

        while (true) {
            iteration++;
//            for (k = 0; k < peakCount; k++) {
//                assignment[k].aminoAcid = peakAssignment[k].aminoAcid;
//                assignment[k].feasibility = peakAssignment[k].feasibility;
//            }

            // PERMUTATIONS:
            //      0: 1->2; 2->3; 3->1;
            //      1: 1->3; 2->1; 3->2;
            //      2: 1->1; 2->3; 3->2;
            //      3: 1->2; 2->1; 3->3;
            //      4: 1->3; 2->2; 3->1;

            pCount = swapList.size();
            for (i = 0; i < pCount; i++) {               // bunu swapliste gore yap
//                System.out.print(i);
                aa1 = peakAssignment[i].aminoAcid;
                for (j = 0; j < peakCount; j++) {
                    aa2 = peakAssignment[j].aminoAcid;
                    for (k = 0; k < peakCount; k++) {
                        aa3 = peakAssignment[k].aminoAcid;

                        if ( (swapList.get(i) == j) && (j == k) ) {
                            delta[j][j][j][0] = Double.POSITIVE_INFINITY;
                            delta[j][j][j][1] = Double.POSITIVE_INFINITY;
                            continue;
                        }

                        // CASE 0: 1->2; 2->3; 3->1;
                        delta[swapList.get(i)][j][k][0] = swapAndCalculate(swapList.get(i), j, k, aa2, aa3, aa1, 0);
                        delta[swapList.get(i)][k][j][0] = delta[swapList.get(i)][j][k][0];
                        delta[j][swapList.get(i)][k][0] = delta[swapList.get(i)][j][k][0];
                        delta[j][k][swapList.get(i)][0] = delta[swapList.get(i)][j][k][0];
                        delta[k][swapList.get(i)][j][0] = delta[swapList.get(i)][j][k][0];
                        delta[k][j][swapList.get(i)][0] = delta[swapList.get(i)][j][k][0];
                        peakAssignment[i].aminoAcid = aa1;  peakAssignment[j].aminoAcid = aa2;  peakAssignment[k].aminoAcid = aa3;

                        // CASE 1: 1->3; 2->1; 3->2;
                        delta[swapList.get(i)][j][k][1] = swapAndCalculate(swapList.get(i), j, k, aa3, aa1, aa2, 1);
                        delta[swapList.get(i)][k][j][1] = delta[swapList.get(i)][j][k][1];
                        delta[j][swapList.get(i)][k][1] = delta[swapList.get(i)][j][k][1];
                        delta[j][k][swapList.get(i)][1] = delta[swapList.get(i)][j][k][1];
                        delta[k][swapList.get(i)][j][1] = delta[swapList.get(i)][j][k][1];
                        delta[k][j][swapList.get(i)][1] = delta[swapList.get(i)][j][k][1];
                        peakAssignment[i].aminoAcid = aa1;  peakAssignment[j].aminoAcid = aa2;  peakAssignment[k].aminoAcid = aa3;

//                        // CASE 2: 1->1; 2->3; 3->2;
//                        if ( delta[i][j][k][2] == 0 ) {
//                            swapAndCalculate(delta, i, j, k, aa1, aa3, aa2, 2, assignment);
//                            if (delta[i][j][k][2] > maxDelta) { maxDelta = delta[i][j][k][2];   maxI = i;   maxJ = j;   maxK = k;   maxC = 2; }
//                            assignment[i].aminoAcid = aa1;  assignment[j].aminoAcid = aa2;  assignment[k].aminoAcid = aa3;
//                        }
//
//                        // CASE 3: 1->2; 2->1; 3->3;
//                        if ( delta[i][j][k][3] == 0 ) {
//                            swapAndCalculate(delta, i, j, k, aa2, aa1, aa3, 3, assignment);
//                            if (delta[i][j][k][3] > maxDelta) { maxDelta = delta[i][j][k][3];   maxI = i;   maxJ = j;   maxK = k;   maxC = 3; }
//                            assignment[i].aminoAcid = aa1;  assignment[j].aminoAcid = aa2;  assignment[k].aminoAcid = aa3;
//                        }
//
//                        // CASE 4: 1->3; 2->2; 3->1;
//                        if ( delta[i][j][k][4] == 0 ) {
//                            swapAndCalculate(delta, i, j, k, aa3, aa2, aa1, 4, assignment);
//                            if (delta[i][j][k][4] > maxDelta) { maxDelta = delta[i][j][k][4];   maxI = i;   maxJ = j;   maxK = k;   maxC = 4; }
//                            assignment[i].aminoAcid = aa1;  assignment[j].aminoAcid = aa2;  assignment[k].aminoAcid = aa3;
//                        }

                    }   // end FOR K
                }   // end FOR J
            }   // end FOR I

            maxI = 0;   maxJ = 0;   maxK = 0;   maxC = -1;
            maxDelta = delta[0][0][0][0];
            for (i=0; i<peakCount; i++)
            for (j=0; j<peakCount; j++)
            for (k=0; k<peakCount; k++)
            for (w=0; w<2; w++)
                if (delta[i][j][k][w] < maxDelta) {
                    maxDelta = delta[i][j][k][w];
                    maxI = i;   maxK = k;
                    maxJ = j;   maxC = w;
                }

            // if we have imporvement
            if (maxDelta < 0) {
                totalScore += maxDelta;
                switch (maxC) {
                    case 0: // 1->2; 2->3; 3->1;
                        aa1 = peakAssignment[maxI].aminoAcid;
                        peakAssignment[maxI].aminoAcid = peakAssignment[maxJ].aminoAcid;
                        peakAssignment[maxJ].aminoAcid = peakAssignment[maxK].aminoAcid;
                        peakAssignment[maxK].aminoAcid = aa1;
//                        System.out.println( maxI+"->"+maxJ+"; "+maxJ+"->"+maxK+"; "+maxK+"->"+maxI+"; " +maxDelta );
                        break;
                    case 1: // 1->3; 2->1; 3->2;
                        aa1 = peakAssignment[maxI].aminoAcid;
                        peakAssignment[maxI].aminoAcid = peakAssignment[maxK].aminoAcid;
                        peakAssignment[maxK].aminoAcid = peakAssignment[maxJ].aminoAcid;
                        peakAssignment[maxJ].aminoAcid = aa1;
//                        System.out.println( maxI+"->"+maxK+"; "+maxJ+"->"+maxI+"; "+maxK+"->"+maxJ+"; " +maxDelta );
                        break;
//                    case 2:  // 1->1; 2->3; 3->2;
//                        aa1 = peakAssignment[maxJ].aminoAcid;
//                        peakAssignment[maxJ].aminoAcid = peakAssignment[maxK].aminoAcid;
//                        peakAssignment[maxK].aminoAcid = aa1;
////                        System.out.println( maxJ+"<->"+maxK+"; " +maxDelta );
//                        break;
//                    case 3:  // 1->2; 2->1; 3->3;
//                        aa1 = peakAssignment[maxI].aminoAcid;
//                        peakAssignment[maxI].aminoAcid = peakAssignment[maxJ].aminoAcid;
//                        peakAssignment[maxJ].aminoAcid = aa1;
////                        System.out.println( maxJ+"<->"+maxI+"; " +maxDelta );
//                        break;
//                    case 4:  // 1->3; 2->2; 3->1;
//                        aa1 = peakAssignment[maxI].aminoAcid;
//                        peakAssignment[maxI].aminoAcid = peakAssignment[maxK].aminoAcid;
//                        peakAssignment[maxK].aminoAcid = aa1;
////                        System.out.println( maxI+"<->"+maxK+"; " +maxDelta );
//                        break;
                    default: System.out.println("AAAAAAAAAAA!!!!!");
                }
                calculateDeltaScore(peakAssignment);

                swapList = new ArrayList<Integer>();         // update swap list
                swapList.add(maxI);
                swapList.add(maxJ);
                swapList.add(maxK);
                noe = noeList[maxI];
                pCount = noe.size();
                for (i=0; i<pCount; i++)
                    if (swapList.indexOf(noe.get(i)) == -1) swapList.add(noe.get(i));
                noe = noeList[maxJ];
                pCount = noe.size();
                for (i=0; i<pCount; i++)
                    if (swapList.indexOf(noe.get(i)) == -1) swapList.add(noe.get(i));
                noe = noeList[maxK];
                pCount = noe.size();
                for (i=0; i<pCount; i++)
                    if (swapList.indexOf(noe.get(i)) == -1) swapList.add(noe.get(i));
                k = swapList.size();
                for (j=3; j<k; j++) {
                    noe = noeList[swapList.get(j)];
                    pCount = noe.size();
                    for (i=0; i<pCount; i++)
                        if (swapList.indexOf(noe.get(i)) == -1) swapList.add(noe.get(i));
                }

            } else break;

//            if (localSearchStop==2 || localSearchStop==3) break;
        }   // end WHILE
        setUpdateScore();
    }

    private double swapAndCalculate( int i, int j, int k, int a1, int a2, int a3, int seq) {

        double tempScore = 0;
        if (!unavailableSelect)
            if ((combinedScore[i][a1] == unavailableValue)
             || (combinedScore[j][a2] == unavailableValue)
             || (combinedScore[k][a3] == unavailableValue))
               return Double.NEGATIVE_INFINITY;

        if (peakAssignment[i].feasibility) {
            if (seq == 0) tempScore -= combinedScore[i][a3];
            if (seq == 1) tempScore -= combinedScore[i][a2];
        } else tempScore -= infeasibleScore;

        if (peakAssignment[j].feasibility) {
            if (seq == 0) tempScore -= combinedScore[j][a1];
            if (seq == 1) tempScore -= combinedScore[j][a3];
        } else tempScore -= infeasibleScore;

        if (peakAssignment[k].feasibility) {
            if (seq == 0) tempScore -= combinedScore[k][a2];
            if (seq == 1) tempScore -= combinedScore[k][a1];
        } else tempScore -= infeasibleScore;

        peakAssignment[i].aminoAcid = a1;
        peakAssignment[j].aminoAcid = a2;
        peakAssignment[k].aminoAcid = a3;
        return deltaScoreThree(peakAssignment, i, j, k, tempScore);
    }

    public double accuracy() {
        if ( !completed() || stopWorking ) return -1;
        
        double res = 0;
        for (int i=0; i<peakCount; i++)
            if ( peakAssignment[i].aminoAcid == i ) res++;
        return (double) Math.round( res / peakCount * 10000 ) / 100;
    }

    public void quicksort(int[] array, int[] index, int left, int right) {
        int pivot, leftIdx = left, rightIdx = right, k;
        if (right - left > 0) {
            pivot = (left + right) / 2;
            while (leftIdx <= pivot && rightIdx >= pivot) {
                while (array[leftIdx] > array[pivot] && leftIdx <= pivot) leftIdx++;
                while (array[rightIdx] < array[pivot] && rightIdx >= pivot) rightIdx--;
                k = array[leftIdx]; array[leftIdx] = array[rightIdx];   array[rightIdx] = k;
                k = index[leftIdx]; index[leftIdx] = index[rightIdx];   index[rightIdx] = k;
                leftIdx++;  rightIdx--;
                if (leftIdx - 1 == pivot) pivot = ++rightIdx;
                else if (rightIdx + 1 == pivot) pivot = --leftIdx;
            }
            quicksort(array, index, left, pivot - 1);
            quicksort(array, index, pivot + 1, right);
        }
    }

    // works only for the last total score calculation agreement. (total += feas ? score : penalty; )
    private double deltaScore( Assignment[] assignment, int peak1, int peak2, double score ) {
        double res = score;
        ArrayList<Integer> noe, listPeak = new ArrayList<Integer>(2), listNoe = new ArrayList<Integer>();
        int i, c, j, aa, curPeak, noePeak;
        boolean feas;

        // changed peaks analysis
        listPeak.add(peak1);
        listPeak.add(peak2);
        for (i=0; i<listPeak.size(); i++) {
            curPeak = listPeak.get(i);
            feas = true;

            // noe analysis
            noe = noeList[curPeak];
            aa = assignment[curPeak].aminoAcid;
            c = noe.size();
            for (j = 0; j < c; j++) {
                noePeak = noe.get(j);
                if (binaryDistance[aa][assignment[noePeak].aminoAcid] == 0) {        // if there is noe violation
                    feas = false;                                                       // set current's peak feas. to false
                    if (assignment[noePeak].feasibility)                             // if it was feas. assignment set it infeasible
                        res = res + infeasibleScore - combinedScore[noePeak][assignment[noePeak].aminoAcid];
                } else if (!assignment[noePeak].feasibility && noePeak!=peak1 && noePeak!=peak2 &&
                          (combinedScore[noePeak][assignment[noePeak].aminoAcid]!=unavailableValue) )
                    listNoe.add(noePeak);       // if it is infeas. assignment but have binary=1 so it is candidate to be feas.
            }

            // assignment analysis
            if (combinedScore[curPeak][aa] == unavailableValue || !feas)  // peak assigned to 10e9
                res += infeasibleScore;
            else                                                          // peak assigned to normal AA but feas. changed
                res += combinedScore[curPeak][aa];
        }

        // candidate peaks analysis. these peaks' assignments are candidate to be feasible
        for (i=0; i<listNoe.size(); i++) {
            curPeak = listNoe.get(i);
            feas = true;
            // only noe analysis
            noe = noeList[curPeak];
            c = noe.size();
            for (j = 0; j < c; j++)
                if (binaryDistance[assignment[curPeak].aminoAcid][assignment[noe.get(j)].aminoAcid] == 0) {
                    feas = false;
                    break;
                }
            if (feas)
                res = res - infeasibleScore + combinedScore[curPeak][assignment[curPeak].aminoAcid];
        }

        return res;
    }

    private double deltaScoreThree( Assignment[] assignment, int peak1, int peak2, int peak3, double score ) {
        double res = score;
        ArrayList<Integer> noe, listPeak = new ArrayList<Integer>(3), listNoe = new ArrayList<Integer>();
        int i, c, j, aa, curPeak, noePeak;
        boolean feas;

        // changed peaks analysis
        listPeak.add(peak1);
        listPeak.add(peak2);
        listPeak.add(peak3);
        for (i=0; i<listPeak.size(); i++) {
            curPeak = listPeak.get(i);
            feas = true;

            // noe analysis
            noe = noeList[curPeak];
            aa = assignment[curPeak].aminoAcid;
            c = noe.size();
            for (j = 0; j < c; j++) {
                noePeak = noe.get(j);
                if (binaryDistance[aa][assignment[noePeak].aminoAcid] == 0) {        // if there is noe violation
                    feas = false;                                                       // set current's peak feas. to false
                    if (assignment[noePeak].feasibility)                             // if it was feas. assignment set it infeasible
                        res = res + infeasibleScore - combinedScore[noePeak][assignment[noePeak].aminoAcid];
                } else if (!assignment[noePeak].feasibility && noePeak!=peak1 && noePeak!=peak2 &&
                           (combinedScore[noePeak][assignment[noePeak].aminoAcid]!=unavailableValue) )
                    listNoe.add(noePeak);       // if it is infeas. assignment but have binary=1 so it is candidate to be feas.
            }

            // assignment analysis
            if (combinedScore[curPeak][aa] == unavailableValue || !feas)  // peak assigned to 10e9
                res += infeasibleScore;
            else                                                          // peak assigned to normal AA but feas. changed
                res += combinedScore[curPeak][aa];
        }

        // candidate peaks analysis. these peaks' assignments are candidate to be feasible
        for (i=0; i<listNoe.size(); i++) {
            curPeak = listNoe.get(i);
            feas = true;
            // only noe analysis
            noe = noeList[curPeak];
            c = noe.size();
            for (j = 0; j < c; j++)
                if (binaryDistance[assignment[curPeak].aminoAcid][assignment[noe.get(j)].aminoAcid] == 0) {
                    feas = false;
                    break;
                }
            if (feas)
                res = res - infeasibleScore + combinedScore[curPeak][assignment[curPeak].aminoAcid];
        }

        return res;
    }

    public double whatScore() {
        FileReader fr;
        BufferedReader br;
        Assignment[] ass = new Assignment[peakCount];
        int i = 0;
        try {
            fr = new FileReader("New.txt");
            br = new BufferedReader(fr);
            String s;
            while ((s = br.readLine()) != null) {
                ass[i++] = new Assignment(Integer.parseInt(s), true);
            }
            return calculateDeltaScore(ass);
//            return setTotalScore(ass);
        } catch (Exception e) {}

        return -1;
    }
}