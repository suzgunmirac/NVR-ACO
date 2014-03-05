package aco_sba;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.*;
import java.lang.*;

public class NewClass {

    public static final int size = 144; // depends on amino acid
    public static final int peakCount = 335;
   
    public static void main(String[] args) {
        FileReader fr;                              // "-1": peak numbering is from 1 in file
        BufferedReader br;                          // so, peak=peak-1, noePeak=noePeak-1
        try {
            fr = new FileReader("Protein/mbp_real_data/NOE_List - Copy.txt");
            br = new BufferedReader(fr);
            String s;
            StringTokenizer st;
            int peak, row = 0, i, count, noePeak;
            
         //   @SuppressWarnings("unchecked")
            ArrayList<Integer> noeList[];
            noeList = new ArrayList[peakCount];
            for (i=0; i<peakCount; i++) {
                noeList[i] = new ArrayList<Integer>();
            }

            int temp1, temp2;
            while ((s = br.readLine()) != null) {
                st = new StringTokenizer(s);
//                count = Integer.parseInt(st.nextToken());                   // number of noe peaks
//                peak = Integer.parseInt(st.nextToken()) - 1;                // number of peak
//                for (i = 0; i < count; i++) {
//                    noePeak = Integer.parseInt(st.nextToken()) - 1;
//                    if (noeList[peak].indexOf(noePeak) == -1) {
//                        noeList[peak].add(noePeak);
//                        noeList[noePeak].add(peak);
//                    }
//                }
                temp1 = Integer.parseInt(st.nextToken()) - 1;
                temp2 = Integer.parseInt(st.nextToken()) - 1;
                if (temp1 == temp2) continue;
                noeList[temp1].add(temp2);
                noeList[temp2].add(temp1);
            }
            
            String temp;
            FileWriter fw = new FileWriter("Noe_List-.txt");
            for (i = 0; i < noeList.length; i++) {
                temp1 = noeList[i].size();
                if (temp1 == 0) continue;
                temp = "";
                count = temp1;
                for (int j = 0; j < temp1; j++) {
                    if (noeList[i].get(j) > i)
                        temp += " " + (noeList[i].get(j) + 1);
                    else 
                        count--;
                }
                if (count != 0)
                    fw.write(count + " " + (i+1) + temp + " \n");
            }
            fw.close();
                    
            System.out.println("noeListFile is succesfully migrated.");
            br.close();
            fr.close();
            //return 1;

        } catch (Exception exc) { exc.printStackTrace(System.out); return; }
    }
}