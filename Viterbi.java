package COMP462;
import java.io.*;
import java.util.*;

public class Viterbi {
	
	//States: H = Hydrophobic, I = Hydrophilic, M = Mixed.
	//Initial probabilities are equal for all three states.
	private static char[] states = {'H', 'I', 'M'};
	private static double[] initialProbs = {0.33333, 0.33333, 0.33333};
	
	/*
	 * row 1 is hydrophobic
	 * row 2 is hydrophilic
	 * row 3 is mixed
	 * column 1 is hydrophobic
	 * column 2 is hydrophilic
	 * column 3 is mixed
	 */	
	private static double[][] transitionProbs = {{0.8, 0.04, 0.16},
												{0.0375, 0.875, 0.0875},
												{0.071428571, 0.071428571, 0.857142857}};
		
	public static void main(String[] args) throws FileNotFoundException {
		
		
		File file = new File("/Users/nikolanisic/Desktop/FASTA.txt");
		ArrayList<String[]> parsedSeqs = FastaParser(file);
		for(int i = 0; i < parsedSeqs.size(); i++) {
			System.out.println("The protein is: " + parsedSeqs.get(i)[0]);
			System.out.println("The sequence is: " + parsedSeqs.get(i)[1]);
			char[] finalSequence = seqToCharArray(parsedSeqs.get(i));
			ViterbiAlg(finalSequence, states, initialProbs, transitionProbs);
	}
}
	
	/**
	 * 
	 * @param file					Fasta file containing the set of observations
	 * @return						Parsed Fasta file. Each species is dedicated a position in the array list.
	 * 								Each array in the array list contains the species description in slot [0] and the amino acid sequence of interest in slot [1].
	 * @throws FileNotFoundException	FileNotFoundException thrown in case of exception.
	 */
	public static ArrayList<String[]> FastaParser(File file) throws FileNotFoundException {
		Scanner sc = new Scanner(file);
		sc.useDelimiter(">");
		ArrayList<String[]> seqs = new ArrayList<String[]>();
		
		while(sc.hasNext()) {
			seqs.add(sc.next().split("\n", 2));
		}
		//Removing all whitespace in sequences
		for(String[] s : seqs) {
			s[1] = s[1].replaceAll("\\s+", "");
		}
		
		sc.close();
		return seqs;
	}
	

	/**
	 * 
	 * @param seq		 An array of Strings. Each String is 
	 * @return		     Same sequences as arrays of characters.
	 * @author 			 Nikola Nisic
	 */
	public static char[] seqToCharArray(String[] seq) {
		char[] sequence = new char[seq[1].length()];
		for(int i = 0; i < seq[1].length(); i++) {
			sequence[i] = seq[1].charAt(i);
		}
		return sequence;
	}
	
	/**
	 * 
	 * @param observation	Current observed character
	 * @param state			Either Hydrophobic, Hydrophilic, or Mixed
	 * @return 				Probability of emitting the current observation from state in question
	 * @author 				Nikola Nisic
	 */
	public static double getEmissionProb(char observation, char state) {
		double probability = 0;
		if(observation == 'A' || observation == 'V' || observation == 'I' || observation == 'L' || observation == 'M' || observation == 'F' || observation == 'Y' || observation == 'W') {
			if(state == 'H') {
				probability = 0.1125;
			}
			else if(state == 'I') {
				probability = 0.025;
			}
			else if(state == 'M') {
				probability = 0.05;
			}
		}
		else {
			if(state == 'H') {
				probability = 0.00833333;
			}
			else if(state == 'I') {
				probability = 0.066666;
			}
			else if(state == 'M') {
				probability = 0.05;
			}
		}
		return probability;
	}
	
	//TrellisNode class will be used for traceback process at end of Viterbi implementation
	//path field will be array of integers which represent the state
	private static class TrellisNode {
		public double pathProbability;
		public int[] path;
		
		public TrellisNode(int[] path, double pathProbability) {
			this.pathProbability = pathProbability;
			this.path = copy(path);
		}
	}
	
	/**
	 * 
	 * @param array		Array of integers representing the current most probable state sequence of the HMM
	 * @return			Same array of integers
	 */
	private static int[] copy(int[] array) {
		int[] newArray = new int[array.length];
		for(int j = 0; j < array.length; j++) {
			newArray[j] = array[j];
		}
		return newArray;
	}
	
	/**
	 * 
	 * @param array		Array of integers representing the current most probable state sequence of the HMM
	 * @param cell		Integer representing the newest most probable state to add to the growing state sequence
	 * @return			Same array as before with the "cell" integer appended
	 * @author 			Nikola Nisic
	 */
	private static int[] copy(int[] array, int cell) {
		int[] newArray = new int[array.length+1];
		for(int j = 0; j < array.length; j++) {
			newArray[j] = array[j];
		}
		newArray[array.length] = cell;
		return newArray;
	}
	
	/**
	 * 
	 * @param obs			Array containing each individual observed amino acid
	 * @param stateSpace		The possible states 'H', 'I', 'M'
	 * @param initProbs		The probability of beginning in each state
	 * @param transProbs		The transition matrix containing the probabilities of transition between any given states
	 * @return				Array of integers representing the 3 states. 0 = H, 1 = I, 2 = M
	 * @author Nikola Nisic
	 * 
	 * Also prints out the results in readable form.
	 */
	public static int[] ViterbiAlg(char[] obs, char[] stateSpace, double[] initProbs, double[][] transProbs) {
		
		TrellisNode[] T = new TrellisNode[stateSpace.length];
		
		for(int i = 0; i < stateSpace.length; i++) {
			int[] array = new int[1];
			array[0] = i;
			T[i] = new TrellisNode(array, initProbs[i]*getEmissionProb(obs[i], stateSpace[i]));
		}
		
		for(int i = 1; i < obs.length; i++) {
			TrellisNode[] V = new TrellisNode[stateSpace.length];
			
			for(int j = 0; j < stateSpace.length; j++) {
					int[] argMax = new int[0];
					double max = 0;
					
					for(int s = 0; s < stateSpace.length; s++) {
							int[] path = copy(T[s].path);
							double probability = T[s].pathProbability;
							double currentProb = getEmissionProb(obs[i], stateSpace[j])*transProbs[s][j];
							probability = probability*currentProb*20;										//scaled by 20 to avoid underflow
							
							if(probability > max) {
								if(path.length == obs.length) {
									argMax = copy(path);
									max = probability;
								}
								else {
									argMax = copy(path, j);
								}
								max = probability;
							}
					}
					V[j] = new TrellisNode(argMax, max);
			}
			T = V;
		}
		
		//find maximum of final states
		int[] argMax = new int[0];
		double max = 0;
		
		for(int i = 0; i < stateSpace.length; i++) {
			int[] path = copy(T[i].path);
			double probability = T[i].pathProbability;
			
			if(probability > max) {
				argMax = copy(path);
				max = probability;
			}	
		}
		
		System.out.print("Most likely sequence of states: [");
		for(int i = 0; i < argMax.length; i++) {
			System.out.print(states[argMax[i]]);
		}
		System.out.println("]. \nProbability of state sequence: " + max + "\n");
		return argMax;
	}
}
