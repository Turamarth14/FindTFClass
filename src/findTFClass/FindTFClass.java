package findTFClass;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

public class FindTFClass {

	static boolean hierarch = false;
	static Path tempDir;
	static Path tftreeDir;
	static Path fastaSeqFile;
	static Path outFile;
	static TFNode tfTreeRoot;
	static List<List<Pair<Path,String>>> tfTree;
	static List<Pair<String,String>> fastaSeqList;
	public static void main(String[] args) {
		if(args.length < 1){
			printHelp();
			return;
		}
		parseArgs(args);
		if(outFile == null) {
			outFile = Paths.get(fastaSeqFile.toString()+".match");
		}
		//tftreeDir = Paths.get("/home/felix/Downloads/tfclass3_trees/dbd");
		//fastaSeqFile = Paths.get("/home/felix/Dokumente/Sequenz-Klassifikation/TFClass-Daten_(Protein)/8.1.1_mammalia.v2.fasta.txt");
		//outFile = Paths.get("/home/felix/Schreibtisch/Ergebnisse.txt");
		//Creating directory for temporary files
		try {
			tempDir = Files.createTempDirectory("temp" + System.nanoTime());
		}
		catch(IOException e) {
			System.out.println("Error while creating temp directory");
			System.out.println(e.getMessage());
			return;
		}
		//Creating tfTree from directory
		try {
			if(hierarch) {
				buildTFClassHierarchy();
			}
			else {
				buildTFClassTree();
			}
		}
		catch(IOException e) {
			System.out.println("Error while reading the tfclass tree dir");
			System.out.println(e.getMessage());
			return;
		}
		//Read sequences from fasta file
		try {
			fastaSeqList = readFastaFile(fastaSeqFile);
			for(Pair<String,String> entry : fastaSeqList) {
				System.out.println(entry.getFirst());
				System.out.println(entry.getSecond());
			}
		}
		catch(IOException e) {
			System.out.println("Error while deleting temp directory");
			System.out.println(e.getMessage());
			return;
		}	
		//Classifying sequences
		List<Pair<String,String>> resultList;
		try {
			if(hierarch) {
				System.out.println("Level 0 = " + tfTree.get(0).size());
				System.out.println("Level 1 = " + tfTree.get(1).size());
				System.out.println("Level 2 = " + tfTree.get(2).size());
				System.out.println("Level 3 = " + tfTree.get(3).size());
				resultList = classifySeqsHierarchy();
			}
			else {
				resultList = classifySeqsTree();
			}
		}
		catch(IOException | InterruptedException e) {
			System.out.println("Error while classifying the sequences");
			System.out.println(e.getMessage());
			return;
		}	
		try {
			printResults(resultList);
		}
		catch(IOException e) {
			System.out.println("Error while writing results");
			System.out.println(e.getMessage());
			return;
		}
		//Deleting temporary files
		try {
			emptyDir(tempDir);
			Files.delete(tempDir);
		}
		catch(IOException e) {
			System.out.println("Error while deleting temp directory");
			System.out.println(e.getMessage());
			return;
		}
	}
	private static List<Pair<String,String>> classifySeqsHierarchy() throws IOException,InterruptedException{
		List<Pair<String,String>> classifyList = new ArrayList<>();
		int level = tfTree.size()-1;
		//Go through all levels or until all sequences are classified
		while(tfTree.get(level).size() > 0 && !fastaSeqList.isEmpty()) {
			System.out.println("Level : " + level + " #Sequences :" + fastaSeqList.size());
			List<Pair<String,Pair<String,Double>>> resultList = classifier(fastaSeqList, tfTree.get(level));
			for(Pair<String,Pair<String,Double>> mapping : resultList) {
				for(int i = 0; i < fastaSeqList.size(); i++) {
					if(fastaSeqList.get(i).getFirst().equals(mapping.getFirst())) {
						fastaSeqList.remove(i);
						break;
					}
				}
				for(int i = 0; i < tfTree.get(level).size(); i++) {
					if(tfTree.get(level).get(i).getFirst().getFileName().toString().equals(mapping.getSecond().getFirst())) {
						classifyList.add(new Pair<>(mapping.getFirst(),tfTree.get(level).get(i).getSecond()+"\t"+mapping.getSecond().getSecond()));
						break;
					}
				}				
			}
			level--;
		}
		System.out.println("Level : " + level + " #Sequences :" + fastaSeqList.size());
		return classifyList;
	}
	private static List<Pair<String,String>> classifySeqsTree() throws IOException,InterruptedException{
		List<Pair<String,String>> classifyList = new ArrayList<>();
		Set<TFNode> tfnodeSet = new HashSet<>(tfTreeRoot.getLeaves());
		//Go through all levels or until all sequences are classified
		while(!tfnodeSet.isEmpty() && !fastaSeqList.isEmpty()) {
			System.out.println("#TFFiles : " + tfnodeSet.size() + " #Sequences :" + fastaSeqList.size());
			List<Pair<Path,String>> tfclassFileList = new ArrayList<>();
			for(TFNode node : tfnodeSet) {
				for(Path file : node.getFiles()) {
					tfclassFileList.add(new Pair<>(file,node.getTfclass()));
				}
			}
			List<Pair<String,Pair<String,Double>>> resultList = classifier(fastaSeqList, tfclassFileList);
			for(Pair<String,Pair<String,Double>> mapping : resultList) {
				for(int i = 0; i < fastaSeqList.size(); i++) {
					if(fastaSeqList.get(i).getFirst().equals(mapping.getFirst())) {
						fastaSeqList.remove(i);
						break;
					}
				}
				for(Pair<Path,String> file : tfclassFileList) {
					if(file.getFirst().getFileName().toString().equals(mapping.getSecond().getFirst())){
						classifyList.add(new Pair<>(mapping.getFirst(),file.getSecond()+"\t"+mapping.getSecond().getSecond()));
						break;
					}
				}				
			}
			Set<TFNode> tempSet = new HashSet<>();
			for(TFNode node : tfnodeSet) {
				tempSet.add(node.getParent());
			}
			tfnodeSet = tempSet;
		}
		System.out.println("#TFFiles : " + tfnodeSet.size() + " #Sequences :" + fastaSeqList.size());
		return classifyList;
	}
	private static List<Pair<String,Pair<String,Double>>> classifier(List<Pair<String,String>> fastaSeqs,List<Pair<Path,String>> tfclassFileList) throws IOException,InterruptedException{
		List<Pair<String,Pair<String,Double>>> classifyList = new ArrayList<>();
		emptyDir(tempDir);
		tfclassFileList.forEach(entry ->{
			try {
				Files.copy(entry.getFirst(), Paths.get(tempDir.toAbsolutePath().toString(),entry.getFirst().getFileName().toString()));
			} catch (IOException e) {
				System.out.println("Error while copying " + entry.getFirst().getFileName().toString() + " to temp directory");
				e.printStackTrace();
			}
		});
		//Executing HMMER
		Files.walk(tempDir).forEach(file ->{
			if(file.toFile().isDirectory()) {
				return;
			}
			ProcessBuilder pb = new ProcessBuilder("hmmbuild", "-n", file.getFileName().toString(), file.getFileName().toString() + ".hmm",file.getFileName().toString());
			pb.directory(tempDir.toFile());
			try {
				execProc(pb,false);
			}
			catch(IOException | InterruptedException e) {
				System.out.println("Error while executing hmmbuild for " + file.toString());
				System.out.println(e.getMessage());
				return;
			}
		});
		List<String> hmmFiles = Files.walk(tempDir).filter(file -> file.toString().endsWith(".hmm")).map( file -> file.toString()).collect(Collectors.toList());
		hmmFiles.add(0, "cat");
		ProcessBuilder pb = new ProcessBuilder(hmmFiles);
		pb.directory(tempDir.toFile());
		pb.redirectOutput(new File(tempDir.toFile(),"hmmDB"));
		try {
			execProc(pb,false);
		}
		catch(IOException e) {
			throw new IOException("Error while concatenating hmms", e);
		}
		catch(InterruptedException e) {
			throw new InterruptedException("Error while concatenating hmms");
		}
		pb = new ProcessBuilder("hmmpress", "hmmDB");
		pb.directory(tempDir.toFile());
		try {
			execProc(pb,false);
		}
		catch(IOException e) {
			throw new IOException("Error while executing hmmpress", e);
		}
		catch(InterruptedException e) {
			throw new InterruptedException("Error while executing hmmpress");
		}
		PrintWriter searchSeqsPW = new PrintWriter(new BufferedWriter(new FileWriter(new File(tempDir.toFile(), "searchSeqs"))));
		fastaSeqs.forEach(fasta ->{
			searchSeqsPW.println(">" + fasta.getFirst());
			searchSeqsPW.println(fasta.getSecond());
		});
		searchSeqsPW.close();
		pb = new ProcessBuilder("hmmscan", "hmmDB", "searchSeqs");
		pb.directory(tempDir.toFile());
		try {
			Process proc = pb.start();
			BufferedReader procIn = new BufferedReader(new InputStreamReader(proc.getInputStream()));
			String line, query;
			while((line = procIn.readLine()) != null) {
				if(line.startsWith("Query:")) {
					query = line.split("\\s+")[1];					
					while((line = procIn.readLine()) != null) {
						if(line.trim().startsWith("E-value")) {
							procIn.readLine();
							line = procIn.readLine();
							while(!line.trim().isEmpty()) {
								String[] strarray = line.trim().split("\\s+");
								if(strarray.length == 9) {
									String model = strarray[8];
									double score = Double.parseDouble(strarray[0]);
									classifyList.add(new Pair<>(query,new Pair<>(model,score)));
								}
								line = procIn.readLine();
							}
							break;
						}
					}
				}
			}
			procIn.close();
			procIn = new BufferedReader(new InputStreamReader(proc.getErrorStream()));
			while(procIn.readLine() != null) {
				procIn.readLine();
			}
			procIn.close();
			proc.waitFor();
		}
		catch(IOException | InterruptedException e) {
			System.out.println("Error while executing hmmscan");
			System.out.println(e.getMessage());
			return null;
		}	
		emptyDir(tempDir);
		return classifyList;
	}
	private static List<Pair<String,String>> readFastaFile(Path file) throws IOException{
		List<Pair<String,String>> seqList = new ArrayList<>();
		List<String> lines = Files.readAllLines(file);
		StringBuilder seq = new StringBuilder();
		String id = null; 
		for(int i = 0; i < lines.size(); i++) {
			if(lines.get(i).startsWith(">")) {
				if(id != null) {
					seqList.add(new Pair<>(id,seq.toString()));
				}
				id = lines.get(i).substring(1);
				seq.setLength(0);
			}
			else {
				seq.append(lines.get(i));
			}
		}
		if(id != null) {
			seqList.add(new Pair<>(id,seq.toString()));
		}
		return seqList;
	}
	private static void buildTFClassHierarchy() throws IOException{
		tfTree = new ArrayList<>();
		Files.walk(tftreeDir, 1).filter(f -> (Files.isDirectory(f) && !f.equals(tftreeDir))).forEach(f ->{
			System.out.println(f.toAbsolutePath().toString());
			try {
				readTFClassDir(f,1);
			}
			catch(IOException e) {
				System.out.println("Error while reading the reading dir " + f.toString());
				System.out.println(e.getMessage());
				return;
			}
		});
	}
	private static void buildTFClassTree() throws IOException{
		tfTreeRoot = new TFNode();
		Files.walk(tftreeDir, 1).filter(f -> (Files.isDirectory(f) && !f.equals(tftreeDir))).forEach(f ->{
			System.out.println(f.toAbsolutePath().toString());
			try {
				readTFClassDirTree(f,tfTreeRoot);
			}
			catch(IOException e) {
				System.out.println("Error while reading the reading dir " + f.toString());
				System.out.println(e.getMessage());
				return;
			}
		});
	}
	private static void readTFClassDirTree(Path dir, TFNode node) throws IOException{
		ArrayList<Path> fileList = new ArrayList<>();
		Files.walk(dir, 1).forEach(f -> {
			if(Files.isDirectory(f) && !f.equals(dir)) {
				try {
					TFNode child = new TFNode(f.getFileName().toString(), node);
					node.addChild(child);
					readTFClassDirTree(f,child);
				}
				catch(IOException e) {
					System.out.println("Error while reading the reading dir " + f.toString());
					System.out.println(e.getMessage());
					return;
				}				
			}
			/**Assumes that there is at most one file with ending "phyml-output.fasta.txt" in each dir**/
			else if(f.toString().endsWith("phyml-output.fasta.txt")) {
				fileList.add(f);
				
			}
		});
		if(fileList.size() > 1) {
			Iterator<Path> it = fileList.iterator();
			while(it.hasNext()) {
				Path temp = it.next();
				if(temp.getFileName().toString().contains("slim")) {
					it.remove();
				}
				else {
					node.addFile(temp);
				}
			}
		}
		else if(fileList.size() == 1){
			node.addFile(fileList.get(0));
		}
	}
	private static void printResults(List<Pair<String,String>> classifyList) throws IOException{
		PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(outFile.toFile())));
		pw.println("FastaID\tTFClass\tE-Value");
		for(Pair<String,String> mapping : classifyList) {
			pw.println(mapping.toString());
		}
		if(!fastaSeqList.isEmpty()) {
			pw.println(fastaSeqList.size() + " sequences could not be classified:");
			for(Pair<String,String> fasta : fastaSeqList) {
				pw.println(fasta.getFirst()+ "\n" + fasta.getSecond());
			}
		}
		pw.close();
	}
	private static void readTFClassDir(Path dir, int level) throws IOException{
		if(tfTree.size() < level) {
			tfTree.add(new ArrayList<>());
		}
		Files.walk(dir, 1).forEach(f -> {
			if(Files.isDirectory(f) && !f.equals(dir)) {
				try {
					readTFClassDir(f,level+1);
				}
				catch(IOException e) {
					System.out.println("Error while reading the reading dir " + f.toString());
					System.out.println(e.getMessage());
					return;
				}				
			}
			/**Assumes that there is at most one file with ending "phyml-output.fasta.txt" in each dir**/
			else if(f.toString().endsWith("phyml-output.fasta.txt")) {
				tfTree.get(level-1).add(new Pair<>(f,dir.getFileName().toString()));
			}
		});
	}
	private static void emptyDir(Path dir) throws IOException {
		Files.walk(dir).forEach( file ->{
			if(file.equals(dir)) {
				return;
			}
			try{
				Files.delete(file);
			} catch (IOException e) {
				System.out.println("Error while deleting "+ file.toString());
				System.out.println(e.getMessage());
			}
		});
	}
	private static void execProc(ProcessBuilder pb, boolean output) throws IOException, InterruptedException{
		String line;
		Process proc = pb.start();
		BufferedReader procIn = new BufferedReader(new InputStreamReader(proc.getInputStream()));
		while((line = procIn.readLine()) != null) {
			if(output) {
				System.out.println(line);
			}
		}
		procIn.close();
		procIn = new BufferedReader(new InputStreamReader(proc.getErrorStream()));
		while((line = procIn.readLine()) != null) {
			if(output) {
				System.out.println(line);
			}
		}
		procIn.close();
		proc.waitFor();
	}
	private static void parseArgs(String[] args){
		for(int i = 0; i < args.length -1; i++){
			switch(args[i]){
			case "-out":
				outFile = Paths.get(args[i+1]);
				i++;
				break;
			case "-hierarch":
				hierarch = true;
				break;
			}		
		}
		tftreeDir = Paths.get(args[args.length-2]);
		fastaSeqFile = Paths.get(args[args.length-1]);
	}
	private static void printHelp(){
		System.out.println("Usage: java -jar FindTFClass.jar {Options} tftreeDir fastaSeqFile");
		System.out.println("Options:");
		System.out.println("-out file\t name of outputfile (default fastaSeqFile.match)");
		System.out.println("-hierarch\t use hierarchical structure of the tftreeDir instead of a tree (default)");
	}

}