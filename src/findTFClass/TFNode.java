package findTFClass;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class TFNode {
	private List<Path> files;
	private String tfclass;
	private List<TFNode> children;
	private TFNode parent;
	
	public TFNode(String tfclass, List<Path> files, TFNode parent) {
		this.tfclass = tfclass;
		this.files = files;
		this.children = new ArrayList<>();
		this.parent = parent;
	}

	public TFNode(String tfclass, TFNode parent) {
		this.tfclass = tfclass;
		this.files = new ArrayList<>();
		this.children = new ArrayList<>();
		this.parent = parent;
	}
	
	public TFNode() {
		this.children = new ArrayList<>();
		this.files = new ArrayList<>();
	}
	
	public String getTfclass() {
		return tfclass;
	}

	public void setTfclass(String tfclass) {
		this.tfclass = tfclass;
	}

	public List<Path> getFiles() {
		return files;
	}

	public void setFiles(List<Path> files) {
		this.files = files;
	}

	public void addFile(Path file) {
		this.files.add(file);
	}
	
	public void addChild(TFNode child) {
		this.children.add(child);
	}
	
	public List<TFNode> getChildren(){
		return this.children;
	}
	
	public TFNode getParent() {
		return parent;
	}

	public void setParent(TFNode parent) {
		this.parent = parent;
	}
	
	public List<TFNode> getLeaves(){
		List<TFNode> leafNodes = new ArrayList<>();
		if(this.getChildren().isEmpty()) {
			leafNodes.add(this);
		}
		else {
			for(TFNode child : this.getChildren()) {
				leafNodes.addAll(child.getLeaves());
			}
		}
		return leafNodes;
	}
}