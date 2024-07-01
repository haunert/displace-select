package BuildingDisplacement;

import java.util.LinkedList;

/**
 * A class which defines a graph with a node for each building and an edge for each
 * pair of adjacent buildings.
 */
public class BuildingGraph {
	
	private LinkedList<Building> vertices;
	private LinkedList<BuildingPair> edges;
	private LinkedList<LinkedList<Building>> components;

	public BuildingGraph() {
		vertices = new LinkedList<Building>();
		edges = new LinkedList<BuildingPair>();
		components = new LinkedList<LinkedList<Building>>();
	}
	
	public void addVertex(Building b) {
		vertices.add(b);
	}
	
	public void addEdge(BuildingPair bp) {
		edges.add(bp);
	}
	
	public boolean containsVertex(Building building) {
		for(Building b : vertices) {
			if(b==building) {
				return true;
			}
		}
		return false;
	}
	
//	public boolean containsEdge(BuildingPair pair) {
//		for(BuildingPair bp : edges) {
//			if(bp==pair) {
//				return true;
//			}
//		}
//		return false;
//	}
	
	public boolean containsEdge(Building b1, Building b2) {
		for(BuildingPair bp : edges) {
			if( (bp.getB1()==b1 && bp.getB2()==b2) || (bp.getB1()==b2 && bp.getB2()==b1)) {
				return true;
			}
		}
		return false;
	}
	
	/**
	 * Determines the connected components.
	 */
	public void depthFirstSearch() {
		for(Building b : vertices) {
			b.setVisited(false);
		}
		for(Building b : vertices) {
			if( ! b.isVisited() ) {
				LinkedList<Building> component = new LinkedList<Building>();
				depthFirstSearchVisit(b, component);
				components.add(component);
			}
		}
	}
	
	public void depthFirstSearchVisit(Building b, LinkedList<Building> component) {
		if( ! b.isVisited()) {
			b.setVisited(true);
			b.setComponent(component);
			component.add(b);
			for(Building bNeighbor : b.getNeighboringBuildings()) {
				depthFirstSearchVisit(bNeighbor, component);
			}
		}
	}

	public LinkedList<LinkedList<Building>> getComponents() {
		return components;
	}

	public LinkedList<Building> getVertices() {
		return vertices;
	}

	public LinkedList<BuildingPair> getEdges() {
		return edges;
	}

}
