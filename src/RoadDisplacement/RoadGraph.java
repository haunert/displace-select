package RoadDisplacement;

import java.util.ArrayList;
import java.util.LinkedList;

import BuildingDisplacement.Building;

/**
 * A class which defines a graph which has a node for each road and an edge for each pair
 * of adjacent roads (= roads which share a node)
 */
public class RoadGraph {
	
	private LinkedList<Road> vertices;
	private LinkedList<RoadPair> edges;
	private int numComponents = 0;  // the number of connected components

	public RoadGraph() {
		vertices = new LinkedList<Road>();
		edges = new LinkedList<RoadPair>();
	}
	
	public void addVertex(Road r) {
		vertices.add(r);
	}
	
	public void addVerticesAndEdgesBetween(ArrayList<Road> roads) {
		for(int i=0; i<roads.size(); i++) {
			Road r_i = roads.get(i);
			vertices.add(r_i);
			for(int j=i+1; j<roads.size(); j++) {
				Road r_j = roads.get(j);
				if(r_i.isAdjacentTo(r_j)) {
					edges.add(new RoadPair(r_i,r_j));
				}
			}
		}
	}
	
	/**
	 * Adds a road r to the graph unless r has only adjacent roads already lying in the graph.
	 * If r is added the connected components of the adjacent roads are merged to one. Returns true
	 * if r has been added.
	 */
	public boolean addVertexIfAdmissibleAndUpdateComponents(Road r) {
		boolean isAdmissible = false;
		int counter=0;
		LinkedList<Road> component = new LinkedList<Road>();
		for(Road rOther : r.getNeighboringRoads()) {
			if(rOther.getComponent()==null) {
				isAdmissible = true;
				break;
			}
			if(counter==0) {
				component = rOther.getComponent();
				counter++;
			}
			else {
				if(rOther.getComponent()!=component) {
					isAdmissible = true;
					break;
				}
			}
		}
		if(isAdmissible) {
			counter=0;
			for(Road rInGraph : vertices) {
				if(r.isAdjacentTo(rInGraph)) {
					edges.add(new RoadPair(r,rInGraph));
					if(counter==0) {
						component = rInGraph.getComponent();
						component.add(r);
						r.setComponent(component);
						counter++;
					}
					else {
						LinkedList<Road> oldComponent = rInGraph.getComponent();
						if(oldComponent!=component) {
							numComponents--;
							for(Road rInOldComponent: oldComponent) {
								rInOldComponent.setComponent(component);
								component.add(rInOldComponent);
							}
						}
					}
				}
				
			}
			addVertex(r);
			if(r.getComponent()==null) {
				r.setComponent(new LinkedList<Road>());
				r.getComponent().add(r);
				numComponents++;
			}
		}
		return isAdmissible;
	}
	
	public void addEdge(RoadPair rp) {
		edges.add(rp);
	}
	
	public boolean containsVertex(Road road) {
		for(Road r : vertices) {
			if(r==road) {
				return true;
			}
		}
		return false;
	}
	
	public boolean containsEdge(Road r1, Road r2) {
		for(RoadPair rp : edges) {
			if( (rp.getR1()==r1 && rp.getR2()==r2) || (rp.getR1()==r2 && rp.getR2()==r1)) {
				return true;
			}
		}
		return false;
	}
	
	/**
	 * Determines the connected components.
	 */
	public void depthFirstSearch() {
		for(Road r : vertices) {
			r.setVisited(false);
		}
		for(Road r : vertices) {
			if( ! r.isVisited() ) {
				LinkedList<Road> component = new LinkedList<Road>();
				depthFirstSearchVisit(r, component);
				numComponents++;
				//components.add(component);
			}
		}
	}
	
	public void depthFirstSearchVisit(Road r, LinkedList<Road> component) {
		if( ! r.isVisited()) {
			r.setVisited(true);
			r.setComponent(component);
			component.add(r);
			for(Road rNeighbor : getNeighborsOf(r)) {
				depthFirstSearchVisit(rNeighbor, component);
			}
		}
	}


	public LinkedList<Road> getVertices() {
		return vertices;
	}

	public LinkedList<RoadPair> getEdges() {
		return edges;
	}
	
	public LinkedList<Road> getNeighborsOf(Road r){
		LinkedList<Road> neighbors = new LinkedList<Road>();
		for(RoadPair rp : edges) {
			if(rp.getR1()==r) {
				neighbors.add(rp.getR2());
			}
			else if(rp.getR2()==r) {
				neighbors.add(rp.getR1());
			}
		}
		return neighbors;
	}

	public int getNumComponents() {
		return numComponents;
	}
	
	public boolean isConnected() {
		return numComponents==1;
	}

}
