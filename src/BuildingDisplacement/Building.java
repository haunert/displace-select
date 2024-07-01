package BuildingDisplacement;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Vector;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.impl.CoordinateArraySequence;

import RoadDisplacement.RoadPoint;
import gurobi.GRBVar;
import general.*;

public class Building {
	
	private Vector<BuildingPoint> shell;   // exterior boundary
	private Vector<Vector<BuildingPoint>> holes;   // interior boundaries (if they exist)
    private String type;
    private int id;
    private GRBVar z;   // encodes the selection status of the building
    private boolean isSelected; // encodes the selection status of the building again (since Gurobi sometimes "forgets" the variable values)
    private double weight;   

    private boolean visited = false;  // for the depth-first search in the class BuildingGraph
    private LinkedList<Building> neighbors;  // the adjacent buildings (for the class BuildingGraph)
    private LinkedList<Building> component;  // the connected component containing this building in the BuildingGraph
    
    private boolean isFixed = false; // for the heuristic: if true, the selection status (isSelected) is not changed anymore
    private double zLP;  // double-value of the Gurobi variable in the relaxation
    
    String selectionWareJones;  // "selected", "unselected" oder "unknown" (selection of the building according to the publication by Ware et al. (2003))
    
	public Building(Vector<BuildingPoint> shell, Vector<Vector<BuildingPoint>> holes, 
			String type) {
		this.shell = shell;
		this.holes = holes;
		this.type = type;
		this.isSelected = true;
		this.weight = 1;   
	}


	public String getType() {
		return type;
	}


	public void setType(String type) {
		this.type = type;
	}


	public int getId() {
		return id;
	}


	public void setId(int id) {
		this.id = id;
	}


	public GRBVar getZ() {
		return z;
	}


	public void setZ(GRBVar z) {
		this.z = z;
	}


	public double getWeight() {
		return weight;
	}


	public void setWeight(double weight) {
		this.weight = weight;
	}


	public Vector<BuildingPoint> getShell() {
		return shell;
	}


	public Vector<Vector<BuildingPoint>> getHoles() {
		return holes;
	}


	public void setHoles(Vector<Vector<BuildingPoint>> holes) {
		this.holes = holes;
	}
	
	
	/**
	 * Returns a list containing all segments (exterior and interior edges) 
	 * of this building.
	 */
	public LinkedList<Segment> getAllSegments(){
		LinkedList<Segment> segments = new LinkedList<Segment>();
		
		// Iterate over the exterior edges
        for(int i = 0; i < shell.size() - 1; i++) {
            segments.add(new Segment(shell.get(i), shell.get(i+1)));
        }
        segments.add(new Segment(shell.get(0), shell.get(shell.size()-1)));
        
        // Iterate over the interior edges
        for(Vector<BuildingPoint> hole : holes){
        	for(int i = 0; i < hole.size() - 1; i++) {
                segments.add(new Segment(hole.get(i), hole.get(i+1)));
            }
            segments.add(new Segment(hole.get(0), hole.get(hole.size()-1)));
        }
        
		return segments;
	}


	public boolean isSelected() {
		return isSelected;
	}


	public void setSelected(boolean isSelected) {
		this.isSelected = isSelected;
	}
	
	
	/**
	 * Removes the Steiner points which have been associated with this building during 
	 * the Delaunay triangulation.
	 */
	public void removeDelaunayPoints() {
		for(int i=0; i<shell.size(); i++) {
			if(shell.get(i).isDelaunay()) {
				shell.remove(i);
				i--;
			}
		}
		for(Vector<BuildingPoint> hole : holes) {
			for(int i=0; i<hole.size(); i++) {
				if(hole.get(i).isDelaunay()) {
					hole.remove(i);
					i--;
				}
			}
		}
	}
	
	
	public Polygon getPolygon() {
		GeometryFactory gf = new GeometryFactory();
		Coordinate[] coords_shell = new Coordinate[shell.size()+1];
		for(int i=0; i<shell.size(); i++) {
			coords_shell[i] = shell.get(i).getCoordinate();
		}
		coords_shell[shell.size()] = shell.get(0).getCoordinate();
		CoordinateArraySequence cas_shell = new CoordinateArraySequence(coords_shell);
		LinearRing shell_ring = new LinearRing(cas_shell, gf);
		LinearRing[] holes_ring = new LinearRing[holes.size()];
		for(int i=0; i<holes.size(); i++) {
			Vector<BuildingPoint> hole = holes.get(i);
			Coordinate[] coords_hole = new Coordinate[hole.size()+1];
			for(int j=0; j<hole.size(); j++) {
				coords_hole[j] = hole.get(j).getCoordinate();
			}
			coords_hole[hole.size()] = hole.get(0).getCoordinate();
			CoordinateArraySequence cas_hole = new CoordinateArraySequence(coords_hole);
			holes_ring[i] = new LinearRing(cas_hole, gf);
		}
		return new Polygon(shell_ring, holes_ring, gf);
	}
	
	
	/**
	 * Returns the number of adjacent buildings.
	 */
	public int getNumNeighbors() {
		int numNeighbors = 0;
		for(BuildingPoint bp : shell) {
			if(bp.getBuildings().size()-1>numNeighbors) {
				numNeighbors = bp.getBuildings().size()-1;
			}
		}
		return numNeighbors;
	}
	
	
	/**
	 * Replaces a point oldPoint on one of the boundaries (ring) by a different
	 * point (newPoint).
	 */
	public void replacePointInRing(BuildingPoint oldPoint, BuildingPoint newPoint, 
			Vector<BuildingPoint> ring) {
		boolean replaced = false;
		for(int i=0; i < ring.size(); i++) {
			BuildingPoint bp = ring.get(i);
			if(bp==oldPoint) {
				if(bp.compareTo(oldPoint) != 0) {
					System.out.println("ERROR in method replacePointInRing in class Building!");
				}
				ring.set(i,newPoint);
				replaced = true;
				break;
			}
		}
		if(replaced==false) {
			System.out.println("ERROR in method replacePointInRing in class Building, the point to be replaced was not found!");
		}
	}
	
	
	/**
	 * Returns a list of the adjacent buildings (unlike getNeighboringBuildings(), recomputes 
	 * this list).
	 */
	public LinkedList<Building> findNeighboringBuildings(){
		neighbors = new LinkedList<Building>();
		for(BuildingPoint bp : shell) {
			ArrayList<Building> buildings = bp.getBuildings();
			for(Building b : buildings) {
				if(b != this) {
					int counter = 0;
					for(Building neighbor : neighbors) {
						if(b==neighbor) { break; }
						counter++;
					}
					if(counter==neighbors.size()) {
						neighbors.add(b);
					}
				}
			}
		}
		return neighbors;
	}

	/**
	 * Returns a list of the adjacent buildings.
	 */
	public LinkedList<Building> getNeighboringBuildings() {
		return neighbors;
	}


	public boolean isVisited() {
		return visited;
	}


	public void setVisited(boolean visited) {
		this.visited = visited;
	}
	
	
	public LinkedList<Building> getComponent() {
		return component;
	}


	public void setComponent(LinkedList<Building> component) {
		this.component = component;
	}


	public boolean isFixed() {
		return isFixed;
	}


	public void setFixed(boolean isFixed) {
		this.isFixed = isFixed;
	}


	public double getzLP() {
		return zLP;
	}


	public void setzLP(double zLP) {
		this.zLP = zLP;
	}


	public String getSelectionWareJones() {
		return selectionWareJones;
	}


	public void setSelectionWareJones(String selectionWareJones) {
		this.selectionWareJones = selectionWareJones;
	}

}
