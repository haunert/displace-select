package RoadDisplacement;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Vector;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.LineString;
import com.vividsolutions.jts.geom.impl.CoordinateArraySequence;

import gurobi.GRBVar;

/**
 * A class which represents a road, geometrically defined by its corners
 */
public class Road {
    private Vector<RoadPoint> vertices;
    private String type;
    private int id;
    private GRBVar z;   // variable for the selection status of the road
    private double weight;   
    private boolean isSelected; // encodes the selection status of the road again (since Gurobi sometimes "forgets" the variable values)
    
    // Variables for the flow formulation => connected road network
    private GRBVar isSink;   // is this road a sink?
    private LinkedList<GRBVar> inFlow = new LinkedList<GRBVar>();  // contains the variables f_qp (incoming flow)
    private LinkedList<GRBVar> outFlow = new LinkedList<GRBVar>(); // contains the variables f_pq (outgoing flow)
    
    private boolean isFixed = false; // for the heuristic: if true, the selection status (isSelected) is not changed anymore
    private double zLP;  // double-value of the Gurobi variable in the relaxation
    private int compID; // auxiliary variable for the method solveViaRelaxationSimpleThreshold
    private LinkedList<Road> component;
    private boolean visited = false;  // for the depth-first search in class RoadGraph
    
    public Road(Vector<RoadPoint> v, String s) {
        vertices = v;
        type = s;
        this.isSelected = true;
    }
    
    public Vector<RoadPoint> getVertices() {
        return vertices;
    }
    
    public String getType() {
        return type;
    }
    
    public LineString getAsLineString() {
        CoordinateArraySequence points = new CoordinateArraySequence(vertices.size());
        for (int i = 0; i < vertices.size(); i++) {
            RoadPoint p = (RoadPoint) vertices.get(i);
            points.setOrdinate(i, 0, p.getX());
            points.setOrdinate(i, 1, p.getY());
        }
        GeometryFactory gf = new GeometryFactory();
        LineString ls = new LineString(points, gf);
        return ls;
    }
        
    public double getWeight() {
//    	double weight = 0;
//    	for (int i = 0; i < vertices.size() - 1; i++) {
//    		RoadPoint p1 = vertices.get(i);
//    		RoadPoint p2 = vertices.get(i + 1);
//    		weight += Math.hypot(p2.getX() - p1.getX(), p2.getY() - p1.getY());
//    	}
    	return weight;
    }

	public int getId() {
		return id;
	}

	public void setId(int id) {
		this.id = id;
	}
	
	/**
	 * Replaces a point in the list vertices by a different point.
	 */
	public void replaceRoadPoint(RoadPoint oldPoint, RoadPoint newPoint) {
		boolean replaced = false;
		for(int i=0; i < vertices.size(); i++) {
			RoadPoint rp = vertices.get(i);
			if(rp==oldPoint) {
				if(rp.compareTo(oldPoint) != 0) {
					System.out.println("ERROR in method replaceRoadPoint in class Road!");
				}
				vertices.set(i,newPoint);
				replaced = true;
				break;
			}
		}
		if(replaced==false) {
			System.out.println("ERROR in method replaceRoadPoint in class Road, the point to be replaced was not found!");
		}
	}

	public GRBVar getZ() {
		return z;
	}

	public void setZ(GRBVar z) {
		this.z = z;
	}

	public void setWeight(double weight) {
		this.weight = weight;
	}
	
	public boolean isSelected() {
		return isSelected;
	}

	public void setSelected(boolean isSelected) {
		this.isSelected = isSelected;
	}
	
	/**
	 * Removes the Steiner points which have been associated with this road during the Delaunay triangulation.
	 */
	public void removeDelaunayPoints() {
		for(int i=0; i<vertices.size(); i++) {
			if(vertices.get(i).isDelaunay()) {
				vertices.remove(i);
				i--;
			}
		}
	}
	
	/**
	 * Returns the number of the original points of this road (without the Steiner
	 * points).
	 */
	public int getNumOriginalPoints() {
		int numPoints = 0;
		for(RoadPoint rp : vertices) {
			if(! rp.isDelaunay()) {
				numPoints++;
			}
		}
		return numPoints;
	}
	
	
	/**
	 * Returns a list of the adjacent roads (= the roads sharing a node with this road).
	 */
	public LinkedList<Road> getNeighboringRoads(){
		LinkedList<Road> neighboringRoads = new LinkedList<Road>();
		for(RoadPoint rp : vertices) {
			ArrayList<Road> roads_rp = rp.getRoads();
			for(Road road_rp : roads_rp) {
				if(road_rp != this) {
					int counter = 0;
					for(Road road : neighboringRoads) {
						if(road==road_rp) {
							break;
						}
						counter++;
					}
					if(counter==neighboringRoads.size()) {
						neighboringRoads.add(road_rp);
					}
				}
			}
		}
		return neighboringRoads;
	}

	public GRBVar getIsSink() {
		return isSink;
	}

	public void setIsSink(GRBVar isSink) {
		this.isSink = isSink;
	}

	public LinkedList<GRBVar> getInFlow() {
		return inFlow;
	}

	public LinkedList<GRBVar> getOutFlow() {
		return outFlow;
	}
	
	public void removeFlowVariables() {
		isSink = null;
		inFlow.clear();
		outFlow.clear();
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
	
	/**
	 * Returns true iff this road and the road rOther share a node.
	 */
	public boolean isAdjacentTo(Road rOther) {
		for(RoadPoint rp : vertices) {
			for(RoadPoint rpOther : rOther.vertices) {
				if(rp == rpOther) {
					return true;
				}
			}
		}
		return false;
	}

	public int getCompID() {
		return compID;
	}

	public void setCompID(int compID) {
		this.compID = compID;
	}

	public boolean isVisited() {
		return visited;
	}

	public void setVisited(boolean visited) {
		this.visited = visited;
	}

	public LinkedList<Road> getComponent() {
		return component;
	}

	public void setComponent(LinkedList<Road> component) {
		this.component = component;
	}
}
