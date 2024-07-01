package general;

import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import org.jgrapht.graph.SimpleWeightedGraph;
import org.locationtech.jts.geom.Coordinate;

import RoadDisplacement.RoadPoint;
import BuildingDisplacement.BuildingPoint;
import gurobi.GRBVar;

public class Point implements Comparable<Object> {
    double x;
    double y;
    int id;
    
    // auxiliary attribute used in the method insertPointsIntoBuildingPolyLine/insertPointsIntoRoadPolyLine in class General Network
    boolean replaced = false;
    BuildingPoint bp;  // only used in the special case that this is a Steiner point lying on the boundary of a building
    RoadPoint rp;  // only used in the special case that this is a Steiner point lying on a road
    
	GRBVar dx;   // displacement in x
	GRBVar dy;   // displacement in y
	double dxDouble=0;   // displacement in x as double (in addition to dx, since Gurobi sometimes seems to "forget" the value of dx!)
	double dyDouble=0;   // displacement in y as double (in addition to dy, since Gurobi sometimes seems to "forget" the value of dy!)
	Coordinate originalCoord;  // the original (input) coordinates of the point
	Coordinate intermediateCoord;  // an inetermediate coordinate of the point (only relevant for the visualization of the bottleneck edges if several iterations are applied)
    
	// Only relevant for a Steiner point: originalBefore and originalAfter are the start and end point 
	// of the segment of the input geometry into which this point has been inserted (meaning they are 
	// not Steiner points).
	Point originalBefore; 
	Point originalAfter;
	
	Collector co;
	boolean incidentToBottleneckEdge=false;
	
	LinkedList<Segment> edgesInDelaunayGraph;  // incident edges in the Delaunay graph

	
    public Point(double myX, double myY) {
        x = myX;
        y = myY;  
        this.id = -1;   // to specify that in fact no id has been given to this point
        originalCoord = new Coordinate(myX,myY);
    }
    
    public Point(double myX, double myY, int id) {
        x = myX;
        y = myY;  
        this.id = id;   
        originalCoord = new Coordinate(myX,myY);
    }
    
    public Point(Coordinate xy) {
        x = xy.getX();
        y = xy.getY();  
        this.id = -1;   // to specify that in fact no id has been given to this point
        originalCoord = xy;
    }
    
    public double getX() {
        return x;
    }
    
    public double getY() {
        return y;
    }
    public void setX(double myX) {
        x = myX;
    }
    
    public void setY(double myY) {
        y = myY;
    }
    
    /**
     * lexicographical order
     */
    public int compareTo(Object o) {
        if(o instanceof Point) {
                Point arg0 = (Point) o;
                if(x > arg0.x || (x == arg0.x && y > arg0.y)) return 1;
                if(x == arg0.x && y == arg0.y) return 0;
                return -1; 
        }
        return 0;
    }
    
    public int hashCode() {
        return Integer.parseInt("" + (int) (x / 1000.0) + (int) (y / 1000.0));
    }

    public boolean equals(Object o) {
        Point p0 = (Point) o;
        return (x == p0.x && y == p0.y);
    }
    
    public String toString() {
    	return "" + this.x + ", " + this.y;
    }
    
    
   	public double getDist(Point p) {
        double dx = p.x - this.x;
        double dy = p.y - this.y; 
        return Math.sqrt(dx * dx + dy * dy);
    }
   	
	/**
	 * Computes the distance from this point to the line going to the points p_i and p_j.
	 */
	public double getDistanceTo(Point p_i, Point p_j) {
		double d_ij = p_i.getDist(p_j);
		double rx = (p_j.getX()-p_i.getX())/d_ij;
		double ry = (p_j.getY()-p_i.getY())/d_ij;
		double dx = this.getX() - p_i.getX();
		double dy = this.getY() - p_i.getY();
		double skp = dx*rx+dy*ry;
		return Math.sqrt((dx-skp*rx)*(dx-skp*rx)+(dy-skp*ry)*(dy-skp*ry));
	}
	
	/**
	 * Computes the distance from this point to the segment which has the points p_i and p_j
	 * as endpoints.
	 */
	public double getDistanceToSegment(Point p_i, Point p_j) {
		double d_ij = p_i.getDist(p_j);
		double rx = (p_j.getX()-p_i.getX())/d_ij;
		double ry = (p_j.getY()-p_i.getY())/d_ij;
		double dx = this.getX() - p_i.getX();
		double dy = this.getY() - p_i.getY();
		double skp = dx*rx+dy*ry;
		if(skp<0) {
			return this.getDist(p_i);
		}
		if(skp>1) {
			return this.getDist(p_j);
		}
		return Math.sqrt((dx-skp*rx)*(dx-skp*rx)+(dy-skp*ry)*(dy-skp*ry));
	}
	
	public Coordinate getCoordinate() {
		return new Coordinate(x,y);
	}

	public boolean isReplaced() {
		return replaced;
	}

	public void setReplaced(boolean replaced) {
		this.replaced = replaced;
	}

	public GRBVar getDx() {
		return dx;
	}

	public void setDx(GRBVar dx) {
		this.dx = dx;
	}

	public GRBVar getDy() {
		return dy;
	}

	public void setDy(GRBVar dy) {
		this.dy = dy;
	}

	public double getDxDouble() {
		return dxDouble;
	}

	public void setDxDouble(double dxDouble) {
		this.dxDouble = dxDouble;
	}

	public double getDyDouble() {
		return dyDouble;
	}

	public void setDyDouble(double dyDouble) {
		this.dyDouble = dyDouble;
	}

	public int getId() {
		return id;
	}

	public void setId(int id) {
		this.id = id;
	}

	public Coordinate getOriginalCoord() {
		return originalCoord;
	}

	public void setOriginalCoord(Coordinate originalCoord) {
		this.originalCoord = originalCoord;
	}
	
	/**
	 * Checks if this point is visible in the map, meaning if at least one of the 
	 * objects to which it belongs is selected.
	 */
	public boolean isVisible() {
		if(this.getClass().toString().equals("class RoadDisplacement.RoadPoint")) {
			return ((RoadPoint) this).isVisibleRoadPoint();
		}
		if(this.getClass().toString().equals("class BuildingDisplacement.BuildingPoint")) {
			return ((BuildingPoint) this).isVisibleBuildingPoint();
		}
		System.out.println("Error in method isVisible(), class Point!");
		return false;
	}

	public Point getOriginalBefore() {
		return originalBefore;
	}

	public void setOriginalBefore(Point originalBefore) {
		this.originalBefore = originalBefore;
	}

	public Point getOriginalAfter() {
		return originalAfter;
	}

	public void setOriginalAfter(Point originalAfter) {
		this.originalAfter = originalAfter;
	}
	
	
	/**
	 * This point and the point p must be connected via a sequence of segments, which all together 
	 * lie on a line in the input. Then this method finds the endpoints of the segment containing the point 
	 * which is at distance d from this point (referring to the original coordinates).
	 */
	public List<Point> findPointsInDirectionOf(Point p, double d,
			             SimpleWeightedGraph<Point,Segment> g) {
		double dist_to_p = this.originalCoord.distance(p.getOriginalCoord());
		if(d>dist_to_p) {
			System.out.println("Error 1 in findPointsInDirectionOf() in class Point!");
			return null;
		}
		Point pBefore = this;
		Point pAfter = this;
		while(pAfter.getOriginalCoord().distance(this.originalCoord) < d) {
			pBefore = pAfter;
			Set<Segment> edges = g.edgesOf(pAfter);
			boolean next_point_found = false;
			for(Segment se : edges) {
				Point pOther = se.getOtherPoint(pAfter);
				double dist_this_pOther = this.originalCoord.distance(pOther.getOriginalCoord());
				double dist_pOther_p = pOther.getOriginalCoord().distance(p.getOriginalCoord());
				if(Math.abs(dist_to_p-dist_this_pOther-dist_pOther_p)<Math.pow(10, -10)
				       && dist_pOther_p < pAfter.getOriginalCoord().distance(p.getOriginalCoord())) {
					pAfter = pOther;
					next_point_found = true;
					break;
				}
			}
			if(!next_point_found) {
				System.out.println("Error 2 in findPointsInDirectionOf() in class Point!");
				return null;
			}
		}
		List<Point> points = new LinkedList<Point>();
		points.add(pBefore);
		points.add(pAfter);
		return points;
	}

	public LinkedList<Segment> getEdgesInDelaunayGraph() {
		return edgesInDelaunayGraph;
	}

	public void setEdgesInDelaunayGraph(LinkedList<Segment> edges) {
		this.edgesInDelaunayGraph = edges;
	}

	public Collector getCollector() {
		return co;
	}

	public void setCollector(Collector co) {
		this.co = co;
	}

	public boolean isIncidentToBottleneckEdge() {
		return incidentToBottleneckEdge;
	}

	public void setIncidentToBottleneckEdge(boolean incidentToBottleneckEdge) {
		this.incidentToBottleneckEdge = incidentToBottleneckEdge;
	}

	public BuildingPoint getBuildingPoint() {
		return bp;
	}

	public void setBuildingPoint(BuildingPoint bp) {
		this.bp = bp;
	}
	
	public RoadPoint getRoadPoint() {
		return rp;
	}

	public void setRoadPoint(RoadPoint rp) {
		this.rp = rp;
	}

	public boolean isDelaunay() {
		if(this.getClass().toString().equals("class RoadDisplacement.RoadPoint")) {
			return ((RoadPoint) this).isDelaunay();
		}
		else return ((BuildingPoint) this).isDelaunay();
	}
	
	public boolean checkIfIsIncidentToBottleneckEdge(SimpleWeightedGraph<Point,Segment> g) {
		if(! g.containsVertex(this)) {
			incidentToBottleneckEdge = false;
			return false;
		}
		for(Segment se : g.edgesOf(this)) {
			if(se.isBottleneckEdge()) {
				incidentToBottleneckEdge = true;
				return true;
			}
		}
		incidentToBottleneckEdge = false;
		return false;
	}

	public Coordinate getIntermediateCoord() {
		return intermediateCoord;
	}

	public void setIntermediateCoord(Coordinate intermediateCoord) {
		this.intermediateCoord = intermediateCoord;
	}
}
