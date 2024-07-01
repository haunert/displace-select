package general;

import java.util.LinkedList;

import org.jgrapht.graph.DefaultWeightedEdge;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.impl.CoordinateArraySequence;

import com.vividsolutions.jts.geom.Envelope;

import BuildingDisplacement.Building;
import BuildingDisplacement.BuildingPoint;
import RoadDisplacement.Road;
import RoadDisplacement.RoadPoint;
import gurobi.GRBVar;

/**
 * A line segment represented by start and end point.
 */

public class Segment extends DefaultWeightedEdge {

	private static final long serialVersionUID = 1L;
	private Point start;
	private Point end;
	private GRBVar deltaX;   // Variable for the distortion in x
	private GRBVar deltaY;   // Variable for the distortion in y
	
	boolean bottleneck_edge = false;  // is this edge a bottleneck edge?
	
	Collector c;  // for each edge that belongs to two buildings
		
	public Segment(Point s, Point e) {
		start = s;
		end = e;
		deltaX = null;
		deltaY = null;
	}
	
	public Point getStart() {
		return start;
	}
	
	public Point getEnd() {
		return end;
	}
	
	public Envelope getEnvelope() {
		double xMin = start.getX();
        double xMax = end.getX();
        if (xMin > xMax) {
            xMin = xMax;
            xMax = start.getX();
        }
        double yMin = start.getY();
        double yMax = end.getY();
        if (yMin > yMax) {
            yMin = yMax;
            yMax = start.getY();
        }
        return new Envelope(xMin, xMax, yMin, yMax);
	}
	
	public String toString() {
		return "(" + start.getX() + " " + start.getY() + ") -> (" + end.getX() + " " + end.getY() + ")";
	}
	
	public String toStringWithOriginalCoordinates() {
		return "(" + start.getOriginalCoord().getX() + " " + 
					start.getOriginalCoord().getY() + ") -> (" + 
					end.getOriginalCoord().getX() + " " + 
					end.getOriginalCoord().getY() + ")";
	}

	public GRBVar getDeltaX() {
		return deltaX;
	}

	public void setDeltaX(GRBVar deltaX) {
		this.deltaX = deltaX;
	}

	public GRBVar getDeltaY() {
		return deltaY;
	}

	public void setDeltaY(GRBVar deltaY) {
		this.deltaY = deltaY;
	}
	
	public Point getOtherPoint(Point p) {
		if(p == start) {
			return end;
		}
		if(p == end) {
			return start;
		}
		System.out.println("Error in method getOtherPoint in class Segment!");
		return null;
	}
	
	public double getLength() {
		return start.getDist(end);
	}

	
	/**
	 * Checks if this edge is visible in the solution graph, i.e. if it belongs to a selected
	 * road or building.
	 */
	public boolean isVisible() {
		if(start.getClass().toString().equals("class RoadDisplacement.RoadPoint")
			&& end.getClass().toString().equals("class RoadDisplacement.RoadPoint")) {
			Road road = RoadPoint.getCommonRoad((RoadPoint) start, (RoadPoint) end);
			if(road==null) return false;
			return road.isSelected();
		}
		if(start.getClass().toString().equals("class BuildingDisplacement.BuildingPoint")
				&& end.getClass().toString().equals("class BuildingDisplacement.BuildingPoint")) {
				LinkedList<Building> buildings = BuildingPoint.getCommonBuildings((BuildingPoint) start, 
																(BuildingPoint) end);
				for(Building b : buildings) {
					if(b.isSelected()) {
						return true;
					}
				}
				return false;
		}
		return false;
	}

	public boolean isBottleneckEdge() {
		return bottleneck_edge;
	}

	public void setBottleneckEdge(boolean bottleneck_edge) {
		this.bottleneck_edge = bottleneck_edge;
	}

	public Collector getCollector() {
		return c;
	}

	public void setCollector(Collector c) {
		this.c = c;
	}
	
	
	/**
	 * Auxiliary method for the class SolverWareJones (for a bottleneck edge). Returns 
	 * true if this edge has as endpoints a road point and a building point and false otherwise.
	 */
	public boolean connectsBuildingWithRoad() {
		if((start.getClass().toString().equals("class RoadDisplacement.RoadPoint")
			&& end.getClass().toString().equals("class BuildingDisplacement.BuildingPoint"))
		|| (start.getClass().toString().equals("class BuildingDisplacement.BuildingPoint")
			 && end.getClass().toString().equals("class RoadDisplacement.RoadPoint"))) {
			return true;
		}
		return false;
	}
}
