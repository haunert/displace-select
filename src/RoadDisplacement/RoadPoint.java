package RoadDisplacement;
import java.util.ArrayList;

import org.locationtech.jts.geom.Coordinate;

import general.Point;
import gurobi.GRBVar;

public class RoadPoint extends Point{
	
	private ArrayList<Road> roads;   // the roads to which this point belongs
	private boolean isCrossing;      // does the point lie on a crossing?
	private boolean isDelaunay;      // is this point a Steiner point?

	public RoadPoint(double x, double y) {
		super(x,y);
		roads = new ArrayList<Road>();
		isCrossing = false;
	}
	
	public RoadPoint(double x, double y, int id) {
		super(x,y,id);
		roads = new ArrayList<Road>();
		isCrossing = false;
	}
	
	public RoadPoint(double x, double y, int id, boolean isDelaunay) {
		super(x,y,id);
		roads = new ArrayList<Road>();
		isCrossing = false;
		this.isDelaunay = isDelaunay;
	}
	
	public void addRoad(Road r) {
		boolean isNewRoad = true;   // first check if this point already contains a reference to the road r
		for(Road rOld : roads) {
			if(r==rOld) {
				isNewRoad = false;
				break;
			}
		}
		if(isNewRoad) {
			roads.add(r);
			if (roads.size() > 1) {
				isCrossing = true;
			}
		}
	}
	
	public void addRoads(ArrayList<Road> otherRoads) {
		for(Road r : otherRoads) {
			addRoad(r);
		}
//		for(Road r : otherRoads) {
//			roads.add(r);	
//		}
//		if (roads.size() > 1) {
//			isCrossing = true;
//		}
	}

	public ArrayList<Road> getRoads() {
		return roads;
	}
	

	public boolean isCrossing() {
		return isCrossing;
	}
	
	/**
	 * Returns the (unique) road to which both rp1 and rp2 belong.
	 */
	public static Road getCommonRoad(RoadPoint rp1, RoadPoint rp2) {
		for(int i=0; i < rp1.getRoads().size(); i++) {
			Road r_i = rp1.getRoads().get(i);
			for(int j=0; j<rp2.getRoads().size(); j++) {
				Road r_j = rp2.getRoads().get(j);
				if (r_i==r_j) {
					return r_i;
				}
			}
		}
		return null;
	}

	public boolean isDelaunay() {
		return isDelaunay;
	}

	public void setDelaunay(boolean isDelaunay) {
		this.isDelaunay = isDelaunay;
	}
	
	
	/**
	 * Checks if this point is visible in the output, i.e. if at least one of the roads to which it belongs 
	 * has been selected.
	 */
	public boolean isVisibleRoadPoint() {
		for(Road r : roads) {
			if(r.isSelected()) {
				return true;
			}
		}
		return false;
	}

}
