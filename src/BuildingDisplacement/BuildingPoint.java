package BuildingDisplacement;

import java.util.ArrayList;
import java.util.LinkedList;

import org.locationtech.jts.geom.Coordinate;

import RoadDisplacement.Road;
import RoadDisplacement.RoadPoint;
import general.Point;
import gurobi.GRBVar;

public class BuildingPoint extends Point{
	
	private ArrayList<Building> buildings;   // the building(s) to which this point belongs
	private boolean isDelaunay;      // is true if this point has been added as a Steiner point, and false otherwise

	public BuildingPoint(double x, double y) {
		super(x,y);
		buildings = new ArrayList<Building>();
	}
	
	public BuildingPoint(double x, double y, int id) {
		super(x,y,id);
		buildings = new ArrayList<Building>();
	}
	
	public BuildingPoint(double x, double y, int id, boolean isDelaunay) {
		super(x,y,id);
		buildings = new ArrayList<Building>();
		this.isDelaunay = isDelaunay;       
	}
	
	
	public void addBuilding(Building b) {
		boolean isNewBuilding = true;   // first check if the point already has a reference to the new building
		for(Building bOld : buildings) {
			if(b==bOld) {
				isNewBuilding = false;
				break;
			}
		}
		if(isNewBuilding) {
			buildings.add(b);
		}
	}

	public ArrayList<Building> getBuildings() {
		return buildings;
	}
	
	/**
	 * Returns a list with the buildings to which both bp1 and bp2 belong. These can be 0, 1 or 
	 * 2 buildings (the last case occurs if both bp1 and bp2 lie on the boundary between these buildings).
	 */
	public static LinkedList<Building> getCommonBuildings(BuildingPoint bp1, BuildingPoint bp2) {
		LinkedList<Building> commonBuildings = new LinkedList<Building>();
		for(int i=0; i < bp1.getBuildings().size(); i++) {
			Building b_i = bp1.getBuildings().get(i);
			for(int j=0; j<bp2.getBuildings().size(); j++) {
				Building b_j = bp2.getBuildings().get(j);
				if (b_i==b_j) {
					commonBuildings.add(b_i);
				}
			}
		}
		return commonBuildings;
	}

	public boolean isDelaunay() {
		return isDelaunay;
	}

	public void setDelaunay(boolean isDelaunay) {
		this.isDelaunay = isDelaunay;
	}
	
	/**
	 * Checks if this point is visible in the output, i.e. if at least one of the 
	 * buildings to which it belongs has been selected.
	 */
	public boolean isVisibleBuildingPoint() {
		for(Building b : buildings) {
			if(b.isSelected()) {
				return true;
			}
		}
		return false;
	}

}
