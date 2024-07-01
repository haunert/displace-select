package general;

import gurobi.GRBVar;

import java.util.LinkedList;

import BuildingDisplacement.*;
import RoadDisplacement.*;

/**
 * A class that groups objects together. This can happen in two cases:
 * 1) A node which is incident to a bottleneck edge belongs to >1 objects.
 * 2) An edge belongs to two buildings.
 */
public class Collector {

	Point p;
	Segment se;
	GRBVar z;

	LinkedList<Building> buildings;
	LinkedList<Road> roads;
	
	public Collector(Point p) {
		this.p = p;
		if(p.getClass().toString().equals("class BuildingDisplacement.BuildingPoint")) {
			buildings = new LinkedList<Building>();
			buildings.addAll(((BuildingPoint) p).getBuildings());
		}
		else if(p.getClass().toString().equals("class RoadDisplacement.RoadPoint")){
			roads = new LinkedList<Road>();
			roads.addAll(((RoadPoint) p).getRoads());
		}
		else {
			System.out.println("Error in class Collector, Constructor!!!");
		}
	}
	
	public Collector(Segment se) {
		this.se = se;
		if(se.getStart().getClass().toString().equals("class BuildingDisplacement.BuildingPoint")
		    && se.getEnd().getClass().toString().equals("class BuildingDisplacement.BuildingPoint")) {
			BuildingPoint bp1 = (BuildingPoint) se.getStart();
			BuildingPoint bp2 = (BuildingPoint) se.getEnd();
			buildings = new LinkedList<Building>();
			buildings.addAll(BuildingPoint.getCommonBuildings(bp1, bp2));
		}
		else {
			System.out.println("Error in class Collector, Constructor!!!");
		}
	}

	public GRBVar getZ() {
		return z;
	}

	public void setZ(GRBVar z) {
		this.z = z;
	}

	public Point getP() {
		return p;
	}

	public LinkedList<Building> getBuildings() {
		return buildings;
	}

	public LinkedList<Road> getRoads() {
		return roads;
	}
	
	

}
