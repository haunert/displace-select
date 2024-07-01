package BuildingDisplacement;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.TreeMap;
import java.util.Vector;

import org.jgrapht.graph.SimpleWeightedGraph;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.MultiLineString;
import org.locationtech.jts.geom.MultiPoint;
import org.locationtech.jts.geom.impl.CoordinateArraySequence;
import org.locationtech.jts.index.kdtree.KdNode;
import org.locationtech.jts.index.kdtree.KdTree;

import com.vividsolutions.jump.feature.Feature;
import com.vividsolutions.jump.feature.FeatureCollection;
import com.vividsolutions.jump.io.DriverProperties;
import com.vividsolutions.jump.io.ShapefileReader;

import DelTri.OwnConformingDelaunayTriangulationBuilder;
import RoadDisplacement.Road;
import RoadDisplacement.RoadNetwork;
import RoadDisplacement.RoadPair;
import RoadDisplacement.RoadPoint;
import general.Point;
import general.Segment;
import main.Main;

/**
 * A class which collects the geometrical information about the buildings
 */
public class BuildingNetwork {
	
	LinkedList<Building> buildings;
    double minX;
    double minY;
    double maxX;
    double maxY;
    KdTree kdTree;   // kd-tree, only containing the original points
    KdTree kdTreeDelPoints;  // kd-tree, only containing the points inserted during the Delaunay triangulation
    static int counterIDs;   // is increased when a new point of class BuildingPoint is created => point gets an ID
    
    // Contains the building points after removing duplicate points. Also contains the Steiner
    // points added during the Delaunay triangulation.
    LinkedList<BuildingPoint> buildingPoints;  
    
    double minWeight;  // the smallest weight of any of the buildings
    
	public BuildingNetwork() {
		buildings = new LinkedList<Building>();
        minX = Double.POSITIVE_INFINITY;
        minY = Double.POSITIVE_INFINITY;
        maxX = Double.NEGATIVE_INFINITY;
        maxY = Double.NEGATIVE_INFINITY;
        counterIDs = 0;
        kdTree = new KdTree(Main.eps);
        kdTreeDelPoints = new KdTree(Main.eps);
        buildingPoints = new LinkedList<BuildingPoint>();
        minWeight = Double.MAX_VALUE;
	}
	
	
	public boolean add(Building b) {
    	if(b.getWeight() < minWeight)  minWeight = b.getWeight();

    	// 1.) Exterior boundary
    	for(BuildingPoint bp : b.getShell()) {  // iterate over all corners
    		KdNode query = kdTree.query(bp.getCoordinate());  // first check if a building point with these coordinates has already been inserted
    		if(query==null) {   // no building point with these coordinates has been inserted yet
				kdTree.insert(bp.getCoordinate(),bp);
				buildingPoints.add(bp);
				bp.setId(counterIDs);  // define the IDs according to the order by which the points are inserted to the kd-tree
				counterIDs++;
				if (bp.getX() < minX) minX = bp.getX();
        		if (bp.getY() < minY) minY = bp.getY();
        		if (bp.getX() > maxX) maxX = bp.getX();
        		if (bp.getY() > maxY) maxY = bp.getY();
			}
    		else {   // a point already exists there (e.g. at the boundary between two buildings)
    			BuildingPoint bpOld = ((BuildingPoint)query.getData());
    			if(bpOld!=bp) {
    				b.replacePointInRing(bp, bpOld, b.getShell());
    				bpOld.addBuilding(b);
    			}
    		}
    	}
    	
    	// 2.) Interior boundaries
    	for(int j=0; j<b.getHoles().size(); j++) {
    		Vector<BuildingPoint> hole_j = b.getHoles().get(j);
    		for(BuildingPoint bp : hole_j) {  // iterate over all corners
        		KdNode query = kdTree.query(bp.getCoordinate());  
        		if(query==null) {   
    				kdTree.insert(bp.getCoordinate(),bp);
    				buildingPoints.add(bp);
    				bp.setId(counterIDs);  
    				counterIDs++;
    				if (bp.getX() < minX) minX = bp.getX();
            		if (bp.getY() < minY) minY = bp.getY();
            		if (bp.getX() > maxX) maxX = bp.getX();
            		if (bp.getY() > maxY) maxY = bp.getY();
    			}
        		else {   
        			BuildingPoint bpOld = ((BuildingPoint)query.getData());
        			if(bpOld!=bp) {
        				b.replacePointInRing(bp, bpOld, hole_j);
        				bpOld.addBuilding(b);
        			}
        		}
        	}
    	}
		return buildings.add(b);
	}
	
	/**
	 * Reads buildings from a shp file. The weights of the buildings are defined as their areas in m^2.
	 */
	public static BuildingNetwork importFromShapefile(String filename) {
        BuildingNetwork myBuildingNetwork = new BuildingNetwork();
        
        if ( filename.endsWith(".shp") ) {
            ShapefileReader shp_input = new ShapefileReader();
            DriverProperties dp = new DriverProperties(filename);
            FeatureCollection myFeatureCollection = null;
            
            try {
                myFeatureCollection = shp_input.read(dp);
                

                // iterate over all buildings
                for (Iterator<?> i = myFeatureCollection.iterator(); i.hasNext(); ) {
            
                  Feature myFeature = (Feature) i.next();
                  String type = "";
                  if (myFeature.getSchema().hasAttribute("type")) {
                     type = (String) myFeature.getAttribute("type");
                  }

                   com.vividsolutions.jts.geom.Geometry myGeometry = myFeature.getGeometry();
                    
                    if (myGeometry.getGeometryType() == "Polygon") {
                    	Building myBuilding = new Building(new Vector<BuildingPoint>(), null, type);
                    	if (myFeature.getSchema().hasAttribute("selection")) {
                    		String selection = (String) myFeature.getAttribute("selection");
                            myBuilding.setSelectionWareJones(selection.replaceAll("\\s+",""));
                         }
                    	myBuilding.setWeight(myGeometry.getArea());
                    	if(myBuilding.getWeight() < myBuildingNetwork.getMinWeight()) {
                    		myBuildingNetwork.setMinWeight(myBuilding.getWeight());
                    	}
                    	myBuildingNetwork.getBuildings().add(myBuilding);
                    	
                    	// 1.) Exterior boundary
                    	com.vividsolutions.jts.geom.LineString shell = 
                    			((com.vividsolutions.jts.geom.Polygon) myGeometry).getExteriorRing();
                    	com.vividsolutions.jts.geom.Coordinate[] xyz = shell.getCoordinates();
                    	for(int j = 0; j < xyz.length-1; j++) {  // iterate over all corners
                    		BuildingPoint bp = new BuildingPoint(xyz[j].x, xyz[j].y);
                    		KdNode query = myBuildingNetwork.getKdTree().query(bp.getCoordinate());  // first check if a building point with these coordinates has already been inserted to the kd-tree
                    		if(query==null) {   // no point with these coordinates has been inserted yet
                				myBuildingNetwork.getKdTree().insert(bp.getCoordinate(),bp);
                				myBuildingNetwork.getBuildingPoints().add(bp);
                				bp.setId(counterIDs);  // define the IDs according to the order by which the building points are inserted to the kd-tree
                				counterIDs++;
                				myBuilding.getShell().add(bp);
                				bp.addBuilding(myBuilding);
                				if (bp.getX() < myBuildingNetwork.getMinX()) myBuildingNetwork.minX = bp.getX();
                        		if (bp.getY() < myBuildingNetwork.getMinY()) myBuildingNetwork.minY = bp.getY();
                        		if (bp.getX() > myBuildingNetwork.getMaxX()) myBuildingNetwork.maxX = bp.getX();
                        		if (bp.getY() > myBuildingNetwork.getMaxY()) myBuildingNetwork.maxY = bp.getY();
                			}
                    		else {   // a point already exists there (e.g. at the boundary between two adjacent buildings)
                    			BuildingPoint bpOld = ((BuildingPoint)query.getData());
                    			myBuilding.getShell().add(bpOld);
                				bpOld.addBuilding(myBuilding);
                    		}
                    	}
                    	
                    	// 2.) Interior boundaries
                    	int numHoles = ((com.vividsolutions.jts.geom.Polygon) myGeometry).getNumInteriorRing();
                    	Vector<Vector<BuildingPoint>> holes = new Vector<Vector<BuildingPoint>>();
                    	for(int j=0; j<numHoles; j++) {
                    		com.vividsolutions.jts.geom.LineString hole = 
                        			((com.vividsolutions.jts.geom.Polygon) myGeometry).getInteriorRingN(j);
                    		com.vividsolutions.jts.geom.Coordinate[] xyz_j = hole.getCoordinates();
                    		Vector<BuildingPoint> hole_j = new Vector<BuildingPoint>();
                    		for(int k = 0; k < xyz_j.length-1; k++) {  // iterate over all corners
                    			BuildingPoint bp = new BuildingPoint(xyz_j[k].x, xyz_j[k].y);
                        		KdNode query = myBuildingNetwork.getKdTree().query(bp.getCoordinate());  
                        		if(query==null) {   
                    				myBuildingNetwork.getKdTree().insert(bp.getCoordinate(),bp);
                    				myBuildingNetwork.getBuildingPoints().add(bp);
                    				bp.setId(counterIDs);  
                    				counterIDs++;
                    				hole_j.add(bp);
                    				bp.addBuilding(myBuilding);
                    			}
                        		else {   
                        			BuildingPoint bpOld = ((BuildingPoint)query.getData());
                        			hole_j.add(bpOld);
                    				bpOld.addBuilding(myBuilding);
                        		}
                        	}
                    		holes.add(hole_j);
                    	}
                    	myBuilding.setHoles(holes);
                    }
                    else {
                        System.out.println("The shapefile does not contain a polygon but a " + myGeometry.getGeometryType());                                 
                    }
                }
                    
            } catch(Exception ex) {
            	System.out.println("shp_read: " + ex);
            }
        }
        return myBuildingNetwork;
    }
	
	
	/**
     * Returns a graph of the geometry defined by the building edges.
     */
    public SimpleWeightedGraph<BuildingPoint, Segment> getGraph() {
        SimpleWeightedGraph<BuildingPoint, Segment> graph = new SimpleWeightedGraph<BuildingPoint, Segment>(Segment.class);

        for(BuildingPoint bp : buildingPoints) {
        	graph.addVertex(bp);
        }

        for(Building building : this.buildings) {
        	
            // Add the edges on the exterior boundary
        	Vector<BuildingPoint> shell = building.getShell();
            for(int i = 0; i < shell.size() - 1; i++) {
//            	if(shell.get(i)==shell.get(i + 1) ) {
//                	shell.remove(i+1);
//                	i--;
//                	continue;
//                }
                if (!graph.containsEdge(shell.get(i), shell.get(i + 1))) {
                    Segment edge = new Segment(shell.get(i), shell.get(i + 1));
                    graph.addEdge(shell.get(i), shell.get(i + 1), edge);
                    if (edge != null)
                        graph.setEdgeWeight(edge, Math.hypot(shell.get(i+1).getX() - shell.get(i).getX() , shell.get(i+1).getY() - shell.get(i).getY()  ));             
                }
            }
            if (!graph.containsEdge(shell.get(0), shell.get(shell.size()-1))) {
                Segment edge = new Segment(shell.get(0), shell.get(shell.size()-1));
                graph.addEdge(shell.get(0), shell.get(shell.size()-1), edge);
                if (edge != null)
                    graph.setEdgeWeight(edge, Math.hypot(shell.get(shell.size()-1).getX() - shell.get(0).getX() , shell.get(shell.size()-1).getY() - shell.get(0).getY()  ));             
            }
            
            // Add the edges on the interior boundary
            Vector<Vector<BuildingPoint>> holes = building.getHoles();
            for(int i=0; i<holes.size(); i++) {
            	Vector<BuildingPoint> hole_i = holes.get(i);
            	for(int j = 0; j < hole_i.size() - 1; j++) {
//            		if(hole_i.get(i)==hole_i.get(i + 1) ) {
//                    	hole_i.remove(i+1);
//                    	i--;
//                    	continue;
//                    }
                    if (!graph.containsEdge(hole_i.get(j), hole_i.get(j + 1)) ) {
                        Segment edge = new Segment(hole_i.get(j), hole_i.get(j + 1));
                        graph.addEdge(hole_i.get(j), hole_i.get(j + 1), edge);
                        if (edge != null)
                            graph.setEdgeWeight(edge, Math.hypot(hole_i.get(j+1).getX() - hole_i.get(j).getX() , hole_i.get(j+1).getY() - hole_i.get(j).getY()  ));             
                    }
                }
                if (!graph.containsEdge(hole_i.get(0), hole_i.get(hole_i.size()-1))) {
                    Segment edge = new Segment(hole_i.get(0), hole_i.get(hole_i.size()-1));
                    graph.addEdge(hole_i.get(0), hole_i.get(hole_i.size()-1), edge);
                    if (edge != null)
                        graph.setEdgeWeight(edge, Math.hypot(hole_i.get(hole_i.size()-1).getX() - hole_i.get(0).getX() , hole_i.get(hole_i.size()-1).getY() - hole_i.get(0).getY()  ));             
                }
            }
        }
        return graph;
    }


    /**
     * Returns the building points as a MultiPoint object (important for the Delaunay triangulation)
     */
    public MultiPoint getBuildingPointsAsMultiPoint(){
    	GeometryFactory gf = new GeometryFactory();
    	org.locationtech.jts.geom.Point[] points = new org.locationtech.jts.geom.Point[buildingPoints.size()];
    	Object[] buildingPointsArr = buildingPoints.toArray();
    	for(int i=0; i < points.length; i++) {
    		BuildingPoint bp = (BuildingPoint) buildingPointsArr[i];
    		Coordinate c = bp.getCoordinate();  
    		Coordinate[] cs = {c};
    		org.locationtech.jts.geom.CoordinateSequence cas = new CoordinateArraySequence(cs);
    		org.locationtech.jts.geom.Point point = new org.locationtech.jts.geom.Point(cas,gf);
    		points[i] = point;
    	}
    	return new MultiPoint(points,gf);
    }
    
    
    /**
     * Returns the segments as a MultiLineString object (important for the Delaunay triangulation)
     */
    public MultiLineString getSegmentsAsMultiLineString(){
    	GeometryFactory gf = new GeometryFactory();
    	ArrayList<LineString> lineStrings = new ArrayList<LineString>();
    	for(Building building : buildings) {
    		LinkedList<Segment> segments = building.getAllSegments();
    		for(Segment se : segments) {
    			Point start = se.getStart();
    			Point end = se.getEnd();
    			Coordinate startCoord = new Coordinate(start.getX(),start.getY());
    			Coordinate endCoord = new Coordinate(end.getX(),end.getY());
    			Coordinate[] cs = {startCoord, endCoord};
    			org.locationtech.jts.geom.CoordinateSequence cas = new CoordinateArraySequence(cs);
    			lineStrings.add(new LineString(cas,gf));
    		}
    	}
    	LineString[] lineStringsArray = new LineString[lineStrings.size()];
    	
    	// Copy the LineStrings from the ArrayList into the normal Array 
    	for(int i=0; i<lineStrings.size(); i++) {  
    		lineStringsArray[i] = lineStrings.get(i);
    	}
    	MultiLineString mls = new MultiLineString(lineStringsArray,gf);
    	return mls;
    }


	public double getMinX() {
		return minX;
	}


	public void setMinX(double minX) {
		this.minX = minX;
	}


	public double getMinY() {
		return minY;
	}


	public void setMinY(double minY) {
		this.minY = minY;
	}


	public double getMaxX() {
		return maxX;
	}


	public void setMaxX(double maxX) {
		this.maxX = maxX;
	}


	public double getMaxY() {
		return maxY;
	}


	public void setMaxY(double maxY) {
		this.maxY = maxY;
	}


	public LinkedList<Building> getBuildings() {
		return buildings;
	}


	public KdTree getKdTree() {
		return kdTree;
	}


	public LinkedList<BuildingPoint> getBuildingPoints() {
		return buildingPoints;
	}


	public static int getCounterIDs() {
		return counterIDs;
	}


	public static void setCounterIDs(int counterIDs) {
		BuildingNetwork.counterIDs = counterIDs;
	}


	public double getMinWeight() {
		return minWeight;
	}


	public void setMinWeight(double minWeight) {
		this.minWeight = minWeight;
	}
	
	
	/**
	 * Returns a graph which contains a node for each building and an edge for each pair of adjacent 
	 * buildings.
	 */
	public BuildingGraph getGraphOfBuildings(){
		BuildingGraph graphOfBuildings = new BuildingGraph();
		for(Building b : buildings) {
			if(!graphOfBuildings.containsVertex(b)) {
				graphOfBuildings.addVertex(b);
			}
			for(Building bOther : b.findNeighboringBuildings()) {
				BuildingPair buildingPair = new BuildingPair(b,bOther);
				if(!graphOfBuildings.containsVertex(bOther)) {
					graphOfBuildings.addVertex(bOther);
				}
				if(!graphOfBuildings.containsEdge(b,bOther)) {
					graphOfBuildings.addEdge(buildingPair);
				}
			}
		}
		graphOfBuildings.depthFirstSearch();
		return graphOfBuildings;
	}
	
	public void setKdTree(KdTree newTree) {
		this.kdTree = newTree;
	}

}
