package RoadDisplacement;
import java.util.LinkedList;
import java.util.List;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
//import java.util.Collections;
import java.util.Vector;
import java.util.Iterator;

import org.jgrapht.graph.SimpleWeightedGraph;

import com.vividsolutions.jump.feature.*;
import com.vividsolutions.jump.io.*;

import DelTri.OwnConformingDelaunayTriangulationBuilder;
import diewald_shapeFile.files.dbf.DBF_File;
import diewald_shapeFile.files.shp.SHP_File;
import diewald_shapeFile.files.shp.shapeTypes.ShpPolyLine;
import diewald_shapeFile.files.shx.SHX_File;
import diewald_shapeFile.shapeFile.ShapeFile;
import general.Segment;
import main.Main;

import com.vividsolutions.jts.geom.*;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.MultiPoint;
import org.locationtech.jts.geom.MultiLineString;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.impl.CoordinateArraySequence;
import org.locationtech.jts.triangulate.ConformingDelaunayTriangulationBuilder;
import org.locationtech.jts.index.kdtree.*;
import org.locationtech.jts.geom.Envelope;

/**
 * A class which represents a network of roads
 */
public class RoadNetwork {
        
    LinkedList<Road> roads;
    double minX;
    double minY;
    double maxX;
    double maxY;
    
    ArrayList<RoadPoint> roadPoints;  // contains points of class RoadPoint ordered by their IDs
    KdTree kdTree;   // kd-tree, only containing the original nodes of the road network
    KdTree kdTreeDelPoints;  // kd-tree, only containing the Steiner nodes created during the Delaunay triangulation
    static int counterIDs;   // is increased when a new point of class RoadPoint is created => point obtains an ID
        
    double minWeight;  // the lowest weight among the weights of the roads
    
    
    public RoadNetwork() {
        roads = new LinkedList<Road>();
        minX = Double.POSITIVE_INFINITY;
        minY = Double.POSITIVE_INFINITY;
        maxX = Double.NEGATIVE_INFINITY;
        maxY = Double.NEGATIVE_INFINITY;
        kdTree = new KdTree(Main.eps);
        kdTreeDelPoints = new KdTree(Main.eps);
        roadPoints = new ArrayList<RoadPoint>();
        counterIDs = 0;
        minWeight = Double.MAX_VALUE;
    }
    

    public boolean add(Road r) {
    	if(r.getWeight() < minWeight) minWeight = r.getWeight();
    	for(RoadPoint rp : r.getVertices()) {
    		KdNode query = kdTree.query(rp.getCoordinate());
    		if(query==null) {   // no point has been inserted yet for these coordinates
    			if (rp.getX() < minX) minX = rp.getX();
    			if (rp.getY() < minY) minY = rp.getY();
    			if (rp.getX() > maxX) maxX = rp.getX();
    			if (rp.getY() > maxY) maxY = rp.getY();
    			kdTree.insert(rp.getCoordinate(),rp);
    			roadPoints.add(rp);
    			rp.setId(counterIDs);  // set the IDs according to the order by which the points are inserted into the kd-tree
    			counterIDs++;
    		}
    		else {     // a point already exists with these coordinates
    			RoadPoint rpOld = ((RoadPoint)query.getData());
    			if(rpOld!=rp) {
    				r.replaceRoadPoint(rp,rpOld);
        			rpOld.addRoad(r);
    			}
    		}
    	}
    	return roads.add(r);
    }
    
    public LinkedList<Road> getRoads() {
        return roads;
    }
      
    /**
     * Reads a road network from a shp file. The weights of the roads are defines as their lengths in m.
     * @param filename the name of the shp file (without .shp)
     * @param filepath the path to the directory containing the file
     */
    public static RoadNetwork importFromShapefile(String filepath, String filename) {
    	RoadNetwork myRoadNetwork = new RoadNetwork();

    	try {
    		DBF_File.LOG_INFO           = false;
    		DBF_File.LOG_ONLOAD_HEADER  = false;
    		DBF_File.LOG_ONLOAD_CONTENT = false;

    		SHX_File.LOG_INFO           = false;
    		SHX_File.LOG_ONLOAD_HEADER  = false;
    		SHX_File.LOG_ONLOAD_CONTENT = false;

    		SHP_File.LOG_INFO           = false;
    		SHP_File.LOG_ONLOAD_HEADER  = false;
    		SHP_File.LOG_ONLOAD_CONTENT = false;

    		ShapeFile shapefile = new ShapeFile(filepath, filename).READ();
    		int number_of_shapes = shapefile.getSHP_shapeCount();
    		for(int i = 0; i < number_of_shapes; i++){
    			ShpPolyLine shape    = shapefile.getSHP_shape(i);
    			int number_of_parts = shape.getNumberOfParts();
    			if(number_of_parts==1) {
    				int number_of_vertices = shape.getNumberOfPoints();
    				double[][] vertices = shape.getPoints();
    				Road myRoad = new Road(new Vector<RoadPoint>(), "");
    				myRoadNetwork.getRoads().add(myRoad);
    				double length = 0;
    				for(int j = 0; j < number_of_vertices; j++){
    					// check if there is already a point at this location and if not, create a new one
    					KdNode query = myRoadNetwork.getKdTree().query(new Coordinate(vertices[j][0], vertices[j][1]));
    					if(query==null) {   // no point inserted yet for these coordinates
    						RoadPoint p = new RoadPoint(vertices[j][0], vertices[j][1]);
    						if (p.getX() < myRoadNetwork.getMinX()) myRoadNetwork.minX = p.getX();
    						if (p.getY() < myRoadNetwork.getMinY()) myRoadNetwork.minY = p.getY();
    						if (p.getX() > myRoadNetwork.getMaxX()) myRoadNetwork.maxX = p.getX();
    						if (p.getY() > myRoadNetwork.getMaxY()) myRoadNetwork.maxY = p.getY();
    						myRoad.getVertices().add(p);
    						p.addRoad(myRoad);
    						myRoadNetwork.getKdTree().insert(p.getCoordinate(),p);
    						myRoadNetwork.getRoadPoints().add(p);
    						p.setId(counterIDs);  // define the IDs according to the order by which the points are inserted to the kd-tree
    						counterIDs++;
    						if(j>0) {
    							RoadPoint previous_point = myRoad.getVertices().get(j-1);
    							length = length + previous_point.getDist(p);
    						}
    					}
    					else {     // a point with these coordinates exists already
    						RoadPoint rpOld = ((RoadPoint)query.getData());
    						myRoad.getVertices().add(rpOld);
    						rpOld.addRoad(myRoad);
    						if(j>0) {
    							RoadPoint previous_point = myRoad.getVertices().get(j-1);
    							length = length + previous_point.getDist(rpOld);
    						}
    					}
    				}
    				myRoad.setWeight(length);
    				if(myRoad.getWeight() < myRoadNetwork.getMinWeight()) {
    					myRoadNetwork.setMinWeight(myRoad.getWeight());
    				}

    			}
    			else {
    				System.out.println("The shapefile contains a MultiLineString!");
    			}
    		}      
    	} catch(Exception ex) {
    		System.out.println("shp_read: " + ex);
    	}
    	return myRoadNetwork;
    }

    public double getMinX() {
    	return minX;
    }

    public double getMinY() {
    	return minY;
    }

    public double getMaxX() {
    	return maxX;
    }

    public double getMaxY() {
    	return maxY;
    }    
    
    /**
     * Returns a graph of the nodes and edges of the road network.
     */
    public SimpleWeightedGraph<RoadPoint, Segment> getGraph() {
        SimpleWeightedGraph<RoadPoint, Segment> graph = new SimpleWeightedGraph<RoadPoint, Segment>(Segment.class);
        
        for(Road road : this.roads) {
            for(RoadPoint RoadPoint : road.getVertices()) {
                graph.addVertex(RoadPoint);
            }
            for(int i = 0; i < road.getVertices().size() - 1; i++) {
//            	if(road.getVertices().get(i)==road.getVertices().get(i + 1) ) {
//                	road.getVertices().remove(i+1);
//                	i--;
//                	continue;
//                }
                if (!graph.containsEdge(road.getVertices().get(i), road.getVertices().get(i + 1))) {
                    Segment edge = new Segment(road.getVertices().get(i), road.getVertices().get(i + 1));
                    graph.addEdge(road.getVertices().get(i), road.getVertices().get(i + 1), edge);
                    if (edge != null)
                        graph.setEdgeWeight(edge, Math.hypot(road.getVertices().get(i + 1).getX() - road.getVertices().get(i).getX() , road.getVertices().get(i + 1).getY() - road.getVertices().get(i).getY()  ));             
                }                
            }
        }
        return graph;
    }
    
    
    /**
     * Returns the closest road point from an input location.
     */
	public RoadPoint getRoadPointAt(Coordinate location) {
		double radius = 2;
		int MAX_SEARCH_POWER = 10;
		for (int i = 0; i <= MAX_SEARCH_POWER; ++i) {
			RoadPoint rp = getRoadPointAt(location, Math.pow(radius, i));
			if (rp != null)
				return rp;
		}
		System.out.println("No vertex found in radius " + Math.pow(radius, MAX_SEARCH_POWER) + " around " + location);
		return null;
	}

	public RoadPoint getRoadPointAt(Coordinate location, double radius) {
		Envelope queryEnv = new Envelope(location.x - radius, location.x + radius, location.y - radius,
				location.y + radius);
		NodeCollector collector = new NodeCollector(location, radius);
		kdTree.query(queryEnv, collector);

		Entry<Double, KdNode> nearestNode = collector.getNearest();
		if (nearestNode == null)
			return null;

		RoadPoint rp = (RoadPoint) nearestNode.getValue().getData();
		return rp;
	}
    
    
    /**
     * Returns all road points as an object of class MultiPoint (important for the Delaunay 
     * triangulation).
     */
    public MultiPoint getRoadPointsAsMultiPoint(){
    	GeometryFactory gf = new GeometryFactory();
    	org.locationtech.jts.geom.Point[] points = new org.locationtech.jts.geom.Point[roadPoints.size()];
    	Object[] roadPointsArr = roadPoints.toArray();
    	for(int i=0; i < points.length; i++) {
    		RoadPoint rp = (RoadPoint) roadPointsArr[i];
    		Coordinate c = rp.getCoordinate();  // new Coordinate(rp.getX(),rp.getY());
    		Coordinate[] cs = {c};
    		org.locationtech.jts.geom.CoordinateSequence cas = new CoordinateArraySequence(cs);
    		org.locationtech.jts.geom.Point point = new org.locationtech.jts.geom.Point(cas,gf);
    		points[i] = point;
    	}
    	return new MultiPoint(points,gf);
    }
    
    /**
     * Returns all segments as an object of class MultiLineString (important for the Delaunay 
     * triangulation).
     */
    public MultiLineString getSegmentsAsMultiLineString(){
    	GeometryFactory gf = new GeometryFactory();
    	ArrayList<LineString> lineStrings = new ArrayList<LineString>();
    	for(Road road : roads) {
    		for(int i = 0; i < road.getVertices().size() - 1; i++) {
    			RoadPoint currentRoadPoint = road.getVertices().get(i);
    			RoadPoint nextRoadPoint = road.getVertices().get(i+1);
    			Coordinate currentCoord = currentRoadPoint.getCoordinate();
    			Coordinate nextCoord = nextRoadPoint.getCoordinate();
    			Coordinate[] cs = {currentCoord, nextCoord};
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

	public KdTree getKdTree() {
		return kdTree;
	}

	public KdTree getKdTreeDelPoints() {
		return kdTreeDelPoints;
	}

	public ArrayList<RoadPoint> getRoadPoints() {
		return roadPoints;
	}

	public static int getCounterIDs() {
		return counterIDs;
	}

	public static void setCounterIDs(int counterIDs) {
		RoadNetwork.counterIDs = counterIDs;
	}

	public double getMinWeight() {
		return minWeight;
	}

	public void setMinWeight(double minWeight) {
		this.minWeight = minWeight;
	}
	
	
	/**
	 * Returns a graph which contains a node for each road and an edge for each pair of 
	 * adjacent roads.
	 */
	public SimpleWeightedGraph<Road,RoadPair> getGraphOfRoads(){
		SimpleWeightedGraph<Road, RoadPair> graphOfRoads = new SimpleWeightedGraph<Road, RoadPair>(RoadPair.class);
		for(Road r : roads) {
			if(!graphOfRoads.containsVertex(r)) {
				graphOfRoads.addVertex(r);
			}
			for(Road rOther : r.getNeighboringRoads()) {
				RoadPair roadPair = new RoadPair(r,rOther);
				if(!graphOfRoads.containsVertex(rOther)) {
					graphOfRoads.addVertex(rOther);
				}
				if(!graphOfRoads.containsEdge(r,rOther)) {
					graphOfRoads.addEdge(r,rOther,roadPair);
				}
			}
		}
		return graphOfRoads;
	}
	
	
	private static class NodeCollector implements KdNodeVisitor {

		private Coordinate origin;
		private double radius;
		private TreeMap<Double, KdNode> foundNodes;

		public NodeCollector(Coordinate origin, double radius) {
			super();
			this.origin = origin;
			this.radius = radius;
			this.foundNodes = new TreeMap<>();
		}

		public Entry<Double, KdNode> getNearest() {
			return foundNodes.firstEntry();
		}

		@Override
		public void visit(KdNode node) {
			double dist = node.getCoordinate().distance(origin);
			if (dist < radius)
				foundNodes.put(dist, node);
		}
	}
	
	public void setKdTree(KdTree newTree) {
		this.kdTree = newTree;
	}


	public void setMinX(double minX) {
		this.minX = minX;
	}


	public void setMinY(double minY) {
		this.minY = minY;
	}


	public void setMaxX(double maxX) {
		this.maxX = maxX;
	}


	public void setMaxY(double maxY) {
		this.maxY = maxY;
	}
}
