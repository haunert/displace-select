package GeneralDisplacement;

import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;
import java.util.Vector;

import org.jgrapht.graph.SimpleWeightedGraph;
import org.locationtech.jts.geom.*;
import org.locationtech.jts.index.kdtree.KdNode;
import org.locationtech.jts.index.kdtree.KdTree;

import BuildingDisplacement.*;
import RoadDisplacement.*;
import general.*;
import general.Point;
import main.Main;
import DelTri.*;

/**
 * A class representing the geometry of the roads and buildings.
 */
public class GeneralNetwork {

	private RoadNetwork rn;
	private BuildingNetwork bn;
	
	private KdTree kdTreeDelPoints;
	
	
	private SimpleWeightedGraph<Point,Segment> T;  		// graph of the Conforming Delaunay Triangulation
	
	// Bounding Box
	private double minX;
    private double minY;
    private double maxX;
    private double maxY;
	
	public GeneralNetwork(RoadNetwork rn, BuildingNetwork bn) {
		this.rn = rn;
		this.bn = bn;
		kdTreeDelPoints = new KdTree(Main.eps);
		
		T = new SimpleWeightedGraph<Point,Segment>(Segment.class);
		
		if(rn.getRoads().size()==0) {
			minX = bn.getMinX();
			maxX = bn.getMaxX();
			minY = bn.getMinY();
			maxY = bn.getMaxY();
		}
		else if(bn.getBuildings().size()==0) {
			minX = rn.getMinX();
			maxX = rn.getMaxX();
			minY = rn.getMinY();
			maxY = rn.getMaxY();
		}
		else {
			minX = Math.min(bn.getMinX(),rn.getMinX());
			minY = Math.min(bn.getMinY(),rn.getMinY());
			maxX = Math.max(bn.getMaxX(),rn.getMaxX());
			maxY = Math.max(bn.getMaxY(),rn.getMaxY());
		}
	}
	
	
	/**
	 * Determines the big constant M for the constraints, with the help of the map extent.
	 */
	public double determineM() {
		return 2*Math.max(maxX-minX, maxY-minY);
	}
	
	/**
	 * Returns a graph of the geometry of the roads and buildings. The nodes belonging 
	 * to buildings are of the class BuildingPoint and the nodes belonging to roads are of the class
	 * RoadPoint.
	 */
	public SimpleWeightedGraph<Point,Segment> getGraph(){
		SimpleWeightedGraph<Point, Segment> graph = 
				new SimpleWeightedGraph<Point, Segment>(Segment.class);
		SimpleWeightedGraph<RoadPoint, Segment> roadGraph = rn.getGraph();
		SimpleWeightedGraph<BuildingPoint, Segment> buildingGraph = bn.getGraph();
		for(RoadPoint rp : roadGraph.vertexSet()) {
			graph.addVertex(rp);
		}
		for(BuildingPoint bp : buildingGraph.vertexSet()) {
			graph.addVertex(bp);
		}
		for(Segment se : roadGraph.edgeSet()) {
			RoadPoint start = (RoadPoint) se.getStart();
			RoadPoint end = (RoadPoint) se.getEnd();
			graph.addEdge(start, end, se);
			graph.setEdgeWeight(se,start.getDist(end));
		}
		for(Segment se : buildingGraph.edgeSet()) {
			BuildingPoint start = (BuildingPoint) se.getStart();
			BuildingPoint end = (BuildingPoint) se.getEnd();
			graph.addEdge(start, end, se);
			graph.setEdgeWeight(se,start.getDist(end));
		}
		return graph;
	}
	
	
	/**
	 * Computes the graph T as a Conforming Delaunay Triangulation of the geometry of the roads and
	 * buildings. The newly created Steiner points are associated with the roads/buildings they 
	 * belong to geometrically.
	 */
	public void computeDelaunayGraph(){
		
		System.out.println("Compute Delaunay triangulation ...");
		
		OwnConformingDelaunayTriangulationBuilder cdtb = new OwnConformingDelaunayTriangulationBuilder();
		cdtb.setSites(this.getPointsAsMultiPoint());
    	cdtb.setConstraints(this.getSegmentsAsMultiLineString());
    	cdtb.setTolerance(1e-2);
    	
    	GeometryFactory gf = new GeometryFactory();
    	Geometry tri = cdtb.getTriangles(gf);
    	
    	// insert the original points into the Delaunay graph T
    	for(BuildingPoint bp : bn.getBuildingPoints()) {
    		T.addVertex(bp);
    	}
    	for(Road road : rn.getRoads()) {
    		for(RoadPoint rp : road.getVertices()) {
    			T.addVertex(rp);
            }
    	}
    	
    	// iterate over the triangles and insert the triangle edges and the Steiner points into the graph.
    	for(int i=0; i < tri.getNumGeometries(); i++) {
    		
    		Geometry triangle = tri.getGeometryN(i);
    		Coordinate[] coordinates = triangle.getCoordinates();
    		if (coordinates.length != 4) {   // the first and last point are the same => 4 points per triangle
    			System.out.println("Error in method getDelaunayGraph() in class GeneralNetwork, a triangle does not have 3 corners!");
    		}
    		else {
    			// Search for the points corresponding to the coordinates of this triangle
    			Point[] pts = new Point[3];
    			for(int j=0; j<3; j++) {
    				KdNode queryNode = rn.getKdTree().query(coordinates[j]);
    				if(queryNode != null) {  // the point already belongs to the road network
    					pts[j] = (RoadPoint) queryNode.getData();
    				}
    				else {
    					queryNode = bn.getKdTree().query(coordinates[j]);
    					if(queryNode != null) {  // the point already belongs to the geometry of the buildings
        					pts[j] = (BuildingPoint) queryNode.getData();
        				}
    					else {  // the point is a Steiner point
    						queryNode = kdTreeDelPoints.query(coordinates[j]);  // check if the point already exists in the kd-tree of the Steiner points
    						if(queryNode != null) {
        						pts[j] = (Point) queryNode.getData();
        					}
        					else {
        						pts[j] = new Point(coordinates[j].getX(),coordinates[j].getY());
        						T.addVertex(pts[j]);        						
            					kdTreeDelPoints.insert(pts[j].getCoordinate(),pts[j]);
        					}
        				} 
    				}
    			}
    			if(pts[0]==null || pts[1]==null || pts[2]==null) {
    				System.out.println("Fehler in getDelaunayGraph() in Klasse GeneralNetwork, mindestens ein Dreieckspunkt ist null!");
    			}
    			Segment edge = new Segment(pts[0], pts[1]);
    			T.addEdge(pts[0], pts[1], edge);
    			T.setEdgeWeight(edge, Math.hypot(pts[0].getX()-pts[1].getX(), pts[0].getY()-pts[1].getY())); 	             
    			edge = new Segment(pts[0], pts[2]);
    			T.addEdge(pts[0], pts[2], edge);
    			T.setEdgeWeight(edge, Math.hypot(pts[0].getX()-pts[2].getX(), pts[0].getY()-pts[2].getY())); 	             
    			edge = new Segment(pts[2], pts[1]);
    			T.addEdge(pts[2], pts[1], edge);
    			T.setEdgeWeight(edge, Math.hypot(pts[2].getX()-pts[1].getX(), pts[2].getY()-pts[1].getY()));  	             
    		}
    	}

        // Associate the Steiner points with the buildings/roads they geometrically belong to
    	for(Building b : bn.getBuildings()) {
    		insertPointsIntoBuildingPolyLine(b.getShell(), b);
    		for(Vector<BuildingPoint> hole : b.getHoles()) {
    			insertPointsIntoBuildingPolyLine(hole, b);
    		}
    	}
    	for(Road r : rn.getRoads()) {
    		insertPointsIntoRoadPolyLine(r);
    	}
    	
    	System.out.println("Delaunay triangulation computed.");
	}
	
	
    /**
     * Finds for a building boundary (interior or exterior) the Steiner points created during 
     * the Delaunay triangulation which lie on this boundary. The Steiner points are added 
     * to this boundary. 
     * @param polyline either one of the interior rings of the building b or the exterior ring
     * @param b the building
     */
    public void insertPointsIntoBuildingPolyLine(Vector<BuildingPoint> polyline, Building b) {
    	for(int i=0; i<polyline.size(); i++) {
    		BuildingPoint p_i = polyline.get(i);
    		BuildingPoint p_iplus1 = null;
    		if(i < polyline.size()-1) {
    			p_iplus1 = polyline.get(i+1);
    		}
    		else {  // the boundary is closed, i.e. the second point of the last edge is the first of the first edge
    			p_iplus1 = polyline.get(0);
    		}
			Envelope queryEnv = new Envelope(p_i.getCoordinate(),p_iplus1.getCoordinate());
			List delPoints = kdTreeDelPoints.query(queryEnv);
			TreeMap<Double,Point> delPointsOrdered = new TreeMap<Double,Point>();
			for(Object obj : delPoints) {  // sort the Steiner points according to their distance from the point p_i ...
				Point p = (Point) ((KdNode) obj).getData();
				if(p.getDistanceTo(p_i,p_iplus1)<Main.eps) {  // ... but only if they are lying on the segment p_i - p_iplus1 
					delPointsOrdered.put(p.getDist(p_i),p);
					//p.setReplaced(true);
				}
			}
			Collection<Point> delPointsColl = delPointsOrdered.values();
			int counter = 1;
			for(Point p : delPointsColl) {  // iterate over the Steiner points according to their distances from p_i
				BuildingPoint bp;
				if(!p.isReplaced()) {  // Case 1: the node has been discovered for the first time.
					// Delete p from the graph and insert a new point of the class BuildingPoint at its place.
					Set<Segment> edges = T.edgesOf(p);
					
					LinkedList<Point> neighbors = new LinkedList<Point>();
					for(Segment se : edges) {
						neighbors.add(se.getOtherPoint(p));
					}
					
					T.removeVertex(p);
					
					bp = new BuildingPoint(p.getX(),p.getY(),BuildingNetwork.getCounterIDs(),true);
					bp.setOriginalBefore(p_i);
					bp.setOriginalAfter(p_iplus1);
					T.addVertex(bp);
					BuildingNetwork.setCounterIDs(BuildingNetwork.getCounterIDs()+1);
					for(Point neighbor : neighbors) {
						Segment seNew = new Segment(bp,neighbor);
						
						T.addEdge(bp,neighbor,seNew);
						T.setEdgeWeight(seNew, Math.hypot(bp.getX()-neighbor.getX(), bp.getY()-neighbor.getY()));
					}
					p.setReplaced(true);
					p.setBuildingPoint(bp);
					bn.getBuildingPoints().add(bp);
				}
				else {  // Case 2: The node already belongs to the boundary of an adjacent building.
					bp = p.getBuildingPoint();
				}
				bp.addBuilding(b);
				if(i+counter < polyline.size()) {
					polyline.add(i+counter,bp);
				}
				else {  // Inserting a point on the segment p_N - p_0
					polyline.add(bp);
				}
				counter++;
			}
			i = i + counter -1;  // index has to be increased since new points have been inserted
    	}
    }
    
    
    /**
     * For a road r, finds the Steiner points created during the Delaunay triangulation which lie on
     * this road. The points are added to the road.
     * @param r the road
     */
    public void insertPointsIntoRoadPolyLine(Road r) {
    	Vector<RoadPoint> polyline = r.getVertices();
    	for(int i=0; i<polyline.size()-1; i++) {
			RoadPoint p_i = polyline.get(i);
			RoadPoint p_iplus1 = polyline.get(i+1);
			Envelope queryEnv = new Envelope(p_i.getCoordinate(),p_iplus1.getCoordinate());
			List delPoints = kdTreeDelPoints.query(queryEnv);
			TreeMap<Double,Point> delPointsOrdered = new TreeMap<Double,Point>();
			for(Object obj : delPoints) {  // sort the Steiner points according to their distance from the point p_i ...
				Point p = (Point) ((KdNode) obj).getData();
				if(p.getDistanceTo(p_i,p_iplus1)<Main.eps) {  // ... but only if they are lying on the segment p_i - p_iplus1
					delPointsOrdered.put(p.getDist(p_i),p);
					//p.setReplaced(true);
				}
			}
			Collection<Point> delPointsColl = delPointsOrdered.values();
			int counter = 1;
			for(Point p : delPointsColl) {  // iterate over the Steiner points according to their distances from p_i
				RoadPoint rp;
				if(!p.isReplaced()) { // Case 1: the node has been discovered for the first time.
					// Delete p from the graph and insert a new point of the class RoadPoint at its place.
					Set<Segment> edges = T.edgesOf(p);
					
					LinkedList<Point> neighbors = new LinkedList<Point>();
					for(Segment se : edges) {
						neighbors.add(se.getOtherPoint(p));
					}
					T.removeVertex(p);
					rp = new RoadPoint(p.getX(),p.getY(),RoadNetwork.getCounterIDs(),true);
					rp.setOriginalBefore(p_i);
					rp.setOriginalAfter(p_iplus1);
					T.addVertex(rp);
					RoadNetwork.setCounterIDs(RoadNetwork.getCounterIDs()+1);
					for(Point neighbor : neighbors) {
						Segment seNew = new Segment(rp,neighbor);
						
						T.addEdge(rp,neighbor,seNew);
						T.setEdgeWeight(seNew, Math.hypot(rp.getX()-neighbor.getX(), rp.getY()-neighbor.getY()));
					}
					p.setReplaced(true);
					p.setRoadPoint(rp);
				}
				else {  // Case 2: The node already belongs to the boundary of a different road (crossing).
					rp = p.getRoadPoint();
				}
				rp.addRoad(r);
				polyline.add(i+counter,rp);
				counter++;
			}
			i = i + counter -1;  // the index has to be increased since new points have been inserted
    	}
    }
    
//    /**
//     * Finds for a road the Steiner points created during the Delaunay triangulation which lie on
//     * this road. The points are added to the road.
//     * @param r the road
//     */
//    public void insertPointsIntoRoadPolyLine(Road r) {
//    	Vector<RoadPoint> polyline = r.getVertices();
//    	for(int i=0; i<polyline.size()-1; i++) {
//			RoadPoint p_i = polyline.get(i);
//			RoadPoint p_iplus1 = polyline.get(i+1);
//			Envelope queryEnv = new Envelope(p_i.getCoordinate(),p_iplus1.getCoordinate());
//			List delPoints = kdTreeDelPoints.query(queryEnv);
//			TreeMap<Double,Point> delPointsOrdered = new TreeMap<Double,Point>();
//			for(Object obj : delPoints) {  // sort the Steiner points according to their distance from the point p_i ...
//				Point p = (Point) ((KdNode) obj).getData();
//				if(p.getDistanceTo(p_i,p_iplus1)<Main.eps && !p.isReplaced()) {  // ... but only if they are lying on the segment p_i - p_iplus1 and have not been associated with another road yet
//					delPointsOrdered.put(p.getDist(p_i),p);
//					p.setReplaced(true);
//				}
//			}
//			Collection<Point> delPointsColl = delPointsOrdered.values();
//			int counter = 1;
//			for(Point p : delPointsColl) {  // iterate over the Steiner points according to their distances from p_i
//				// Delete p from the graph and insert a new point of the class RoadPoint at its place.
//				Set<Segment> edges = T.edgesOf(p);
//				
//				LinkedList<Point> neighbors = new LinkedList<Point>();
//				for(Segment se : edges) {
//					neighbors.add(se.getOtherPoint(p));
//				}
//				T.removeVertex(p);
//				RoadPoint rp = new RoadPoint(p.getX(),p.getY(),RoadNetwork.getCounterIDs(),true);
//				rp.setOriginalBefore(p_i);
//				rp.setOriginalAfter(p_iplus1);
//				T.addVertex(rp);
//				RoadNetwork.setCounterIDs(RoadNetwork.getCounterIDs()+1);
//				for(Point neighbor : neighbors) {
//					Segment seNew = new Segment(rp,neighbor);
//					
//					T.addEdge(rp,neighbor,seNew);
//					T.setEdgeWeight(seNew, Math.hypot(rp.getX()-neighbor.getX(), rp.getY()-neighbor.getY()));
//				}
//				rp.addRoad(r);
//				polyline.add(i+counter,rp);
//				counter++;
//			}
//			i = i + counter -1;  // the index has to be increased since new points have been inserted
//    	}
//    }
	
	public MultiPoint getPointsAsMultiPoint() {
		MultiPoint rpts = rn.getRoadPointsAsMultiPoint();
		MultiPoint bpts = bn.getBuildingPointsAsMultiPoint();
		return (MultiPoint) rpts.union(bpts);
	}
	
	public MultiLineString getSegmentsAsMultiLineString() {
		MultiLineString mls_rn = rn.getSegmentsAsMultiLineString();
		MultiLineString mls_bn = bn.getSegmentsAsMultiLineString();
		return (MultiLineString) mls_rn.union(mls_bn);
	}

	public SimpleWeightedGraph<Point, Segment> getDelaunayGraph() {
		return T;
	}
	
	
	
	/**
	 * Normalizes the weights, so that those of the roads and those of the buildings are comparable.
	 * @param road_factor the factor defining the ratio between the weights of the buildings and
	 * those of the roads (the smallest building obtains the weight 1 and the shortest road obtains the weight
	 * road_factor)
	 */
	public void normalizeWeights(double road_factor) {
		for(Building b : bn.getBuildings()) {
			b.setWeight(b.getWeight()/bn.getMinWeight());
		}
		bn.setMinWeight(1);
		for(Road r : rn.getRoads()) {
			r.setWeight(road_factor*r.getWeight()/rn.getMinWeight());
		}
		rn.setMinWeight(road_factor);
	}
	
	
	/**
	 * Returns the road being the closest road to the centroid of the building b.
	 */
	public Road getClosestRoadOfBuilding(Building b) {
		Polygon poly = b.getPolygon();
		org.locationtech.jts.geom.Coordinate centroid_coord = poly.getCentroid().getCoordinate();
		Point centroid = new Point(centroid_coord.getX(),centroid_coord.getY());
		RoadPoint rp_closest = rn.getRoadPointAt(centroid_coord);
		ArrayList<Road> closest_roads = rp_closest.getRoads();
		if(closest_roads.size()==1) {
			return closest_roads.get(0);
		}
		double minDiff = Double.MAX_VALUE;
		Road rClosest = null;
		for(Road r : closest_roads) {
			for(int i=0; i<r.getVertices().size()-1;i++) {
				double d_se = centroid.getDistanceToSegment(r.getVertices().get(i),
																r.getVertices().get(i+1));
				if(d_se < minDiff) {
					rClosest = r;
					minDiff = d_se;
				}
			}
		}
		return rClosest;
	}
	
	
	/**
	 * Removes the Steiner points from the roads and buildings and clears the kd-tree kdTreeDelPoints.
	 */
	public void removeDelaunayPoints() {
		if(kdTreeDelPoints.size()>0) {
			for(Building b : bn.getBuildings()) {
				b.removeDelaunayPoints();
			}
			for(Road r : rn.getRoads()) {
				r.removeDelaunayPoints();
			}
			kdTreeDelPoints = new KdTree(Main.eps);
		}
	}
	
	/**
	 * Removes the Steiner points from the roads and buildings and clears the kd-tree kdTreeDelPoints.
	 */
	public void removeDelaunayPointsWithoutChecking() {
		for(Building b : bn.getBuildings()) {
			b.removeDelaunayPoints();
		}
		for(Road r : rn.getRoads()) {
			r.removeDelaunayPoints();
		}
		kdTreeDelPoints = new KdTree(Main.eps);
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

	public RoadNetwork getRn() {
		return rn;
	}

	public BuildingNetwork getBn() {
		return bn;
	}


	public void setKdTreeDelPoints(KdTree kdTreeDelPoints) {
		this.kdTreeDelPoints = kdTreeDelPoints;
	}
}
