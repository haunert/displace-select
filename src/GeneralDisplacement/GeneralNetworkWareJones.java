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

public class GeneralNetworkWareJones {

	private RoadNetwork rn;
	private BuildingNetwork bn;
	
	private KdTree kdTreeDelPoints;
	
	
	private SimpleWeightedGraph<Point,Segment> delGraph;
	
	// Bounding Box
	private double minX;
    private double minY;
    private double maxX;
    private double maxY;
	
	public GeneralNetworkWareJones(RoadNetwork rn, BuildingNetwork bn) {
		this.rn = rn;
		this.bn = bn;
		kdTreeDelPoints = new KdTree(Main.eps);
		
		delGraph = new SimpleWeightedGraph<Point,Segment>(Segment.class);
		
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
	 * Returns a graph of the geometry of the roads and buildings. The nodes belonging 
	 * to buildings are of the class BuildingPoint and the nodes belonging to roads are of the class
	 * RoadPoint.
	 */
	public SimpleWeightedGraph<Point,Segment> getGraph(){

		SimpleWeightedGraph<Point, Segment> graph = 
				new SimpleWeightedGraph<Point, Segment>(Segment.class);
		for(Road road : rn.getRoads()) {
			// Insert the nodes of this road 
			for(RoadPoint RoadPoint : road.getVertices()) {
				graph.addVertex(RoadPoint);
			}
			// Insert the segments of this road
			for(int i = 0; i < road.getVertices().size() - 1; i++) {
				Segment edge = new Segment(road.getVertices().get(i), road.getVertices().get(i + 1));
				graph.addEdge(road.getVertices().get(i), road.getVertices().get(i + 1), edge);
				if (edge != null) {
					graph.setEdgeWeight(edge, Math.hypot(road.getVertices().get(i + 1).getX() - road.getVertices().get(i).getX() , road.getVertices().get(i + 1).getY() - road.getVertices().get(i).getY()  ));             
				}
			}
		}
		for(BuildingPoint bp : bn.getBuildingPoints()) {
			graph.addVertex(bp);
		}
		for(Building building : bn.getBuildings()) {
			// Insert the exterior edges of this building
			Vector<BuildingPoint> shell = building.getShell();
			for(int i = 0; i < shell.size() - 1; i++) {
				Segment edge = new Segment(shell.get(i), shell.get(i + 1));
				graph.addEdge(shell.get(i), shell.get(i + 1), edge);
				if (edge != null) {
					graph.setEdgeWeight(edge, Math.hypot(shell.get(i+1).getX() - shell.get(i).getX() , shell.get(i+1).getY() - shell.get(i).getY()  ));             
				}
			}
			Segment edge = new Segment(shell.get(0), shell.get(shell.size()-1));
			graph.addEdge(shell.get(0), shell.get(shell.size()-1), edge);
			if (edge != null) {
				graph.setEdgeWeight(edge, Math.hypot(shell.get(shell.size()-1).getX() - shell.get(0).getX() , shell.get(shell.size()-1).getY() - shell.get(0).getY()  ));             
			}
			// Insert the interior edges of this building
			Vector<Vector<BuildingPoint>> holes = building.getHoles();
			for(int i=0; i<holes.size(); i++) {
				Vector<BuildingPoint> hole_i = holes.get(i);
				for(int j = 0; j < hole_i.size() - 1; j++) {
					Segment edge_hole = new Segment(hole_i.get(j), hole_i.get(j + 1));
					graph.addEdge(hole_i.get(j), hole_i.get(j + 1), edge_hole);
					if (edge_hole != null) {
						graph.setEdgeWeight(edge_hole, Math.hypot(hole_i.get(j+1).getX() - hole_i.get(j).getX() , hole_i.get(j+1).getY() - hole_i.get(j).getY()  ));             
					}
				}
				Segment edge_last = new Segment(hole_i.get(0), hole_i.get(hole_i.size()-1));
				graph.addEdge(hole_i.get(0), hole_i.get(hole_i.size()-1), edge_last);
				if (edge_last != null) {
					graph.setEdgeWeight(edge_last, Math.hypot(hole_i.get(hole_i.size()-1).getX() - hole_i.get(0).getX() , hole_i.get(hole_i.size()-1).getY() - hole_i.get(0).getY()  ));             
				}
			}
		}
		return graph;
	}
	
	
	public double determineM() {
		return 2*Math.max(maxX-minX, maxY-minY);
	}
	
	
	public void computeDelaunayGraph(){
		
		OwnConformingDelaunayTriangulationBuilder cdtb = new OwnConformingDelaunayTriangulationBuilder();
		cdtb.setSites(this.getPointsAsMultiPoint());
    	cdtb.setConstraints(this.getSegmentsAsMultiLineString());
    	cdtb.setTolerance(1e-2);
    	
    	GeometryFactory gf = new GeometryFactory();
    	Geometry tri = cdtb.getTriangles(gf);
    	
    	// insert the original points of the roads and buildings into the graph
    	for(BuildingPoint bp : bn.getBuildingPoints()) {
    		delGraph.addVertex(bp);
    	}
    	for(Road road : rn.getRoads()) {
    		for(RoadPoint rp : road.getVertices()) {
    			delGraph.addVertex(rp);
            }
    	}
    	
    	// iterate over the triangles and add insert the triangle edges and the new Steiner points to the graph.
    	for(int i=0; i < tri.getNumGeometries(); i++) {
    		
    		Geometry triangle = tri.getGeometryN(i);
    		Coordinate[] coordinates = triangle.getCoordinates();
    		if (coordinates.length != 4) {   // the first and last point are the same => 4 points per triangle
    			System.out.println("Error in method getDelaunayGraph() in class GeneralNetworkWareJones, a triangle does not have 3 corners!");
    		}
    		else {
    			// Find the points belonging to the coordinates of this triangle
    			Point[] pts = new Point[3];
    			for(int j=0; j<3; j++) {
    				KdNode queryNode = rn.getKdTree().query(coordinates[j]);
    				if(queryNode != null) {  // the point already exists in the road network
    					pts[j] = (RoadPoint) queryNode.getData();
    				}
    				else {
    					queryNode = bn.getKdTree().query(coordinates[j]);
    					if(queryNode != null) {  // the point already belongs to the geometry formed by the buildings
        					pts[j] = (BuildingPoint) queryNode.getData();
        				}
    					else {  // the point is a newly created Steiner point
    						queryNode = kdTreeDelPoints.query(coordinates[j]);  // first check if the point already exists in the kd-tree of the Steiner points
    						if(queryNode != null) {
        						pts[j] = (Point) queryNode.getData();
        					}
        					else {
        						pts[j] = new Point(coordinates[j].getX(),coordinates[j].getY());
            					
        						delGraph.addVertex(pts[j]);        						
            					kdTreeDelPoints.insert(pts[j].getCoordinate(),pts[j]);
        					}
        				} 
    				}
    			}
    			if(pts[0]==null || pts[1]==null || pts[2]==null) {
    				System.out.println("Fehler in getDelaunayGraph() in Klasse GeneralNetworkWareJones, mindestens ein Dreieckspunkt ist null!");
    			}
    			Segment edge = new Segment(pts[0], pts[1]);
    			delGraph.addEdge(pts[0], pts[1], edge);
    			delGraph.setEdgeWeight(edge, Math.hypot(pts[0].getX()-pts[1].getX(), pts[0].getY()-pts[1].getY())); 	             
    			edge = new Segment(pts[0], pts[2]);
    			delGraph.addEdge(pts[0], pts[2], edge);
    			delGraph.setEdgeWeight(edge, Math.hypot(pts[0].getX()-pts[2].getX(), pts[0].getY()-pts[2].getY())); 	             
    			edge = new Segment(pts[2], pts[1]);
    			delGraph.addEdge(pts[2], pts[1], edge);
    			delGraph.setEdgeWeight(edge, Math.hypot(pts[2].getX()-pts[1].getX(), pts[2].getY()-pts[1].getY()));  	             
    		}
    	}

        // Associate the new Steiner points with the roads/buildings they belong to geometrically
    	for(Building b : bn.getBuildings()) {
    		insertPointsIntoBuildingPolyLine(b.getShell(), b);
    		for(Vector<BuildingPoint> hole : b.getHoles()) {
    			insertPointsIntoBuildingPolyLine(hole, b);
    		}
    	}
    	for(Road r : rn.getRoads()) {
    		insertPointsIntoRoadPolyLine(r);
    	}
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
				if(p.getDistanceTo(p_i,p_iplus1)<Main.eps && !p.isReplaced()) {  // ... but only if they are lying on the segment p_i - p_iplus1 and have not been associated with a different building yet
					delPointsOrdered.put(p.getDist(p_i),p);
					p.setReplaced(true);
				}
			}
			Collection<Point> delPointsColl = delPointsOrdered.values();
			int counter = 1;
			for(Point p : delPointsColl) {  // iterate over the Steiner points according to their distance from p_i
				// Delete p from the graph and insert a new point of the class BuildingPoint at its place.
				
				Set<Segment> edges = delGraph.edgesOf(p);
				//LinkedList<Segment> edges = p.getEdgesInDelaunayGraph();
				
				LinkedList<Point> neighbors = new LinkedList<Point>();
				for(Segment se : edges) {
					neighbors.add(se.getOtherPoint(p));
				}
				
				delGraph.removeVertex(p);
				
				BuildingPoint bp = new BuildingPoint(p.getX(),p.getY(),BuildingNetwork.getCounterIDs(),true);
				bp.setOriginalBefore(p_i);
				bp.setOriginalAfter(p_iplus1);
				delGraph.addVertex(bp);
				BuildingNetwork.setCounterIDs(BuildingNetwork.getCounterIDs()+1);
				for(Point neighbor : neighbors) {
					Segment seNew = new Segment(bp,neighbor);
					
//					delGraph.addEdge(seNew);
					delGraph.addEdge(bp,neighbor,seNew);
					delGraph.setEdgeWeight(seNew, Math.hypot(bp.getX()-neighbor.getX(), bp.getY()-neighbor.getY()));
				}
				
				bp.addBuilding(b);
				if(i+counter < polyline.size()) {
					polyline.add(i+counter,bp);
				}
				else {  // Inserting a point within the segments p_N - p_0
					polyline.add(bp);
				}
				bn.getBuildingPoints().add(bp);
				counter++;
			}
			i = i + counter -1;  // index has to be increased since new points have been inserted
    	}
    }
    
    
    /**
     * Finds for a road the Steiner points created during the Delaunay triangulation which lie on
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
				if(p.getDistanceTo(p_i,p_iplus1)<Main.eps && !p.isReplaced()) {  // ... but only if they are lying on the segment p_i - p_iplus1 and have not been associated with another road yet
					delPointsOrdered.put(p.getDist(p_i),p);
					p.setReplaced(true);
				}
			}
			Collection<Point> delPointsColl = delPointsOrdered.values();
			int counter = 1;
			for(Point p : delPointsColl) {  // iterate over the Steiner points according to their distance from p_i
				// Delete p from the graph and insert a new point of the class RoadPoint at its place.
				Set<Segment> edges = delGraph.edgesOf(p);
				//LinkedList<Segment> edges = p.getEdgesInDelaunayGraph();
				
				LinkedList<Point> neighbors = new LinkedList<Point>();
				for(Segment se : edges) {
					neighbors.add(se.getOtherPoint(p));
				}
				delGraph.removeVertex(p);
				RoadPoint rp = new RoadPoint(p.getX(),p.getY(),RoadNetwork.getCounterIDs(),true);
				rp.setOriginalBefore(p_i);
				rp.setOriginalAfter(p_iplus1);
				delGraph.addVertex(rp);
				RoadNetwork.setCounterIDs(RoadNetwork.getCounterIDs()+1);
				for(Point neighbor : neighbors) {
					Segment seNew = new Segment(rp,neighbor);
					
//					delGraph.addEdge(seNew);
					delGraph.addEdge(rp,neighbor,seNew);
					delGraph.setEdgeWeight(seNew, Math.hypot(rp.getX()-neighbor.getX(), rp.getY()-neighbor.getY()));
				}
				rp.addRoad(r);
				polyline.add(i+counter,rp);
				counter++;
			}
			i = i + counter -1;  // the index has to be increased since new points have been inserted
    	}
    }
	
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
		return delGraph;
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
}
