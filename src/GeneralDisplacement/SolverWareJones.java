package GeneralDisplacement;
import gurobi .*;
import gurobi.GRB.DoubleAttr;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.Vector;

import org.jgrapht.graph.SimpleGraph;
import org.jgrapht.graph.SimpleWeightedGraph;

import org.jgrapht.alg.shortestpath.DijkstraShortestPath;
import org.jheaps.*;
import org.locationtech.jts.geom.Coordinate;

import BuildingDisplacement.Building;
import BuildingDisplacement.BuildingGraph;
import BuildingDisplacement.BuildingNetwork;
import BuildingDisplacement.BuildingPair;
import BuildingDisplacement.BuildingPoint;
import RoadDisplacement.*;
//import general.Graph;
import general.Point;
import general.Collector;
import general.Segment;

public class SolverWareJones {

	private SimpleWeightedGraph<Point, Segment> G;   // graph representing the input geometry of the roads and buildings
	private SimpleWeightedGraph<Point, Segment> G_prime 
	= new SimpleWeightedGraph<Point, Segment>(Segment.class);  // enriched graph with the input geometry plus the Steiner nodes which are incident to bottleneck edges

	private SimpleWeightedGraph<Point, Segment> T;   // graph of the Delaunay triangulation

	private BuildingNetwork bn;     // the geometry of the buildings
	private RoadNetwork rn;         // the geometry of the roads
	private double eps;  // tolerance for two points being regarded as one [m]
	private double M;    // big constant for the constraints
	private double w_pos;   // weight for point displacements 
	private double w_edge;   // weight for edge distortions
	private double w_select;   // weight for unselecting objects

	private double t; // threshold for the maximally allowed stretch factor of an edge

	private LinkedList<Segment> bottleneckEdges;

	// Bounding Box of the input coordinates
	private double minX;
	private double minY;
	private double maxX;
	private double maxY;


	//private SimpleWeightedGraph<Point,Segment> graphWithDelHeuristic;  // Hilfsgraph zur Auswertung des relaxierten Verfahrens

	double road_factor = 1.0; // Weighting of road edges => by this factor higher than bottleneck edges
	
	// ratio between the weights of bottleneck edges connecting a building with a road and those 
	// connecting two buildings => the first ones are by the factor conflict-factor higher weighted
	// than the second ones (see paper by Ware et al., p. 750)
	double conflict_factor = 10.0; 
	
	public SolverWareJones(SimpleWeightedGraph<Point,Segment> G, SimpleWeightedGraph<Point,Segment> G_prime,
			SimpleWeightedGraph<Point,Segment> T, GeneralNetworkWareJones gn, 
			double eps, double t,
			double m, double w_pos, double w_edge, double w_select, double road_factor) {
		this.G = G;
		this.G_prime = G_prime;
		this.T = T;
		this.bn = gn.getBn();
		this.rn = gn.getRn();
		this.eps = eps;
		this.M = m;
		this.w_pos = w_pos;
		this.w_edge = w_edge;
		this.w_select = w_select;
		this.bottleneckEdges = new LinkedList<Segment>();
		this.minX = gn.getMinX();
		this.maxX = gn.getMaxX();
		this.minY = gn.getMinY();
		this.maxY = gn.getMaxY();
		this.t = t;
		this.road_factor = road_factor;
	}


	/**
	 * Solves the exact program.
	 * @param weightObjects if true: objects obtain individual weights (here: proportional to length/area)
	 * @param withSelectionBuildings if true: selection of buildings is possible, if false: only displacement
	 * @param stayWithinMap if true: the objects have to remain within the bounding box of the input points
	 * @param downweightLongBottleneckEdges if true: bottleneck edges longer than epsilon obtain smaller weights (see method downweightFactor)
	 * @param enforceSelectionWJ if true, the selection of buildings is done according to the selection from the paper by Ware et al.
	 * @param withDisplacementRoads if true, road nodes may be displaced
	 */
	public void solve(boolean weightObjects, boolean withSelectionBuildings, 
			boolean withDisplacementRoads, boolean stayWithinMap, boolean enforceSelectionWJ,
			boolean downweightLongBottleneckEdges) {

		if(bottleneckEdges.size()==0) {
			addBottleneckEdges();
			removeUnusedSteinerPoints();
		}
		if(bottleneckEdges.size() > 0) {
			try {
				GRBEnv env = new GRBEnv ("GeneralNetworkWareJones.log ");
				GRBModel model = new GRBModel ( env );

				// Define the objective
				GRBQuadExpr obj = new GRBQuadExpr ();

				// Define the variables for the points (for road points only if roads are allowed 
				// to be displaced)
				for(Point p : G_prime.vertexSet()) {
					if(p.getClass().toString().equals("class BuildingDisplacement.BuildingPoint")
							|| (p.getClass().toString().equals("class RoadDisplacement.RoadPoint")
									&& withDisplacementRoads)){
						GRBVar dx = model.addVar(-M , M , 0.0 , GRB.CONTINUOUS , "");
						GRBVar dy = model.addVar(-M , M , 0.0 , GRB.CONTINUOUS , "");
						p.setDx(dx);
						p.setDy(dy);
						obj.addTerm(w_pos, dx, dx );
						obj.addTerm(w_pos, dy, dy );

						if(stayWithinMap) {
							GRBLinExpr expr1 = new GRBLinExpr ();
							expr1.addTerm(1.0 , dx ); 
							double rhs1 = -p.getX() + minX;
							model . addConstr ( expr1 , GRB.GREATER_EQUAL , rhs1 , "");
							GRBLinExpr expr2 = new GRBLinExpr ();
							expr2.addTerm(1.0 , dx ); 
							double rhs2 = -p.getX() + maxX;
							model . addConstr ( expr2 , GRB.LESS_EQUAL , rhs2 , "");
							GRBLinExpr expr3 = new GRBLinExpr ();
							expr3.addTerm(1.0 , dy ); 
							double rhs3 = -p.getY() + minY;
							model . addConstr ( expr3 , GRB.GREATER_EQUAL , rhs3 , "");
							GRBLinExpr expr4 = new GRBLinExpr ();
							expr4.addTerm(1.0 , dy ); 
							double rhs4 = -p.getY() + maxY;
							model . addConstr ( expr4 , GRB.LESS_EQUAL , rhs4 , "");
						}
					}
				}

				if(withSelectionBuildings) {
					// Define the variables for the buildings
					for(Building b : bn.getBuildings()) {
						// Variable for the selection of a building
						GRBVar z;
						if(!enforceSelectionWJ || b.getSelectionWareJones().equals("unknown")) {
							z = model.addVar(0.0 , 1.0 , 0.0 , GRB.BINARY , "");
						}
						else if(b.getSelectionWareJones().equals("selected")) {
							z = model.addVar(1.0 , 1.0 , 0.0 , GRB.BINARY , "");
						}
						else {
							z = model.addVar(0.0 , 0.0 , 0.0 , GRB.BINARY , "");
						}
						if(!weightObjects) {  // case 1: all buildings obtain the same weight
							obj.addTerm(-w_select,z);   
						}
						else {   // case 2: buildings obtain individual weights
							obj.addTerm(-w_select*b.getWeight(),z);
						}
						b.setZ(z);

					}

					// define the collector for each node which belongs to > 1 buildings
					for(Point p : G_prime.vertexSet()) {
						if(p.getClass().toString().equals("class BuildingDisplacement.BuildingPoint")) {
							BuildingPoint bp = (BuildingPoint) p;
							if(bp.getBuildings().size()>1) {
								Collector co = new Collector(p);
								bp.setCollector(co);
								GRBVar z_co = model.addVar(0.0 , 1.0 , 0.0 , GRB.BINARY , "");
								co.setZ(z_co);
								for(Building b : bp.getBuildings()) {
									// Constraint z_Building <= z_CO (select building only if the collector is selected too)
									GRBLinExpr exprCO = new GRBLinExpr ();
									exprCO.addTerm(1.0 , b.getZ() );
									exprCO.addTerm(-1.0 , z_co );
									model.addConstr(exprCO, GRB.LESS_EQUAL, 0, "");
								}
							}
						}
					}
				}

				// Setting up the variables and constraints for the bottleneck edges
				for(Segment se : bottleneckEdges) {
					Point p_i = se.getStart();
					Point p_j = se.getEnd();
					GRBVar deltaX = model.addVar(-M , M , 0.0 , GRB.CONTINUOUS , "");
					GRBVar deltaY = model.addVar(-M , M , 0.0 , GRB.CONTINUOUS , "");
					double alpha = 1.0;
					if(downweightLongBottleneckEdges) {
						alpha = downweightFactor(p_i, p_j);
					}
					if(se.connectsBuildingWithRoad()) {
						alpha = alpha*conflict_factor;
					}
					obj.addTerm(alpha*w_edge, deltaX, deltaX );
					obj.addTerm(alpha*w_edge, deltaY, deltaY );
					se.setDeltaX(deltaX);
					se.setDeltaY(deltaY);

					double d_ij = p_i.getDist(p_j);
					double s = getS(d_ij);

					if(p_i.getClass().toString().equals("class RoadDisplacement.RoadPoint")
							&& p_j.getClass().toString().equals("class BuildingDisplacement.BuildingPoint")) {
						GRBVar dx_j = p_j.getDx();
						GRBVar dy_j = p_j.getDy();
						GRBVar dx_i = null;
						GRBVar dy_i = null;
						if(withDisplacementRoads) {
							dx_i = p_i.getDx();
							dy_i = p_i.getDy();
						}
						if(withSelectionBuildings) {
							GRBVar z_j = ((BuildingPoint) p_j).getBuildings().get(0).getZ();
							if(((BuildingPoint) p_j).getBuildings().size()>1) {
								z_j = p_j.getCollector().getZ();
							}
							addConstraints(model, deltaX, deltaY, p_i,
									p_j, dx_i, dy_i, dx_j, dy_j, null, z_j, s, true);
						}
						else {
							addConstraints(model, deltaX, deltaY, p_i,
									p_j, dx_i, dy_i, dx_j, dy_j, null, null, s, false);
						}
					}
					else if(p_i.getClass().toString().equals("class BuildingDisplacement.BuildingPoint")
							&& p_j.getClass().toString().equals("class RoadDisplacement.RoadPoint")) {
						GRBVar dx_i = p_i.getDx();
						GRBVar dy_i = p_i.getDy();
						GRBVar dx_j = null;
						GRBVar dy_j = null;
						if(withDisplacementRoads) {
							dx_j = p_j.getDx();
							dy_j = p_j.getDy();
						}
						if(withSelectionBuildings) {
							GRBVar z_i = ((BuildingPoint) p_i).getBuildings().get(0).getZ();
							if(((BuildingPoint) p_i).getBuildings().size()>1) {
								z_i = p_i.getCollector().getZ();
							}
							addConstraints(model, deltaX, deltaY, p_i,
									p_j, dx_i, dy_i, dx_j, dy_j, z_i, null, s, true);
						}
						else {
							addConstraints(model, deltaX, deltaY, p_i,
									p_j, dx_i, dy_i, dx_j, dy_j, null, null, s, false);
						}
					}
					else if(p_i.getClass().toString().equals("class BuildingDisplacement.BuildingPoint")
							&& p_j.getClass().toString().equals("class BuildingDisplacement.BuildingPoint")) {
						GRBVar dx_i = p_i.getDx();
						GRBVar dy_i = p_i.getDy();
						GRBVar dx_j = p_j.getDx();
						GRBVar dy_j = p_j.getDy();
						LinkedList<Building> commonBuildings = 
								BuildingPoint.getCommonBuildings((BuildingPoint) p_i, (BuildingPoint) p_j);
						if(commonBuildings.size()==0) { // if the two points belong to the same building the constraint is dropped
							if(withSelectionBuildings) {
								GRBVar z_i = ((BuildingPoint) p_i).getBuildings().get(0).getZ();
								if(((BuildingPoint) p_i).getBuildings().size()>1) {
									z_i = p_i.getCollector().getZ();
								}
								GRBVar z_j = ((BuildingPoint) p_j).getBuildings().get(0).getZ();
								if(((BuildingPoint) p_j).getBuildings().size()>1) {
									z_j = p_j.getCollector().getZ();
								}
								addConstraints(model, deltaX, deltaY, p_i,
										p_j, dx_i, dy_i, dx_j, dy_j, z_i, z_j, s, true);
							}
							else {
								addConstraints(model, deltaX, deltaY, p_i,
										p_j, dx_i, dy_i, dx_j, dy_j, null, null, s, false);
							}
						}
					}
					else if(p_i.getClass().toString().equals("class RoadDisplacement.RoadPoint")
							&& p_j.getClass().toString().equals("class RoadDisplacement.RoadPoint")
							&& withDisplacementRoads) {
						GRBVar dx_i = p_i.getDx();
						GRBVar dy_i = p_i.getDy();
						GRBVar dx_j = p_j.getDx();
						GRBVar dy_j = p_j.getDy();
						addConstraints(model, deltaX, deltaY, p_i,
								p_j, dx_i, dy_i, dx_j, dy_j, null, null, s, false);
					}
				}

				// Define the edge variables also for the remaining edges. For the buildings they are forced to
				// 0 to prevent them from being distorted. For the roads too, if withDisplacementRoads is set to false.
				for(Segment se : G_prime.edgeSet()) {
					Point p_i = se.getStart();
					Point p_j = se.getEnd();
					if(se.getDeltaX()==null) {
						GRBVar deltaX=null, deltaY=null;
						if(p_i.getClass().toString().equals("class BuildingDisplacement.BuildingPoint")) {
							deltaX = model.addVar(0.0 , 0.0 , 0.0 , GRB.CONTINUOUS , "");
							deltaY = model.addVar(0.0 , 0.0 , 0.0 , GRB.CONTINUOUS , "");
						}
						else if(withDisplacementRoads) {
							deltaX = model.addVar(-M , M , 0.0 , GRB.CONTINUOUS , "");
							deltaY = model.addVar(-M , M , 0.0 , GRB.CONTINUOUS , "");
						}
						if(deltaX!=null && deltaY!=null) {
							if(p_i.getClass().toString().equals("class RoadDisplacement.RoadPoint")
									&& p_j.getClass().toString().equals("class RoadDisplacement.RoadPoint")) {
								obj.addTerm(road_factor*w_edge, deltaX, deltaX );  // distortion of a road edge is stronger penalized than that of a bottleneck edge
								obj.addTerm(road_factor*w_edge, deltaY, deltaY );
							}
							else {
								obj.addTerm(w_edge, deltaX, deltaX );  
								obj.addTerm(w_edge, deltaY, deltaY );
							}
							se.setDeltaX(deltaX);
							se.setDeltaY(deltaY);
						}

						if(p_i.getClass().toString().equals("class BuildingDisplacement.BuildingPoint")
								&& p_j.getClass().toString().equals("class BuildingDisplacement.BuildingPoint")) {
							GRBVar dx_i = p_i.getDx();
							GRBVar dx_j = p_j.getDx();
							GRBVar dy_i = p_i.getDy();
							GRBVar dy_j = p_j.getDy();
							if(withSelectionBuildings) {
								LinkedList<Building> buildings = BuildingPoint.getCommonBuildings((BuildingPoint) p_i, (BuildingPoint) p_j);
								if(buildings.size()==1) {
									addConstraints(model, deltaX, deltaY, p_i,
											p_j, dx_i, dy_i, dx_j, dy_j, buildings.get(0).getZ(), buildings.get(0).getZ(), 1, true);
								}
								else {
									Collector c = new Collector(se);
									se.setCollector(c);
									GRBVar z_co = model.addVar(0.0 , 1.0 , 0.0 , GRB.BINARY , "");
									c.setZ(z_co);
									for(Building b : buildings) {
										// Constraint z_Building <= z_co (select building only if collector is selected too)
										GRBLinExpr exprCO = new GRBLinExpr ();
										exprCO.addTerm(1.0 , b.getZ() );
										exprCO.addTerm(-1.0 , z_co );
										model.addConstr(exprCO, GRB.LESS_EQUAL, 0, "");
									}
									addConstraints(model, deltaX, deltaY, p_i,
											p_j, dx_i, dy_i, dx_j, dy_j, z_co, z_co, 1, true);
								}
							}
							else {
								addConstraints(model, deltaX, deltaY, p_i,
										p_j, dx_i, dy_i, dx_j, dy_j, null, null, 1, false);
							}
						}
						else if(withDisplacementRoads) {  // road nodes may be displaced
							GRBVar dx_i = p_i.getDx();
							GRBVar dx_j = p_j.getDx();
							GRBVar dy_i = p_i.getDy();
							GRBVar dy_j = p_j.getDy();
							addConstraints(model, deltaX, deltaY, p_i,
									p_j, dx_i, dy_i, dx_j, dy_j, null, null, 1, false);
						}
					}
				}

				model.setObjective (obj, GRB.MINIMIZE);			
				model.optimize ();

				for(Point p : G_prime.vertexSet()) {
					if(p.getDx()!=null && p.getDy()!= null) {
						p.setDxDouble(p.getDx().get(GRB.DoubleAttr.X));
						p.setDyDouble(p.getDy().get(GRB.DoubleAttr.X));
					}
				}	
				for(Building b : bn.getBuildings()) {
					if(b.getZ()!=null) {
						if(Math.abs(b.getZ().get(GRB.DoubleAttr.X)-1)<0.01) {
							b.setSelected(true);
						}
						else {
							b.setSelected(false);
						}
					}
				}
				for(Road r : rn.getRoads()) {
					r.setSelected(true);
				}

				//removeBottleneckEdges();
				model.dispose();

			}catch ( GRBException e ) {
				e.printStackTrace();
				System . out . println (" Error code : " + e . getErrorCode () + ". " +
						e . getMessage ());
			}
		}
		else {
			// no optimization necessary since constraints are already met
			for(Point p : G_prime.vertexSet()) {
				p.setDxDouble(0);
				p.setDyDouble(0);
			}	
		}

	}


	/**
	* Sets up the constraints for an edge (5-8 in the paper). 
	* If the edge is associated to one object/collector only (|o(e)|=1), the
	* same variable must be given for z_i and z_j as input.
	 */
	public void addConstraints(GRBModel model, GRBVar deltaX, GRBVar deltaY, Point p_i,
			Point p_j, GRBVar dx_i, GRBVar dy_i, GRBVar dx_j, GRBVar dy_j, GRBVar z_i,
			GRBVar z_j, double s, boolean withRemoval) throws GRBException {
		if(withRemoval) {
			// 1st constraint
			GRBLinExpr expr1 = new GRBLinExpr ();
			expr1.addTerm(1.0 , deltaX ); 
			if(dx_j!=null) expr1.addTerm(-1.0 , dx_j );
			if(dx_i!=null) expr1.addTerm(1.0, dx_i);
			if(z_i!=null) expr1.addTerm(-M, z_i);
			if(z_j!=null) expr1.addTerm(-M, z_j);
			double rhs1 = (1-s)*(p_j.getX()-p_i.getX()) - M;
			if(z_i!=null && z_j!=null) {
				rhs1 = (1-s)*(p_j.getX()-p_i.getX()) - 2*M;
			}
			else if(z_i==null && z_j==null) {
				rhs1 = (1-s)*(p_j.getX()-p_i.getX());
			}
			model . addConstr ( expr1 , GRB.GREATER_EQUAL , rhs1 , "");

			// 2nd constraint
			GRBLinExpr expr2 = new GRBLinExpr ();
			expr2.addTerm(1.0 , deltaY ); 
			if(dy_j!=null) expr2.addTerm(-1.0 , dy_j );
			if(dy_i!=null) expr2.addTerm(1.0, dy_i);
			if(z_i!=null) expr2.addTerm(-M, z_i);
			if(z_j!=null) expr2.addTerm(-M, z_j);
			double rhs2 = (1-s)*(p_j.getY()-p_i.getY()) - M;
			if(z_i!=null && z_j!=null) {
				rhs2 = (1-s)*(p_j.getY()-p_i.getY()) - 2*M;
			}
			else if(z_i==null && z_j==null) {
				rhs2 = (1-s)*(p_j.getY()-p_i.getY());
			}
			model . addConstr ( expr2 , GRB.GREATER_EQUAL , rhs2 , "");

			// 3rd constraint
			GRBLinExpr expr3 = new GRBLinExpr ();
			expr3.addTerm(1.0 , deltaX ); 
			if(dx_j!=null) expr3.addTerm(1.0 , dx_j );
			if(dx_i!=null) expr3.addTerm(-1.0, dx_i);
			if(z_i!=null) expr3.addTerm(-M, z_i);
			if(z_j!=null) expr3.addTerm(-M, z_j);
			double rhs3 = (s-1)*(p_j.getX()-p_i.getX()) - M;
			if(z_i!=null && z_j!=null) {
				rhs3 = (s-1)*(p_j.getX()-p_i.getX()) - 2*M;
			}
			else if(z_i==null && z_j==null) {
				rhs3 = (s-1)*(p_j.getX()-p_i.getX());
			}
			model . addConstr ( expr3 , GRB.GREATER_EQUAL , rhs3 , "");

			// 4th constraint
			GRBLinExpr expr4 = new GRBLinExpr ();
			expr4.addTerm(1.0 , deltaY ); 
			if(dy_j!=null) expr4.addTerm(1.0 , dy_j );
			if(dy_i!=null) expr4.addTerm(-1.0, dy_i);
			if(z_i!=null) expr4.addTerm(-M, z_i);
			if(z_j!=null) expr4.addTerm(-M, z_j);
			double rhs4 = (s-1)*(p_j.getY()-p_i.getY()) - M;
			if(z_i!=null && z_j!=null) {
				rhs4 = (s-1)*(p_j.getY()-p_i.getY()) - 2*M;
			}
			else if(z_i==null && z_j==null) {
				rhs4 = (s-1)*(p_j.getY()-p_i.getY());
			}
			model . addConstr ( expr4 , GRB.GREATER_EQUAL , rhs4 , "");
		}
		else {

			// 1st constraint
			GRBLinExpr expr1 = new GRBLinExpr ();
			expr1.addTerm(1.0 , deltaX ); 
			if(dx_j!=null) expr1.addTerm(-1.0 , dx_j );
			if(dx_i!=null) expr1.addTerm(1.0, dx_i);
			double rhs1 = (1-s)*(p_j.getX()-p_i.getX());
			model . addConstr ( expr1 , GRB.EQUAL , rhs1 , "");

			// 2nd constraint
			GRBLinExpr expr2 = new GRBLinExpr ();
			expr2.addTerm(1.0 , deltaY ); 
			if(dy_j!=null) expr2.addTerm(-1.0 , dy_j );
			if(dy_i!=null) expr2.addTerm(1.0, dy_i);
			double rhs2 = (1-s)*(p_j.getY()-p_i.getY());
			model . addConstr ( expr2 , GRB.EQUAL , rhs2 , "");

		}

	}

	/**
	 * Builds the graph containing the nodes at their final positions (without the Steiner
	 * nodes).
	 */
	public SimpleWeightedGraph<Point, Segment> getSolutionGraph() throws GRBException{
		SimpleWeightedGraph<Point, Segment> solGraph = 
				new SimpleWeightedGraph<Point, Segment>(Segment.class);
		ArrayList<BuildingPoint> BuildingPointsSol = new ArrayList<BuildingPoint>();  // BuildingPoints of the solution graph, ordered by their ids
		ArrayList<RoadPoint> RoadPointsSol = new ArrayList<RoadPoint>();  // RoadPoints of the solution graph, ordered by their ids
		for(Point p : G.vertexSet()) {
			if(p.isDelaunay())	System.out.println("ERROR in getSolutionGraphFrom() !!!!!!!!!!!!!!!!!!");
			double dx = p.getDxDouble();
			double dy = p.getDyDouble();
			if(p.getClass().toString().equals("class BuildingDisplacement.BuildingPoint")) {
				BuildingPoint bpSol = new BuildingPoint(p.getX()+dx, p.getY()+dy, p.getId());
				bpSol.setOriginalCoord(p.getOriginalCoord());
				if(p.getId() >= BuildingPointsSol.size()) {
					while (p.getId()>BuildingPointsSol.size()) {
						BuildingPointsSol.add(null);
					}
					BuildingPointsSol.add(bpSol);
				}
				else {
					BuildingPointsSol.set(p.getId(),bpSol);
				}
				solGraph.addVertex(bpSol);
			}
			else if(p.getClass().toString().equals("class RoadDisplacement.RoadPoint")) {
				RoadPoint rpSol = new RoadPoint(p.getX()+dx, p.getY()+dy, p.getId());
				rpSol.setOriginalCoord(p.getOriginalCoord());
				if(p.getId() >= RoadPointsSol.size()) {
					while (p.getId()>RoadPointsSol.size()) {
						RoadPointsSol.add(null);
					}
					RoadPointsSol.add(rpSol);
				}
				else {
					RoadPointsSol.set(p.getId(),rpSol);
				}
				solGraph.addVertex(rpSol);
			}
			else {
				System.out.println("Error in getSolutionGraphFrom() in class Solver"
						+ " (Setting up the nodes)!");
			}

		}
		for(Segment se : G.edgeSet()) {
			Point start = se.getStart();
			Point end = se.getEnd();
			if(start.getClass().toString().equals("class BuildingDisplacement.BuildingPoint")
					&& end.getClass().toString().equals("class BuildingDisplacement.BuildingPoint")) {
				LinkedList<Building> buildings = 
						BuildingPoint.getCommonBuildings((BuildingPoint) start,(BuildingPoint) end);
				for(Building b : buildings) {  // edge is drawn if at least one of the common buildings is selected
					if(b.isSelected()){  
						BuildingPoint startSol = BuildingPointsSol.get(start.getId());
						BuildingPoint endSol = BuildingPointsSol.get(end.getId());
						if(startSol!=null && endSol!=null) {
							Segment edge = new Segment(startSol, endSol);
							solGraph.addEdge(startSol,endSol,edge);
						}
						else {
							System.out.println("Warning, a BuildingPoint is missing (method "
									+ "getSolutionGraphFrom, class Solver)!");
						}
						break;
					}
				}
			}
			else if(start.getClass().toString().equals("class RoadDisplacement.RoadPoint")
					&& end.getClass().toString().equals("class RoadDisplacement.RoadPoint")) {
				Road road = RoadPoint.getCommonRoad((RoadPoint) start,(RoadPoint) end);					
				if(road.isSelected()){  
					RoadPoint startSol = RoadPointsSol.get(start.getId());
					RoadPoint endSol = RoadPointsSol.get(end.getId());
					if(startSol!=null && endSol!=null) {
						Segment edge = new Segment(startSol, endSol);
						solGraph.addEdge(startSol,endSol,edge);
					}
					else {
						System.out.println("Warning, a RoadPoint is missing (method "
								+ "getSolutionGraphFrom, class Solver)!");
					}
				}
			}
			else {
				System.out.println("Error in getSolutionGraphFrom() in class Solver"
						+ " (Inserting the edges)!");
			}
		}
		return solGraph;
	}
	
//	/**
//	 * Stellt den Graphen auf, der die Knoten an den finalen Positionen enth�lt. Falls
//	 * withDelaunayPoints auf true gesetzt ist, enth�lt der L�sungsgraph auch die Punkte, die in
//	 * der Delaunay-Triangulierung neu dazugekommen sind.
//	 */
//	public SimpleWeightedGraph<Point, Segment> getSolutionGraph(boolean withDelaunayPoints) 
//			throws GRBException{
//		if(withDelaunayPoints) {
//			return getSolutionGraphFrom(G_prime);
//		}
//		else {
//			return getSolutionGraphFrom(G);
//		}
//	}
//
//	public SimpleWeightedGraph<Point, Segment> getSolutionGraphFrom(
//			SimpleWeightedGraph<Point,Segment> g) throws GRBException{
//		SimpleWeightedGraph<Point, Segment> solGraph = 
//				new SimpleWeightedGraph<Point, Segment>(Segment.class);
//		ArrayList<BuildingPoint> BuildingPointsSol = new ArrayList<BuildingPoint>();  // BuildingPoints of the solution graph, ordered by their ids
//		ArrayList<RoadPoint> RoadPointsSol = new ArrayList<RoadPoint>();  // RoadPoints of the solution graph, ordered by their ids
//		for(Point p : g.vertexSet()) {
//			double dx = p.getDxDouble();
//			double dy = p.getDyDouble();
//			if(p.getClass().toString().equals("class BuildingDisplacement.BuildingPoint")) {
//				BuildingPoint bpSol = new BuildingPoint(p.getX()+dx, p.getY()+dy, p.getId());
//				bpSol.setOriginalCoord(p.getOriginalCoord());
//				if(p.getId() >= BuildingPointsSol.size()) {
//					while (p.getId()>BuildingPointsSol.size()) {
//						BuildingPointsSol.add(null);
//					}
//					BuildingPointsSol.add(bpSol);
//				}
//				else {
//					BuildingPointsSol.set(p.getId(),bpSol);
//				}
//				solGraph.addVertex(bpSol);
//			}
//			else if(p.getClass().toString().equals("class RoadDisplacement.RoadPoint")) {
//				RoadPoint rpSol = new RoadPoint(p.getX()+dx, p.getY()+dy, p.getId());
//				rpSol.setOriginalCoord(p.getOriginalCoord());
//				if(p.getId() >= RoadPointsSol.size()) {
//					while (p.getId()>RoadPointsSol.size()) {
//						RoadPointsSol.add(null);
//					}
//					RoadPointsSol.add(rpSol);
//				}
//				else {
//					RoadPointsSol.set(p.getId(),rpSol);
//				}
//				solGraph.addVertex(rpSol);
//			}
//			else {
//				System.out.println("Fehler in getSolutionGraphFrom() in Klasse SolverWareJones"
//						+ " (Erstellen der Knoten)!");
//			}
//		}
//		for(Segment se : g.edgeSet()) {
//			Point start = se.getStart();
//			Point end = se.getEnd();
//			if(start.getClass().toString().equals("class BuildingDisplacement.BuildingPoint")
//					&& end.getClass().toString().equals("class BuildingDisplacement.BuildingPoint")) {
//				LinkedList<Building> buildings = 
//						BuildingPoint.getCommonBuildings((BuildingPoint) start,(BuildingPoint) end);
//				for(Building b : buildings) {  // Kante wird gezeichnet, wenn mindestens eines der gemeinsamen Geb�ude ausgew�hlt ist
//					if(b.isSelected()){  
//						BuildingPoint startSol = BuildingPointsSol.get(start.getId());
//						BuildingPoint endSol = BuildingPointsSol.get(end.getId());
//						if(startSol!=null && endSol!=null) {
//							Segment edge = new Segment(startSol, endSol);
//							solGraph.addEdge(startSol,endSol,edge);
//						}
//						else {
//							System.out.println("Warnung, ein BuildingPoint fehlt (Methode "
//									+ "getSolutionGraphFrom, Klasse SolverWareJones)!");
//						}
//						break;
//					}
//				}
//			}
//			else if(start.getClass().toString().equals("class RoadDisplacement.RoadPoint")
//					&& end.getClass().toString().equals("class RoadDisplacement.RoadPoint")) {
//				Road road = RoadPoint.getCommonRoad((RoadPoint) start,(RoadPoint) end);					
//				// pr�fe, ob die Stra�e, zu der diese Kante geh�rt, ausgew�hlt ist
//				if(road.isSelected()){  
//					RoadPoint startSol = RoadPointsSol.get(start.getId());
//					RoadPoint endSol = RoadPointsSol.get(end.getId());
//					if(startSol!=null && endSol!=null) {
//						Segment edge = new Segment(startSol, endSol);
//						solGraph.addEdge(startSol,endSol,edge);
//					}
//					else {
//						System.out.println("Warnung, ein RoadPoint fehlt (Methode "
//								+ "getSolutionGraphFrom, Klasse SolverWareJones)!");
//					}
//				}
//			}
//			else {
//				System.out.println("Fehler in getSolutionGraphFrom() in Klasse SolverWareJones"
//						+ " (Einf�gen der Kanten)!");
//			}
//		}
//		return solGraph;
//	}


	/**
	 * Updates the coordinates of the buildings and roads according to the solution.
	 */
	public void updateNetworks() {
		for(Road r : rn.getRoads()) {
			for(RoadPoint rp : r.getVertices()) {
				rp.setX(rp.getX()+rp.getDxDouble());
				rp.setY(rp.getY()+rp.getDyDouble());
				rp.setDxDouble(0);
				rp.setDyDouble(0);
			}
		}
		for(Building b : bn.getBuildings()) {
			for(BuildingPoint bp : b.getShell()) {
				bp.setX(bp.getX()+bp.getDxDouble());
				bp.setY(bp.getY()+bp.getDyDouble());
				bp.setDxDouble(0);
				bp.setDyDouble(0);
			}
			for(Vector<BuildingPoint> hole : b.getHoles()) {
				for(BuildingPoint bp : hole) {
					bp.setX(bp.getX()+bp.getDxDouble());
					bp.setY(bp.getY()+bp.getDyDouble());
					bp.setDxDouble(0);
					bp.setDyDouble(0);
				}
			}
		}
	}


	/**
	 * Finds the bottleneck edges and adds them to the graph G_prime.
	 * They are also stored in the list bottleneckEdges.
	 */
	public void addBottleneckEdges() {

		System.out.println("Add bottleneck edges ...");
		LinkedList<Segment> listAll = new LinkedList<Segment>();

		for (Segment se : T.edgeSet()) {   // for (Segment se : T.getAllEdges())
			if( !G_prime.containsEdge(se) ) {
				listAll.add(se);
			}
		}

		Comparator<Segment> c = new Comparator<Segment>() {
			@Override
			public int compare(Segment s1, Segment s2) { 
				Double d1 = s1.getLength(); //T.getEdgeWeight(s1);
				Double d2 = s2.getLength(); //T.getEdgeWeight(s2);
				return d1.compareTo(d2); }
		};
		Collections.sort(listAll, c);

		for ( Segment se : listAll) {
			Point s = se.getStart();  //T.getEdgeSource(se);
			Point u = se.getEnd(); //T.getEdgeTarget(se);
			double edgeWeight = se.getLength(); //T.getEdgeWeight(se);
			double radius = edgeWeight*t;
			DijkstraShortestPath<Point,Segment> dsp = new DijkstraShortestPath<Point,Segment>(G_prime,radius);
			double pathWeight = dsp.getPathWeight(s,u);
			//System.out.println("pathWeight: " + pathWeight + ", edgeWeight: " + edgeWeight);
			double dilation = pathWeight / edgeWeight;
			if (dilation > t) {
				G_prime.addEdge(s,u,se);
				G_prime.setEdgeWeight(se, se.getLength());
				bottleneckEdges.add(se);
				se.setBottleneckEdge(true);
				se.getStart().setIncidentToBottleneckEdge(true);
				se.getEnd().setIncidentToBottleneckEdge(true);
			}
		}
		System.out.println("Number of bottleneck edges added: " + bottleneckEdges.size());
	}


	/**
	 * Removes the bottleneck edges from G_prime and clears the list bottleneckEdges.
	 */
	public void removeBottleneckEdges() {
		for(Segment se : bottleneckEdges) {
			G_prime.removeEdge(se);
			se.getStart().setIncidentToBottleneckEdge(false);
			se.getEnd().setIncidentToBottleneckEdge(false);
		}
		bottleneckEdges.clear();
	}

	/**
	 * Determines the constant s for a point pair p_i, p_j with distance d_ij.
	 */
	public double getS(double d_ij) {
		if(d_ij >= eps) {
			return 1;
		}
		else {
			if(d_ij < Math.pow(10,-4)) {
				System.out.println("Problem in getS(), class SolverWareJones: point distance is very small!");
				return 1;
			}
			return eps/d_ij;
		}
	}


	/**
	 * Computes the factor alpha<1 by which the weight of a long bottleneck edge is 
	 * decreased (if downweightLongBottleneckEdges is set to true in the solve/solveViaRelaxation
	 * method).
	 * @param start 1. Punkt der bottleneck edge
	 * @param end 2. Punkt der bottleneck edge
	 */
	public double downweightFactor(Point start, Point end) {
		double d = start.getDist(end);
		double alpha = 1;
		if(d>eps)	alpha = eps*eps/(d*d);
		return alpha;
	}


	public SimpleWeightedGraph<Point,Segment> getBottleneckGraph(){
		SimpleWeightedGraph<Point,Segment> bottleneckGraph = 
				new SimpleWeightedGraph<Point, Segment>(Segment.class);
		for(Segment se : bottleneckEdges) {
			Point start = se.getStart();
			Point end = se.getEnd();
			if(!bottleneckGraph.containsVertex(start)) {
				bottleneckGraph.addVertex(start);
			}
			if(!bottleneckGraph.containsVertex(end)) {
				bottleneckGraph.addVertex(end);
			}
			bottleneckGraph.addEdge(start,end,se);
		}
		return bottleneckGraph;
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
	 * Evaluates the solution according to the objective.
	 * @param gWithDel graph of the original geometry, enriched with the Steiner points
	 */
	public double evaluate(SimpleWeightedGraph<Point,Segment> gWithDel,
			Set<Segment> bottleneckEdges,LinkedList<Building> buildings, 
			boolean downweightLongBEdges, boolean withDisplacementRoads) {
		System.out.println("\n----------------------- Evaluation "
				+ "-------------------------");
		System.out.println("# Points: " + gWithDel.vertexSet().size());
		System.out.println("# original edges: " + gWithDel.edgeSet().size());
		System.out.println("# initial bottleneck edges: " + bottleneckEdges.size());
		System.out.println("# buildings: " + buildings.size());
		double obj = 0;
		double term = 0;
		double largest_displacement = 0;  // largest displacement of a point
		for(Point p : gWithDel.vertexSet()) {
			if(p.isVisible()) {
				double dx = p.getX()+p.getDxDouble()-p.getOriginalCoord().getX();
				double dy = p.getY()+p.getDyDouble()-p.getOriginalCoord().getY();
				obj = obj + w_pos*(dx*dx+dy*dy);
				term = term + w_pos*(dx*dx+dy*dy);
				double displacement = Math.sqrt(dx*dx+dy*dy);
				if(displacement>largest_displacement) {
					largest_displacement = displacement;
				}
			}
		}
		System.out.println("Contribution of point movement: " + term);
		System.out.println("Largest displacement of a point: " + largest_displacement);
		term = 0;
		//counter=0;
		for(Segment se : gWithDel.edgeSet()) {
			if(se.isVisible()) {
				//counter++;
				Point start = se.getStart();
				Point end = se.getEnd();
				double dxAfter = start.getX()+start.getDxDouble()-end.getX()-end.getDxDouble();
				double dyAfter = start.getY()+start.getDyDouble()-end.getY()-end.getDyDouble();
				double dxBefore = start.getOriginalCoord().getX()-end.getOriginalCoord().getX();
				double dyBefore = start.getOriginalCoord().getY()-end.getOriginalCoord().getY();
				double deltaX = Math.abs(dxAfter-dxBefore);
				double deltaY = Math.abs(dyAfter-dyBefore);
				obj = obj + w_edge*(deltaX*deltaX+deltaY*deltaY);
				term = term + w_edge*(deltaX*deltaX+deltaY*deltaY);
			}
		}
		System.out.println("Contribution of the distortion of input edges: " + term);
		term = 0;
		for(Segment se : bottleneckEdges) {
			Point start = se.getStart();
			Point end = se.getEnd();
			if(start.isVisible() && end.isVisible()) {
				
				//--------- No weighting of bottleneck edges whose endpoints belong to the same object, as in the paper by Ware et al.
				if(start.getClass().toString().equals("class BuildingDisplacement.BuildingPoint")
					&& end.getClass().toString().equals("class BuildingDisplacement.BuildingPoint")) {
					if(BuildingPoint.getCommonBuildings((BuildingPoint) start, (BuildingPoint) end).size()>0) {
						continue;
					}
				}
				if(start.getClass().toString().equals("class RoadDisplacement.RoadPoint")
					&& end.getClass().toString().equals("class RoadDisplacement.RoadPoint")
					&& !withDisplacementRoads) {
					continue;
				}
				//---------
				double dxAfter = start.getX()+start.getDxDouble()-end.getX()-end.getDxDouble();
				double dyAfter = start.getY()+start.getDyDouble()-end.getY()-end.getDyDouble();
				double dxBefore = start.getOriginalCoord().getX()-end.getOriginalCoord().getX();
				double dyBefore = start.getOriginalCoord().getY()-end.getOriginalCoord().getY();
				double s = getS(start.getOriginalCoord().distance(end.getOriginalCoord()));
				double deltaX = Math.abs(dxAfter-s*dxBefore);
				double deltaY = Math.abs(dyAfter-s*dyBefore);
				double w_BE = 1;
				if(downweightLongBEdges) {
					w_BE = downweightFactor(new Point(start.getOriginalCoord()),
							new Point(end.getOriginalCoord()));
				}
				if(se.connectsBuildingWithRoad()) {
					w_BE = w_BE*conflict_factor;
				}
				obj = obj + w_edge*w_BE*(deltaX*deltaX+deltaY*deltaY);
				term = term + w_edge*w_BE*(deltaX*deltaX+deltaY*deltaY);
			}
		}
		System.out.println("Contribution of the initial bottleneck edges: " + term);
		term = 0;
		double sum_of_weights = 0;
		int numBuildingsRemoved = 0;
		for(Building b : buildings) {
			sum_of_weights = sum_of_weights + b.getWeight();
			if(!b.isSelected()) {
				obj = obj + w_select*b.getWeight();
				term = term + w_select*b.getWeight();
				numBuildingsRemoved++;
			}
		}
		System.out.println("# buildings removed: " + numBuildingsRemoved);
		System.out.println("Contribution of buildings removed: " + term);
		return obj;
	}


	/**
	 * Removes the Steiner points which are not incident to bottleneck edges.
	 */
	public void removeUnusedSteinerPoints() {
		System.out.println("Removing unused Steiner points...");
		for(Road r : rn.getRoads()) {
			for(int i=0; i<r.getVertices().size(); i++) {
				RoadPoint rp = r.getVertices().get(i);
				if(rp.isDelaunay() && !rp.isIncidentToBottleneckEdge()) {
					Object[] edges = G_prime.edgesOf(rp).toArray();
					if(edges.length==2) {
						Segment s1 = (Segment) edges[0];
						Segment s2 = (Segment) edges[1];
						RoadPoint rp1 = (RoadPoint) s1.getOtherPoint(rp);
						RoadPoint rp2 = (RoadPoint) s2.getOtherPoint(rp);
						G_prime.addEdge(rp1, rp2, new Segment(rp1,rp2));
						G_prime.removeVertex(rp);
						r.getVertices().remove(rp);
						i--;
					}
					else {
						System.out.println("Warning in method removeUnusedSteinerPoints in class SolverWareJones, a Steiner node could not be removed.");
					}
				}
			}
		}
		for(Building b : bn.getBuildings()) {
			Vector<BuildingPoint> shell = b.getShell();
			for(int i=0; i<shell.size(); i++) {
				BuildingPoint bp = shell.get(i);
				if(bp.isDelaunay() && !bp.isIncidentToBottleneckEdge()) {
					Object[] edges = G_prime.edgesOf(bp).toArray();
					if(edges.length==2) {
						Segment s1 = (Segment) edges[0];
						Segment s2 = (Segment) edges[1];
						BuildingPoint bp1 = (BuildingPoint) s1.getOtherPoint(bp);
						BuildingPoint bp2 = (BuildingPoint) s2.getOtherPoint(bp);
						G_prime.addEdge(bp1, bp2, new Segment(bp1, bp2));
						G_prime.removeVertex(bp);
						shell.remove(bp);
						i--;
					}
					else {
						System.out.println("Warning in method removeUnusedSteinerPoints in class SolverWareJones, a Steiner node could not be removed.");
					}
				}
			}
			Vector<Vector<BuildingPoint>> holes = b.getHoles();
			for(Vector<BuildingPoint> hole : holes) {
				for(int i=0; i<hole.size(); i++) {
					BuildingPoint bp = hole.get(i);
					if(bp.isDelaunay() && !bp.isIncidentToBottleneckEdge()) {
						Object[] edges = G_prime.edgesOf(bp).toArray();
						if(edges.length==2) {
							Segment s1 = (Segment) edges[0];
							Segment s2 = (Segment) edges[1];
							BuildingPoint bp1 = (BuildingPoint) s1.getOtherPoint(bp);
							BuildingPoint bp2 = (BuildingPoint) s2.getOtherPoint(bp);
							G_prime.addEdge(bp1, bp2, new Segment(bp1, bp2));
							G_prime.removeVertex(bp);
							hole.remove(bp);
							i--;
						}
						else {
							System.out.println("Error in method removeUnusedSteinerPoints in class Solver!!!");
						}
					}
				}
			}
		}
		System.out.println("Unused Steiner points removed.");
	}
	

}
