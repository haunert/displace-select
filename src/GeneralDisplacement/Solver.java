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
import general.Point;
import general.Collector;
import general.Segment;

public class Solver {
	
	SimpleWeightedGraph<Point, Segment> G;   // graph representing the input geometry of the roads and buildings
	SimpleWeightedGraph<Point, Segment> G_prime 
		    = new SimpleWeightedGraph<Point, Segment>(Segment.class);  // enriched graph with the input geometry plus the Steiner nodes which are incident to bottleneck edges
	SimpleWeightedGraph<Point, Segment> T;   // graph of the Delaunay triangulation
	BuildingNetwork bn;     // the geometry of the buildings
	RoadNetwork rn;         // the geometry of the roads
	GeneralNetwork gn;      // the geometry of buildings and roads together, with additional methods
	double eps;  // tolerance for two points being regarded as one [m]
	double M;    // big constant for the constraints
	double w_pos;   // weight for point displacements 
	double w_edge;   // weight for edge distortions
	double w_select;   // weight for unselecting objects
	double w_depend;  // weight for treating adjacent buildings the same (in a group of >= 3 buildings), if coupledSelectionBuildings is true
	
	double t;  // threshold for the maximally allowed stretch factor of an edge
	
	LinkedList<Segment> bottleneckEdges;
	
	// Bounding Box of the input coordinates
	double minX;
    double minY;
    double maxX;
    double maxY;
    
	double theta;  // threshold for the selection of objects in the heuristic

    SimpleWeightedGraph<Point,Segment> bottleneckGraph;  // graph for the bottleneck edges
    
	public Solver(SimpleWeightedGraph<Point, Segment> G, SimpleWeightedGraph<Point, Segment> G_prime,
			SimpleWeightedGraph<Point, Segment> T, GeneralNetwork gn, 
			double eps, double t, 
			double m, double w_pos, double w_edge, double w_select, double w_depend, 
			double theta) {
		this.G = G;
		this.G_prime = G_prime;
		this.T = T;
		this.gn = gn;
		this.bn = gn.getBn();
		this.rn = gn.getRn();
		this.eps = eps;
		this.M = m;
		this.w_pos = w_pos;
		this.w_edge = w_edge;
		this.w_select = w_select;
		this.w_depend = w_depend;
		this.bottleneckEdges = new LinkedList<Segment>();
		this.minX = gn.getMinX();
		this.maxX = gn.getMaxX();
		this.minY = gn.getMinY();
		this.maxY = gn.getMaxY();
		this.t = t;
		this.theta = theta;   
	}
	
	
	/**
	 * Solves the exact program.
	 * @param weightObjects if true: objects obtain individual weights (here: proportional to length/area)
	 * @param withSelection if true: selection is possible, if false: only displacement
	 * @param stayWithinMap if true: the objects have to remain within the bounding box of the input points
	 * @param downweightLongBottleneckEdges if true: bottleneck edges longer than epsilon obtain smaller weights (see method downweightFactor) 
	 * @param enforceConnectivityRoads if false: road network may be split into several components, if true: road network remains connected
	 * @param buildingOnlyIfRoad if true: additional constraints for the coupled selection of buildings and roads are applied
	 * @param coupledSelectionBuildings if true: additional soft constraint that adjacent buildings should be treated the same (in a group of >= 3 buildings)
	 */
	public void solve(boolean weightObjects, boolean withSelection, boolean stayWithinMap, 
				boolean downweightLongBottleneckEdges, 
				boolean enforceConnectivityRoads, boolean buildingOnlyIfRoad, 
				boolean coupledSelectionBuildings) {
		
		if(bottleneckEdges.size()==0) {
			addBottleneckEdges();
			removeUnusedSteinerPoints();
		}
				
		if(bottleneckEdges.size() > 0) {
			try {
				GRBEnv env = new GRBEnv ("GeneralNetwork.log ");
				GRBModel model = new GRBModel ( env );
				
				// Define the objective
				GRBQuadExpr obj = new GRBQuadExpr ();
				
				// Define the variables for the points
				for(Point p : G_prime.vertexSet()) {
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

				if(withSelection) {
					// Define the variables for the roads ...		
					SimpleWeightedGraph<Road,RoadPair> graphOfRoads = rn.getGraphOfRoads();
					for(Road r : graphOfRoads.vertexSet()) {
						// Variable for the selection of a road
						GRBVar z = model.addVar(0.0 , 1.0 , 0.0 , GRB.BINARY , "");
						if(!weightObjects) {  // case 1: all roads obtain the same weight
							obj.addTerm(-w_select,z);   // constant term of the objective is omitted here
						}
						else {   // case 2: roads obtain individual weights
							obj.addTerm(-w_select*r.getWeight(),z);
						}
						r.setZ(z);
						if(enforceConnectivityRoads) {
							GRBVar w_r = model.addVar(0.0 , 1.0 , 0.0 , GRB.BINARY , "");
							r.setIsSink(w_r);
						}
					}
					if(enforceConnectivityRoads){
						System.out.println("Flussbedingungen werden aufgestellt ...");
						for(RoadPair roadPair : graphOfRoads.edgeSet()) {
							Road r1 = roadPair.getR1();
							Road r2 = roadPair.getR2();
							GRBVar f_r1_r2 = model.addVar(0.0 , (double) rn.getRoads().size() , 0.0 , GRB.CONTINUOUS , "");
							GRBVar f_r2_r1 = model.addVar(0.0 , (double) rn.getRoads().size() , 0.0 , GRB.CONTINUOUS , "");
							r1.getOutFlow().add(f_r1_r2);
							r2.getOutFlow().add(f_r2_r1);
							r1.getInFlow().add(f_r2_r1);
							r2.getInFlow().add(f_r1_r2);
						}
						addConstraintsConnectivity(model);
						System.out.println("Flussbedingungen sind aufgestellt.");
					}
					// ... and for the buildings
					for(Building b : bn.getBuildings()) {
						// Variable for the selection of a building
						GRBVar z = model.addVar(0.0 , 1.0 , 0.0 , GRB.BINARY , "");
						if(!weightObjects) {  // case 1: all buildings obtain the same weight
							obj.addTerm(-w_select,z);   // constant term of the objective is omitted here
						}
						else {   // case 2: buildings obtain individual weights
							obj.addTerm(-w_select*b.getWeight(),z);
						}
						b.setZ(z);

						if(buildingOnlyIfRoad) {
							// building may only be selected if the closest road is selected too
							Road r = gn.getClosestRoadOfBuilding(b);
							GRBLinExpr expr = new GRBLinExpr ();
							expr.addTerm(1.0 , r.getZ() ); 
							expr.addTerm(-1.0 , z );
							model . addConstr ( expr , GRB.GREATER_EQUAL , 0 , "");
						}
					}
					if(coupledSelectionBuildings) {
						BuildingGraph graphOfBuildings = bn.getGraphOfBuildings();
						for(BuildingPair bPair : graphOfBuildings.getEdges()) {
							Building b1 = bPair.getB1();
							Building b2 = bPair.getB2();
							if(b1.getComponent().size()>2) {
								GRBVar z_b1_b2 = model.addVar(0.0 , 1.0 , 0.0 , GRB.BINARY , "");
								obj.addTerm(w_depend, z_b1_b2);
								GRBLinExpr expr1 = new GRBLinExpr (); // z_b1_b2 >= z_b1 - z_b2
								expr1.addTerm(1.0, z_b1_b2);
								expr1.addTerm(-1.0, b1.getZ());
								expr1.addTerm(1.0, b2.getZ());
								model.addConstr(expr1, GRB.GREATER_EQUAL, 0, "");
								GRBLinExpr expr2 = new GRBLinExpr (); // z_b1_b2 >= z_b2 - z_b1
								expr2.addTerm(1.0, z_b1_b2);
								expr2.addTerm(1.0, b1.getZ());
								expr2.addTerm(-1.0, b2.getZ());
								model.addConstr(expr2, GRB.GREATER_EQUAL, 0, "");
							}
						}
					}
					
					// define the collector for each node which belongs to > 1 objects and
					// which is incident to a bottleneck edge
					for(Point p : G_prime.vertexSet()) {
						if(p.isIncidentToBottleneckEdge()) {
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
							else {
								RoadPoint rp = (RoadPoint) p;
								if(rp.getRoads().size()>1) {
									Collector co = new Collector(p);
									rp.setCollector(co);
									GRBVar z_co = model.addVar(0.0 , 1.0 , 0.0 , GRB.BINARY , "");
									co.setZ(z_co);
									for(Road r : rp.getRoads()) {
										// Constraint z_Road <= z_CO (select road only if the collector is selected too)
										GRBLinExpr exprCO = new GRBLinExpr ();
										exprCO.addTerm(1.0 , r.getZ() );
										exprCO.addTerm(-1.0 , z_co );
										model.addConstr(exprCO, GRB.LESS_EQUAL, 0, "");
									}
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
					if(downweightLongBottleneckEdges) {
						double alpha = downweightFactor(p_i, p_j);
						obj.addTerm(alpha*w_edge, deltaX, deltaX );
						obj.addTerm(alpha*w_edge, deltaY, deltaY );
					}
					else {
						obj.addTerm(w_edge, deltaX, deltaX );
						obj.addTerm(w_edge, deltaY, deltaY );
					}
					se.setDeltaX(deltaX);
					se.setDeltaY(deltaY);

					GRBVar dx_i = p_i.getDx();
					GRBVar dx_j = p_j.getDx();
					GRBVar dy_i = p_i.getDy();
					GRBVar dy_j = p_j.getDy();
					double d_ij = p_i.getDist(p_j);
					double s = getS(d_ij);

					if(withSelection) {
						if(p_i.getClass().toString().equals("class RoadDisplacement.RoadPoint")
								&& p_j.getClass().toString().equals("class RoadDisplacement.RoadPoint")) {
							GRBVar z_i = ((RoadPoint) p_i).getRoads().get(0).getZ();
							if(((RoadPoint) p_i).getRoads().size()>1) {
								z_i = p_i.getCollector().getZ();
							}
							GRBVar z_j = ((RoadPoint) p_j).getRoads().get(0).getZ();
							if(((RoadPoint) p_j).getRoads().size()>1) {
								z_j = p_j.getCollector().getZ();
							}
							addConstraints(model, deltaX, deltaY, p_i,
									p_j, dx_i, dy_i, dx_j, dy_j, z_i, z_j, s, true);
							
						}
						else if(p_i.getClass().toString().equals("class RoadDisplacement.RoadPoint")
								&& p_j.getClass().toString().equals("class BuildingDisplacement.BuildingPoint")) {
							GRBVar z_i = ((RoadPoint) p_i).getRoads().get(0).getZ();
							if(((RoadPoint) p_i).getRoads().size()>1) {
								z_i = p_i.getCollector().getZ();
							}
							GRBVar z_j = ((BuildingPoint) p_j).getBuildings().get(0).getZ();
							if(((BuildingPoint) p_j).getBuildings().size()>1) {
								z_j = p_j.getCollector().getZ();
							}
							addConstraints(model, deltaX, deltaY, p_i,
										p_j, dx_i, dy_i, dx_j, dy_j, z_i, z_j, s, true);
						}
						else if(p_i.getClass().toString().equals("class BuildingDisplacement.BuildingPoint")
								&& p_j.getClass().toString().equals("class RoadDisplacement.RoadPoint")) {
							GRBVar z_i = ((BuildingPoint) p_i).getBuildings().get(0).getZ();
							if(((BuildingPoint) p_i).getBuildings().size()>1) {
								z_i = p_i.getCollector().getZ();
							}
							GRBVar z_j = ((RoadPoint) p_j).getRoads().get(0).getZ();
							if(((RoadPoint) p_j).getRoads().size()>1) {
								z_j = p_j.getCollector().getZ();
							}
							addConstraints(model, deltaX, deltaY, p_i,
									p_j, dx_i, dy_i, dx_j, dy_j, z_i, z_j, s, true);
						}
						else if(p_i.getClass().toString().equals("class BuildingDisplacement.BuildingPoint")
								&& p_j.getClass().toString().equals("class BuildingDisplacement.BuildingPoint")) {
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
							System.out.println("Fehler in solveWithBottleneckEdges() in Klasse Solver, "
									+ "ein Punkt ist nicht klassifizierbar!");
						}
					}
					else {
						addConstraints(model, deltaX, deltaY, p_i,
									p_j, dx_i, dy_i, dx_j, dy_j, null, null, s, false);
					}	
				}
				
				// Setting up the variables and constraints for the remaining edges of G_prime (with s = 1)
				for(Segment se : G_prime.edgeSet()) {
					if(se.getDeltaX()==null) {
						GRBVar deltaX = model.addVar(-M , M , 0.0 , GRB.CONTINUOUS , "");
						GRBVar deltaY = model.addVar(-M , M , 0.0 , GRB.CONTINUOUS , "");
						obj.addTerm(w_edge, deltaX, deltaX );  
						obj.addTerm(w_edge, deltaY, deltaY );
						se.setDeltaX(deltaX);
						se.setDeltaY(deltaY);
						
						Point p_i = se.getStart();
						Point p_j = se.getEnd();
						GRBVar dx_i = p_i.getDx();
						GRBVar dx_j = p_j.getDx();
						GRBVar dy_i = p_i.getDy();
						GRBVar dy_j = p_j.getDy();
						if(dx_i==null || dx_j==null || dy_i==null || dy_j==null) {
							System.out.println("p_i: " + p_i.toString());
							System.out.println("p_j: " + p_j.toString());
						}
						
						if(withSelection) {
							if(p_i.getClass().toString().equals("class RoadDisplacement.RoadPoint")
									&& p_j.getClass().toString().equals("class RoadDisplacement.RoadPoint")) {
								Road road = RoadPoint.getCommonRoad((RoadPoint) p_i, (RoadPoint) p_j);
								addConstraints(model, deltaX, deltaY, p_i,
										p_j, dx_i, dy_i, dx_j, dy_j, road.getZ(), road.getZ(), 1, true);
							}
							else if(p_i.getClass().toString().equals("class BuildingDisplacement.BuildingPoint")
									&& p_j.getClass().toString().equals("class BuildingDisplacement.BuildingPoint")) {
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
								System.out.println("Fehler in Methode solve in Klasse Solver, Punkte geh�ren nicht zur selben Klasse!");
							}
						}
						else {
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
					if(r.getZ()!=null) {
						if(Math.abs(r.getZ().get(GRB.DoubleAttr.X)-1)<0.01) {
							r.setSelected(true);
						}
						else {
							r.setSelected(false);
						}
					}
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
			// no optimization neccessary since no conflicts exist
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
							GRBVar z_j, double s, boolean withSelection) throws GRBException {
			if(withSelection) {
				// 1st Constraint
				GRBLinExpr expr1 = new GRBLinExpr ();
				expr1.addTerm(1.0 , deltaX ); 
				expr1.addTerm(-1.0 , dx_j );
				expr1.addTerm(1.0, dx_i);
				expr1.addTerm(-M, z_i);
				double rhs1 = (1-s)*(p_j.getX()-p_i.getX()) - M;
				if( z_j != z_i) {
					expr1.addTerm(-M, z_j);
					rhs1 = rhs1 - M;
				}
				model . addConstr ( expr1 , GRB.GREATER_EQUAL , rhs1 , "");
				
				// 2nd Constraint
				GRBLinExpr expr2 = new GRBLinExpr ();
				expr2.addTerm(1.0 , deltaY ); 
				expr2.addTerm(-1.0 , dy_j );
				expr2.addTerm(1.0, dy_i);
				expr2.addTerm(-M, z_i);
				double rhs2 = (1-s)*(p_j.getY()-p_i.getY()) - M;
				if( z_j != z_i) {
					expr2.addTerm(-M, z_j);
					rhs2 = rhs2 - M;
				}
				model . addConstr ( expr2 , GRB.GREATER_EQUAL , rhs2 , "");
				
				// 3rd Constraint
				GRBLinExpr expr3 = new GRBLinExpr ();
				expr3.addTerm(1.0 , deltaX ); 
				expr3.addTerm(1.0 , dx_j );
				expr3.addTerm(-1.0, dx_i);
				expr3.addTerm(-M, z_i);
				double rhs3 = (s-1)*(p_j.getX()-p_i.getX()) - M;
				if( z_j != z_i) {
					expr3.addTerm(-M, z_j);
					rhs3 = rhs3 - M;
				}
				model . addConstr ( expr3 , GRB.GREATER_EQUAL , rhs3 , "");
				
				// 4th Constraint
				GRBLinExpr expr4 = new GRBLinExpr ();
				expr4.addTerm(1.0 , deltaY ); 
				expr4.addTerm(1.0 , dy_j );
				expr4.addTerm(-1.0, dy_i);
				expr4.addTerm(-M, z_i);
				double rhs4 = (s-1)*(p_j.getY()-p_i.getY()) - M;
				if( z_j != z_i) {
					expr4.addTerm(-M, z_j);
					rhs4 = rhs4 - M;
				}
				model . addConstr ( expr4 , GRB.GREATER_EQUAL , rhs4 , "");
			}
			else {
				
				// 1st Constraint
				GRBLinExpr expr1 = new GRBLinExpr ();
				expr1.addTerm(1.0 , deltaX ); 
				expr1.addTerm(-1.0 , dx_j );
				expr1.addTerm(1.0, dx_i);
				double rhs1 = (1-s)*(p_j.getX()-p_i.getX());
				model . addConstr ( expr1 , GRB.EQUAL , rhs1 , "");
				
				// 2nd Constraint
				GRBLinExpr expr2 = new GRBLinExpr ();
				expr2.addTerm(1.0 , deltaY ); 
				expr2.addTerm(-1.0 , dy_j );
				expr2.addTerm(1.0, dy_i);
				double rhs2 = (1-s)*(p_j.getY()-p_i.getY());
				model . addConstr ( expr2 , GRB.EQUAL , rhs2 , "");
			}
			
		}
		
		
		
		
		/**
		 * Adds the additional flow constraints to enforce a connected road network.
		 */
		public void addConstraintsConnectivity(GRBModel model) throws GRBException {
			double n = (double)  rn.getRoads().size();
			
			// 1. Constraint: There is exactly one sink
			GRBLinExpr expr1 = new GRBLinExpr();
			
			for(Road r : rn.getRoads()) {
				
				expr1.addTerm(1.0,r.getIsSink());
				
				// 2nd Constraint: Only a selected road may be a sink
				GRBLinExpr expr2 = new GRBLinExpr ();
				expr2.addTerm(1.0,r.getZ());
				expr2.addTerm(-1.0,r.getIsSink());
				model.addConstr(expr2, GRB.GREATER_EQUAL, 0, "");
				
				// 3rd Constraint: Each source contributes a positive amount of flow.
				// 4th Constraint: No flow abundance
				GRBLinExpr expr3 = new GRBLinExpr();
				GRBLinExpr expr4 = new GRBLinExpr();
				for(GRBVar f_r_rOther : r.getOutFlow()) {
					expr3.addTerm(1.0,f_r_rOther);
					expr4.addTerm(1.0,f_r_rOther);
				}
				for(GRBVar f_rOther_r : r.getInFlow()) {
					expr3.addTerm(-1.0,f_rOther_r);
				}
				expr3.addTerm(n, r.getIsSink());
				expr3.addTerm(-1.0, r.getZ());
				expr4.addTerm(1.0-n, r.getZ());
				model.addConstr(expr3, GRB.GREATER_EQUAL, 0, "");
				model.addConstr(expr4, GRB.LESS_EQUAL, 0, "");
			}
			model.addConstr(expr1, GRB.EQUAL, 1, "");
		}
		
		

		/**
		 * Solving in several iterations (either the exact or the heuristic method), specified by numIter. 
		 * Each iteration takes as input the roads and buildings which have been selected in the previous iteration.
		 * The Delaunay triangulation is recomputed in each iteration. Within this method, you can specify 
		 * by removing the comments if the method uses the original coordinates of the points or the 
		 * coordinates computed in the previous iteration.
		 */
		public void solveWithIterations(int numIter, boolean weightObjects, boolean withSelection,
				boolean stayWithinMap, boolean downweightLongBottleneckEdges,
				boolean enforceConnectivityRoads, boolean buildingOnlyIfRoad,
				boolean coupledSelectionBuildings, boolean relaxation) throws GRBException {
			int i = 0;
			while(i < numIter) {
				System.out.println("\n----------------------------------------- "
						+ "Iteration " + String.valueOf(i+1) + " ------------------------------------------");
				if(!relaxation){
					solve(weightObjects, withSelection, stayWithinMap, 
							downweightLongBottleneckEdges, 
							enforceConnectivityRoads, buildingOnlyIfRoad, 
							coupledSelectionBuildings);
				}
				else{
					solveViaRelaxationSimpleThreshold(weightObjects, stayWithinMap, 
							downweightLongBottleneckEdges,enforceConnectivityRoads, buildingOnlyIfRoad, 
							coupledSelectionBuildings);									
				}
			
				if(bottleneckEdges.size()>0) {
					buildBottleneckGraph(); // ONLY RELEVANT FOR THE VISUALISATION!!!
				}
				
				removeBottleneckEdges();
				
				// ----- Case 1: Subsequent iterations shall use the new coordinates.
				// updateNetworks();   
				// -------------------------------------------------------
				
				// ----- Case 2: Subsequent iterations shall still use the original coordinates.
				if(i==numIter-1) {  
					updateNetworks();  
				}
				else {
					for(Road r : rn.getRoads()) {  
						for(RoadPoint rp : r.getVertices()) {
							rp.setDxDouble(0);
							rp.setDyDouble(0);
						}
					}
					for(Building b : bn.getBuildings()) {
						for(BuildingPoint bp : b.getShell()) {
							bp.setDxDouble(0);
							bp.setDyDouble(0);
						}
						for(Vector<BuildingPoint> hole : b.getHoles()) {
							for(BuildingPoint bp : hole) {
								bp.setDxDouble(0);
								bp.setDyDouble(0);
							}
						}
					}
				}
				// ---------------------------------------------------------
				
				RoadNetwork rnNew = new RoadNetwork();
				BuildingNetwork bnNew = new BuildingNetwork();
				gn.removeDelaunayPoints();
				for(Road r : rn.getRoads()) {
					r.removeFlowVariables();
					if(r.isSelected()) {
						rnNew.add(r);
					}
					else {
						for(RoadPoint rp : r.getVertices()) {
							for(int j=0; j<rp.getRoads().size(); j++) {
								if(rp.getRoads().get(j) == r) {
									rp.getRoads().remove(j);
									break;
								}
							}
						}
					}
				}
				for(Building b : bn.getBuildings()) {
					//b.removeFlowVariables();
					if(b.isSelected()) {
						bnNew.add(b);
					}
					else {
						for(BuildingPoint bp : b.getShell()) {
							for(int j=0; j<bp.getBuildings().size(); j++) {
								if(bp.getBuildings().get(j) == b) {
									bp.getBuildings().remove(j);
									break;
								}
							}
						}
						for(Vector<BuildingPoint> hole : b.getHoles()) {
							for(BuildingPoint bp : hole) {
								for(int j=0; j<bp.getBuildings().size(); j++) {
									if(bp.getBuildings().get(j) == b) {
										bp.getBuildings().remove(j);
										break;
									}
								}
							}
						}
					}
				}
				GeneralNetwork gnNew = new GeneralNetwork(rnNew, bnNew);
				this.rn = rnNew;
				this.bn = bnNew;
				this.gn = gnNew;
				G = gn.getGraph();
				if(i<numIter-1) {
					gn.computeDelaunayGraph();
					T = gn.getDelaunayGraph();
					G_prime = gn.getGraph();
					addBottleneckEdges();
					removeUnusedSteinerPoints();
					double minBottleneckLength = Double.MAX_VALUE;
					for(Segment se : bottleneckEdges) {
						double length = se.getStart().getCoordinate().distance(se.getEnd().getCoordinate());
						if(length<minBottleneckLength)  minBottleneckLength = length;
					}
					if(minBottleneckLength>=eps) break;
					System.out.println("Shortest bottleneck edge after iteration " + 
							String.valueOf(i+1) + ": " + minBottleneckLength);	
				}
				i++;
			}
		}


		/**
		 * Solves the relaxed program as part of the heuristic algorithm.
		 * @param weightObjects if true: objects obtain individual weights (here: proportional to length/area)
		 * @param withSelection if true: selection is possible, if false: only displacement
		 * @param stayWithinMap if true: the objects have to remain within the bounding box of the input points
		 * @param downweightLongBottleneckEdges if true: bottleneck edges longer than epsilon obtain smaller weights (see method downweightFactor) 
		 * @param enforceConnectivityRoads if false: road network may be split into several components, if true: road network remains connected
		 * @param buildingOnlyIfRoad if true: additional constraints for the coupled selection of buildings and roads are applied
		 * @param coupledSelectionBuildings if true: additional soft constraint that adjacent buildings should be treated the same (in a group of >= 3 buildings)
		 */
		public void solveViaRelaxation(boolean weightObjects, boolean withSelection, boolean stayWithinMap, 
				boolean downweightLongBottleneckEdges, boolean buildingOnlyIfRoad, 
				boolean coupledSelectionBuildings) {

			if(bottleneckEdges.size()==0) {
				addBottleneckEdges();
				removeUnusedSteinerPoints();
			}
			if(bottleneckEdges.size() > 0) {
				try {
					GRBEnv env = new GRBEnv ("GeneralNetwork.log ");
					GRBModel model = new GRBModel ( env );

					// Define the objective
					GRBQuadExpr obj = new GRBQuadExpr ();

					// Define the variables for the graph nodes
					for(Point p : G_prime.vertexSet()) {
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

					if(withSelection) {
						// Define the variables for the roads ...		
						SimpleWeightedGraph<Road,RoadPair> graphOfRoads = rn.getGraphOfRoads();
						for(Road r : graphOfRoads.vertexSet()) {
							GRBVar z;
							if(!r.isFixed()) {
								z = model.addVar(0.0 , 1.0 , 0.0 , GRB.CONTINUOUS , "");
							}
							else {
								if(r.isSelected()) {
									z = model.addVar(1.0 , 1.0 , 0.0 , GRB.CONTINUOUS , "");
								}
								else {
									z = model.addVar(0.0 , 0.0 , 0.0 , GRB.CONTINUOUS , "");
								}
							}
							if(!weightObjects) {  // Case 1: all roads obtain the same weight
								obj.addTerm(-w_select,z);   
							}
							else {   // Case 2 (standard): the roads obtain individual weights
								obj.addTerm(-w_select*r.getWeight(),z);
							}
							r.setZ(z);
						}
						// ... and for the buildings
						for(Building b : bn.getBuildings()) {
							GRBVar z = null;
							if(!b.isFixed()) {
								z = model.addVar(0.0 , 1.0 , 0.0 , GRB.CONTINUOUS , "");
							}
							else {
								if(b.isSelected()) {
									z = model.addVar(1.0 , 1.0 , 0.0 , GRB.CONTINUOUS , "");
								}
								else {
									z = model.addVar(0.0 , 0.0 , 0.0 , GRB.CONTINUOUS , "");
								}
							}
							if(!weightObjects) {  // Case 1: all buildings obtain the same weight
								obj.addTerm(-w_select,z);   
							}
							else {   // Case 2 (standard): individual weights for the buildings
								obj.addTerm(-w_select*b.getWeight(),z);
							}
							b.setZ(z);

							if(buildingOnlyIfRoad) {
								Road r = gn.getClosestRoadOfBuilding(b);
								GRBLinExpr expr = new GRBLinExpr ();
								expr.addTerm(1.0 , r.getZ() ); 
								expr.addTerm(-1.0 , z );
								model . addConstr ( expr , GRB.GREATER_EQUAL , 0 , "");
							}
						}
						if(coupledSelectionBuildings) {
							BuildingGraph graphOfBuildings = bn.getGraphOfBuildings();
							for(BuildingPair bPair : graphOfBuildings.getEdges()) {
								Building b1 = bPair.getB1();
								Building b2 = bPair.getB2();
								if(b1.getComponent().size()>2) {
									GRBVar z_b1_b2 = model.addVar(0.0 , 1.0 , 0.0 , GRB.CONTINUOUS , "");
									obj.addTerm(w_depend, z_b1_b2);
									GRBLinExpr expr1 = new GRBLinExpr (); // z_b1_b2 >= z_b1 - z_b2
									expr1.addTerm(1.0, z_b1_b2);
									expr1.addTerm(-1.0, b1.getZ());
									expr1.addTerm(1.0, b2.getZ());
									model.addConstr(expr1, GRB.GREATER_EQUAL, 0, "");
									GRBLinExpr expr2 = new GRBLinExpr (); // z_b1_b2 >= z_b2 - z_b1
									expr2.addTerm(1.0, z_b1_b2);
									expr2.addTerm(1.0, b1.getZ());
									expr2.addTerm(-1.0, b2.getZ());
									model.addConstr(expr2, GRB.GREATER_EQUAL, 0, "");
								}
							}
						}
						
						// if the node belongs to > 1 objects and is incident to a bottleneck edge,
						// define the collector
						for(Point p : G_prime.vertexSet()) {
							if(p.isIncidentToBottleneckEdge()) {
								if(p.getClass().toString().equals("class BuildingDisplacement.BuildingPoint")) {
									BuildingPoint bp = (BuildingPoint) p;
									if(bp.getBuildings().size()>1) {
										Collector c = new Collector(p);
										bp.setCollector(c);
										GRBVar z_c = model.addVar(0.0 , 1.0 , 0.0 , GRB.CONTINUOUS , "");
										c.setZ(z_c);
										for(Building b : bp.getBuildings()) {
											// Constraint z_Building <= z_c (select building only if collector is selected too)
											GRBLinExpr exprC = new GRBLinExpr ();
											exprC.addTerm(1.0 , b.getZ() );
											exprC.addTerm(-1.0 , z_c );
											model.addConstr(exprC, GRB.LESS_EQUAL, 0, "");
										}
									}
								}
								else {
									RoadPoint rp = (RoadPoint) p;
									if(rp.getRoads().size()>1) {
										Collector c = new Collector(p);
										rp.setCollector(c);
										GRBVar z_c = model.addVar(0.0 , 1.0 , 0.0 , GRB.CONTINUOUS , "");
										c.setZ(z_c);
										for(Road r : rp.getRoads()) {
											// Constraint z_Road <= z_c (select road only if collector is selected too)
											GRBLinExpr exprC = new GRBLinExpr ();
											exprC.addTerm(1.0 , r.getZ() );
											exprC.addTerm(-1.0 , z_c );
											model.addConstr(exprC, GRB.LESS_EQUAL, 0, "");
										}
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
						if(downweightLongBottleneckEdges) {
							double alpha = downweightFactor(p_i, p_j);
							obj.addTerm(alpha*w_edge, deltaX, deltaX );
							obj.addTerm(alpha*w_edge, deltaY, deltaY );
						}
						else {
							obj.addTerm(w_edge, deltaX, deltaX );
							obj.addTerm(w_edge, deltaY, deltaY );
						}
						se.setDeltaX(deltaX);
						se.setDeltaY(deltaY);

						GRBVar dx_i = p_i.getDx();
						GRBVar dx_j = p_j.getDx();
						GRBVar dy_i = p_i.getDy();
						GRBVar dy_j = p_j.getDy();
						double d_ij = p_i.getDist(p_j);
						double s = getS(d_ij);

						if(withSelection) {
							if(p_i.getClass().toString().equals("class RoadDisplacement.RoadPoint")
									&& p_j.getClass().toString().equals("class RoadDisplacement.RoadPoint")) {
								GRBVar z_i = ((RoadPoint) p_i).getRoads().get(0).getZ();
								if(((RoadPoint) p_i).getRoads().size()>1) {
									z_i = p_i.getCollector().getZ();
								}
								GRBVar z_j = ((RoadPoint) p_j).getRoads().get(0).getZ();
								if(((RoadPoint) p_j).getRoads().size()>1) {
									z_j = p_j.getCollector().getZ();
								}
								addConstraints(model, deltaX, deltaY, p_i,
										p_j, dx_i, dy_i, dx_j, dy_j, z_i, z_j, s, true);
							}
							else if(p_i.getClass().toString().equals("class RoadDisplacement.RoadPoint")
									&& p_j.getClass().toString().equals("class BuildingDisplacement.BuildingPoint")) {
								GRBVar z_i = ((RoadPoint) p_i).getRoads().get(0).getZ();
								if(((RoadPoint) p_i).getRoads().size()>1) {
									z_i = p_i.getCollector().getZ();
								}
								GRBVar z_j = ((BuildingPoint) p_j).getBuildings().get(0).getZ();
								if(((BuildingPoint) p_j).getBuildings().size()>1) {
									z_j = p_j.getCollector().getZ();
								}
								addConstraints(model, deltaX, deltaY, p_i,
										p_j, dx_i, dy_i, dx_j, dy_j, z_i, z_j, s, true);
							}
							else if(p_i.getClass().toString().equals("class BuildingDisplacement.BuildingPoint")
									&& p_j.getClass().toString().equals("class RoadDisplacement.RoadPoint")) {
								GRBVar z_i = ((BuildingPoint) p_i).getBuildings().get(0).getZ();
								if(((BuildingPoint) p_i).getBuildings().size()>1) {
									z_i = p_i.getCollector().getZ();
								}
								GRBVar z_j = ((RoadPoint) p_j).getRoads().get(0).getZ();
								if(((RoadPoint) p_j).getRoads().size()>1) {
									z_j = p_j.getCollector().getZ();
								}
								addConstraints(model, deltaX, deltaY, p_i,
										p_j, dx_i, dy_i, dx_j, dy_j, z_i, z_j, s, true);
							}
							else if(p_i.getClass().toString().equals("class BuildingDisplacement.BuildingPoint")
									&& p_j.getClass().toString().equals("class BuildingDisplacement.BuildingPoint")) {
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
								System.out.println("Fehler in solveWithBottleneckEdges() in Klasse Solver, "
										+ "ein Punkt ist nicht klassifizierbar!");
							}
						}
						else {
							addConstraints(model, deltaX, deltaY, p_i,
										p_j, dx_i, dy_i, dx_j, dy_j, null, null, s, false);
						}	
				}
					
					// Define the distortion variables also for the remaining edges of the enriched graph
					// (with s=1) to keep their lengths as much unchanged as possible.
					for(Segment se : G_prime.edgeSet()) {
						if(!se.isBottleneckEdge()) {
							GRBVar deltaX = model.addVar(-M , M , 0.0 , GRB.CONTINUOUS , "");
							GRBVar deltaY = model.addVar(-M , M , 0.0 , GRB.CONTINUOUS , "");
							obj.addTerm(w_edge, deltaX, deltaX );  
							obj.addTerm(w_edge, deltaY, deltaY );
							se.setDeltaX(deltaX);
							se.setDeltaY(deltaY);
							
							Point p_i = se.getStart();
							Point p_j = se.getEnd();
							GRBVar dx_i = p_i.getDx();
							GRBVar dx_j = p_j.getDx();
							GRBVar dy_i = p_i.getDy();
							GRBVar dy_j = p_j.getDy();
							if(dx_i==null || dx_j==null || dy_i==null || dy_j==null) {
								System.out.println("p_i: " + p_i.toString());
								System.out.println("p_j: " + p_j.toString());
							}
							
							if(withSelection) {
								if(p_i.getClass().toString().equals("class RoadDisplacement.RoadPoint")
										&& p_j.getClass().toString().equals("class RoadDisplacement.RoadPoint")) {
									Road road = RoadPoint.getCommonRoad((RoadPoint) p_i, (RoadPoint) p_j);
									addConstraints(model, deltaX, deltaY, p_i,
											p_j, dx_i, dy_i, dx_j, dy_j, road.getZ(), road.getZ(), 1, true);
								}
								else if(p_i.getClass().toString().equals("class BuildingDisplacement.BuildingPoint")
										&& p_j.getClass().toString().equals("class BuildingDisplacement.BuildingPoint")) {
									LinkedList<Building> buildings = BuildingPoint.getCommonBuildings((BuildingPoint) p_i, (BuildingPoint) p_j);
									if(buildings.size()==1) {
										addConstraints(model, deltaX, deltaY, p_i,
												p_j, dx_i, dy_i, dx_j, dy_j, buildings.get(0).getZ(), buildings.get(0).getZ(), 1, true);
									}
									else {
										Collector c = new Collector(se);
										se.setCollector(c);
										GRBVar z_co = model.addVar(0.0 , 1.0 , 0.0 , GRB.CONTINUOUS , "");
										c.setZ(z_co);
										for(Building b : buildings) {
											// Constraint z_Building <= z_c (select building only if collector is selected too)
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
									System.out.println("Error in method solveViaRelaxation in class Solver, points do not belong to the same class!");
								}
							}
							else {
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
					if(withSelection) {
						for(Building b : bn.getBuildings()) {
							if(b.getZ()!=null) {
								b.setzLP(b.getZ().get(GRB.DoubleAttr.X));
							}
						}
						for(Road r : rn.getRoads()) {
							if(r.getZ()!=null) {
								r.setzLP(r.getZ().get(GRB.DoubleAttr.X));
							}
						}
					}
					
					model.dispose();
					
				}catch ( GRBException e ) {
					e.printStackTrace();
					System . out . println (" Error code : " + e . getErrorCode () + ". " +
							e . getMessage ());
				}
			}
			else {
				// constraints are already met => no optimization needed
				for(Point p : G_prime.vertexSet()) {
					p.setDxDouble(0);
					p.setDyDouble(0);
				}	
			}
			
		}

		
		
		/**
		 * Implementation of the heuristic method.
		 * Important: If you call this method directly (not indirectly within solveWithIterations),
		 * you have to call solver.updateNetworks() afterwards to associate the points with the 
		 * computed coordinates!
		 */
		public void solveViaRelaxationSimpleThreshold(boolean weightObjects,
				boolean stayWithinMap, boolean downweightLongBottleneckEdges,
				boolean enforceConnectivityRoads, boolean buildingOnlyIfRoad,
				boolean coupledSelectionBuildings) throws GRBException {
			solveViaRelaxation(weightObjects, true, stayWithinMap, 
					downweightLongBottleneckEdges, buildingOnlyIfRoad, coupledSelectionBuildings);  
			for(Building b : bn.getBuildings()) {
				b.setSelected(b.getzLP()>=theta);
				if(b.isSelected())	b.setFixed(true);
			}
			for(Road r : rn.getRoads()) {
				r.setSelected(r.getzLP()>=theta);
				if(r.isSelected())	r.setFixed(true);
			}
			
			if(enforceConnectivityRoads) { 
				TreeMap<Double,Road> remaining_roads = new TreeMap<Double,Road>();
				ArrayList<Road> selected_roads = new ArrayList<Road>();
				for(Road r : rn.getRoads()) {
					if(r.isSelected()) {
						selected_roads.add(r);
					}
					else {
						remaining_roads.put(r.getzLP(), r);
					}
					r.setComponent(null);
				}
				RoadGraph graphOfRoads = new RoadGraph();
				graphOfRoads.addVerticesAndEdgesBetween(selected_roads);
				graphOfRoads.depthFirstSearch();
				System.out.println("# nodes in the graph: " + graphOfRoads.getVertices().size());
				System.out.println("road network is connected: " + graphOfRoads.isConnected());
				while(remaining_roads.size()>0 && !graphOfRoads.isConnected()) {
					Map.Entry<Double, Road> entry = remaining_roads.pollLastEntry();
					Road rCandidate = entry.getValue();
					boolean added = graphOfRoads.addVertexIfAdmissibleAndUpdateComponents(
																rCandidate);
					if(added) {
						rCandidate.setSelected(true);
						rCandidate.setFixed(true);
					}
				}
				System.out.println("road network is connected: " + graphOfRoads.isConnected());
			}

			RoadNetwork rnNew = new RoadNetwork();
			BuildingNetwork bnNew = new BuildingNetwork();
			for(Road r : rn.getRoads()) {
				r.removeFlowVariables();
				if(r.isSelected()) rnNew.add(r);
				else {
					for(RoadPoint rp : r.getVertices()) {
						for(int j=0; j<rp.getRoads().size(); j++) {
							if(rp.getRoads().get(j) == r) {
								rp.getRoads().remove(j);
								break;
							}
						}
					}
				}
			}
			for(Building b : bn.getBuildings()) {
				//b.removeFlowVariables();
				if(b.isSelected()) bnNew.add(b);
				else {
					for(BuildingPoint bp : b.getShell()) {
						for(int j=0; j<bp.getBuildings().size(); j++) {
							if(bp.getBuildings().get(j) == b) {
								bp.getBuildings().remove(j);
								break;
							}
						}
					}
					for(Vector<BuildingPoint> hole : b.getHoles()) {
						for(BuildingPoint bp : hole) {
							for(int j=0; j<bp.getBuildings().size(); j++) {
								if(bp.getBuildings().get(j) == b) {
									bp.getBuildings().remove(j);
									break;
								}
							}
						}
					}
				}
			}
			GeneralNetwork gnNew = new GeneralNetwork(rnNew, bnNew);
			this.rn = rnNew;
			this.bn = bnNew;
			this.gn = gnNew;
			Object[] array_of_points = G_prime.vertexSet().toArray();
			for(int i=0; i<array_of_points.length; i++) {
				Point p = (Point) array_of_points[i];
				if(! p.isVisible()) {
					G_prime.removeVertex(p);
					G.removeVertex(p);
				}
			}
			LinkedList<Segment> edgesToRemove = new LinkedList<Segment>();
			for(Segment se : bottleneckEdges) {
				if(!se.getStart().isVisible() || !se.getEnd().isVisible()) {
					edgesToRemove.add(se);
				}
			}
			for(Segment se : edgesToRemove) {
				bottleneckEdges.remove(se);
				G_prime.removeEdge(se);
				se.getStart().checkIfIsIncidentToBottleneckEdge(G_prime);
				se.getEnd().checkIfIsIncidentToBottleneckEdge(G_prime);
			}

			solveViaRelaxation(weightObjects, false, stayWithinMap, downweightLongBottleneckEdges, 
						buildingOnlyIfRoad, coupledSelectionBuildings);
			
			buildBottleneckGraph();
			
			removeBottleneckEdges();
			
			//updateNetworks();
			gn.removeDelaunayPointsWithoutChecking();
			
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

		
		/**
		 * Updates the coordinates of the buildings and roads according to the solution.
		 */
		public void updateNetworks() {
			
			//-------- NUR RELEVANT FUER DIE VISUALISIERUNG DER BOTTLENECK
			//-------- EDGES NACH MEHREREN ITERATIONEN 
//			for(Road r : rn.getRoads()) {
//				for(RoadPoint rp : r.getVertices()) {
//					rp.setIntermediateCoord(rp.getCoordinate());
//				}
//			}
//			for(Building b : bn.getBuildings()) {
//				for(BuildingPoint bp : b.getShell()) {
//					bp.setIntermediateCoord(bp.getCoordinate());
//				}
//				for(Vector<BuildingPoint> hole : b.getHoles()) {
//					for(BuildingPoint bp : hole) {
//						bp.setIntermediateCoord(bp.getCoordinate());;
//					}
//				}
//			}
			//------------------------------------------------
			
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
			
			for (Segment se : T.edgeSet()) {
				if( !G_prime.containsEdge(se) ) {
					listAll.add(se);
				}
			}
			
			Comparator<Segment> c = new Comparator<Segment>() {
				@Override
				public int compare(Segment s1, Segment s2) { 
					Double d1 = T.getEdgeWeight(s1);
					Double d2 = T.getEdgeWeight(s2);
					return d1.compareTo(d2); }
			};
			Collections.sort(listAll, c);

			for ( Segment se : listAll) {
				Point s = T.getEdgeSource(se);
				Point u = T.getEdgeTarget(se);
				double edgeWeight = T.getEdgeWeight(se);
				double radius = edgeWeight*t;
				DijkstraShortestPath<Point,Segment> dsp = new DijkstraShortestPath<Point,Segment>(G_prime,radius);
				double pathWeight = dsp.getPathWeight(s,u);
				double dilation = pathWeight / edgeWeight;
				if (dilation > t) {
					G_prime.addEdge(s,u,se);
					G_prime.setEdgeWeight(se, T.getEdgeWeight(se));
					bottleneckEdges.add(se);
					se.setBottleneckEdge(true);

					se.getStart().setIncidentToBottleneckEdge(true);
					se.getEnd().setIncidentToBottleneckEdge(true);
				}
			}
			System.out.println("Number of bottleneck edges added: " + bottleneckEdges.size());
			System.out.println("Number of nodes in augmented graph: " + G_prime.vertexSet().size());
			System.out.println("Number of edges in augmented graph: " + G_prime.edgeSet().size());
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
					System.out.println("Problem in getS(), class Solver: point distance is very small!");
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


		/**
		 * Sets up the graph consisting of the bottleneck edges as edges.
		 */
		public void buildBottleneckGraph(){
			bottleneckGraph = new SimpleWeightedGraph<Point, Segment>(Segment.class);
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
		}
		
		/**
		 * Returns the graph consisting of the bottleneck edges as edges.
		 */
		public SimpleWeightedGraph<Point,Segment> getBottleneckGraph(){
			if(bottleneckGraph==null) {
				buildBottleneckGraph();
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
				Set<Segment> bottleneckEdges,List<Building> buildings, 
				Set<Road> roads,boolean downweightLongBEdges,
        		boolean coupledSelectionBuildings, LinkedList<BuildingPair> buildingAdjacencies, 
				double weightConnBuildings) {
			System.out.println("\n----------------------- Evaluation "
					+ "-------------------------");
			System.out.println("# points: " + gWithDel.vertexSet().size());
			System.out.println("# original edges: " + gWithDel.edgeSet().size());
			System.out.println("# initial bottleneck edges: " + bottleneckEdges.size());
			System.out.println("# buildings: " + buildings.size());
			System.out.println("# roads: " + roads.size());
			double obj = 0;
			double term = 0;
			for(Point p : gWithDel.vertexSet()) {
				if(p.isVisible()) {
					double dx = p.getX()+p.getDxDouble()-p.getOriginalCoord().getX();
					double dy = p.getY()+p.getDyDouble()-p.getOriginalCoord().getY();
					obj = obj + w_pos*(dx*dx+dy*dy);
					term = term + w_pos*(dx*dx+dy*dy);
				}
			}
			System.out.println("Contribution of point movement: " + term);
			term = 0;
			for(Segment se : gWithDel.edgeSet()) {
				if(se.isVisible()) {
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
					obj = obj + w_edge*w_BE*(deltaX*deltaX+deltaY*deltaY);
					term = term + w_edge*w_BE*(deltaX*deltaX+deltaY*deltaY);
				}
			}
			System.out.println("Contribution of the initial bottleneck edges: " + term);
			term = 0;
			double tension_measure = 0;
			double weighted_tension_measure = 0;
			double sum_of_weights = 0;
			for(Building b : buildings) {
				sum_of_weights = sum_of_weights + b.getWeight();
				if(!b.isSelected()) {
					obj = obj + w_select*b.getWeight();
					term = term + w_select*b.getWeight();
					tension_measure = tension_measure + 1.0;
					weighted_tension_measure = weighted_tension_measure + b.getWeight();
				}
			}
			System.out.println("Contribution of buildings removed: " + term);
			term = 0;
			for(Road r : roads) {
				sum_of_weights = sum_of_weights + r.getWeight();
				if(!r.isSelected()) {
					obj = obj + w_select*r.getWeight();
					term = term + w_select*r.getWeight();
					tension_measure = tension_measure + 1.0;
					weighted_tension_measure = weighted_tension_measure + r.getWeight();
				}
			}
			System.out.println("Contribution of roads removed: " + term);
			term = 0;
			if(coupledSelectionBuildings) {
				for(BuildingPair bp : buildingAdjacencies) {
					if((bp.getB1().isSelected() && !bp.getB2().isSelected())
						|| (bp.getB2().isSelected() && !bp.getB1().isSelected())) {
						obj = obj + weightConnBuildings;
						term = term + weightConnBuildings;
					}
				}
			}
			System.out.println("Contribution of building connectivity constraints: " + term);
			tension_measure = 100* tension_measure/((double) buildings.size() + (double) roads.size());
			weighted_tension_measure = 100*weighted_tension_measure/sum_of_weights;
			System.out.println("Weighted Percentage of removed objects [%]: " + weighted_tension_measure);
			return obj;
		}

		
		
		/**
		 * Removes the Steiner points which are not incident to bottleneck edges.
		 */
		public void removeUnusedSteinerPoints() {
			System.out.println("Removing unused Steiner points...");
			int counter = 0;
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
							counter++;
							r.getVertices().remove(rp);
							rp.getRoads().clear();
							i--;
						}
						else {
							System.out.println("Warning in method removeUnusedSteinerPoints in class Solver, a Steiner point could not be removed.");
						}
					}
				}
			}
			for(Building b : bn.getBuildings()) {
				Vector<BuildingPoint> shell = b.getShell();
				for(int i=0; i<shell.size(); i++) {
					BuildingPoint bp = shell.get(i);
					if(bp.isDelaunay() && !bp.isIncidentToBottleneckEdge()) {
						if(G_prime.containsVertex(bp)) {
							Object[] edges = G_prime.edgesOf(bp).toArray();
							if(edges.length==2) {
								Segment s1 = (Segment) edges[0];
								Segment s2 = (Segment) edges[1];
								BuildingPoint bp1 = (BuildingPoint) s1.getOtherPoint(bp);
								BuildingPoint bp2 = (BuildingPoint) s2.getOtherPoint(bp);
								G_prime.addEdge(bp1, bp2, new Segment(bp1, bp2));
								G_prime.removeVertex(bp);
								
								for(Building building : bp.getBuildings()) {
									building.getShell().remove(bp);
								}
								//shell.remove(bp);
								
								counter++;
								bp.getBuildings().clear();
								i--;
							}
							else {
								System.out.println("Warnung in Methode removeUnusedSteinerPoints in Klasse Solver, ein Steiner-Punkt konnte nicht entfernt werden.");
							}
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
								counter++;
								bp.getBuildings().clear();
								i--;
							}
							else {
								System.out.println("Error in method removeUnusedSteinerPoints in class Solver!!!");
							}
						}
					}
				}
			}
			System.out.println("# unused Steiner points removed: " + counter);
		}


		public LinkedList<Segment> getBottleneckEdges() {
			return bottleneckEdges;
		}
		
		

}

