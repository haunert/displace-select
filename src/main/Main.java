package main;
import java.awt.Color;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.Vector;

import javax.swing.JFrame;

import org.apache.batik.dom.GenericDOMImplementation;
import org.apache.batik.svggen.SVGGraphics2D;
import org.jgrapht.graph.SimpleWeightedGraph;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Polygon;
import org.w3c.dom.DOMImplementation;
import org.w3c.dom.Document;

import BuildingDisplacement.*;
import GeneralDisplacement.GeneralNetwork;
import GeneralDisplacement.Solver;
import RoadDisplacement.*;
//import general.Graph;
import general.Point;
import general.Segment;
import gurobi.GRBException;
import io.structures.Feature;
import viewer.base.ListLayer;
import viewer.base.MapPanel;
import viewer.symbols.BasicSymbolFactory;
import visualization.MapFrame;
import visualization.OwnBasicSymbolFactory;
import visualization.OwnMapPanel;

public class Main {

	/**
	 * ------------ Parameters set by the user --------------------------------------
	 */
	static String buildingPath = "data\\geb-goetheallee.shp";
	static String roadPath = "data\\";
	static String roadName = "goetheallee";
	
	static double w_pos = 0.0001; // weight for point displacements 
	static double w_edge = 0.8;    // weight for edge distortions
	static double w_select = 1-w_pos-w_edge;   // weight for unselecting objects
	static double w_depend = 0.5;  // weight for treating adjacent buildings the same (in a group of >= 3 buildings), if coupledSelectionBuildings is true

	static double epsilon = 7.5;  // [m]: threshold for a bottleneck edge representing a conflict
	static double t = 5;  // threshold for the maximally allowed stretch factor of an edge
	
	static boolean weightObjects = true;  // if true: objects obtain individual weights (here: proportional to length/area)
	static boolean withSelection = true;   // if true: selection is possible, if false: only displacement
	static boolean stayWithinMap = true; // if true: the objects have to remain within the bounding box of the input points
	static boolean downweightLongBottleneckEdges = true; // if true: bottleneck edges longer than epsilon obtain smaller weights
	static boolean enforceConnectivityRoads = false; // if false: road network may be split into several components, if true: road network remains connected (exact method -> flow model, heuristic method -> postprocessing)
	static boolean buildingOnlyIfRoad = true; // if true: additional constraints for the coupled selection of buildings and roads are applied
	static boolean coupledSelectionBuildings = true; // if true: additional soft constraint that adjacent buildings should be treated the same (in a group of >= 3 buildings)
	static boolean heuristic = false;  // if false: exact method, if true: heuristic
	static double theta = 0.998;  // threshold for the selection of objects in the heuristic
	
	static int numIter = 1;   // number of iterations 
	static double road_factor = 10;  // specifies the relation between the weights of the roads and those of the buildings (shortest road has the weight road_factor, whereas the smallest building has the weight 1)
	public static final double eps = Math.pow(10,-4); // tolerance [m] for point queries within the program (kd tree, point to line etc.)
	public static double maxDistortion = 0.6;  // the value obtaining the darkest red in the distortion plot
	public static double maxDisplacement = 1.5; // the value obtaining the darkest red in the displacement plot
	
	/**
	 * ----------------------------------------------------------------------------------------------------
	 */
	
	// auxiliary variables for some iterations below
    static int counter = 0, counterDisplacementMap = 0;
    
    // auxiliary variables for measuring the computation time
    static double start, end, duration;
    static double startExcluded, endExcluded, durationExcluded=0; // for processes not to be included in the time measurement


	public static void main(String[] args) throws GRBException, IOException {

		start = System.nanoTime();
		System.out.println("Reading shapefiles ...");

		GeometryFactory gf = new GeometryFactory();
		ArrayList<Color> colorGradient = getColorGradient();
		OwnMapPanel myMapPanel = new OwnMapPanel();
		OwnMapPanel myMapPanelBefore = new OwnMapPanel();
		OwnMapPanel myMapPanelAfter = new OwnMapPanel();
		OwnMapPanel myMapPanelWithDelaunay = new OwnMapPanel();
		OwnMapPanel myMapPanelEdges = new OwnMapPanel(colorGradient, maxDistortion, "Distortion");
		OwnMapPanel myMapPanelDisplacements = new OwnMapPanel(colorGradient, maxDisplacement, "Displacement");

		BuildingNetwork bn = BuildingNetwork.importFromShapefile(buildingPath);
		RoadNetwork rn = RoadNetwork.importFromShapefile(roadPath, roadName);
		
		startExcluded = System.nanoTime();
		BuildingGraph bg = bn.getGraphOfBuildings();
		SimpleWeightedGraph<Road,RoadPair> roadGraph = rn.getGraphOfRoads();
		LinkedList<BuildingPair> buildingAdjacencies = new LinkedList<BuildingPair>();
		for(BuildingPair bp : bg.getEdges()) {
			if(bp.getB1().getComponent().size() >= 3) {
				buildingAdjacencies.add(bp);
			}
		}
		endExcluded = System.nanoTime();
		durationExcluded = durationExcluded + (endExcluded-startExcluded)/1000000000;
		
		GeneralNetwork gn = new GeneralNetwork(rn, bn);
		SimpleWeightedGraph<Point, Segment> G = gn.getGraph();
		
		startExcluded = System.nanoTime();
		ListLayer polygon_layer_before = new ListLayer(new BasicSymbolFactory(Color.BLACK,Color.GRAY));
		for (Building b : bn.getBuildings()) {
			Polygon poly = b.getPolygon();
			Feature f = new Feature(poly);
			polygon_layer_before.add(f);
		}
		myMapPanelBefore.getMap().addLayer(polygon_layer_before, 1);
		endExcluded = System.nanoTime();
		durationExcluded = durationExcluded + (endExcluded-startExcluded)/1000000000;
		ListLayer input_layer = new ListLayer(new BasicSymbolFactory(Color.LIGHT_GRAY,Color.LIGHT_GRAY,(float) 3));
		ListLayer input_layer_polygon = new ListLayer(new BasicSymbolFactory(Color.BLACK,Color.BLACK,(float) 2));
		for (Segment s : G.edgeSet()) {
			Coordinate[] coords = new Coordinate[2];
			coords[0] = s.getStart().getCoordinate();
			coords[1] = s.getEnd().getCoordinate();
			LineString ls = gf.createLineString(coords);
			Feature f = new Feature(ls);
			input_layer.add(f);
			input_layer_polygon.add(f);
		}     
		myMapPanel.getMap().addLayer(input_layer, 1);
		myMapPanelWithDelaunay.getMap().addLayer(input_layer,1);
		myMapPanelDisplacements.getMap().addLayer(input_layer, 1);
		myMapPanelBefore.getMap().addLayer(input_layer_polygon, 2);
		System.out.println("Reading shapefiles done." + "\n");
		
		gn.computeDelaunayGraph();
		SimpleWeightedGraph<Point, Segment> T = gn.getDelaunayGraph();
		SimpleWeightedGraph<Point, Segment> G_prime = gn.getGraph();
		System.out.println("Additional Steiner points added to the graph." + "\n");

		double M = gn.determineM();
		gn.normalizeWeights(road_factor);
		Solver solver = new Solver(G, G_prime, T, gn, epsilon, t, M, 
									w_pos, w_edge, w_select, w_depend, theta);
		solver.addBottleneckEdges();
		solver.removeUnusedSteinerPoints();
		SimpleWeightedGraph<Point, Segment> G_prime_copy = 
						new SimpleWeightedGraph<Point, Segment>(Segment.class);
		for(Point p : G_prime.vertexSet()) G_prime_copy.addVertex(p);
		for(Segment se : G_prime.edgeSet()) {
			if(!se.isBottleneckEdge())  G_prime_copy.addEdge(se.getStart(), se.getEnd(), se);
		}

		SimpleWeightedGraph<Point, Segment> bEdges_initial = solver.getBottleneckGraph();
		solver.solveWithIterations(numIter,weightObjects,withSelection,stayWithinMap,downweightLongBottleneckEdges, 
						enforceConnectivityRoads,buildingOnlyIfRoad,
						coupledSelectionBuildings,heuristic);
		
        if(numIter == 1) {
        	double obj = solver.evaluate(G_prime_copy,bEdges_initial.edgeSet(),bg.getVertices(), 
    				roadGraph.vertexSet(),downweightLongBottleneckEdges,coupledSelectionBuildings,
    				buildingAdjacencies, w_depend);
            System.out.println("Objective value: " + obj);
        }
        

		// solution graph as an additional layer (without the Steiner nodes)
        counter=0;
		SimpleWeightedGraph<Point, Segment> solGraph = solver.getSolutionGraph();      
		ListLayer sol_layer = new ListLayer(new BasicSymbolFactory(Color.BLACK,null,(float) 2));
		for (Segment s : solGraph.edgeSet()) {
			Coordinate[] coords = new Coordinate[2];
			coords[0] = s.getStart().getCoordinate();
			coords[1] = s.getEnd().getCoordinate();
			LineString ls = gf.createLineString(coords);
			Feature f = new Feature(ls);
			sol_layer.add(f);
			counter++;
		}
		myMapPanel.getMap().addLayer(sol_layer, 3);

		// bottleneck edges as an additional layer
		ListLayer bEdges_short_layer = new ListLayer(new OwnBasicSymbolFactory(Color.RED,null,(float) 2));
		ListLayer bEdges_long_layer = new ListLayer(new OwnBasicSymbolFactory(Color.BLUE,null,(float) 2));
		for (Segment s : bEdges_initial.edgeSet())  {     
			Coordinate[] coords = new Coordinate[2];
			coords[0] = s.getStart().getOriginalCoord();  // bottleneck edges displayed at the original coordinates
			coords[1] = s.getEnd().getOriginalCoord();
			LineString ls = gf.createLineString(coords);
			Feature f = new Feature(ls);
			if(coords[0].distance(coords[1])<epsilon) {
				bEdges_short_layer.add(f);
			}
			else {
				bEdges_long_layer.add(f);
			}
		}
		myMapPanel.getMap().addLayer(bEdges_short_layer, 4);
		myMapPanel.getMap().addLayer(bEdges_long_layer, 5);

		// solution graph with the additional Steiner nodes (not meaningful in case of several iterations,
		// since Steiner nodes are recomputed then)
		ListLayer sol_layer_del = new ListLayer(new BasicSymbolFactory(Color.BLACK,null,(float) 2));
		counter = 0;
		for (Segment s : G_prime_copy.edgeSet()) {  
			if(s.isVisible()) {
				Coordinate[] coords = new Coordinate[2];
				coords[0] = new Coordinate(s.getStart().getX()+s.getStart().getDxDouble(),
									s.getStart().getY()+s.getStart().getDyDouble());
				coords[1] = new Coordinate(s.getEnd().getX()+s.getEnd().getDxDouble(),
									s.getEnd().getY()+s.getEnd().getDyDouble());
				LineString ls = gf.createLineString(coords);
				Feature f = new Feature(ls);
				sol_layer_del.add(f);
				counter++;
			}
		}
		myMapPanelWithDelaunay.getMap().addLayer(sol_layer_del, 2);
		myMapPanelWithDelaunay.getMap().addLayer(bEdges_short_layer, 3);
		myMapPanelWithDelaunay.getMap().addLayer(bEdges_long_layer, 4);

		startExcluded = System.nanoTime();
		// color the edges according to their distortions
		double[] distortions = new double[G_prime_copy.edgeSet().size()];
		counter = 0;
		for(Segment s : G_prime_copy.edgeSet()) {
			if(s.isVisible()) {
				double dx_j = s.getEnd().getX()+s.getEnd().getDxDouble()-s.getEnd().getOriginalCoord().getX();
				double dx_i = s.getStart().getX()+s.getStart().getDxDouble()-s.getStart().getOriginalCoord().getX();
				double dy_j = s.getEnd().getY()+s.getEnd().getDyDouble()-s.getEnd().getOriginalCoord().getY();
				double dy_i = s.getStart().getY()+s.getStart().getDyDouble()-s.getStart().getOriginalCoord().getY();
				double delta_x_ij = Math.abs(dx_j-dx_i);
				double delta_y_ij = Math.abs(dy_j-dy_i);
				distortions[counter] = Math.sqrt(delta_x_ij*delta_x_ij+delta_y_ij*delta_y_ij);
				Coordinate[] coords = new Coordinate[2];
				coords[0] = s.getStart().getOriginalCoord();
				coords[1] = s.getEnd().getOriginalCoord();
				double dist_rel = Math.min(distortions[counter]/maxDistortion,1);
				int colorIdx = (int) Math.floor(dist_rel*(colorGradient.size()-1));
				Color color = colorGradient.get(colorIdx);
				ListLayer current_Edge = new ListLayer(new BasicSymbolFactory(color,null,3));
				LineString ls = gf.createLineString(coords);
				Feature f = new Feature(ls);
				current_Edge.add(f);
				myMapPanelEdges.getMap().addLayer(current_Edge, counter);
				
				// preparation for the plot with the point displacements
				ListLayer current_Edge_grey = new ListLayer(new BasicSymbolFactory(Color.LIGHT_GRAY,null,3));
				current_Edge_grey.add(f);
				myMapPanelDisplacements.getMap().addLayer(current_Edge_grey, counter);
				
				counter++;
			}
		}
		counterDisplacementMap = counter;
		
		// color the bottleneck edges according to their distortions
		for (Segment s : bEdges_initial.edgeSet())  { 
			Point start = s.getStart();
			Point end = s.getEnd();
			if(start.isVisible() && end.isVisible()) {
				double dx_ij_before = end.getOriginalCoord().getX()-start.getOriginalCoord().getX();
				double dy_ij_before = end.getOriginalCoord().getY()-start.getOriginalCoord().getY();
				double dx_ij_after = end.getX()+end.getDxDouble()-start.getX()-start.getDxDouble();
				double dy_ij_after = end.getY()+end.getDyDouble()-start.getY()-start.getDyDouble();
				double scale = solver.getS(start.getOriginalCoord().distance(end.getOriginalCoord()));
				double deltaX = dx_ij_after - scale*dx_ij_before;
				double deltaY = dy_ij_after - scale*dy_ij_before;
				double distortion = Math.sqrt(deltaX*deltaX+deltaY*deltaY);
				Coordinate[] coords = new Coordinate[2];
				coords[0] = start.getOriginalCoord();
				coords[1] = end.getOriginalCoord();
				double dist_rel = Math.min(distortion/maxDistortion,1);
				int colorIdx = (int) Math.floor(dist_rel*(colorGradient.size()-1));
				Color color = colorGradient.get(colorIdx);
				ListLayer current_Edge = new ListLayer(new OwnBasicSymbolFactory(color,null,(float) 3));
				LineString ls = gf.createLineString(coords);
				Feature f = new Feature(ls);
				current_Edge.add(f);
				myMapPanelEdges.getMap().addLayer(current_Edge, counter++);
			}
		}
		
		// display the point displacements
		double[] displacements = new double[G_prime_copy.vertexSet().size()];
		counter = 0;
		for(Point p : G_prime_copy.vertexSet()) {
			if(p.isVisible()) {
				double dx_i = p.getX()+p.getDxDouble()-p.getOriginalCoord().getX();
				double dy_i = p.getY()+p.getDyDouble()-p.getOriginalCoord().getY();
				displacements[counter] = Math.sqrt(dx_i*dx_i + dy_i*dy_i);
				double displ_rel = Math.min(displacements[counter]/maxDisplacement,1);
				int colorIdxDispl = (int) Math.floor(displ_rel*(colorGradient.size()-1));
				Color color = colorGradient.get(colorIdxDispl);
				Feature fPoint = new Feature(gf.createPoint(p.getOriginalCoord()));
				ListLayer ll_Point = new ListLayer(new OwnBasicSymbolFactory(color,null,10));
				ll_Point.add(fPoint);
				myMapPanelDisplacements.getMap().addLayer(ll_Point, counterDisplacementMap);
				counter++;
				counterDisplacementMap++;
			}
		}

		// display the buildings as filled polygons 
		gn.removeDelaunayPoints();
		//solver.updateNetworks();
		ListLayer polygon_layer = new ListLayer(new BasicSymbolFactory(Color.BLACK,Color.GRAY));
		for (Building b : bn.getBuildings()) {
			if(b.isSelected()) {
				Polygon poly = b.getPolygon();
				Feature f = new Feature(poly);
				polygon_layer.add(f);
			} 
		}
		myMapPanelAfter.getMap().addLayer(polygon_layer, 1);
		myMapPanelAfter.getMap().addLayer(sol_layer, 2);

		MapFrame myFrameBefore = new MapFrame("Before optimization",myMapPanelBefore);
		myMapPanelBefore.getMap().setFrameRatio(0.01);
		myFrameBefore.setVisible(true);
		myFrameBefore.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		endExcluded = System.nanoTime();
		durationExcluded = durationExcluded + (endExcluded-startExcluded)/1000000000;
		
		MapFrame myFrame = new MapFrame("Before and after optimization and bottleneck edges (red/blue)",myMapPanel);
		myMapPanel.getMap().setFrameRatio(0.01);
		myFrame.setVisible(true);
		myFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		
		MapFrame myFrameWithDelaunay = new MapFrame("With additional Steiner points",myMapPanelWithDelaunay);
		myMapPanelWithDelaunay.getMap().setFrameRatio(0.01);
		myFrameWithDelaunay.setVisible(true);
		myFrameWithDelaunay.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		
		startExcluded = System.nanoTime();
		MapFrame myFrameEdges = new MapFrame("Distortions of the remaining edges",myMapPanelEdges);
		myMapPanelEdges.getMap().setFrameRatio(0.01);
		myFrameEdges.setVisible(true);
		myFrameEdges.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		
		MapFrame myFrameDisplacements = new MapFrame("Displacements",myMapPanelDisplacements);
		myMapPanelDisplacements.getMap().setFrameRatio(0.01);
		myFrameDisplacements.setVisible(true);
		myFrameDisplacements.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

		MapFrame myFrameAfter = new MapFrame("After optimization",myMapPanelAfter);
		myMapPanelAfter.getMap().setFrameRatio(0.01);
		myFrameAfter.setVisible(true);
		myFrameAfter.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		
		endExcluded = System.nanoTime();
		durationExcluded = durationExcluded + (endExcluded-startExcluded)/1000000000;
		
		end = System.nanoTime();
		duration = (end-start)/1000000000;
		System.out.println("\ncomputation time [s]: " + String.valueOf(duration-durationExcluded));
		
		
		
		/**
		 * ------- Save as SVG -----------------------------
		 */
		// Get a DOMImplementation.
        DOMImplementation domImpl = GenericDOMImplementation.getDOMImplementation();

        // Create an instance of org.w3c.dom.Document.
        String svgNS = "http://www.w3.org/2000/svg";
        Document document = domImpl.createDocument(svgNS, "svg", null);

        // Create an instance of the SVG Generator.
        SVGGraphics2D svgGenerator = new SVGGraphics2D(document);

        // Ask the test to render into the SVG Graphics2D implementation.
        myMapPanelAfter.getMap().paint(svgGenerator);

        // Finally, stream out SVG to the standard output using
        // UTF-8 encoding.
        boolean useCSS = true; // we want to use CSS style attributes
        Writer out;
        try {
            out = new FileWriter("output.svg");
            svgGenerator.stream(out, useCSS);   
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println("\nSaved to svg.");
		
	}

	
	
	
	
	public static ArrayList<Color> getColorGradient(){
		ArrayList<Color> colorGradient = new ArrayList<Color>();
		colorGradient.add(new Color(192,192,192));
		colorGradient.add(new Color(255,152,152));
		colorGradient.add(new Color(255,0,0));
		colorGradient.add(new Color(161,0,0));
		return colorGradient;
	}

	public static double getMax(double[] values) {
		double max = -Double.MAX_VALUE;
		for(int i=0; i<values.length; i++) {
			if(values[i]>max) {
				max = values[i];
			}
		}
		return max;
	}


	public static double[] getNormalizedValues(double[] values) {
		double max = getMax(values);
		for(int i=0; i<values.length; i++) {
			values[i] = values[i]/max;
		}
		return values;
	}
	

}
