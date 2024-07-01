package main;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.geom.AffineTransform;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import javax.imageio.ImageIO;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

import org.apache.batik.dom.GenericDOMImplementation;
import org.apache.batik.svggen.SVGGraphics2D;
import org.jgrapht.graph.SimpleWeightedGraph;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.MultiLineString;
import org.locationtech.jts.geom.MultiPoint;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.index.kdtree.KdNode;
import org.locationtech.jts.index.kdtree.KdTree;
import org.w3c.dom.DOMImplementation;
import org.w3c.dom.Document;

import com.vividsolutions.jump.feature.AttributeType;
import com.vividsolutions.jump.feature.BasicFeature;
import com.vividsolutions.jump.feature.FeatureCollection;
import com.vividsolutions.jump.feature.FeatureDataset;
import com.vividsolutions.jump.feature.FeatureSchema;
import com.vividsolutions.jump.io.DriverProperties;
import com.vividsolutions.jump.io.ShapefileWriter;

import BuildingDisplacement.*;
import GeneralDisplacement.GeneralNetworkWareJones;
import GeneralDisplacement.Solver;
import GeneralDisplacement.SolverWareJones;
import RoadDisplacement.Road;
import RoadDisplacement.RoadNetwork;
import RoadDisplacement.RoadPair;
import RoadDisplacement.RoadPoint;
//import RoadDisplacement.SolverRoads;
//import general.Graph;
import general.Point;
import general.Segment;
import gurobi.GRB;
import gurobi.GRBException;
import io.shp.FeatureReader;
import io.structures.Feature;
import viewer.base.Layer;
import viewer.base.ListLayer;
import viewer.base.MapPanel;
import viewer.symbols.BasicSymbolFactory;
import viewer.symbols.Symbol;
import viewer.symbols.SymbolFactory;
import visualization.MapFrame;
import visualization.OwnBasicSymbolFactory;
import visualization.OwnMapPanel;

/**
 * Main class for the comparison with the method by Ware et al. (2003): Automated map generalization with 
 * multiple operators: a simulated annealing approach.
 * - roads are not unselected
 * - buildings are not distorted (and roads only if the variable withDisplacementRoads is set to true below)
 * - buildings are all isolated
 * - we use our exact algorithm (not the heuristic)
 */
public class MainWareJones {

	/**
	 * ------------------- Parameters set by the user --------------------------------
	 */
	static String buildingPath = "data\\Ware and Jones\\buildings_input_labelled_fig6d.shp";
	static String roadPath = "data\\Ware and Jones\\";
	static String roadName = "roads_input";
	
	static double w_pos = 0.0001;     // weight for point displacements
	static double w_edge = 0.8;     // weight for edge distortions
	static double w_select = 0.1999;  // weight for unselecting objects
	
	static double epsilon = 7.5;  // [m]: threshold for a bottleneck edge representing a conflict
	static double t = 5;  // threshold for the maximally allowed stretch factor of an edge
	
	static boolean weightObjects = true;  // if true: objects obtain individual weights (here: proportional to length/area)
	static boolean withSelectionBuildings = true;   // if true: selection (of buildings) is possible, if false: only displacement
	static boolean withDisplacementRoads = false;  // if true: road nodes may be displaced
	static boolean stayWithinMap = true; // if true: the objects have to remain within the bounding box of the input points
	static boolean downweightLongBottleneckEdges = false; // if true: bottleneck edges longer than epsilon obtain smaller weights
	static boolean enforceSelectionWJ = true; // if true: the buildings are selected/unselected according to the results from the paper by Ware and Jones
	
	static double road_factor = 10;  // only meaningful if withDisplacementRoads is set to true: Specifies the factor by which the road edges are higher weighted than the bottleneck edges
	public static final double eps = Math.pow(10,-4); // tolerance [m] for point queries within the program (kd tree, point to line etc.)
	
	/**
	 * ----------------------------------------------------------------------------------------------------
	 */


	public static void main(String[] args) throws GRBException, IOException {

		BuildingNetwork bn = BuildingNetwork.importFromShapefile(buildingPath);
		RoadNetwork rn = RoadNetwork.importFromShapefile(roadPath, roadName);
		
		double scale_factor = 400/375.2;   
		rn.setKdTree(new KdTree(Main.eps));
		bn.setKdTree(new KdTree(Main.eps));
		rn.setMinX(scale_factor*rn.getMinX()+20000);
		rn.setMinY(scale_factor*rn.getMinY()+20000);
		rn.setMaxX(scale_factor*rn.getMaxX()+20000);
		rn.setMaxY(scale_factor*rn.getMaxY()+20000);
		bn.setMinX(scale_factor*bn.getMinX()+20000);
		bn.setMinY(scale_factor*bn.getMinY()+20000);
		bn.setMaxX(scale_factor*bn.getMaxX()+20000);
		bn.setMaxY(scale_factor*bn.getMaxY()+20000);
		for(RoadPoint rp : rn.getRoadPoints()) {
			rp.setX(scale_factor*rp.getX()+20000);
			rp.setY(scale_factor*rp.getY()+20000);
			rp.setOriginalCoord(new Coordinate(rp.getX(),rp.getY()));
			KdNode query = rn.getKdTree().query(rp.getOriginalCoord());
			if(query==null) {
				rn.getKdTree().insert(rp.getCoordinate(),rp);
			}
		}
		for(BuildingPoint bp : bn.getBuildingPoints()) {
			bp.setX(scale_factor*bp.getX()+20000);
			bp.setY(scale_factor*bp.getY()+20000);
			bp.setOriginalCoord(new Coordinate(bp.getX(),bp.getY()));
			KdNode query = bn.getKdTree().query(bp.getOriginalCoord());
			if(query==null) {
				bn.getKdTree().insert(bp.getCoordinate(),bp);
			}
		}
		
		GeneralNetworkWareJones gn = new GeneralNetworkWareJones(rn, bn);
		SimpleWeightedGraph<Point, Segment> G = gn.getGraph();
		System.out.println("Reading input done.");
		
		//visualization
		GeometryFactory gf = new GeometryFactory();
		OwnMapPanel myMapPanel = new OwnMapPanel();
		OwnMapPanel myMapPanelBefore = new OwnMapPanel();
		OwnMapPanel myMapPanelAfter = new OwnMapPanel();
		OwnMapPanel myMapPanelWithDelaunay = new OwnMapPanel();

		ListLayer polygon_layer_before = new ListLayer(new BasicSymbolFactory(Color.BLACK,Color.GRAY));
		for (Building b : bn.getBuildings()) {
			Polygon poly = b.getPolygon();
			Feature f = new Feature(poly);
			polygon_layer_before.add(f);
		}
		myMapPanelBefore.getMap().addLayer(polygon_layer_before, 1);

		ListLayer input_layer = new ListLayer(new BasicSymbolFactory(Color.LIGHT_GRAY,Color.LIGHT_GRAY,(float) 2.5));
		for (Segment s : G.edgeSet()) {
			Coordinate[] coords = new Coordinate[2];
			coords[0] = new Coordinate(s.getStart().getX(), s.getStart().getY());
			coords[1] = new Coordinate(s.getEnd().getX(), s.getEnd().getY());
			LineString ls = gf.createLineString(coords);
			Feature f = new Feature(ls);
			input_layer.add(f);
		}     
		myMapPanel.getMap().addLayer(input_layer, 1);
		myMapPanelBefore.getMap().addLayer(input_layer, 2);
		myMapPanelWithDelaunay.getMap().addLayer(input_layer,1);
		
		gn.computeDelaunayGraph();
		SimpleWeightedGraph<Point, Segment> T = gn.getDelaunayGraph();
		SimpleWeightedGraph<Point, Segment> G_prime = gn.getGraph();
		System.out.println("Additional Steiner points added to the graph." + "\n");	

		double M = gn.determineM();
		SolverWareJones solver = new SolverWareJones(G, G_prime, T, gn, epsilon, t, 
				M, w_pos, w_edge, w_select, road_factor);
		solver.addBottleneckEdges();
		solver.removeUnusedSteinerPoints();
		SimpleWeightedGraph<Point, Segment> G_prime_copy = new 
				SimpleWeightedGraph<Point, Segment>(Segment.class);
		for(Point p : G_prime.vertexSet()) G_prime_copy.addVertex(p);
		for(Segment se : G_prime.edgeSet()) {
			if(!se.isBottleneckEdge())  G_prime_copy.addEdge(se.getStart(), se.getEnd(), se);
		}
		SimpleWeightedGraph<Point, Segment> bEdges_initial = solver.getBottleneckGraph();
		solver.solve(weightObjects, withSelectionBuildings, withDisplacementRoads, 
			         stayWithinMap, enforceSelectionWJ, downweightLongBottleneckEdges);
		double obj = solver.evaluate(G_prime_copy, bEdges_initial.edgeSet(),bn.getBuildings(), 
				downweightLongBottleneckEdges, withDisplacementRoads);
		System.out.println("Objective value: " + obj);
		
		// Solution graph as an additional layer
		SimpleWeightedGraph<Point, Segment> solGraph = solver.getSolutionGraph();      
		ListLayer sol_layer = new ListLayer(new BasicSymbolFactory(Color.BLACK,null,(float) 1.5));
		ListLayer sol_layer_2 = new ListLayer(new BasicSymbolFactory(Color.BLACK,null,1));
		for (Segment s : solGraph.edgeSet()) {
			Coordinate[] coords = new Coordinate[2];
			coords[0] = new Coordinate(s.getStart().getX(), s.getStart().getY());
			coords[1] = new Coordinate(s.getEnd().getX(), s.getEnd().getY());
			LineString ls = gf.createLineString(coords);
			Feature f = new Feature(ls);
			sol_layer.add(f);
			sol_layer_2.add(f);
		}
		myMapPanel.getMap().addLayer(sol_layer, 3);

		// bottleneck edges as an additional layer
		ListLayer bEdges_short_layer = new ListLayer(new OwnBasicSymbolFactory(Color.RED,null,(float) 1));
		ListLayer bEdges_long_layer = new ListLayer(new OwnBasicSymbolFactory(Color.BLUE,null,(float) 1 ));
		for (Segment s : bEdges_initial.edgeSet()) {
			// ---- bottleneck edges connecting two nodes of the same building are not displayed since
			// they are ignored in the optimization
			if(s.getStart().getClass().toString().equals("class BuildingDisplacement.BuildingPoint")
				&& s.getEnd().getClass().toString().equals("class BuildingDisplacement.BuildingPoint")) {
				if(BuildingPoint.getCommonBuildings((BuildingPoint) s.getStart(), (BuildingPoint) s.getEnd()).size()>0) {
					continue;
				}
			}
			// ---------------
			Coordinate[] coords = new Coordinate[2];
			coords[0] = new Coordinate(s.getStart().getX(), s.getStart().getY());
			coords[1] = new Coordinate(s.getEnd().getX(), s.getEnd().getY());
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
		
		// solution graph with the additional Steiner nodes     
		ListLayer sol_layer_del = new ListLayer(new BasicSymbolFactory(Color.BLACK,null,(float) 1.5));
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
			}
		}
		myMapPanelWithDelaunay.getMap().addLayer(sol_layer_del, 2);
		myMapPanelWithDelaunay.getMap().addLayer(bEdges_short_layer, 3);
		myMapPanelWithDelaunay.getMap().addLayer(bEdges_long_layer, 4);


		// buildings displayed as filled polygons 
		gn.removeDelaunayPoints();
		solver.updateNetworks();
		ListLayer polygon_layer = new ListLayer(new BasicSymbolFactory(Color.BLACK,Color.GRAY));
		for (Building b : bn.getBuildings()) {
			if(b.isSelected()) {
				Polygon poly = b.getPolygon();
				Feature f = new Feature(poly);
				polygon_layer.add(f);
			} 
		}
		myMapPanelAfter.getMap().addLayer(polygon_layer, 1);
		myMapPanelAfter.getMap().addLayer(sol_layer_2, 2);

		MapFrame myFrameBefore = new MapFrame("Before optimization",myMapPanelBefore);
		myMapPanelBefore.getMap().setFrameRatio(0.01);
		myFrameBefore.setVisible(true);
		myFrameBefore.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

		MapFrame myFrame = new MapFrame("Before and after optimization and bottleneck edges (red/blue)",myMapPanel);
		myMapPanel.getMap().setFrameRatio(0.01);
		myFrame.setVisible(true);
		myFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		
		MapFrame myFrameWithDelaunay = new MapFrame("With additional Delaunay points",myMapPanelWithDelaunay);
		myMapPanelWithDelaunay.getMap().setFrameRatio(0.01);
		myFrameWithDelaunay.setVisible(true);
		myFrameWithDelaunay.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

		MapFrame myFrameAfter = new MapFrame("After optimization",myMapPanelAfter);
		myMapPanelAfter.getMap().setFrameRatio(0.01);
		myFrameAfter.setVisible(true);
		myFrameAfter.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		
		
		
		
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

	}


	public static ArrayList<Color> getColorGradient(){
		ArrayList<Color> colorGradient = new ArrayList<Color>();
		colorGradient.add(new Color(142,142,142));
		colorGradient.add(new Color(219,136,136));
		colorGradient.add(new Color(247,74,74));
		colorGradient.add(new Color(255,0,0));
		colorGradient.add(new Color(161,5,5));
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

