package visualization;

import viewer.base.*;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridLayout;
import java.util.ArrayList;

import javax.swing.BorderFactory;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.border.EtchedBorder;

import java.awt.Color;
import java.awt.ComponentOrientation;

/**
 * panel that contains a OwnMap and a panel with the coordinates of the mouse
 * cursor
 * 
 * @author haunert
 */
public class OwnMapPanel extends MapPanel {

	private static final long serialVersionUID = 6063108022820064493L;

	/**
	 * the OwnMap displayed in this panel
	 */
	private OwnMap myMap;

	/**
	 * the label displaying the x-coordinate
	 */
	private JLabel xLabel;

	/**
	 * the label displaying the y-coordinate
	 */
	private JLabel yLabel;

	/**
	 * constructor for crating an empty OwnMapPanel
	 */
	public OwnMapPanel() {

		setLayout(new BorderLayout());

		// labels showing the world coordinates of the mouse pointer
		xLabel = new JLabel();
		yLabel = new JLabel();

		myMap = new OwnMap(this);

		// when fitting the content to the OwnMap extend, use some empty space at the
		// boundary
		myMap.setFrameRatio(0.1);

		JPanel myPanel = new JPanel();
		myPanel.setLayout(new GridLayout(2, 1));
		myPanel.add(xLabel);
		myPanel.add(yLabel);
		myPanel.setBorder(BorderFactory.createTitledBorder(new EtchedBorder(), "Map Coordinates", 0, 0));

		// some GUI attributes
		myMap.fitMapToDisplay();

		myPanel.setMinimumSize(new Dimension(200, 80));
		myPanel.setPreferredSize(new Dimension(200, 80));

		this.add(BorderLayout.NORTH, myPanel);
		this.add(BorderLayout.CENTER, myMap);
	}
	
	/**
	 * constructor for creating an empty OwnMapPanel
	 */
	public OwnMapPanel(ArrayList<Color> colorGradient, double maxValue, String legendText) {

		setLayout(new BorderLayout());

		// labels showing the world coordinates of the mouse pointer
		xLabel = new JLabel();
		yLabel = new JLabel();

		myMap = new OwnMap(this);

		// when fitting the content to the OwnMap extend, use some empty space at the
		// boundary
		myMap.setFrameRatio(0.1);

		JPanel myPanel = new JPanel();
		myPanel.setLayout(new GridLayout(1, 2));
		
		JPanel coordinatePanel = new JPanel();
		coordinatePanel.setLayout(new GridLayout(2,1));
		coordinatePanel.add(xLabel);
		coordinatePanel.add(yLabel);
		coordinatePanel.setBorder(BorderFactory.createTitledBorder(new EtchedBorder(), "Map Coordinates", 0, 0));

		JPanel legendPanel = new JPanel();
		int numRows;
		if(colorGradient.size()%2==0) {
			numRows = colorGradient.size()/2;
		}
		else {
			numRows = colorGradient.size()/2+1;
		}
		legendPanel.setLayout(new GridLayout(numRows,2));
		double value1 = 0;
		for(int i = 1; i < colorGradient.size(); i++) {
			double value2 = ((double) i)/(colorGradient.size()-1) * maxValue;
			value2 = Math.round(value2*10.0)/10.0;
			String text = String.valueOf(value1) + " m to " + String.valueOf(value2) + " m";
			JLabel textLabel = new JLabel(text);
			JLabel colorLabel = new JLabel();
			colorLabel.setBackground(colorGradient.get(i-1));
			colorLabel.setOpaque(true);
			JPanel textPlusColor = new JPanel();
			GridLayout gl = new GridLayout(1,2);
			textPlusColor.setLayout(gl);
			textPlusColor.add(textLabel);			
			JPanel colorPanel = new JPanel();
			colorPanel.setLayout(new GridLayout(1,4));
			colorPanel.add(new JLabel());
			colorPanel.add(colorLabel);
			colorPanel.add(new JLabel());
			colorPanel.add(new JLabel());
			textPlusColor.add(colorPanel);
			legendPanel.add(textPlusColor);
			value1 = value2;
		}
		JLabel textLabel = new JLabel("> " + String.valueOf(value1) + " m");
		JLabel colorLabel = new JLabel();
		colorLabel.setBackground(colorGradient.get(colorGradient.size()-1));
		colorLabel.setOpaque(true);
		colorLabel.setPreferredSize(new Dimension(15,15));
		JPanel textPlusColor = new JPanel();
		GridLayout gl = new GridLayout(1,2);
		textPlusColor.setLayout(gl);
		textPlusColor.add(textLabel);			
		JPanel colorPanel = new JPanel();
		colorPanel.setLayout(new GridLayout(1,4));
		colorPanel.add(new JLabel());
		colorPanel.add(colorLabel);
		colorPanel.add(new JLabel());
		colorPanel.add(new JLabel());
		textPlusColor.add(colorPanel);
		legendPanel.add(textPlusColor);
		legendPanel.setBorder(BorderFactory.createTitledBorder(new EtchedBorder(), legendText, 0, 0));
		
		myPanel.add(coordinatePanel);
		myPanel.add(legendPanel);
		
		// some GUI attributes
		myMap.fitMapToDisplay();

		myPanel.setMinimumSize(new Dimension(200, 80));
		myPanel.setPreferredSize(new Dimension(200, 80));

		this.add(BorderLayout.NORTH, myPanel);
		this.add(BorderLayout.CENTER, myMap);
	}
	

	/**
	 * setter method used to update the coordinates displayed as text
	 * 
	 * @param x: the x-coordinate that is displayed
	 * @param y: the y-coordinate that is displayed
	 */
	public void setXYLabelText(double x, double y) {
		xLabel.setText(" x = " + x);
		xLabel.repaint();
		yLabel.setText(" y = " + y);
		yLabel.repaint();
	}

	/**
	 * the OwnMap displayed in this panel
	 * 
	 * @return the OwnMap
	 */
	public OwnMap getMap() {
		return myMap;
	}
}
