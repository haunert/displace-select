package visualization;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.GridLayout;
import java.awt.RenderingHints;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.geom.AffineTransform;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

import viewer.base.MapPanel;

public class MapFrame extends JFrame{

	private JPanel screenshotPanel;
	private JPanel bbox;
	private JButton printBBox;
	private JLabel bboxLabel;
	private JTextField bboxField;
	private JButton bboxButton;
	private JButton screenshotButton;
	
	public MapFrame(String title, OwnMapPanel myMapPanel) {
		super(title);
		screenshotPanel = new JPanel();
		bbox = new JPanel();
		printBBox = new JButton("Print current map extent");
		printBBox.addMouseListener(new MouseAdapter(){
			@Override
			public void mouseClicked(MouseEvent e) {
				double xMin = myMapPanel.getMap().getMyTransformation().getX(0);
				double xMax = myMapPanel.getMap().getMyTransformation().getX(myMapPanel.getMap().getWidth());
				double yMin = myMapPanel.getMap().getMyTransformation().getY(myMapPanel.getMap().getHeight());
				double yMax = myMapPanel.getMap().getMyTransformation().getY(0);
//				double xMin = myMapPanel.getMap().getxMin();
//				double yMin = myMapPanel.getMap().getyMin();
//				double xMax = myMapPanel.getMap().getxMax();
//				double yMax = myMapPanel.getMap().getyMax();
				System.out.println("Map extent (xMin,yMin,xMax,yMax): " + xMin + "," 
						+ yMin + "," + xMax + "," + yMax);
			}
		});
		bboxLabel = new JLabel("Set map extent:");
		bboxField = new JTextField("xMin,yMin,xMax,yMax");
		bboxField.setPreferredSize(new Dimension(100, 20));
		bboxButton = new JButton("OK");
		bboxButton.addMouseListener(new MouseAdapter(){
			@Override
			public void mouseClicked(MouseEvent e) {
				String input = bboxField.getText();
				if(!input.equals("xMin,yMin,xMax,yMax")) {
					String[] coordinates = input.split(",");
					double xMin = Double.parseDouble(coordinates[0]);
					double yMin = Double.parseDouble(coordinates[1]);
					double xMax = Double.parseDouble(coordinates[2]);
					double yMax = Double.parseDouble(coordinates[3]);
					myMapPanel.getMap().fitBoxToDisplay(xMin,yMin,xMax,yMax);
				}
			}
		});
		//bbox.setLayout(new GridLayout(1,4));
		JPanel printBBoxPanel = new JPanel();
		JLabel spaceholder = new JLabel();
		spaceholder.setPreferredSize(new Dimension(40, 20));
		printBBoxPanel.add(printBBox, BorderLayout.WEST);
		printBBoxPanel.add(spaceholder, BorderLayout.EAST);
		bbox.add(printBBoxPanel,BorderLayout.WEST);
		bbox.add(bboxLabel,BorderLayout.CENTER);
		bbox.add(bboxField,BorderLayout.CENTER);
		bbox.add(bboxButton,BorderLayout.EAST);
		screenshotPanel.add(bbox,BorderLayout.WEST);
		screenshotButton = new JButton("Screenshot");
		screenshotPanel.add(spaceholder, BorderLayout.EAST);
	    screenshotPanel.add(screenshotButton, BorderLayout.EAST);
		setSize(640, 480);
		add(screenshotPanel,BorderLayout.NORTH);
		screenshotButton.addMouseListener(new MouseAdapter() {
	      	@Override
	      	public void mouseClicked(MouseEvent e) {
	      		BufferedImage img = new BufferedImage(4*myMapPanel.getMap().getWidth(), 4*myMapPanel.getMap().getHeight(), BufferedImage.TYPE_INT_RGB);
	      		Graphics2D g2d = img.createGraphics();
	      		g2d.setBackground(Color.WHITE);
	      		g2d.clearRect(0, 0, 4*myMapPanel.getMap().getWidth(), 4*myMapPanel.getMap().getHeight());
	      		g2d.setRenderingHint(RenderingHints.KEY_ALPHA_INTERPOLATION, RenderingHints.VALUE_ALPHA_INTERPOLATION_QUALITY);
	      		g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
	      		g2d.setRenderingHint(RenderingHints.KEY_COLOR_RENDERING, RenderingHints.VALUE_COLOR_RENDER_QUALITY);
	      		g2d.setRenderingHint(RenderingHints.KEY_DITHERING, RenderingHints.VALUE_DITHER_ENABLE);
	      		g2d.setRenderingHint(RenderingHints.KEY_FRACTIONALMETRICS, RenderingHints.VALUE_FRACTIONALMETRICS_ON);
	      		g2d.setRenderingHint(RenderingHints.KEY_INTERPOLATION, RenderingHints.VALUE_INTERPOLATION_BILINEAR);
	      		g2d.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
	      		g2d.setRenderingHint(RenderingHints.KEY_STROKE_CONTROL, RenderingHints.VALUE_STROKE_PURE);
	      		g2d.setTransform(AffineTransform.getScaleInstance(4, 4));
	      		myMapPanel.getMap().printAll(g2d);
	      		g2d.dispose();
	      		try {
	      			ImageIO.write(img, "png", new File("screenshot.png"));
	      		} catch (IOException ex) {
	      			ex.printStackTrace();
	      		}
	      	}
	      });
		add(myMapPanel);
		
	}

}
