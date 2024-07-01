package visualization;

import java.awt.BasicStroke;
import java.awt.Color;

import io.structures.Feature;
import viewer.symbols.LineSymbol;
import viewer.symbols.PointSymbol;
import viewer.symbols.PolygonSymbol;
import viewer.symbols.Symbol;
import viewer.symbols.SymbolFactory;

/***
 * For dashed lines und points of a user-defined size (see method
 * createSymbol)
 */
public class OwnBasicSymbolFactory implements SymbolFactory {

	private Color strokeColor;
	private Color fillColor;
	private float strokeWidth;

	public OwnBasicSymbolFactory(Color strokeColor, Color fillColor) {
		this.strokeColor = strokeColor;
		this.fillColor = fillColor;
		this.strokeWidth = 1;
	}
	
	public OwnBasicSymbolFactory(Color strokeColor, Color fillColor, float strokeWidth) {
		this.strokeColor = strokeColor;
		this.fillColor = fillColor;
		this.strokeWidth = strokeWidth;
	}

	@Override
	public Symbol createSymbol(Feature feature) {
//		if (feature.getGeometryType().startsWith("Multi"))
//			System.err.println("WARNING! Feature has geometry of Multi-type (" + feature.getGeometryType()
//					+ ") only first feature is drawn!");

		switch (feature.getGeometryType()) {
		case "Point": return new PointSymbol(feature, strokeColor, (int) strokeWidth);
		case "MultiPoint":
			return new PointSymbol(feature, strokeColor, (int) strokeWidth);

		case "LineString":
		case "MultiLineString":
			return new LineSymbol(feature, strokeColor, new BasicStroke(strokeWidth,
					BasicStroke.CAP_BUTT, BasicStroke.JOIN_BEVEL, 1,
					new float[] { 4, 4 }, 0));

		case "Polygon":
		case "MultiPolygon":
			return new PolygonSymbol(feature, strokeColor, new BasicStroke(strokeWidth), fillColor);

		default:
			System.err.println("Unknown geometry type! (" + feature.getGeometryType() + ")");
			return null;
		}
	}
}
