package DelTri;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.triangulate.Segment;

/**
 * An interface for factories which create a {@link ConstraintVertex}
 * 
 * @author Martin Davis
 */
public interface OwnConstraintVertexFactory {
    OwnConstraintVertex createVertex(Coordinate p, Segment constraintSeg);
}
