package RoadDisplacement;


/**
 * A class which represents a pair of adjacent roads (= two roads which have a common node)
 */
public class RoadPair {
	
	private Road r1;
	private Road r2;

	public RoadPair(Road r1, Road r2) {
		this.r1 = r1;
		this.r2 = r2;
	}

	public Road getR1() {
		return r1;
	}

	public Road getR2() {
		return r2;
	}
	
	public Road getOtherRoad(Road r) {
		if(r==r1) {
			return r2;
		}
		if(r==r2) {
			return r1;
		}
		System.out.println("Error in class RoadPair, method getOtherRoad!");
		return null;
	}

}
