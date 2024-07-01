package BuildingDisplacement;

/**
 * A class representing a pair of two adjacent buildings (i.e. which share at least one edge).
 */
public class BuildingPair {
	
	private Building b1;
	private Building b2;

	public BuildingPair(Building b1, Building b2) {
		this.b1 = b1;
		this.b2 = b2;
	}

	public Building getB1() {
		return b1;
	}

	public Building getB2() {
		return b2;
	}
	
	public Building getOtherBuilding(Building b) {
		if(b==b1) {
			return b2;
		}
		if(b==b2) {
			return b1;
		}
		System.out.println("Error in class BuildingPair, method getOtherBuilding!");
		return null;
	}

}
