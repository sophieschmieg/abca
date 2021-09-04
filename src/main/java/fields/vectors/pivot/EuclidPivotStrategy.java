package fields.vectors.pivot;

import java.math.BigInteger;
import java.util.List;

import fields.interfaces.Element;
import fields.interfaces.Ring;

public class EuclidPivotStrategy<T extends Element<T>> implements PivotStrategy<T> {
	private Ring<T> ring;

	public EuclidPivotStrategy(Ring<T> ring) {
		this.ring = ring;
	}
	
	@Override
	public int pivot(List<T> elements) {
		BigInteger minValue = null;
		int pivot = -1;
		for (int i = 0; i < elements.size(); i++) {
			if (elements.get(i).equals(ring.zero())) {
				continue;
			}
			BigInteger value = ring.euclidMeasure(elements.get(i));
			if (minValue == null || value.compareTo(minValue) < 0) {
				pivot = i;
				minValue = value;
			}
		}
		return pivot;
	}

}
